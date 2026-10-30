#= /*
* This file is part of OpenModelica.
*
* Copyright (c) 1998-2026, Open Source Modelica Consortium (OSMC),
* c/o Linköpings universitet, Department of Computer and Information Science,
* SE-58183 Linköping, Sweden.
*
* All rights reserved.
*
* THIS PROGRAM IS PROVIDED UNDER THE TERMS OF AGPL VERSION 3 LICENSE OR
* THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.8.
* ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
* RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GNU AGPL
* VERSION 3, ACCORDING TO RECIPIENTS CHOICE.
*
* The OpenModelica software and the OSMC (Open Source Modelica Consortium)
* Public License (OSMC-PL) are obtained from OSMC, either from the above
* address, from the URLs:
* http://www.openmodelica.org or
* https://github.com/OpenModelica/ or
* http://www.ida.liu.se/projects/OpenModelica,
* and in the OpenModelica distribution.
*
* GNU AGPL version 3 is obtained from:
* https://www.gnu.org/licenses/licenses.html#GPL
*
* This program is distributed WITHOUT ANY WARRANTY; without
* even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
* IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
*
* See the full OSMC Public License conditions for more details.
*
*/ =#

#= Traversal of SimCode expressions, statements and equations, by `Traversal` (IR/Traversal.jl):
   what a node's field holds is declared here, per type; every traversal follows from it. =#

"Expressions and statements: a traversal's children, handed to its `recur`."
const TraversalChild = Union{Exp, Statement}

"Iterators, branches, when-clauses and equations: they hold children, traversed field by field."
const TraversalContainer = Union{Iterator, Branch, WHEN_STMTS, Equation, Construct}

#= Anything else (types, references, sources, attributes) is a leaf, kept. =#
Traversal.traverseField(x::TraversalChild, recur::R, arg::A) where {R, A} = recur(x, arg)
Traversal.traverseField(x::TraversalContainer, recur::R, arg::A) where {R, A} = traverseFields(x, recur, arg)
Traversal.traverseField(xs::Vector{<:Union{TraversalChild, TraversalContainer}}, recur::R, arg::A) where {R, A} =
  traverseElements(xs, recur, arg)
Traversal.traverseField(xss::Vector{<:Vector{<:Union{TraversalChild, TraversalContainer}}}, recur::R, arg::A) where {R, A} =
  traverseElements(xss, recur, arg)

"""
    children(e)

The direct subexpressions of `e`, in field order (a reduction's iterator ranges after its body).
"""
children(e::Exp)::Vector{Exp} = last(traverseFields(e, _pushChild, Exp[]))

_pushChild(x::Exp, out::Vector{Exp})::Tuple{Exp, Vector{Exp}} = (push!(out, x); (x, out))

"""
    mapChildren(g, e)

`e` with `g` applied to each of its direct subexpressions, one level; `e` itself if none changed.
"""
mapChildren(g::G, e::Exp) where {G} = first(traverseFields(e, _Apply(g), nothing))

"`g` as a traversal's `recur`: a child replaced by `g(child)`."
struct _Apply{G}
  g::G
end
(m::_Apply)(x::Exp, a::Nothing) = (m.g(x), a)

"`mapExp` of `f` as a callable, for `mapChildren`."
struct _BottomUp{F}
  f::F
end
(m::_BottomUp)(e::Exp) = mapExp(m.f, e)

"""
    mapExp(f, e)

`e` with every subexpression `s` replaced by `f(s)`, bottom-up: `f` sees a node whose children
are already mapped.
"""
mapExp(f::F, e::Exp) where {F} = f(mapChildren(_BottomUp(f), e))

"""
    foreachExp(f, e)

Call `f` on `e` and on every subexpression, top-down.
"""
function foreachExp(f::F, e::Exp)::Nothing where {F}
  f(e)
  traverseFields(e, _Visit(f), nothing)
  return nothing
end

"`foreachExp` of `f` as a traversal's `recur`."
struct _Visit{F}
  f::F
end
(v::_Visit)(x::Exp, a::Nothing) = (foreachExp(v.f, x); (x, a))

"Whether `p` holds for `e` or one of its subexpressions."
function anyExp(p::F, e::Exp)::Bool where {F}
  p(e) && return true
  return last(traverseFields(e, _Any(p), false))
end

"`anyExp` of `p` as a traversal's `recur`: the children after the first that holds are not visited."
struct _Any{F}
  p::F
end
(s::_Any)(x::Exp, found::Bool) = (x, found || anyExp(s.p, x))

"""
    mapStatementExps(g, s)

`s` with each of its expressions `e` (of nested statements too) replaced by `g(e)`, whole;
the targets of `reinit` and of OM.jl's structural changes (references, not expressions) kept.
"""
mapStatementExps(g::G, s::Statement) where {G} = first(traverseFields(s, _StatementExps(g), nothing))

mapStatementExps(g::G, body::Vector{Statement}) where {G} = Statement[mapStatementExps(g, s) for s in body]

"`g` as a traversal's `recur` over a statement: an expression replaced by `g` of it, a nested statement mapped in turn."
struct _StatementExps{G}
  g::G
end
(m::_StatementExps)(e::Exp, a::Nothing) = (m.g(e), a)
(m::_StatementExps)(s::Statement, a::Nothing) = (mapStatementExps(m.g, s), a)

"`s` with `f` mapped (as by `mapExp`) over its expressions; the targets of `reinit` and of OM.jl's structural changes kept."
mapStatement(f::F, s::Statement) where {F} = mapStatementExps(e -> mapExp(f, e), s)

mapStatements(f::F, body::Vector{Statement}) where {F} = Statement[mapStatement(f, s) for s in body]

"""
    crefNames!(names, e)

The names of the variables `e` reads, in order: a reference's name (with its subscripts,
`x[1]` of an array element), a derivative's state's, a `pre(x)`'s variable;
an array read by literal subscripts also as its element (`R_T[1][2]`: `crefName` of its
`elementRead`). An array read by a subscript that reads a variable is the array's read only; with
`dynamic`, that array (its `arrayBase`) recorded there.
"""
function crefNames!(names::OrderedSet{String}, e::Exp; dynamic::Union{Nothing, Set{Symbol}} = nothing)::OrderedSet{String}
  foreachRead(e) do c::SimCref, index::Union{Nothing, Vector{Int}}, isDynamic::Bool
    if isDynamic
      (dynamic !== nothing && c.ref !== nothing) && push!(dynamic, arrayBase(c.ref))
    else
      push!(names, index === nothing ? crefName(c) : crefName(c, index))
    end
  end
  return names
end

"""
    foreachRead(f, e)

`f(c, index, dynamic)` for each reference `e` reads, in order, the one traversal of `crefNames!`, `crefReads!` and
`dynamicArrays!`: a reference, a derivative's state, a `pre(x)`'s variable (`index` nothing); an array read by literal
subscripts, `c` the array and `index` the subscripts (its element's name `crefName(c, index)`, its reference
`elementRead(c, index)`, made only where asked); by a subscript that reads a variable, the array, `dynamic`.
"""
function foreachRead(f::F, e::Exp)::Nothing where {F}
  if e isa Union{EXP_CREF, DER, PRE}
    f(e.cref, nothing, false)
  elseif e isa ASUB && e.exp isa EXP_CREF
    if all(s::Exp -> s isa ICONST, e.subs)
      f(e.exp.cref, Int[(s::ICONST).value for s in e.subs], false)
    elseif any(s::Exp -> anyExp(x::Exp -> x isa Union{EXP_CREF, DER, PRE}, s), e.subs)
      f(e.exp.cref, nothing, true)
    end
  end
  for c in children(e)
    foreachRead(f, c)
  end
  return nothing
end

"""
    Read(cref, dynamic)

A reference an expression reads (`crefReads!`): a variable's, an array element's by literal subscripts (`elementRead`),
or, `dynamic`, an array read by a subscript that reads a variable, whose element no subscript names.
"""
struct Read
  cref::SimCref
  dynamic::Bool
end
Base.:(==)(a::Read, b::Read)::Bool = a.cref == b.cref && a.dynamic == b.dynamic
Base.isequal(a::Read, b::Read)::Bool = a == b
Base.hash(r::Read, h::UInt)::UInt = hash(r.dynamic, hash(r.cref, h))

"""
    crefReads!(reads, e) -> Vector{Read}

The references `e` reads, in order (`foreachRead`), an element's by literal subscripts as its reference (`elementRead`).
"""
function crefReads!(reads::Vector{Read}, e::Exp)::Vector{Read}
  foreachRead(e) do c::SimCref, index::Union{Nothing, Vector{Int}}, isDynamic::Bool
    push!(reads, Read(index === nothing ? c : elementRead(c, index), isDynamic))
  end
  return reads
end

"The arrays `e` reads by a subscript that reads a variable (`foreachRead`'s dynamic reads, by `arrayBase`), pushed onto `dynamic`."
function dynamicArrays!(dynamic::Set{Symbol}, e::Exp)::Set{Symbol}
  foreachRead(e) do c::SimCref, _::Union{Nothing, Vector{Int}}, isDynamic::Bool
    (isDynamic && c.ref !== nothing) && push!(dynamic, arrayBase(c.ref))
  end
  return dynamic
end

"""
    elementRead(c, index) -> SimCref

The element `index` (literal subscripts) of the array reference `c`: its key and subscripts
followed by `index`, its IR reference the array's subscripted on its last part (`x[1,2]` of `x`).
"""
function elementRead(c::SimCref, index::Vector{Int})::SimCref
  local ref::Union{Nothing, IR.ComponentRef} = c.ref === nothing ? nothing :
    IR.appendSubscripts(c.ref, IR.Subscript[IR.INDEX(IR.INTEGER(k)) for k in index])
  return SimCref(c.sym, Int[c.subs; index], ref)
end

"""
    statementExps!(out, body) -> Vector{Exp}

The expressions of the statements `body`: their conditions, sides and values, a `reinit`'s state as a reference, a
loop's range; nested bodies' too; pushed onto `out`.
"""
statementExps!(out::Vector{Exp}, body::Vector{Statement})::Vector{Exp} = (foreach(s::Statement -> statementExps!(out, s), body); out)
statementExps!(out::Vector{Exp}, s::ASSIGN)::Vector{Exp} = push!(out, s.left, s.right)
statementExps!(out::Vector{Exp}, s::REINIT)::Vector{Exp} = push!(out, EXP_CREF(s.stateVar, TYPE_REAL()), s.value)
statementExps!(out::Vector{Exp}, s::NORETCALL)::Vector{Exp} = push!(out, s.exp)
statementExps!(out::Vector{Exp}, s::ASSERT)::Vector{Exp} = push!(out, s.condition, s.message, s.level)
statementExps!(out::Vector{Exp}, s::TERMINATE)::Vector{Exp} = push!(out, s.message)
function statementExps!(out::Vector{Exp}, s::Union{IF, WHEN})::Vector{Exp}
  for b in s.branches
    push!(out, b.condition)
    statementExps!(out, b.body)
  end
  s isa IF && statementExps!(out, s.elseBody)
  return out
end
statementExps!(out::Vector{Exp}, s::FOR)::Vector{Exp} = (push!(out, s.iterator.range); statementExps!(out, s.body))
statementExps!(out::Vector{Exp}, s::WHILE)::Vector{Exp} = (push!(out, s.condition); statementExps!(out, s.body))
statementExps!(out::Vector{Exp}, s::RECOMPILATION)::Vector{Exp} = push!(out, EXP_CREF(s.componentToChange, TYPE_REAL()), s.newValue)
statementExps!(out::Vector{Exp}, ::Statement)::Vector{Exp} = out

"A reference's name: its symbol, with its Integer subscripts (`A[1][2]`)."
crefName(c::SimCref)::String = isempty(c.subs) ? String(c.sym) : string(c.sym, join(("[$(i)]" for i in c.subs)))

"The name of the element `index` of the array reference `c` (`crefName` of `elementRead`), without making its reference."
crefName(c::SimCref, index::Vector{Int})::String = string(crefName(c), join(("[$(i)]" for i in index)))

#= ==================================================================================== =#
#= Traversals with an accumulator (the passes')                                         =#
#= ==================================================================================== =#

"""
    traverseExpTopDown(e, visitor, arg)

`e` traversed top-down by `visitor(x, arg) -> (newX, descend, newArg)`: `newX` replaces `x`,
its children are traversed in turn if `descend`; `arg` is threaded through. A node whose
children did not change is returned itself (`===`).
"""
function traverseExpTopDown(e::Exp, visitor::F, arg::T)::Tuple{Exp, T} where {F, T}
  local (out, descend, a) = visitor(e, arg)
  descend || return (out, a)
  return traverseFields(out, (x, acc) -> traverseExpTopDown(x, visitor, acc), a)
end

"""
    traverseExpBottomUp(e, visitor, arg)

`e` traversed bottom-up by `visitor(x, arg) -> (newX, newArg)`: `x`'s children are traversed
first.
"""
function traverseExpBottomUp(e::Exp, visitor::F, arg::T)::Tuple{Exp, T} where {F, T}
  local (out, a) = traverseFields(e, (x, acc) -> traverseExpBottomUp(x, visitor, acc), arg)
  return visitor(out, a)
end
