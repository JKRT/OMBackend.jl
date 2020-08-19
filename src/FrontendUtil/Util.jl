module Util

import ..DAE
import DoubleEnded

using MetaModelica

const Type_a = Any
const Argument = Any

function traverseExpTopDown(inExp::DAE.Exp, func::Function, ext_arg::Type_a) ::Tuple{DAE.Exp, Type_a}
  local outArg::Type_a
  local outExp::DAE.Exp
  local cont::Bool
  (outExp, cont, outArg) = func(inExp, ext_arg)
  (outExp, outArg) = traverseExpTopDown1(cont, outExp, func, outArg)
  (outExp, outArg)
end

function traverseExpTopDown1(continueTraversal::Bool, inExp::DAE.Exp, func::Function, inArg::Type_a) ::Tuple{DAE.Exp, Type_a}
  local outArg::Type_a
  local outExp::DAE.Exp
  (outExp, outArg) = begin
    local aliases::List{List{String}}
    local attr::DAE.CallAttributes
    local cr::DAE.ComponentRef
    local cr_1::DAE.ComponentRef
    local dim::ModelicaInteger
    local e1::DAE.Exp
    local e1_1::DAE.Exp
    local e2::DAE.Exp
    local e2_1::DAE.Exp
    local e3::DAE.Exp
    local e3_1::DAE.Exp
    local e4::DAE.Exp
    local e4_1::DAE.Exp
    local e::DAE.Exp
    local et::Type
    local expl::List{DAE.Exp}
    local expl_1::List{DAE.Exp}
    local ext_arg::Type_a
    local ext_arg_1::Type_a
    local ext_arg_2::Type_a
    local ext_arg_3::Type_a
    local fieldNames::List{String}
    local fn::Absyn.Path
    local i::ModelicaInteger
    local index_::ModelicaInteger
    local isExpisASUB::Option{Tuple{DAE.Exp, ModelicaInteger, ModelicaInteger}}
    local lstexpl::List{List{DAE.Exp}}
    local lstexpl_1::List{List{DAE.Exp}}
    local op::Operator
    local reductionInfo::DAE.ReductionInfo
    local rel::Function
    local riters::DAE.ReductionIterators
    local scalar::Bool
    local t::Type
    local tp::Type

    if !continueTraversal
      return (inExp, inArg)
    end

    @match (inExp, func, inArg) begin
      (DAE.ICONST(_), _, ext_arg)  => begin
        (inExp, ext_arg)
      end

      (DAE.RCONST(_), _, ext_arg)  => begin
        (inExp, ext_arg)
      end

      (DAE.SCONST(_), _, ext_arg)  => begin
        (inExp, ext_arg)
      end

      (DAE.BCONST(_), _, ext_arg)  => begin
        (inExp, ext_arg)
      end

      (DAE.ENUM_LITERAL(__), _, ext_arg)  => begin
        (inExp, ext_arg)
      end

      (DAE.CREF(cr, tp), rel, ext_arg)  => begin
        (cr_1, ext_arg_1) = traverseExpTopDownCrefHelper(cr, rel, ext_arg)
        (if referenceEq(cr, cr_1)
           inExp
         else
           DAE.CREF(cr_1, tp)
         end, ext_arg_1)
      end

      (DAE.UNARY(operator = op, exp = e1), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (if referenceEq(e1, e1_1)
           inExp
         else
           DAE.UNARY(op, e1_1)
         end, ext_arg_1)
      end

      (DAE.BINARY(exp1 = e1, operator = op, exp2 = e2), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
        (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
         inExp
         else
         DAE.BINARY(e1_1, op, e2_1)
         end, ext_arg_2)
      end

      (DAE.LUNARY(operator = op, exp = e1), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (if referenceEq(e1, e1_1)
         inExp
         else
         DAE.LUNARY(op, e1_1)
         end, ext_arg_1)
      end

      (DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
        (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
         inExp
         else
         DAE.LBINARY(e1_1, op, e2_1)
         end, ext_arg_2)
      end

      (DAE.RELATION(exp1 = e1, operator = op, exp2 = e2, index = index_, optionExpisASUB = isExpisASUB), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
        (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
         inExp
         else
         DAE.RELATION(e1_1, op, e2_1, index_, isExpisASUB)
         end, ext_arg_2)
      end

      (DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
        (e3_1, ext_arg_3) = traverseExpTopDown(e3, rel, ext_arg_2)
        (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1) && referenceEq(e3, e3_1)
         inExp
         else
         DAE.IFEXP(e1_1, e2_1, e3_1)
         end, ext_arg_3)
      end

      (DAE.CALL(path = fn, expLst = expl, attr = attr), rel, ext_arg)  => begin
        (expl_1, ext_arg_1) = traverseExpListTopDown(expl, rel, ext_arg)
        (DAE.CALL(fn, expl_1, attr), ext_arg_1)
      end

      (DAE.RECORD(path = fn, exps = expl, comp = fieldNames, ty = tp), rel, ext_arg)  => begin
        (expl_1, ext_arg_1) = traverseExpListTopDown(expl, rel, ext_arg)
        (DAE.RECORD(fn, expl_1, fieldNames, tp), ext_arg_1)
      end

      (DAE.PARTEVALFUNCTION(fn, expl, tp, t), rel, ext_arg)  => begin
        (expl_1, ext_arg_1) = traverseExpListTopDown(expl, rel, ext_arg)
        (DAE.PARTEVALFUNCTION(fn, expl_1, tp, t), ext_arg_1)
      end

      (DAE.ARRAY(ty = tp, scalar = scalar, array = expl), rel, ext_arg)  => begin
        (expl_1, ext_arg_1) = traverseExpListTopDown(expl, rel, ext_arg)
        (DAE.ARRAY(tp, scalar, expl_1), ext_arg_1)
      end

      (DAE.MATRIX(ty = tp, integer = dim, matrix = lstexpl), rel, ext_arg)  => begin
        (lstexpl_1, ext_arg_1) = traverseExpMatrixTopDown(lstexpl, rel, ext_arg)
        (DAE.MATRIX(tp, dim, lstexpl_1), ext_arg_1)
      end

      (DAE.RANGE(ty = tp, start = e1, step = NONE(), stop = e2), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
        (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
         inExp
         else
         DAE.RANGE(tp, e1_1, NONE(), e2_1)
         end, ext_arg_2)
      end

      (DAE.RANGE(ty = tp, start = e1, step = SOME(e2), stop = e3), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
        (e3_1, ext_arg_3) = traverseExpTopDown(e3, rel, ext_arg_2)
        (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1) && referenceEq(e3, e3_1)
         inExp
         else
         DAE.RANGE(tp, e1_1, SOME(e2_1), e3_1)
         end, ext_arg_3)
      end

      (DAE.TUPLE(PR = expl), rel, ext_arg)  => begin
        (expl_1, ext_arg_1) = traverseExpListTopDown(expl, rel, ext_arg)
        (DAE.TUPLE(expl_1), ext_arg_1)
      end

      (DAE.CAST(ty = tp, exp = e1), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (DAE.CAST(tp, e1_1), ext_arg_1)
      end

      (DAE.ASUB(exp = e1, sub = expl_1), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (expl_1, ext_arg_2) = traverseExpListTopDown(expl_1, rel, ext_arg_1)
        (makeASUB(e1_1, expl_1), ext_arg_2)
      end

      (DAE.TSUB(e1, i, tp), rel, ext_arg)  => begin
        (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
        (DAE.TSUB(e1_1, i, tp), ext_arg_1)
      end

     (e1 && DAE.RSUB(__), rel, ext_arg)  => begin
       (e1_1, ext_arg_1) = traverseExpTopDown(e1.exp, rel, ext_arg)
       if ! referenceEq(e1.exp, e1_1)
         e1.exp = e1_1
       end
     (e1, ext_arg_1)
    end

     (DAE.SIZE(exp = e1, sz = NONE()), rel, ext_arg)  => begin
       (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
       (DAE.SIZE(e1_1, NONE()), ext_arg_1)
     end

     (DAE.SIZE(exp = e1, sz = SOME(e2)), rel, ext_arg)  => begin
       (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
       (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
       (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
         inExp
        else
          DAE.SIZE(e1_1, SOME(e2_1))
       end, ext_arg_2)
     end

     (DAE.CODE(__), _, ext_arg)  => begin
       (inExp, ext_arg)
     end

     (DAE.REDUCTION(reductionInfo = reductionInfo, expr = e1, iterators = riters), rel, ext_arg)  => begin
       (e1, ext_arg) = traverseExpTopDown(e1, rel, ext_arg)
       (riters, ext_arg) = traverseReductionIteratorsTopDown(riters, rel, ext_arg)
       (DAE.REDUCTION(reductionInfo, e1, riters), ext_arg)
     end

     (DAE.EMPTY(__), _, _)  => begin
       (inExp, inArg)
     end

     (DAE.CONS(e1, e2), rel, ext_arg)  => begin
       (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
       (e2_1, ext_arg_2) = traverseExpTopDown(e2, rel, ext_arg_1)
       (if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
          inExp
        else
         DAE.CONS(e1_1, e2_1)
        end, ext_arg_2)
     end

     (DAE.LIST(expl), rel, ext_arg)  => begin
       (expl_1, ext_arg_1) = traverseExpListTopDown(expl, rel, ext_arg)
       (DAE.LIST(expl_1), ext_arg_1)
     end

     (DAE.UNBOX(e1, tp), rel, ext_arg)  => begin
       (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
       (DAE.UNBOX(e1_1, tp), ext_arg_1)
     end

     (DAE.BOX(e1), rel, ext_arg)  => begin
       (e1_1, ext_arg_1) = traverseExpTopDown(e1, rel, ext_arg)
       (DAE.BOX(e1_1), ext_arg_1)
     end

     (DAE.PATTERN(__), _, ext_arg)  => begin
       (inExp, ext_arg)
     end

     (DAE.SHARED_LITERAL(__), _, ext_arg)  => begin
      (inExp, ext_arg)
     end

     _  => begin
       throw("Error: traverseExpTopDown1 failed")
     end

   end
  end
  (outExp, outArg)
end

function traverseExpListTopDown(expLst::List{DAE.Exp}, func::Function, inArg::Type_a)
  outArg::Type_a = inArg
  for e in expLst
    (_, outArg) = traverseExpTopDown1(true, e, func, outArg)
  end
  return (expLst, outArg)
end

function traverseExpTopDownCrefHelper(inCref::DAE.ComponentRef, rel::Function, iarg::Argument) ::Tuple{DAE.ComponentRef, Argument}
  local outArg::Argument
  local outCref::DAE.ComponentRef
  (outCref, outArg) = begin
    local arg::Argument
    local cr::DAE.ComponentRef
    local cr_1::DAE.ComponentRef
    local name::String
    local subs::List{DAE.Subscript}
    local subs_1::List{DAE.Subscript}
    local ty::Type
    @match (inCref, rel, iarg) begin
      (DAE.CREF_QUAL(ident = name, identType = ty, subscriptLst = subs, componentRef = cr), _, arg)  => begin
        (subs_1, arg) = traverseExpTopDownSubs(subs, rel, arg)
        (cr_1, arg) = traverseExpTopDownCrefHelper(cr, rel, arg)
        (if referenceEq(subs, subs_1) && referenceEq(cr, cr_1)
         inCref
         else
         DAE.CREF_QUAL(name, ty, subs_1, cr_1)
         end, arg)
      end

      (DAE.CREF_IDENT(ident = name, identType = ty, subscriptLst = subs), _, arg)  => begin
        (subs_1, arg) = traverseExpTopDownSubs(subs, rel, arg)
        (if referenceEq(subs, subs_1)
         inCref
         else
         DAE.CREF_IDENT(name, ty, subs_1)
         end, arg)
      end

      (DAE.WILD(__), _, arg)  => begin
        (inCref, arg)
      end
    end
  end
  (outCref, outArg)
end

function traverseExpTopDownSubs(inSubscript::List{<:DAE.Subscript}, rel::Function, iarg::Argument) ::Tuple{List{DAE.Subscript}, Argument}
  local allEq::Bool = true
  local arg::Argument = iarg
  local delst::DoubleEnded.MutableList{DAE.Subscript}
  local exp::DAE.Exp
  local nEq::ModelicaInteger = 0
  local nsub::DAE.Subscript
  local outSubscript::List{DAE.Subscript}

  for sub in inSubscript
    nsub = begin
      @match sub begin
        DAE.WHOLEDIM(__)  => begin
          sub
        end

        DAE.SLICE(__)  => begin
          (exp, arg) = traverseExpTopDown(sub.exp, rel, arg)
          if referenceEq(sub.exp, exp)
            sub
          else
            DAE.SLICE(exp)
          end
        end

        DAE.INDEX(__)  => begin
          (exp, arg) = traverseExpTopDown(sub.exp, rel, arg)
          if referenceEq(sub.exp, exp)
            sub
          else
            DAE.INDEX(exp)
          end
        end

        DAE.WHOLE_NONEXP(__)  => begin
          (exp, arg) = traverseExpTopDown(sub.exp, rel, arg)
          if referenceEq(sub.exp, exp)
            sub
          else
            DAE.WHOLE_NONEXP(exp)
          end
        end
      end
    end
    if if allEq
      ! referenceEq(nsub, sub)
    else
      false
    end
      allEq = false
      delst = DoubleEnded.empty(nsub)
      for elt in inSubscript
        if nEq < 1
          break
        end
        DoubleEnded.push_back(delst, elt)
        nEq = nEq - 1
      end
    end
    if allEq
      nEq = nEq + 1
    else
      DoubleEnded.push_back(delst, nsub)
    end
  end
  #=  Preserve reference equality without any allocation if nothing changed =#
  outSubscript = if allEq
    inSubscript
  else
    DoubleEnded.toListAndClear(delst)
  end
  (outSubscript, arg)
end

function transposeNestedList(lstlst::List{List{T}})::List{List{T}} where{T}
  transposeNestedListAccumulator(lstlst, nil)
end

function transposeNestedListAccumulator(lstlst::List{List{T}}, acc::List{List{T}})::List{List{T}} where{T}
  local rest::List{List{T}}=nil
  local tmpLst::List{T}
  local tmp::T
  @match tmpLst <| _ = lstlst
  if listLength(tmpLst) == 0
    for lst in lstlst
      tmp <| lst = lst
      tmpLst = tmp <| tmpLst
      rest = lst <| rest
    end
    transposeNestedListAccumulator(listReverse(rest), tmpLst <| acc)
  end
  return acc
end


#= author: lochel
This function extracts all crefs from the input expression, except 'time'. =#
function getAllCrefs(inExp::DAE.Exp)::List{DAE.ComponentRef}
  local outCrefs::List{DAE.ComponentRef}
  (_, outCrefs) = traverseExpBottomUp(inExp, getAllCrefs2, nil)
  outCrefs
end

function getAllCrefs2(inExp::DAE.Exp, inCrefList::List{<:DAE.ComponentRef})::Tuple{DAE.Exp, List{DAE.ComponentRef}}
  local outCrefList::List{DAE.ComponentRef} = inCrefList
  local outExp::DAE.Exp = inExp
  local cr::DAE.ComponentRef
  if isCref(inExp)
    @match DAE.CREF(componentRef = cr) = inExp
    if ! listMember(cr, inCrefList)
      outCrefList = _cons(cr, outCrefList)
    end
  end
  (outExp, outCrefList)
end


end #=End Util=#
