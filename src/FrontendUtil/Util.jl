module Util

import DAE
import DoubleEnded

using MetaModelica

const Type_a = Any
const Argument = Any

"""
  Traverses an expression top down.
  The traversal function is expected to be on the following format
  function <name>(exp, dictonary).
  The function is expected to return a tuple of three elements.
  The first returning an expression, the second returning a boolean indicating if the traversal should continue
  and the last is the out argument.
"""
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

"""
  Traverses all subexpressions of an expression.
  Takes a function and an extra argument passed through the traversal.
  The function can potentially change the expression. In such cases,
  the changes are made bottom-up, i.e. a subexpression is traversed
  and changed before the complete expression is traversed.
  NOTE: The user-provided function is not allowed to fail! If you want to
  detect a failure, return NONE() in your user-provided datatype.
"""
function traverseExpBottomUp(inExp::DAE.Exp, inFunc::Function, inExtArg::T)  where {T}
  local outExtArg::T
  local outExp::DAE.Exp
  (outExp, outExtArg) = begin
    local e1_1::DAE.Exp
    local e::DAE.Exp
    local e1::DAE.Exp
    local e2_1::DAE.Exp
    local e2::DAE.Exp
    local e3_1::DAE.Exp
    local e3::DAE.Exp
    local e4::DAE.Exp
    local e4_1::DAE.Exp
    local ext_arg::T
    local op::Operator
    local rel::FuncExpType
    local expl_1::List{DAE.Exp}
    local expl::List{DAE.Exp}
    local fn::Absyn.Path
    local scalar::Bool
    local tp::Type
    local t::Type
    local i::Integer
    local lstexpl_1::List{List{DAE.Exp}}
    local lstexpl::List{List{DAE.Exp}}
    local dim::Integer
    local str::String
    local localDecls::List{DAE.Element}
    local fieldNames::List{String}
    local attr::DAE.CallAttributes
    local cases::List{DAE.MatchCase}
    local cases_1::List{DAE.MatchCase}
    local matchTy::DAE.MatchType
    local index_::Integer
    local isExpisASUB::Option{Tuple{DAE.Exp, Integer, Integer}}
    local reductionInfo::DAE.ReductionInfo
    local riters::DAE.ReductionIterators
    local riters_1::DAE.ReductionIterators
    local cr::DAE.ComponentRef
    local cr_1::DAE.ComponentRef
    local aliases::List{List{String}}
    local clk::DAE.ClockKind
    local clk1::DAE.ClockKind
    local typeVars::List{DAE.Type}
    @match inExp begin
      DAE.EMPTY(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.ICONST(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.RCONST(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.SCONST(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.BCONST(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.CLKCONST(clk)  => begin
        (clk1, ext_arg) = traverseExpClk(clk, inFunc, inExtArg)
        e = if referenceEq(clk1, clk)
          inExp
        else
          DAE.CLKCONST(clk1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.ENUM_LITERAL(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.CREF(cr, tp)  => begin
        (cr_1, ext_arg) = traverseExpCref(cr, inFunc, inExtArg)
        e = if referenceEq(cr, cr_1)
          inExp
        else
          DAE.CREF(cr_1, tp)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.UNARY(operator = op, exp = e1)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.UNARY(op, e1_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
          inExp
        else
          DAE.BINARY(e1_1, op, e2_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.LUNARY(operator = op, exp = e1)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.LUNARY(op, e1_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
          inExp
        else
          DAE.LBINARY(e1_1, op, e2_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2, index = index_, optionExpisASUB = isExpisASUB)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
          inExp
        else
          DAE.RELATION(e1_1, op, e2_1, index_, isExpisASUB)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        (e3_1, ext_arg) = traverseExpBottomUp(e3, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1) && referenceEq(e3, e3_1)
          inExp
        else
          DAE.IFEXP(e1_1, e2_1, e3_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.CALL(path = fn, expLst = expl, attr = attr)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.CALL(fn, expl_1, attr)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.RECORD(path = fn, exps = expl, comp = fieldNames, ty = tp)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.RECORD(fn, expl_1, fieldNames, tp)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.PARTEVALFUNCTION(fn, expl, tp, t)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.PARTEVALFUNCTION(fn, expl_1, tp, t)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.ARRAY(ty = tp, scalar = scalar, array = expl)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.ARRAY(tp, scalar, expl_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.MATRIX(ty = tp, integer = dim, matrix = lstexpl)  => begin
        (lstexpl_1, ext_arg) = traverseExpMatrix(lstexpl, inFunc, inExtArg)
        e = if referenceEq(lstexpl, lstexpl_1)
          inExp
        else
          DAE.MATRIX(tp, dim, lstexpl_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.RANGE(ty = tp, start = e1, step = NONE(), stop = e2)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
          inExp
        else
          DAE.RANGE(tp, e1_1, NONE(), e2_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.RANGE(ty = tp, start = e1, step = SOME(e2), stop = e3)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        (e3_1, ext_arg) = traverseExpBottomUp(e3, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1) && referenceEq(e3, e3_1)
          inExp
        else
          DAE.RANGE(tp, e1_1, SOME(e2_1), e3_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.TUPLE(PR = expl)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.TUPLE(expl_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.CAST(ty = tp, exp = e1)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.CAST(tp, e1_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.ASUB(exp = e1, sub = expl)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(expl, expl_1)
          inExp
        else
          makeASUB(e1_1, expl_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.TSUB(e1, i, tp)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.TSUB(e1_1, i, tp)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      e1 && DAE.RSUB(__)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1.exp, inFunc, inExtArg)
        if ! referenceEq(e1.exp, e1_1)
          e1.exp = e1_1
        end
        (e1, ext_arg) = inFunc(e1, ext_arg)
        (e1, ext_arg)
      end
      
      DAE.SIZE(exp = e1, sz = NONE())  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.SIZE(e1_1, NONE())
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.SIZE(exp = e1, sz = SOME(e2))  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
          inExp
        else
          DAE.SIZE(e1_1, SOME(e2_1))
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.REDUCTION(reductionInfo = reductionInfo, expr = e1, iterators = riters)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (riters_1, ext_arg) = traverseReductionIterators(riters, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(riters, riters_1)
          inExp
        else
          DAE.REDUCTION(reductionInfo, e1_1, riters_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.CONS(e1, e2)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        (e2_1, ext_arg) = traverseExpBottomUp(e2, inFunc, ext_arg)
        e = if referenceEq(e1, e1_1) && referenceEq(e2, e2_1)
          inExp
        else
          DAE.CONS(e1_1, e2_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.LIST(expl)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.LIST(expl_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.META_TUPLE(expl)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.META_TUPLE(expl_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.META_OPTION(NONE())  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.META_OPTION(SOME(e1))  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.META_OPTION(SOME(e1_1))
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.BOX(e1)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.BOX(e1_1)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.UNBOX(e1, tp)  => begin
        (e1_1, ext_arg) = traverseExpBottomUp(e1, inFunc, inExtArg)
        e = if referenceEq(e1, e1_1)
          inExp
        else
          DAE.UNBOX(e1_1, tp)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.METARECORDCALL(fn, expl, fieldNames, i, typeVars)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        e = if referenceEq(expl, expl_1)
          inExp
        else
          DAE.METARECORDCALL(fn, expl_1, fieldNames, i, typeVars)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.MATCHEXPRESSION(matchTy, expl, aliases, localDecls, cases, tp)  => begin
        (expl_1, ext_arg) = traverseExpList(expl, inFunc, inExtArg)
        (cases_1, ext_arg) = Patternm.traverseCases(cases, inFunc, ext_arg)
        e = if referenceEq(expl, expl_1) && referenceEq(cases, cases_1)
          inExp
        else
          DAE.MATCHEXPRESSION(matchTy, expl_1, aliases, localDecls, cases_1, tp)
        end
        (e, ext_arg) = inFunc(e, ext_arg)
        (e, ext_arg)
      end
      
      DAE.SHARED_LITERAL(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.PATTERN(__)  => begin
        (e, ext_arg) = inFunc(inExp, inExtArg)
        (e, ext_arg)
      end
      
      DAE.CODE(__)  => begin
        (inExp, inExtArg)
      end
      
      _  => begin
        str = string(inExp)
        str = "Expression.traverseExpBottomUp or one of the user-defined functions using it is not implemented correctly: " + str
        @error str
        fail()
      end
    end
  end
  (outExp, outExtArg)
end


function traverseExpSubs(inSubscript::List{DAE.Subscript}, rel::Function, iarg::Type_a) ::Tuple{List{DAE.Subscript}, Type_a} 
  local outArg::Type_a
  local outSubscript::List{DAE.Subscript}

  (outSubscript, outArg) = begin
    local sub_exp::DAE.Exp
    local sub_exp_1::DAE.Exp
    local rest::List{DAE.Subscript}
    local res::List{DAE.Subscript}
    local arg::Type_a
    @match (inSubscript, rel, iarg) begin
      ( nil(), _, arg)  => begin
        (inSubscript, arg)
      end
      
      (DAE.WHOLEDIM(__) <| rest, _, arg)  => begin
        (res, arg) = traverseExpSubs(rest, rel, arg)
        res = if referenceEq(rest, res)
          inSubscript
        else
          _cons(DAE.WHOLEDIM(), res)
        end
        (res, arg)
      end
      
      (DAE.SLICE(exp = sub_exp) <| rest, _, arg)  => begin
        (sub_exp_1, arg) = traverseExpBottomUp(sub_exp, rel, arg)
        (res, arg) = traverseExpSubs(rest, rel, arg)
        res = if referenceEq(sub_exp, sub_exp_1) && referenceEq(rest, res)
          inSubscript
        else
          _cons(DAE.SLICE(sub_exp_1), res)
        end
        (res, arg)
      end
      
      (DAE.INDEX(exp = sub_exp) <| rest, _, arg)  => begin
        (sub_exp_1, arg) = traverseExpBottomUp(sub_exp, rel, arg)
        (res, arg) = traverseExpSubs(rest, rel, arg)
        res = if referenceEq(sub_exp, sub_exp_1) && referenceEq(rest, res)
          inSubscript
        else
          _cons(DAE.INDEX(sub_exp_1), res)
        end
        (res, arg)
      end
      
      (DAE.WHOLE_NONEXP(exp = sub_exp) <| rest, _, arg)  => begin
        (sub_exp_1, arg) = traverseExpBottomUp(sub_exp, rel, arg)
        (res, arg) = traverseExpSubs(rest, rel, arg)
        res = if referenceEq(sub_exp, sub_exp_1) && referenceEq(rest, res)
          inSubscript
        else
          _cons(DAE.WHOLE_NONEXP(sub_exp_1), res)
        end
        (res, arg)
      end
    end
  end
  (outSubscript, outArg)
end


function traverseExpCref(inCref::DAE.ComponentRef, rel::Function, iarg::Type_a) ::Tuple{DAE.ComponentRef, Type_a} 
  local outArg::Type_a
  local outCref::DAE.ComponentRef
  (outCref, outArg) = begin
    local name::String
    local cr::DAE.ComponentRef
    local cr_1::DAE.ComponentRef
    local ty::Type
    local subs::List{DAE.Subscript}
    local subs_1::List{DAE.Subscript}
    local arg::Type_a
    local ix::Integer
    local instant::String
    @match (inCref, rel, iarg) begin
      (DAE.CREF_QUAL(ident = name, identType = ty, subscriptLst = subs, componentRef = cr), _, arg)  => begin
        (subs_1, arg) = traverseExpSubs(subs, rel, arg)
        (cr_1, arg) = traverseExpCref(cr, rel, arg)
        cr = if referenceEq(cr, cr_1) && referenceEq(subs, subs_1)
          inCref
        else
          DAE.CREF_QUAL(name, ty, subs_1, cr_1)
        end
        (cr, arg)
      end
      
      (DAE.CREF_IDENT(ident = name, identType = ty, subscriptLst = subs), _, arg)  => begin
        (subs_1, arg) = traverseExpSubs(subs, rel, arg)
        cr = if referenceEq(subs, subs_1)
          inCref
        else
          DAE.CREF_IDENT(name, ty, subs_1)
        end
        (cr, arg)
      end
      
      (DAE.CREF_ITER(ident = name, index = ix, identType = ty, subscriptLst = subs), _, arg)  => begin
        (subs_1, arg) = traverseExpSubs(subs, rel, arg)
        cr = if referenceEq(subs, subs_1)
          inCref
        else
          DAE.CREF_ITER(name, ix, ty, subs_1)
        end
        (cr, arg)
      end
      
      (DAE.OPTIMICA_ATTR_INST_CREF(componentRef = cr, instant = instant), _, arg)  => begin
        (cr_1, arg) = traverseExpCref(cr, rel, arg)
        cr = if referenceEq(cr, cr_1)
          inCref
        else
          DAE.OPTIMICA_ATTR_INST_CREF(cr_1, instant)
        end
        (cr, arg)
      end
      
      (DAE.WILD(__), _, arg)  => begin
        (inCref, arg)
      end
      
      _  => begin
        @error "Expression.traverseExpCref: Unknown cref " * string(inCref)
        fail()
      end
    end
  end
  (outCref, outArg)
end


function evaluateCref(icr::DAE.ComponentRef, iels::List{<:DAE.Element})::Option{DAE.Exp} 
  local oexp::Option{DAE.Exp}
  local e::DAE.Exp
  local ee::DAE.Exp
  local crefs::List{DAE.ComponentRef}
  local oexps::List{Option{DAE.Exp}}
  local o::Option{DAE.Exp}
  oexp = getVarBinding(iels, icr)
  if isSome(oexp)
    @match SOME(e) = oexp
#    (e, _) = ExpressionSimplify.simplify(e)
    if isConst(e)
      oexp = SOME(e)
      return oexp
    end
    crefs = getAllCrefs(e)
    oexps = ListUtil.map1(crefs, evaluateCref, iels)
    for c in crefs
      @match _cons(SOME(ee), oexps) = oexps
      e = replaceCref(e, (c, ee))
      (e, _) = ExpressionSimplify.simplify(e)
    end
    oexp = SOME(e)
  end
  oexp
end

"""
  Replaces a component reference with an expression
"""
function replaceCref(inExp::DAE.Exp, inTpl::Tuple{<:DAE.ComponentRef, DAE.Exp})::Tuple{DAE.Exp, Tuple{DAE.ComponentRef, DAE.Exp}} 
  local otpl::Tuple{DAE.ComponentRef, DAE.Exp}
  local outExp::DAE.Exp
  (outExp, otpl) = begin
    local target::DAE.Exp
    local cr::DAE.ComponentRef
    local cr1::DAE.ComponentRef
    @match (inExp, inTpl) begin
      (DAE.CREF(componentRef = cr), (cr1, target)) where (string(cr) == string(cr1))  => begin
        (target, inTpl)
      end
      _  => begin
        (inExp, inTpl)
      end
    end
  end
  (outExp, otpl)
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


" author: lochel
  This function extracts all crefs from the input expression, except 'time'.
"
function getAllCrefs(inExp::DAE.Exp)::List{DAE.ComponentRef}
  local outCrefs::List{DAE.ComponentRef}
  (_, outCrefs) = traverseExpTopDown(inExp, getAllCrefs2, nil)
  outCrefs
end

function getAllCrefs2(inExp::DAE.Exp, inCrefList::List{<:DAE.ComponentRef})::Tuple{DAE.Exp, Bool, List{DAE.ComponentRef}}
  local outCrefList::List{DAE.ComponentRef} = inCrefList
  local outExp::DAE.Exp = inExp
  local cr::DAE.ComponentRef
  if isCref(inExp)
    @match DAE.CREF(componentRef = cr) = inExp
    if ! listMember(cr, inCrefList)
      outCrefList = _cons(cr, outCrefList)
    end
  end
  (outExp, true, outCrefList)
end

function isCref(inExp::DAE.Exp)
  res = @match inExp begin
    DAE.CREF(__) => true
    _ => false
  end
  return res
end


function isEvaluatedConst(inExp::DAE.Exp) ::Bool 
  local outBoolean::Bool
  outBoolean = begin
    @match inExp begin
      DAE.ICONST(__)  => begin
        true
      end
      DAE.RCONST(__)  => begin
        true
      end
      DAE.BCONST(__)  => begin
        true
      end
      DAE.SCONST(__)  => begin
        true
      end      
      DAE.ENUM_LITERAL(__)  => begin
        true
      end      
      _  => begin
        false
      end
    end
  end
  outBoolean
end

end #=End Util=#
