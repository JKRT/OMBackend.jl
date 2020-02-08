#= /*
* This file is part of OpenModelica.
*
* Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
* c/o Linköpings universitet, Department of Computer and Information Science,
* SE-58183 Linköping, Sweden.
*
* All rights reserved.
*
* THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
* THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
* ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
* RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
* ACCORDING TO RECIPIENTS CHOICE.
*
* The OpenModelica software and the Open Source Modelica
* Consortium (OSMC) Public License (OSMC-PL) are obtained
* from OSMC, either from the above address,
* from the URLs: http:www.ida.liu.se/projects/OpenModelica or
* http:www.openmodelica.org, and in the OpenModelica distribution.
* GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
*
* This program is distributed WITHOUT ANY WARRANTY; without
* even the implied warranty of  MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
* IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
*
* See the full OSMC Public License conditions for more details.
*
=#


module BackendDump

using MetaModelica

#= ExportAll is not good practice but it makes it so that we do not have to write export after each function :( =#
using ExportAll

import Absyn
import DAE

import BackendDAE
import BackendDAEUtil

const DOUBLE_LINE = "============================================"::String
const LINE = "---------------------------------------------"::String


function dumpBackendDAEStructure(dae::BackendDAE.BackendDAEStructure, heading::String)
  print(DOUBLE_LINE + "\n")
  print("BackendDAE: " + heading + "\n")
  print(DOUBLE_LINE + "\n")

  for eq in dae.eqs
    print("\nEqs:\n")
    print(LINE + "\n")
    BackendDAEUtil.mapEqSystemEquationsNoUpdate(eq, printEqTraverse, 0)
    print("\nVars:\n")
    print(LINE + "\n")
    BackendDAEUtil.mapEqSystemVariablesNoUpdate(eq, printVarTraverse, 0)
  end
end

function printAnyTraverse(any, extArg)
  print(any)
  print("\n")
  (extArg)
end

function printVarTraverse(var::BackendDAE.Var, extArg)
  print(var.varName.ident)
  print(" | ")
  print(var.varKind)
  print("\n")
  (extArg)
end

# kabdelhak: Very ugly, this needs improvement!
# We need a reasonable printExp() function
function printEqTraverse(eq::BackendDAE.Equation, extArg)
  _ = begin
    local lhs::DAE.Exp
    local rhs::DAE.Exp
    local cref::DAE.ComponentRef
    @match eq begin
      BackendDAE.EQUATION(lhs = lhs, rhs = rhs) => begin
        print(expStringify(lhs) + " = " + expStringify(rhs) + "\n")
      end

      BackendDAE.SOLVED_EQUATION(componentRef = cref, exp = rhs) => begin
        print(crefStr(cref) + " = " + expStringify(rhs) + "\n")
      end

      BackendDAE.RESIDUAL_EQUATION(exp = rhs) => begin
        print("0 = " + expStringify(rhs) + "\n")
      end
    end
  end
end

"need to add subscripts and other cases!"
function crefStr(cr::DAE.ComponentRef)::String
  str = begin
    local ident::String
    local cref::DAE.ComponentRef
    @match cr begin
      DAE.CREF_QUAL(ident = ident, componentRef = cref) => begin
        (ident + "." + crefStr(cref))
      end
      DAE.CREF_IDENT(ident = ident) => begin
        (ident)
      end
      DAE.CREF_ITER(ident = ident) => begin
        (ident)
      end
    end
  end
end

function opStr(op::DAE.Operator)::String
  str = begin
    @match op begin
      DAE.ADD() => begin
        ("+")
      end

      DAE.SUB() => begin
        ("-")
      end

      DAE.MUL() => begin
        ("*")
      end

      DAE.DIV() => begin
        ("/")
      end

      DAE.POW() => begin
        ("^")
      end

      DAE.UMINUS() => begin
        ("-")
      end

      DAE.UMINUS_ARR() => begin
        ("-")
      end

      DAE.ADD_ARR() => begin
        ("+")
      end

      DAE.SUB_ARR() => begin
        ("-")
      end

      DAE.MUL_ARR() => begin
        ("*")
      end

      DAE.DIV_ARR() => begin
        ("/")
      end

      DAE.MUL_ARRAY_SCALAR() => begin
        ("*")
      end

      DAE.ADD_ARRAY_SCALAR() => begin
        ("+")
      end

      DAE.SUB_SCALAR_ARRAY() => begin
        ("-")
      end

      DAE.MUL_SCALAR_PRODUCT() => begin
        ("*")
      end

      DAE.MUL_MATRIX_PRODUCT() => begin
        ("*")
      end

      DAE.DIV_ARRAY_SCALAR() => begin
        ("/")
      end

      DAE.DIV_SCALAR_ARRAY() => begin
        ("/")
      end

      DAE.POW_ARRAY_SCALAR() => begin
        ("^")
      end

      DAE.POW_SCALAR_ARRAY() => begin
        ("^")
      end

      DAE.POW_ARR() => begin
        ("^")
      end

      DAE.POW_ARR2() => begin
        ("^")
      end

      DAE.AND() => begin
        ("and")
      end

      DAE.OR() => begin
        ("or")
      end

      DAE.NOT() => begin
        ("not")
      end

      DAE.LESS() => begin
        ("<")
      end

      DAE.LESSEQ() => begin
        ("<=")
      end

      DAE.GREATER() => begin
        (">")
      end

      DAE.GREATEREQ() => begin
        (">=")
      end

      DAE.EQUAL() => begin
        ("=")
      end

      DAE.NEQUAL() => begin
        ("<>")
      end

      DAE.USERDEFINED() => begin
        ("[UNDEF OP]")
      end

    end
  end
end

function expStringify(exp::DAE.Exp)::String
  str = begin
    local int::ModelicaInteger
    local real::ModelicaReal
    local bool::Bool
    local tmpStr::String
    local cr::DAE.ComponentRef
    local e1::DAE.Exp
    local e2::DAE.Exp
    local e3::DAE.Exp
    local expl::List{DAE.Exp}
    @match exp begin
      DAE.ICONST(int) => begin
        string(int)
      end

      DAE.RCONST(real)  => begin
        string(real)
      end

      DAE.SCONST(tmpStr)  => begin
        (tmpStr)
      end

      DAE.BCONST(bool)  => begin
        string(bool)
      end

      DAE.ENUM_LITERAL((Absyn.IDENT(str), int))  => begin
        (str + "()" + string(int) + ")")
      end

      DAE.CREF(cr, _)  => begin
        crefStr(cr)
      end

      DAE.UNARY(operator = op, exp = e1) => begin
        ("(" + opStr(op) + " " + expStringify(e1) + ")")
      end

      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expStringify(e1) + " " + opStr(op) + " " + expStringify(e2))
      end

      DAE.LUNARY(operator = op, exp = e1)  => begin
        ("(" + opStr(op) + " " + expStringify(e1) + ")")
      end

      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expStringify(e1) + " " + opStr(op) + " " + expStringify(e2))
      end

      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        (expStringify(e1) + " " + opStr(op) + " " + expStringify(e2))
      end

      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        (expStringify(e1) + " " + expStringify(e2) + " " + expStringify(e3))
      end

      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = expl)  => begin
        tmpStr = tmpStr + "(" + expLstStringify(expl, ", ") + ")"
      end

      DAE.RECORD(path = Absyn.IDENT(tmpStr), exps = expl)  => begin
        tmpStr = tmpStr + "(" + expLstStringify(expl, ", ") + ")"
      end

      DAE.PARTEVALFUNCTION(path = Absyn.IDENT(tmpStr), expList = expl)  => begin
        tmpStr = tmpStr + "(" + expLstStringify(expl, ", ") + ")"
      end

      DAE.ARRAY(scalar = e1)  => begin
        expStringify(eq)
      end

      _ => begin
        str = ""
      end

"""
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
     """
    end
  end
end

function expLstStringify(expLst::List{DAE.Exp}, seperator::String)::String
  str = begin
    local e::DAE.Exp
    local rest::List{DAE.Exp}
    @match expLst begin
      (e <| rest) => begin
        str = expStringify(e)
        for r in rest
          str = str + seperator + expStringify(r)
        end
        (str)
      end
      _ => begin
        str = ""
      end
    end
  end
end

@exportAll()
end
