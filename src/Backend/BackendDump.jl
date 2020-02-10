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
    print("\nVars:\n")
    print(LINE + "\n")
    BackendDAEUtil.mapEqSystemVariablesNoUpdate(eq, printVarTraverse, 0)
    print("\nEqs:\n")
    print(LINE + "\n")
    BackendDAEUtil.mapEqSystemEquationsNoUpdate(eq, printEqTraverse, 0)
  end

  print("\n")
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

function printEqTraverse(eq::BackendDAE.Equation, extArg)
  _ = begin
    local lhs::DAE.Exp
    local rhs::DAE.Exp
    local cref::DAE.ComponentRef
    local whenEquation::BackendDAE.WhenEquation
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

      BackendDAE.WHEN_EQUATION(whenEquation = whenEquation) => begin
        printWhenEquation(whenEquation)
      end
    end
  end
end


function printWhenEquation(whenEq::BackendDAE.WhenEquation)
  local elseWhen::BackendDAE.WhenEquation
  print("when " + expStringify(whenEq.condition) + " then\n")

  for op in whenEq.whenStmtLst
    print("  " + whenOperatorStr(op) + "\n")
  end

  if isSome(whenEq.elsewhenPart)
    SOME(elseWhen) = whenEq.elseWhenPart
    print("else \n")
    printWhenEquation(elseWhen)
  end
  print("end;\n")
end

function whenOperatorStr(whenOp::BackendDAE.WhenOperator)::String
  str = begin
    local e1::DAE.Exp
    local e2::DAE.Exp
    local cref::DAE.ComponentRef
    @match whenOp begin
      BackendDAE.ASSIGN(left = e1, right = e2) => begin
        (expStringify(e1) + " := " + expStringify(e2))
      end
      BackendDAE.REINIT(stateVar = cref, value = e1) => begin
        ("reinit(" + crefStr(cref) + ", " + expStringify(e1) + ")")
      end
      BackendDAE.ASSERT(condition = e1, message = e2) => begin
        ("assert(" + expStringify(e1) + ", " + expStringify(e2) + ")")
      end
      BackendDAE.TERMINATE(message = e1) => begin
        ("[TERMINATE]" + expStringify(e1))
      end
      BackendDAE.NORETCALL(exp = e1) => begin
        ("[NORET]" + expStringify(e1))
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
    local lstexpl::List{List{DAE.Exp}}
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
        tmpStr = tmpStr + "[REC(" + expLstStringify(expl, ", ") + ")"
      end

      DAE.PARTEVALFUNCTION(path = Absyn.IDENT(tmpStr), expList = expl)  => begin
        tmpStr = tmpStr + "[PARTEVAL](" + expLstStringify(expl, ", ") + ")"
      end

      DAE.ARRAY(array = expl)  => begin
        "[ARR]" + expLstStringify(expl, ", ")
      end

      DAE.MATRIX(matrix = lstexpl)  => begin
        str = "[MAT]"
        for lst in lstexp
          str = str + "{" + expLstStringify(lst, ", ") + "}"
        end
        (str)
      end

      DAE.RANGE(start = e1, step = NONE(), stop = e2)  => begin
         expStringify(e1) + ":" + expStringify(e2)
      end

      DAE.RANGE(start = e1, step = SOME(e2), stop = e3)  => begin
         expStringify(e1) + ":" + expStringify(e2) + ":" + expStringify(e3)
      end

      DAE.TUPLE(PR = expl) => begin
         "[TPL](" + expLstStringify(expl, ", ") + ")"
      end

      DAE.CAST(exp = e1)  => begin
         "[CAST]" + expStringify(e1)
      end

      DAE.ASUB(exp = e1, sub = expl)  => begin
         "[ASUB]" + expStringify(e1) + "{" + expLstStringify(expl, ", ") + "}"
      end

      DAE.TSUB(exp = e1, ix = int) => begin
         "[TSUB]" + expStringify(e1) + "(" + string(int) + ")"
      end

      DAE.RSUB(exp = e1)  => begin
        "[RSUB]" + expStringify(e1)
      end

      DAE.SIZE(exp = e1, sz = NONE())  => begin
        "[SIZE]" + expStringify(e1)
      end

      DAE.SIZE(exp = e1, sz = SOME(e2))  => begin
         "[SIZE]" + expStringify(e1) + "(" + expStringify(e2) + ")"
      end

     DAE.CODE(__) => begin
       "[CODE]"
     end

     DAE.REDUCTION(expr = e1) => begin
       "[REDUCTION]" + expStringify(e1)
     end

     DAE.EMPTY(__)  => begin
       "[EMPTY]"
     end

     DAE.CONS(e1, e2)  => begin
       "[CONS]" + "{" + expStringify(e1) + ", " + expStringify(e2) + "}"
     end

     DAE.LIST(expl)  => begin
       "[LST]" + "{" + expLstStringify(expl, ", ") + " }"
     end

    _ => begin
      str = ""
    end

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

function dictPrettyPrint(d::Dict, pre=1)
    for (k,v) in d
        if typeof(v) <: Dict
            s = "$(repr(k)) => "
            println(join(fill(" ", pre)) * s)
            pretty_print(v, pre+1+length(s))
        else
            println(join(fill(" ", pre)) * "$(repr(k)) => $(repr(v))")
        end
    end
    nothing
end

@exportAll()
end
