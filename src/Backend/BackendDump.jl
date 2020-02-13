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
import SimulationCode

const HEAD_LINE = "#############################################"::String
const DOUBLE_LINE = "============================================"::String
const LINE = "---------------------------------------------"::String


function stringHeading1(in::Any, heading::String)::String
  str = heading1(heading) + "\n" + string(in)
end

function stringHeading2(in::Any, heading::String)::String
  str = heading2(heading) + "\n" + string(in)
end

function stringHeading3(in::Any, heading::String)::String
  str = heading3(heading) + "\n" + string(in)
end


function heading1(heading::String)::String
  str = HEAD_LINE + "\n" + heading + "\n" + HEAD_LINE + "\n\n"
end

function heading2(heading::String)::String
  str = DOUBLE_LINE + "\n" + heading + "\n" + DOUBLE_LINE + "\n\n"
end

function heading3(heading::String)::String
  str = heading + ":\n" + LINE + "\n"
end

function string(dae::BackendDAE.BackendDAEStructure)::String
  str::String = ""
  for i in 1:arrayLength(dae.eqs)
    str = str + stringHeading2(dae.eqs[i], "EqSystem " + Base.string(i)) + "\n"
  end
  return str
end

function string(eq::BackendDAE.EqSystem)::String
  str::String = ""
  str = str + heading3("Variables") + BackendDAEUtil.mapEqSystemVariablesNoUpdate(eq, stringTraverse, "") + "\n"
  str = str + heading3("Equations") + BackendDAEUtil.mapEqSystemEquationsNoUpdate(eq, stringTraverse, "") + "\n"
end

function stringTraverse(in, str)::String
  str = str + string(in)
end

function string(var::BackendDAE.Var)::String
  str = var.varName.ident + " | " + string(var.varKind) + "\n"
end

function string(varKind::BackendDAE.VarKind)::String
  str = begin
    @match varKind begin
      BackendDAE.VARIABLE() => begin "VARIABLE" end

      BackendDAE.STATE() => begin "STATE" end

      BackendDAE.STATE_DER() => begin "STATE_DER" end

      BackendDAE.DUMMY_DER() => begin "DUMMY_DER" end

      BackendDAE.DUMMY_STATE() => begin "DUMMY_STATE" end

      BackendDAE.CLOCKED_STATE() => begin "CLOCKED_STATE" end

      BackendDAE.DISCRETE() => begin "DISCRETE" end

      BackendDAE.PARAM() => begin "PARAM" end

      BackendDAE.CONST() => begin "CONST" end

      BackendDAE.EXTOBJ() => begin "EXTOBJ" end

      BackendDAE.JAC_VAR() => begin "JAC_VAR" end

      BackendDAE.JAC_DIFF_VAR() => begin "JAC_DIFF_VAR" end

      BackendDAE.SEED_VAR() => begin "SEED_VAR" end

      BackendDAE.OPT_CONSTR() => begin "OPT_CONSTR" end

      BackendDAE.OPT_FCONSTR() => begin "OPT_FCONSTR" end

      BackendDAE.OPT_INPUT_WITH_DER() => begin "OPT_INPUT_WITH_DER" end

      BackendDAE.OPT_INPUT_DER() => begin "OPT_INPUT_DER" end

      BackendDAE.OPT_TGRID() => begin "OPT_TGRID" end

      BackendDAE.OPT_LOOP_INPUT() => begin "OPT_LOOP_INPUT" end

      BackendDAE.ALG_STATE() => begin "ALG_STATE" end

      BackendDAE.ALG_STATE_OLD() => begin "ALG_STATE_OLD" end

      BackendDAE.DAE_RESIDUAL_VAR() => begin "DAE_RESIDUAL_VAR" end

      BackendDAE.DAE_AUX_VAR() => begin "DAE_AUX_VAR" end

      BackendDAE.LOOP_ITERATION() => begin "LOOP_ITERATION" end

      BackendDAE.LOOP_SOLVED() => begin "LOOP_SOLVED" end

    end
  end
end

function string(eq::BackendDAE.Equation)::String
  str = begin
    local lhs::DAE.Exp
    local rhs::DAE.Exp
    local cref::DAE.ComponentRef
    local whenEquation::BackendDAE.WhenEquation
    @match eq begin
      BackendDAE.EQUATION(lhs = lhs, rhs = rhs) => begin
        (string(lhs) + " = " + string(rhs))
      end

      BackendDAE.SOLVED_EQUATION(componentRef = cref, exp = rhs) => begin
        (string(cref) + " = " + string(rhs))
      end

      BackendDAE.RESIDUAL_EQUATION(exp = rhs) => begin
        ("0 = " + string(rhs))
      end

      BackendDAE.WHEN_EQUATION(whenEquation = whenEquation) => begin
        string(whenEquation)
      end

      BackendDAE.IF_EQUATION() => begin
        local strTmp::String
        local conditions::List{DAE.Exp}
        local condition::DAE.Exp
        local trueEquations::List{List{BackendDAE.Equation}}
        local trueEquation::List{BackendDAE.Equation}

        condition <| conditions = ifEq.conditions
        trueEquation <| trueEquations = ifEq.eqnstrue

        strTmp = "if " + string(condition) + " then\n"
        for eq in trueEquation
          strTmp = strTmp + "  " + string(eq) + "\n"
        end

        while listLength(conditions) != 0
          condition <| conditions = conditions
          trueEquation <| trueEquations = trueEquations

          strTmp = strTmp + "else if " + string(condition) + " then\n"
          for eq in trueEquation
            strTmp = strTmp + "  " + string(eq) + "\n"
          end
        end

        strTmp = strTmp + "else\n"
        for eq in eq.eqnsfalse
          strTmp = strTmp + "  " + string(eq) + "\n"
        end
        strTmp = strTmp + "end"
      end
    end
  end
  return str + "\n"
end

function string(whenEq::BackendDAE.WhenEquation)::String
  local elseWhen::BackendDAE.WhenEquation
  str = "when " + string(whenEq.condition) + " then\n"

  for op in whenEq.whenStmtLst
    str = str + "  " + string(op) + "\n"
  end

  if isSome(whenEq.elsewhenPart)
    SOME(elseWhen) = whenEq.elseWhenPart
    str = str + "else \n" + string(elseWhen)
  end
  return str + "end;\n"
end

function string(whenOp::BackendDAE.WhenOperator)::String
  str = begin
    local e1::DAE.Exp
    local e2::DAE.Exp
    local cref::DAE.ComponentRef
    @match whenOp begin
      BackendDAE.ASSIGN(left = e1, right = e2) => begin
        (string(e1) + " := " + string(e2))
      end
      BackendDAE.REINIT(stateVar = cref, value = e1) => begin
        ("reinit(" + string(cref) + ", " + string(e1) + ")")
      end
      BackendDAE.ASSERT(condition = e1, message = e2) => begin
        ("assert(" + string(e1) + ", " + string(e2) + ")")
      end
      BackendDAE.TERMINATE(message = e1) => begin
        ("[TERMINATE]" + string(e1))
      end
      BackendDAE.NORETCALL(exp = e1) => begin
        ("[NORET]" + string(e1))
      end
    end
  end
end

"need to add subscripts and other cases!"
function string(cr::DAE.ComponentRef)::String
  str = begin
    local ident::String
    local cref::DAE.ComponentRef
    @match cr begin
      DAE.CREF_QUAL(ident = ident, componentRef = cref) => begin
        (ident + "." + string(cref))
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

function string(op::DAE.Operator)::String
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

function string(exp::DAE.Exp)::String
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
        Base.string(int)
      end

      DAE.RCONST(real)  => begin
        Base.string(real)
      end

      DAE.SCONST(tmpStr)  => begin
        (tmpStr)
      end

      DAE.BCONST(bool)  => begin
        Base.string(bool)
      end

      DAE.ENUM_LITERAL((Absyn.IDENT(str), int))  => begin
        (str + "()" + string(int) + ")")
      end

      DAE.CREF(cr, _)  => begin
        string(cr)
      end

      DAE.UNARY(operator = op, exp = e1) => begin
        ("(" + string(op) + " " + string(e1) + ")")
      end

      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (string(e1) + " " + string(op) + " " + string(e2))
      end

      DAE.LUNARY(operator = op, exp = e1)  => begin
        ("(" + string(op) + " " + string(e1) + ")")
      end

      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (string(e1) + " " + string(op) + " " + string(e2))
      end

      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        (string(e1) + " " + string(op) + " " + string(e2))
      end

      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        (string(e1) + " " + string(e2) + " " + string(e3))
      end

      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = expl)  => begin
        tmpStr = tmpStr + "(" + lstString(expl, ", ") + ")"
      end

      DAE.RECORD(path = Absyn.IDENT(tmpStr), exps = expl)  => begin
        tmpStr = tmpStr + "[REC(" + lstString(expl, ", ") + ")"
      end

      DAE.PARTEVALFUNCTION(path = Absyn.IDENT(tmpStr), expList = expl)  => begin
        tmpStr = tmpStr + "[PARTEVAL](" + lstString(expl, ", ") + ")"
      end

      DAE.ARRAY(array = expl)  => begin
        "[ARR]" + lstString(expl, ", ")
      end

      DAE.MATRIX(matrix = lstexpl)  => begin
        str = "[MAT]"
        for lst in lstexp
          str = str + "{" + lstString(lst, ", ") + "}"
        end
        (str)
      end

      DAE.RANGE(start = e1, step = NONE(), stop = e2)  => begin
         string(e1) + ":" + v(e2)
      end

      DAE.RANGE(start = e1, step = SOME(e2), stop = e3)  => begin
         string(e1) + ":" + string(e2) + ":" + string(e3)
      end

      DAE.TUPLE(PR = expl) => begin
         "[TPL](" + lstString(expl, ", ") + ")"
      end

      DAE.CAST(exp = e1)  => begin
         "[CAST]" + string(e1)
      end

      DAE.ASUB(exp = e1, sub = expl)  => begin
         "[ASUB]" + string(e1) + "{" + lstString(expl, ", ") + "}"
      end

      DAE.TSUB(exp = e1, ix = int) => begin
         "[TSUB]" + string(e1) + "(" + string(int) + ")"
      end

      DAE.RSUB(exp = e1)  => begin
        "[RSUB]" + string(e1)
      end

      DAE.SIZE(exp = e1, sz = NONE())  => begin
        "[SIZE]" + string(e1)
      end

      DAE.SIZE(exp = e1, sz = SOME(e2))  => begin
         "[SIZE]" + string(e1) + "(" + string(e2) + ")"
      end

     DAE.CODE(__) => begin
       "[CODE]"
     end

     DAE.REDUCTION(expr = e1) => begin
       "[REDUCTION]" + string(e1)
     end

     DAE.EMPTY(__)  => begin
       "[EMPTY]"
     end

     DAE.CONS(e1, e2)  => begin
       "[CONS]" + "{" + string(e1) + ", " + string(e2) + "}"
     end

     DAE.LIST(expl)  => begin
       "[LST]" + "{" + lstString(expl, ", ") + " }"
     end

    _ => begin
      str = ""
    end

    end
  end
end

function lstString(expLst::List{T}, seperator::String)::String where{T}
  str = begin
    local e::T
    local rest::List{T}
    @match expLst begin
      (e <| rest) => begin
        str = string(e)
        for r in rest
          str = str + seperator + string(r)
        end
        (str)
      end
      _ => begin
        str = ""
      end
    end
  end
end

function string(simCode::SimulationCode.SIM_CODE)::String
  str = stringHeading3(simCode.crefToSimVarHT, "SimCodeVars")
  str = str + heading3("SimCodeEquations")
  for eq in simCode.equations
    str = str + string(eq)
  end
  return str + "\n"
end

function string(d::Dict)::String
  str = ""
  for (k,v) in d
    str = str + "  $(repr(k)) => $(repr(v))\n"
  end
  return str + "\n"
end

@exportAll()
end
