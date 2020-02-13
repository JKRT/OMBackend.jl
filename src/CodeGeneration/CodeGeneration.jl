#= /*
* This file is part of OpenModelica.
*
* Copyright (c) 1998-2020, Open Source Modelica Consortium (OSMC),
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

module CodeGeneration

import Absyn
import BackendDump

using BackendDAE
using SimulationCode

include("CodeGenerationUtil.jl")
include("simulationCodeTransformation.jl")

"""
The header string with the necessary imports
"""
const HEADER_STRING ="
$(copyRightString())

using DifferentialEquations
using Sundials
using Plots
"

function generate_state_boolean_vector()
end


"""
Write  modelName_DAE_equations() to file.
"""
function writeDAE_equationsToFile(fileName::String, contents::String)
  local fdesc = open(fileName, "w")
  write(fdesc, contents)
  close(fdesc)
end


"""
Generate a julia file containing functions to simulate the DAE
"""
function generateCode(simCode::SimulationCode.SIM_CODE)
  local stateVariables::Array = []
  local algVariables::Array = []
  local stateDerivatives::Array = []
  local parameters::Array = []
  local stateMarkings::Array = []
  local parameterEquations = ""
  local crefToSimVarHT = simCode.crefToSimVarHT
  local modelName::String = simCode.name
  local exp::DAE.Exp
  #= An array of 0:s=#
  local residuals::Array = [0 for i in 1:length(simCode.equations)]
  for varName in keys(crefToSimVarHT)
    ixAndTy = crefToSimVarHT[varName]
    local varType = ixAndTy[2]
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algVariables, varName)
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end

  for i in stateVariables
    push!(stateMarkings, true)
  end
  for i in algVariables
    push!(stateMarkings, false)
  end
  local differentialVarsFunction ="
function $(modelName)DifferentialVars()
  return $stateMarkings
end
"
  local startCondtions ="
function $(modelName)StartConditions(p, t0)
  x0 = [1.0]
  dx0 = [p[1]*x0[1]]
  return x0, dx0
end
"
  #= Generate $modelName_DAE_equations=#
  local DAE_EQUATIONS = let
    local eqStr = ""
    for (equationCounter,eq) in enumerate(simCode.equations)
      eqStr *= eqTraverseAppendToString(eq, simCode, equationCounter)
    end
    eqStr[1:end-1]
  end

  local dae_equation_function ="
function $(modelName)DAE_equations(res, dx #=The state derivatives =#, x #= State & alg variables =#, p, t #=time=#)
$DAE_EQUATIONS
end
"
  for param in parameters
    (index, simVarType) = crefToSimVarHT[param]
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => begin exp
    end
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
    parameterEquations *= "  p[$index] #= $param =# = $(expStringify(bindExp, simCode))\n"
  end
  local parameterVars ="
function $(modelName)ParameterVars()
  p = Array{Float64}(undef, $(arrayLength(parameters)))
$(parameterEquations)  return p
end
"

local runnable ="
function $(modelName)Simulate(tspan = (0.0, 1.0))
  # Define problem
  p_is = $(modelName)ParameterVars()
  (x0, dx0) =$(modelName)StartConditions(p_is, tspan[1])
  differential_vars = $(modelName)DifferentialVars()
  #= Pass the residual equations =#
  problem = DAEProblem($(modelName)DAE_equations, dx0, x0, tspan, p_is, differential_vars=differential_vars)
  # Solve with IDA:)
  solution = solve(problem, IDA())
  return solution
end
"
  # Return file content
  return ("$(modelName)",
   HEADER_STRING * startCondtions * differentialVarsFunction
   * dae_equation_function * parameterVars * runnable)
end

"""
TODO: Make less messy
"""
function eqTraverseAppendToString(eq::BackendDAE.Equation, simCode::SimulationCode.SIM_CODE, resNumber)
  _ = begin
    local lhs::DAE.Exp
    local rhs::DAE.Exp
    local cref::DAE.ComponentRef
    local whenEquation::BackendDAE.WhenEquation
    local result::String = ""
    @match eq begin
      BackendDAE.RESIDUAL_EQUATION(exp = rhs) => begin
        result = result * ("  res[$resNumber] = " + expStringify(rhs, simCode) + "\n")
      end
      BackendDAE.WHEN_EQUATION(whenEquation = whenEquation) => begin
        ErrorException("When equations not yet supported")

      end
      _ =>
        ErrorException("traversalError for $eq")
    end
  end
  return result
end

function expStringify(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE)::String
  hashTable = simCode.crefToSimVarHT
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
        varName = BackendDump.string(cr)
        indexAndType = hashTable[varName]
        @match indexAndType[2] begin
          SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
          SimulationCode.STATE(__) => "x[$(indexAndType[1])] #= $varName =#"
          SimulationCode.PARAMETER(__) => "p[$(indexAndType[1])] #= $varName =#"
          SimulationCode.ALG_VARIABLE(__) => "x[$(indexAndType[1])] #= $varName =#"
          SimulationCode.STATE_DERIVATIVE(__) => "dx[$(indexAndType[1])] #= der($varName) =#"
        end
      end

      DAE.UNARY(operator = op, exp = e1) => begin
        ("(" + BackendDump.string(op) + " " + expStringify(e1, simCode) + ")")
      end

      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expStringify(e1, simCode) + " " + BackendDump.string(op) + " " + expStringify(e2, simCode))
      end

      DAE.LUNARY(operator = op, exp = e1)  => begin
        ("(" + BackendDump.string(op) + " " + expStringify(e1, simCode) + ")")
      end

      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expStringify(e1, simCode) + " " + BackendDump.string(op) + " " + expStringify(e2, simCode))
      end

      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        (expStringify(e1, simCode) + " " + BackendDump.string(op) + " " + expStringify(e2, simCode))
      end

      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        (expStringify(e1, simCode) + " " + expStringify(e2, simCode) + " " + expStringify(e3, simCode))
      end

      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = expl)  => begin
        #=
          TODO: Keeping it simple for now=, we assume we only have one argument in the call
          We handle derivitives seperatly
        =#
        varName = BackendDump.string(listHead(expl))
        (index, type) = hashTable[varName]
        @match tmpStr begin
          "der" => "dx[$index]  #= der($varName) =#"
          _  =>  begin
            tmpStr = tmpStr + "(" + BackendDump.lstStr(expl, ", ") + ")"
          end
        end
      end

      DAE.RECORD(path = Absyn.IDENT(tmpStr), exps = expl)  => begin
        tmpStr = tmpStr + "[REC(" + BackendDump.lstStr(expl, ", ") + ")"
      end

      DAE.PARTEVALFUNCTION(path = Absyn.IDENT(tmpStr), expList = expl)  => begin
        tmpStr = tmpStr + "[PARTEVAL](" + BackendDump.lstStr(expl, ", ") + ")"
      end

      DAE.ARRAY(array = expl)  => begin
        "[ARR]" + BackendDump.lstStr(expl, ", ")
      end

      DAE.MATRIX(matrix = lstexpl)  => begin
        str = "[MAT]"
        for lst in lstexp
          str = str + "{" + BackendDump.lstStr(lst, ", ") + "}"
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
         "[TPL](" + BackendDump.lstStr(expl, ", ") + ")"
      end

      DAE.CAST(exp = e1)  => begin
         expStringify(e1, simCode)
      end

      DAE.ASUB(exp = e1, sub = expl)  => begin
         "[ASUB]" + expStringify(e1, simCode) + "{" + BackendDump.lstStr(expl, ", ") + "}"
      end

      DAE.TSUB(exp = e1, ix = int) => begin
         "[TSUB]" + expStringify(e1, simCode) + "(" + string(int) + ")"
      end

      DAE.RSUB(exp = e1)  => begin
        "[RSUB]" + expStringify(e1, simCode)
      end

      DAE.SIZE(exp = e1, sz = NONE())  => begin
        "[SIZE]" + expStringify(e1, simCode)
      end

      DAE.SIZE(exp = e1, sz = SOME(e2))  => begin
         "[SIZE]" + expStringify(e1, simCode) + "(" + expStringify(e2, simCode) + ")"
      end

     DAE.CODE(__) => begin
       "[CODE]"
     end

     DAE.REDUCTION(expr = e1) => begin
       "[REDUCTION]" + expStringify(e1, simCode)
     end

     DAE.EMPTY(__)  => begin
       "[EMPTY]"
     end

     DAE.CONS(e1, e2)  => begin
       "[CONS]" + "{" + expStringify(e1, simCode) + ", " + expStringify(e2, simCode) + "}"
     end

     DAE.LIST(expl)  => begin
       "[LST]" + "{" + BackendDump.lstStr(expl, ", ") + " }"
     end

    _ => begin
      str = ""
    end
    end
  end
str = "(" * str * ")"
end

end #= End CodeGeneration=#
