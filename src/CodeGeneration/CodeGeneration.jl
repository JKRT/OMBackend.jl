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
  local startEquations = ""
  local crefToSimVarHT = simCode.crefToSimVarHT
  local modelName::String = simCode.name
  local exp::DAE.Exp
  local index::Integer
  #= An array of 0:s=#
  local residuals::Array = [0 for i in 1:length(simCode.equations)]
  for varName in keys(crefToSimVarHT)
    simVar = crefToSimVarHT[varName]
    @match simVar.varKind  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algVariables, varName)
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end

  # Generate state variable marking and start equations for each state
  for var in stateVariables
    push!(stateMarkings, true)
    simVar = crefToSimVarHT[var]
    (startExp, index) = @match simVar begin
      SimulationCode.SIMVAR(attributes=SOME(DAE.VAR_ATTR_REAL(start=SOME(exp))), index=SOME(index)) => (exp, index)
      #= if there is no start expression take default 0 =#
      SimulationCode.SIMVAR(index=SOME(index)) => (DAE.RCONST(0), index)
      _ => begin
        ErrorException("SimVar has no index: $(var)")
        (DAE.RCONST(0), -2)
      end
    end
    startEquations *= "  x0[$(index)] #= $var =# = $(expStringify(startExp, simCode))\n"
  end

  for var in algVariables
    push!(stateMarkings, false)
  end
  local differentialVarsFunction ="
function $(modelName)DifferentialVars()
  return $stateMarkings
end
"

  # Generate start values
  #for var in hcat(stateVariables, algVariables, stateDerivatives)
    # TODO match VariableAttributes to get start
    #@match var.attributes begin
      #println(var.attributes)
    #end
    #startEquations *= "$var\n"
  #end
  local startCondtions ="
function $(modelName)StartConditions(p, t0)
  x0 = Array{Float64}(undef, $(length(stateVariables)+length(algVariables)))
  dx0 = Array{Float64}(undef, $(length(stateDerivatives)))

$startEquations  return x0, dx0
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
    simVar = crefToSimVarHT[param]
    _ = @match simVar begin
      SimulationCode.SIMVAR(varKind=SimulationCode.PARAMETER(bindExp = SOME(exp)), index = SOME(index)) => begin
        parameterEquations *= "  p[$(index)] #= $param =# = $(expStringify(exp, simCode))\n"
      end
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
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
    local index::Integer
    @match exp begin
      DAE.ICONST(int) => string(int)

      DAE.RCONST(real) => string(real)

      DAE.SCONST(tmpStr) => tmpStr

      DAE.BCONST(bool)  => string(bool)

      DAE.ENUM_LITERAL((Absyn.IDENT(str), int)) => (str + "(" + string(int) + ")")

      DAE.CREF(cr, _)  => begin
        varName = BackendDump.string(cr)
        simVar = hashTable[varName]
        @match simVar begin
          SimulationCode.SIMVAR(varKind=SimulationCode.INPUT(__), index = SOME(index)) => @error "INPUT not supported in CodeGen"
          SimulationCode.SIMVAR(varKind=SimulationCode.STATE(__), index = SOME(index)) => "x[$(index)] #= $varName =#"
          SimulationCode.SIMVAR(varKind=SimulationCode.PARAMETER(__), index = SOME(index)) => "p[$(index)] #= $varName =#"
          SimulationCode.SIMVAR(varKind=SimulationCode.ALG_VARIABLE(__), index = SOME(index)) => "x[$(index)] #= $varName =#"
          SimulationCode.SIMVAR(varKind=SimulationCode.STATE_DERIVATIVE(__), index = SOME(index)) => "dx[$(index)] #= der($varName) =#"
        end
      end

      DAE.UNARY(operator = op, exp = e1) => ("(" + BackendDump.string(op) + " " + expStringify(e1, simCode) + ")")

      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => (expStringify(e1, simCode) + " " + BackendDump.string(op) + " " + expStringify(e2, simCode))

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
        simVar = hashTable[varName]
        @match (tmpStr, simVar.index) begin
          ("der", SOME(index)) => "dx[$(index)]  #= der($varName) =#"
          _  =>  begin
            tmpStr = tmpStr + "(" + BackendDump.lstStr(expl, ", ") + ")"
          end
        end
      end

      DAE.RECORD(path = Absyn.IDENT(tmpStr), exps = expl)  => begin
        tmpStr = tmpStr + "[REC](" + BackendDump.lstStr(expl, ", ") + ")"
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
