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

using DataStructures
using MetaModelica
using Setfield

using ..FrontendUtil
using ..Backend
using ..SimulationCode

import ..Backend.BDAE
import DAE
import Absyn

include("codeGenerationUtil.jl")

"""
  The header string with the necessary imports
"""
const HEADER_STRING ="
$(copyRightString())

using DiffEqBase
using DifferentialEquations
using Plots
using Sundials

"

global callbackCounter = 0
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

function generateCode(simCode::SimulationCode.EXPLICIT_SIM_CODE)
  local stateVariables::Array = []
  local algVariables::Array = []
  local stateDerivatives::Array = []
  local parameters::Array = []
  local stateMarkings::Array = []
  local crefToSimVarHT = simCode.nameToVar
  local modelName::String = simCode.name
  local exp::DAE.Exp
  #= An array of 0:s=#
  local residuals::Array = [0 for _ in 1:length(simCode.residualEquations)]
  @debug "Residuals:" simCode.residualEquations
  @debug "Mapping:" simCode.eqVariableMapping
end


"""
  Generate Julia code from SimCode.
"""
function generateCode(simCode::SimulationCode.SIM_CODE)
  local stateVariables::Array = []
  local algVariables::Array = []
  local stateDerivatives::Array = []
  local parameters::Array = []
  local stateMarkings::Array = []
  local crefToSimVarHT = simCode.crefToSimVarHT
  local modelName::String = simCode.name
  local exp::DAE.Exp
  #= An array of 0:s=#
  local residuals::Array = [0 for _ in 1:length(simCode.residualEquations)]
  for varName in keys(crefToSimVarHT)
    ixAndVar = crefToSimVarHT[varName]
    local varType = ixAndVar[2].varKind
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algVariables, varName)
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end
  local stateMarkings = vcat([true for _ in stateVariables], [false for _ in algVariables])
  #=Start conditions of algebraic and state variables=#
  local DAE_EQUATIONS = createResidualEquations(simCode.residualEquations, simCode)
  #= Discrete events=#
  local WHEN_EQUATIONS = createEquations(simCode.whenEquations, simCode)
  #=Continuous events=#
  local IF_EQUATIONS = createEquations(simCode.ifEquations, simCode)
  local PARAMETER_EQUATIONS = createParameterEquations(parameters, simCode)
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquations(algVariables, stateVariables, simCode)
  @debug("Generating start conditions")
  local START_CONDTIONS ="
  function $(modelName)StartConditions(p, t0)
    local x0 = Array{Float64}(undef, $(length(stateVariables) + length(algVariables)))
    local dx0 = Array{Float64}(undef, $(length(stateVariables) + length(algVariables)))
    $START_CONDTIONS_EQUATIONS
    return x0, dx0
  end
  "
  local DAE_EQUATION_FUNCTION ="
  function $(modelName)DAE_equations(res, dx #=The state derivatives =#, x #= State & alg variables =#, p, t #=time=#)
  $DAE_EQUATIONS
  end
  "
  local DIFFERENTIAL_VARS_FUNCTION ="
  function $(modelName)DifferentialVars()
    return $stateMarkings
  end
  "
  local PARAMETER_VARS = "
  function $(modelName)ParameterVars()
    p = Array{Float64}(undef, $(arrayLength(parameters)))
  $(PARAMETER_EQUATIONS)  return p
  end
  "
  local WHEN_EQUATIONS_FN ="
  function $(modelName)CallbackSet(p#=parameters=#)
    #= WHEN EQUATIONS. Discrete callbacks =#
    $WHEN_EQUATIONS
    #= IF EQUATIONS =#
    $IF_EQUATIONS
    return CallbackSet($(evaluateCallBackset()))
  end
  "
  local RUNNABLE ="
  function $(modelName)Simulate(tspan = (0.0, 1.0))
    # Define problem
    p = $(modelName)ParameterVars()
    (x0, dx0) =$(modelName)StartConditions(p, tspan[1])
    differential_vars = $(modelName)DifferentialVars()
    #= Pass the residual equations =#
    problem = DAEProblem($(modelName)DAE_equations, dx0, x0, 
                         tspan, p, differential_vars=differential_vars, 
                         callback=$(modelName)CallbackSet(p))
    # Solve with IDA:)
    solution = solve(problem, IDA())
    return solution
  end
  "
  @debug("Code-generation done")
  global callbackCounter = 0
  # Return file content
  return ("$(modelName)",
          HEADER_STRING
          * START_CONDTIONS
          * DIFFERENTIAL_VARS_FUNCTION
          * DAE_EQUATION_FUNCTION
          * PARAMETER_VARS
          * WHEN_EQUATIONS_FN
          * RUNNABLE)
end


function evaluateCallBackset()::String
    local tmpStr::String = ""
    for t in 1:callbackCounter
      tmpStr += "cb$(t),"
    end
    return tmpStr
end

"""
    TODO: We use state for everything now. I am unsure how to differentiate between the two
"""
function createStartConditionsEquations(algVariables::Array,
                                        stateVariables::Array,
                                        simCode::SimulationCode.SIM_CODE)
  return getStartConditions(algVariables, "x0", simCode) *
         getStartConditions(stateVariables, "x0", simCode)
end

function createResidualEquations(equations::Array, simCode::SimulationCode.SIM_CODE)
  local eqStr = ""
  for (equationCounter, eq) in enumerate(equations)
    eqStr *= residualEqtoJulia(eq, simCode, equationCounter)
  end
  eqStr[1:end-1]
end

function createEquations(equations::Array, simCode::SimulationCode.SIM_CODE)
  local eqStr = ""
  for (equationCounter, eq) in enumerate(equations)
    eqStr *= eqToJulia(eq, simCode, equationCounter)
  end
  eqStr[1:end-1]
end

function createParameterEquations(parameters::Array, simCode::SimulationCode.SIM_CODE)
  local parameterEquations::String = ""
  local hT = simCode.crefToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    local simVarType::SimulationCode.SimVarType = simVar.varKind
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => exp
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
    parameterEquations *= "  p[$index] #= $param =# = $(expToJL(bindExp, simCode))\n"
  end
  return parameterEquations
end

"""
  Transforms a given equation into Julia code
"""
function residualEqtoJulia(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, resNumber)::String
  local rhs::DAE.Exp
  local result::String = ""
  result  *= @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
      "res[$resNumber] = " + expToJL(rhs, simCode) + "\n"
    end
    _ => begin
      "traversalError for $eq"
    end
  end
end

"""
  Transforms a given equation into Julia code
The optional parameter prefix specifices the prefix to which the result is written. 
"""
function eqToJulia(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, resNumber; prefix="res")::String
  local rhs::DAE.Exp
  local result::String = ""
  result  *= @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
      @info "Value of eqToJulia $prefix"
      "  $(prefix)[$resNumber] = " + expToJL(rhs, simCode; varPrefix = prefix) + "\n"
    end
    BDAE.WHEN_EQUATION(_, wEq, _) => begin
      global callbackCounter += 1
      local whenStmts = createWhenStatements(wEq.whenStmtLst, simCode, prefix="integerator.u")
      local cond = wEq.condition
      "condition(x,t,integrator) = bool(floor($(expToJL(cond, simCode))))
       affect!(integrator) = let
        $whenStmts
       end
       cb$(callbackCounter) = DiscreteCallback(condition, affect!)
      "
    end
    BDAE.IF_EQUATION(conds, trueEqs, falseEqs) => begin
      global callbackCounter += 1
      f(x) = eqToJulia(x, simCode, resNumber, prefix = "integerator.u")
      local condsJL = tuple(map((x) -> expToJL(x, simCode; varPrefix="x"), conds)...)
      local trueEqJL = map(f, trueEqs)
      local falseEqJL = map(f, falseEqs)
      truePart = "condition(x,t,integrator) = bool(floor($(condsJL...)))
       affect$(callbackCounter)!(integrator) = let
        $(trueEqJL...)
       end
       cb$(callbackCounter)$(callbackCounter) = ContiniuousCallback(condition, affect!)
      "
      global callbackCounter += 1
      falsePart =
       "condition$(callbackCounter)(x,t,integrator) = ! bool(floor($(condsJL...)))
       affect$(callbackCounter)!(integrator) = let
        $(falseEqJL...)
       end
       cb$(callbackCounter) = ContinuousCallback(condition, affect!)
      "
      truePart + falsePart
    end
    _ => begin
      "traversalError for $eq"
    end
  end
end


"""
  Creates Julia code for  when equation.
  E.g a set of Discrete callbacks
  TODO: Handle the other cases
"""
function createWhenStatements(whenStatements::List, simCode::SimulationCode.SIM_CODE; prefix="x")
  local res::String = ""
  @debug "Calling createWhenStatements with: $whenStatements"
  for wStmt in  whenStatements
    res *= @match wStmt begin
      BDAE.ASSIGN(__) => begin
        "$(expToJL(wStmt.left, simCode, varPrefix=prefix)) = $(expToJL(wStmt.right, simCode, varPrefix=prefix)) \n"
      end
      BDAE.REINIT(__) => begin
        "$(BDAE.string(wStmt.stateVar)) = $(expToJL(wStmt.value, simCode, varPrefix=prefix))\n"
      end
      _ => ErrorException("$whenStatements in @__FUNCTION__ not supported")
    end
  end
  return res
end

"
 Conver DAE.Exp into a Julia string. 
"
function expToJL(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE; varPrefix="x")::String
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
      DAE.ICONST(int) => string(int)
      DAE.RCONST(real)  => string(real)
      DAE.SCONST(tmpStr)  => (tmpStr)
      DAE.BCONST(bool)  => string(bool)
      DAE.ENUM_LITERAL((Absyn.IDENT(str), int))  => str + "()" + string(int) + ")"
      DAE.CREF(cr, _)  => begin
        varName = BDAE.string(cr)
        indexAndVar = hashTable[varName]
        varKind::SimulationCode.SimVarType = indexAndVar[2].varKind
        @match varKind begin
          SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
          SimulationCode.STATE(__) => "$varPrefix[$(indexAndVar[1])] #= $varName =#"
          SimulationCode.PARAMETER(__) => "p[$(indexAndVar[1])] #= $varName =#"
          SimulationCode.ALG_VARIABLE(__) => "$varPrefix[$(indexAndVar[1])] #= $varName =#"
          SimulationCode.STATE_DERIVATIVE(__) => "dx[$(indexAndVar[1])] #= der($varName) =#"
        end
      end
      DAE.UNARY(operator = op, exp = e1) => begin
        ("(" + BDAE.string(op) + " " + expToJL(e1, simCode) + ")")
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expToJL(e1, simCode, varPrefix=varPrefix) + " " + BDAE.string(op) + " " + expToJL(e2, simCode, varPrefix=varPrefix))
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        ("(" + BDAE.string(op) + " " + expToJL(e1, simCode, varPrefix=varPrefix) + ")")
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expToJL(e1, simCode, varPrefix=varPrefix) + " " + BDAE.string(op) + " " + expToJL(e2, simCode, varPrefix=varPrefix))
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        (expToJL(e1, simCode, varPrefix=varPrefix) + " " + BDAE.string(op) + " " + expToJL(e2, simCode,varPrefix=varPrefix))
      end
      #=TODO?=#
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        "if" + expToJL(e1, simCode, varPrefix=varPrefix) + "\n" + expToJL(e2, simCode,varPrefix=varPrefix) + "else\n" + expToJL(e3, simCode,varPrefix=varPrefix) + "\nend"
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = expl)  => begin
        #=
          TODO: Keeping it simple for now, we assume we only have one argument in the call
          We handle derivitives seperatly
        =#
        varName = BDAE.string(listHead(expl))
        (index, type) = hashTable[varName]
        @match tmpStr begin
          "der" => "dx[$index]  #= der($varName) =#"
          "pre" => begin
            indexForVar = hashTable[varName][1]
            "(integrator.uprev[$(indexForVar)])"
          end
          "edge" =>  begin
             indexForVar = hashTable[varName][1]
             string(tuple(map((x) -> expToJL(x, simCode, varPrefix=varPrefix), expl)...)...) + " && ! integrator.uprev[$(indexForVar)]"
          end
          _  =>  begin
            tmpStr *= string(tuple(map((x) -> expToJL(x, simCode, varPrefix=varPrefix), expl)...)...)
          end
        end
      end
      DAE.CAST(exp = e1)  => begin
         expToJL(e1, simCode)
      end
      DAE.ARRAY(DAE.T_ARRAY(DAE.T_BOOL(__)), scalar, array) => begin
        local arrayExp = "#= Array exp=# reduce(|, ["
        for e in array
          arrayExp *= expToJL(e, simCode) + ","
        end
        arrayExp *= "])"
      end
      _ =>  throw(ErrorException("$exp not yet supported"))
    end
  end
  return "(" + str + ")"
end

"""
  Generates the starting conditions
"""
function getStartConditions(vars::Array, condName::String, simCode::SimulationCode.SIM_CODE)
  local startCondStr::String = ""
  local ht::Dict = simCode.crefToSimVarHT
  if length(vars) == 0
    return startCondStr
  end
  for var in vars
    (index, simVar) = ht[var]
    simVarType = simVar.varKind
    local optAttributes::Option{DAE.VariableAttributes} = simVar.attributes
    @debug "simVar.attributes: $simVar.attributes"
    if simVar.attributes == nothing
      continue
    end
    _ = @match optAttributes begin
        SOME(attributes) => begin
          startCondStr *= @match attributes.start begin
            SOME(start) => begin
              @debug "Start value is:" start
              "$condName[$index] #= $var =# = $(expToJL(start, simCode))\n"
            end
            NONE() => ""
          end
          NONE() => ""
        end
    end
  end
  @debug "Returning $startCondStr"
  return startCondStr
end

end #= End CodeGeneration=#
