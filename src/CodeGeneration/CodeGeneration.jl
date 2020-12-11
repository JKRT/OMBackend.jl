#=
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
  Author: John Tinnerholm, john.tinnerholm@liu.se
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
$(copyRightString())"

"Counter for generating callbacks for a model"
global CALLBACK_COUNTER = 0

"""
  Write  modelName_DAE_equations() to file.
"""
function writeDAE_equationsToFile(fileName::String, contents::String)
  local fdesc = open(fileName, "w")
  write(fdesc, contents)
  close(fdesc)
end

"
  ODE-mode code generation
"
function generateCode(simCode::SimulationCode.EXPLICIT_SIM_CODE)
  throw("Not implemented..")
end


"""
  Generate Julia code from SimCode.
Returns an expression in memory representing the program and a string
"""
function generateCode(simCode::SimulationCode.SIM_CODE)::Tuple{String, Expr}
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
  local START_CONDTIONS = quote
      function $(Symbol("$(modelName)StartConditions"))(p, t0)
        local x0 = zeros($(length(stateVariables) + length(algVariables)))
        local dx0 = zeros($(length(stateVariables) + length(algVariables)))
        $START_CONDTIONS_EQUATIONS
        return x0, dx0
      end
    end
    local DAE_EQUATION_FUNCTION = quote
        function $(Symbol("$(modelName)DAE_equations"))(res, dx, x, p, t)
          $(DAE_EQUATIONS...)
        end
      end
      local DIFFERENTIAL_VARS_FUNCTION = quote
        begin
          function $(Symbol("$(modelName)DifferentialVars"))()
            return $stateMarkings
          end
        end
      end
  local PARAMETER_VARS = quote
    function $(Symbol("$(modelName)ParameterVars"))()
      p = Array{Float64}(undef, $(arrayLength(parameters)))
      $(PARAMETER_EQUATIONS...)
      return p
    end
  end
  local WHEN_EQUATIONS_FN = quote
    function $(Symbol("$(modelName)CallbackSet"))(p)
      $(LineNumberNode((@__LINE__), "WHEN EQUATIONS"))
      $(WHEN_EQUATIONS...)
      $(LineNumberNode((@__LINE__), "IF EQUATIONS"))
      $(IF_EQUATIONS...)
      return CallbackSet($(evaluateCallBackset()))
    end
  end
  local RUNNABLE = quote
    function $(Symbol("$(modelName)Simulate"))(tspan = (0.0, 1.0))
      # Define problem    
      p = $(Symbol("$(modelName)ParameterVars"))()
      (x0, dx0) =$(Symbol("$(modelName)StartConditions"))(p, tspan[1])
      differential_vars = $(Symbol("$(modelName)DifferentialVars"))()
      #= Pass the residual equations =#
      problem = DAEProblem($(Symbol("$(modelName)DAE_equations")), dx0, x0, 
                           tspan, p, differential_vars=differential_vars, 
                           callback=$(Symbol("$(modelName)CallbackSet"))(p))
      # Solve with IDA:)
      solution = solve(problem, IDA())
      return solution
    end
  end
  @debug("Code-generation done")
  global CALLBACK_COUNTER = 0
  program = quote    
    $(HEADER_STRING)
    using DiffEqBase
    using DifferentialEquations
    using Plots
    using Sundials
    $(START_CONDTIONS)
    $(DIFFERENTIAL_VARS_FUNCTION)
    $(DAE_EQUATION_FUNCTION)
    $(PARAMETER_VARS)
    $(WHEN_EQUATIONS_FN)
    $(RUNNABLE)
  end
  # Return file content
  return ("$(modelName)", program)
end

function evaluateCallBackset()::Expr
  local cbs::Array{Expr} = []
  for t in 1:CALLBACK_COUNTER
    push!(cbs,
          quote
          $(Symbol("cb$(t)"))
          end)
  end
  expr::Expr = if length(cbs) < 1
    quote end
  else
    quote
      $((cbs...))
    end
  end
  return expr 
end

"""
    TODO: We use state for everything now. 
    I am unsure how to differentiate between the two
"""
function createStartConditionsEquations(algVariables::Array,
                                        stateVariables::Array,
                                        simCode::SimulationCode.SIM_CODE)::Expr
  return quote
    $(getStartConditions(algVariables, "x0", simCode))
    $(getStartConditions(stateVariables, "x0", simCode)) #Should be DX0?
  end
end

function createResidualEquations(equations::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  local eqs::Array{Expr} = []
  for (equationCounter, eq) in enumerate(equations)
    push!(eqs, residualEqtoJulia(eq, simCode, equationCounter))
  end
  return eqs
end

function createEquations(equations::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  local eqs::Array{Expr} = []
  for (equationCounter, eq) in enumerate(equations)
    push!(eqs, eqToJulia(eq, simCode, equationCounter))
  end
  return eqs
end

function createParameterEquations(parameters::Array, simCode::SimulationCode.SIM_CODE)
  local parameterEquations::Array = []
  local hT = simCode.crefToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    local simVarType::SimulationCode.SimVarType = simVar.varKind
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => exp
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
    push!(parameterEquations,
          quote
          $(LineNumberNode(@__LINE__, "$param"))
          p[$index] = $(expToJuliaExp(bindExp, simCode))
          end
          )
  end
  return parameterEquations
end

"""
  Transforms a given equation into Julia code
"""
function residualEqtoJulia(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, resNumber)::Expr
  local rhs::DAE.Exp
  local result::Expr = @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
      quote
        res[$(resNumber)] = $(expToJuliaExp(rhs, simCode))
      end
    end
    _ => begin
      throw("traversalError for $eq")
    end
  end
  return result
end

"""
    Transforms a given equation into Julia code
    The optional parameter prefix specifices the prefix.
    The prefix is used for...
"""
function eqToJulia(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, resNumber; prefix="res")::Expr
  local rhs::DAE.Exp
  @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
      @info "Value of eqToJulia $prefix"
      "  $(prefix)[$resNumber] = " + expToJuliaExp(rhs, simCode; varPrefix = prefix) + "\n"
    end
    BDAE.WHEN_EQUATION(_, wEq, _) => begin
      local whenStmts = createWhenStatements(wEq.whenStmtLst, simCode, prefix="integrator.u")
      local cond = prepForZeroCrossing(wEq.condition)
      global CALLBACK_COUNTER += 1
      quote
        function $(Symbol("condition$(CALLBACK_COUNTER)"))(x,t,integrator) 
          $(expToJuliaExp(cond, simCode))
        end
        function $(Symbol("affect$(CALLBACK_COUNTER)!"))(integrator)
          $(whenStmts...)
        end
        $(Symbol("cb$(CALLBACK_COUNTER)")) = ContinuousCallback($(Symbol("condition$(CALLBACK_COUNTER)")), 
                                                                $(Symbol("affect$(CALLBACK_COUNTER)!")),
                                                                rootfind=true, save_positions=(false,true),
                                                                affect_neg! = $(Symbol("affect$(CALLBACK_COUNTER)!")),)
      end
    end
    BDAE.IF_EQUATION(conds, trueEqs, falseEqs) => begin
      global CALLBACK_COUNTER += 1
      f(x) = eqToJulia(x, simCode, resNumber, prefix = "integrator.u")
      local condsJL = tuple(map((x) -> expToJuliaExp(x, simCode; varPrefix="x"), conds)...)
      local trueEqJL = map(f, trueEqs)
      local falseEqJL = map(f, falseEqs)
      quote 
        function $(Symbol("condition$(CALLBACK_COUNTER)"))(x,t,integrator)
          Bool(floor($(condsJL...)))
        end
        function $(Symbol("affect$(CALLBACK_COUNTER)!"))(integrator)
          $(trueEqJL...)
        end
        $(Symbol("cb$(CALLBACK_COUNTER)")) = ContinuousCallback(Symbol("condition$(CALLBACK_COUNTER)"),
                                                            Symbol("affect$(CALLBACK_COUNTER)!"))
        $(global CALLBACK_COUNTER += 1)
        function $(Symbol("condition$(CALLBACK_COUNTER)"))(x,t,integrator)
           ! Bool(floor($(condsJL...)))
         end
         function $(Symbol("affect$(CALLBACK_COUNTER)!"))(integrator)
           $(falseEqJL...)
         end
         createSymbol("cb$(CALLBACK_COUNTER)") = ContinuousCallback(
           Symbol("condition$(CALLBACK_COUNTER)"),
           Symbol("affect$(CALLBACK_COUNTER)!"))
      end
    end
    _ => begin
      throw("Error when generating code for $eq")
    end
  end
end


"""
  Creates Julia code for  when equation.
  E.g a set of Discrete callbacks
  TODO: Handle the other cases
"""
function createWhenStatements(whenStatements::List, simCode::SimulationCode.SIM_CODE; prefix="x")::Array{Expr}
  local res::Array{Expr} = []
  @debug "Calling createWhenStatements with: $whenStatements"
  for wStmt in  whenStatements
    @match wStmt begin
      BDAE.ASSIGN(__) => begin
        push!(res, quote
              $(expToJuliaExp(wStmt.left, simCode, varPrefix=prefix)) = $(expToJuliaExp(wStmt.right, simCode, varPrefix=prefix))
              end)
      end
      BDAE.REINIT(__) => begin
        (index, type) = simCode.crefToSimVarHT[BDAE.string(wStmt.stateVar)]
        push!(res, quote 
              integrator.u[$(index)] = $(expToJuliaExp(wStmt.value, simCode, varPrefix=prefix))
              end)
      end
      _ => ErrorException("$whenStatements in @__FUNCTION__ not supported")
    end
  end
  return res
end


"Converts a DAE expression into a Julia expression"
function expToJuliaExp(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE; varPrefix="x")::Expr
  hashTable = simCode.crefToSimVarHT
  local expr::Expr = begin
    local int::Int64
    local real::Float64
    local bool::Bool
    local tmpStr::String
    local cr::DAE.ComponentRef
    local e1::DAE.Exp
    local e2::DAE.Exp
    local e3::DAE.Exp
    local expl::List{DAE.Exp}
    local lstexpl::List{List{DAE.Exp}}
    @match exp begin
      DAE.BCONST(bool) => quote $bool end
      DAE.ICONST(int) => quote $int end
      DAE.RCONST(real) => quote $real end
      DAE.SCONST(tmpStr) => quote $tmpStr end
      DAE.CREF(cr, _)  => begin
        varName = BDAE.string(cr)
        builtin = if varName == "time"
          true
        else
          false
        end
        if ! builtin
          #= If we refeer to time, we simply return t instead of a concrete variable =#
          indexAndVar = hashTable[varName]
          varKind::SimulationCode.SimVarType = indexAndVar[2].varKind
          @match varKind begin
            SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
            SimulationCode.STATE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName state"))
              $(Symbol(varPrefix))[$(indexAndVar[1])]
            end
            SimulationCode.PARAMETER(__) => quote
              $(LineNumberNode(@__LINE__, "$varName parameter"))
              p[$(indexAndVar[1])]
              end
            SimulationCode.ALG_VARIABLE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName, algebraic"))
              $(Symbol(varPrefix))[$(indexAndVar[1])]
            end
            SimulationCode.STATE_DERIVATIVE(__) => :(dx[$(indexAndVar[1])] #= der($varName) =#)
          end
        else #= Currently only time is a builtin variabe. Time is represented as t in the generated code =#
         quote
            t
         end
        end
      end
      DAE.UNARY(operator = op, exp = e1) => begin
        o = DAE_OP_toJuliaOperator(op)
        quote
          $(o)($(expToJuliaExp(e1, simCode)))
        end
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        a = expToJuliaExp(e1, simCode, varPrefix=varPrefix)
        b = expToJuliaExp(e2, simCode, varPrefix=varPrefix)
        o = DAE_OP_toJuliaOperator(op)
        quote
          $o($(a), $(b))
        end
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        quote
          $("(" + BDAE.string(op) + " " + expToJuliaExp(e1, simCode, varPrefix=varPrefix) + ")")
        end
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExp(e1, simCode, varPrefix=varPrefix) + " " + BDAE.string(op) + " " + expToJuliaExp(e2, simCode, varPrefix=varPrefix))
        end
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExp(e1, simCode,
                          varPrefix=varPrefix) + " "
            + BDAE.string(op) + " " + expToJuliaExp(e2, simCode,varPrefix=varPrefix))
        end
      end
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        throw(ErrorException("If expressions not allowed in backend code"))
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst)  => begin
        DAECallExpressionToJuliaCallExpression(tmpStr, explst, simCode, hashTable, varPrefix=varPrefix)
      end
      _ =>  throw(ErrorException("$exp not yet supported"))
    end
  end
  return expr
end


"""
    Generates the start conditions
  TODO:
    Fix generation here further.. intermediate strings are probably not necessary.. 
"""
function getStartConditions(vars::Array, condName::String, simCode::SimulationCode.SIM_CODE)::Expr
  local startExprs::Array{Expr} = []
  local ht::Dict = simCode.crefToSimVarHT
  if length(vars) == 0
    return quote
    end
  end
  for var in vars
    (index, simVar) = ht[var]
    local simVarType = simVar.varKind
    local optAttributes::Option{DAE.VariableAttributes} = simVar.attributes
    if simVar.attributes == nothing
      continue
    end
    () = @match optAttributes begin
      SOME(attributes) => begin
        () = @match (attributes.start, attributes.fixed) begin
          (SOME(start), SOME(fixed)) || (SOME(start), _)  => begin
            @debug "Start value is:" start
            push!(startExprs,
                  quote
                  $(LineNumberNode(@__LINE__, "$var"))
                  $(Symbol("$condName"))[$index] = $(expToJuliaExp(start, simCode))
                  end)
            ()
          end
          (NONE(), SOME(fixed)) => begin
            push!(startExprs, :($(condName)[$(index)] = 0.0))
            ()
          end
          (_, _) => ()
        end
        NONE() => ()
      end
    end
  end
  return quote
    $(startExprs...)
  end
end

end #= End CodeGeneration=#
