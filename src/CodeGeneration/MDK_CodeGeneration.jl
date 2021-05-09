#=
  Author: John Tinnerholm

TODO: Remember the state derivative scheme. What the heck did I mean  with that? 
TODO: Make duplicate code better...
=#
using ModelingToolkit
import OMBackend
"""
  Generates simulation code targetting modeling toolkit
"""
function generateMDKCode(simCode::SimulationCode.SIM_CODE)
  if OMBackend.MODE == MODELING_TOOLKIT_DAE_MODE
    DAE_MODE_MDK(simCode::SimulationCode.SIM_CODE)
  else
    ODE_MODE_MDK(simCode::SimulationCode.SIM_CODE)
  end
end

function ODE_MODE_MDK(simCode::SimulationCode.SIM_CODE)
  @debug "Runnning: generateMDKCode"
  local crefToSimVarHT = simCode.crefToSimVarHT
  local equations::Array = []
  local exp::DAE.Exp
  local modelName::String = simCode.name
  local parameters::Array = []
  local stateDerivatives::Array = []
  local stateMarkings::Array = []
  local stateVariables::Array = []
  local algebraicVariables::Array = []
  for varName in keys(crefToSimVarHT)
    ixAndVar = crefToSimVarHT[varName]
    local varType = ixAndVar[2].varKind
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algebraicVariables #= Should be alg if more specialised scheme is used.=#, varName)
      #=TODO: Do I need to modify this?=#
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end
  equations = simCode.residualEquations
  equations = createResidualEquationsMDK(stateVariables, algebraicVariables, equations, simCode::SimulationCode.SIM_CODE)
  #= Symbolic names =#
  local stateVariablesSym = [:($(Symbol(v))(t)) for v in stateVariables]
  local algebraicVariablesSym = [Symbol(v) for v in algebraicVariables]
  local parVariablesSym = [Symbol(p) for p in parameters]
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquationsMDK(stateVariables, algebraicVariables, simCode)
  local PARAMETER_EQUATIONS = createParameterEquationsMDK(parameters, simCode)
  #= Formulate the problem as a DAE Problem=#
  program = quote
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations
    function $(Symbol("$(modelName)Model"))(tspan = (0.0, 1.0))
      pars = ModelingToolkit.@parameters begin
        ($(parVariablesSym...), t)
      end
      vars = ModelingToolkit.@variables begin
        ($(stateVariablesSym...), $(algebraicVariablesSym...))
      end
      der = Differential(t)
      eqs = [$(equations...)]
      nonLinearSystem = ModelingToolkit.ODESystem(eqs, t, vars, pars, name=:($(Symbol($modelName))))
      pars = Dict($(PARAMETER_EQUATIONS...), t => tspan[1])
      initialValues = [$(START_CONDTIONS_EQUATIONS...)]
      firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
      reducedSystem = ModelingToolkit.dae_index_lowering(firstOrderSystem)
      problem = ModelingToolkit.ODEProblem(reducedSystem, initialValues, tspan, pars)
      return problem
    end
    $(Symbol("$(modelName)Model_problem")) = $(Symbol("$(modelName)Model"))()
    function $(Symbol("$(modelName)Simulate"))(tspan = (0.0, 1.0))     
      solve($(Symbol("$(modelName)Model_problem"), radau5()))
    end
  end
  return (modelName, program)
end

function DAE_MODE_MDK(simCode::SimulationCode.SIM_CODE)
  @debug "Runnning: generateMDKCode"
  local crefToSimVarHT = simCode.crefToSimVarHT
  local equations::Array = []
  local exp::DAE.Exp
  local modelName::String = simCode.name
  local parameters::Array = []
  local stateDerivatives::Array = []
  local stateMarkings::Array = []
  local stateVariables::Array = []
  local algebraicVariables::Array = []
  for varName in keys(crefToSimVarHT)
    ixAndVar = crefToSimVarHT[varName]
    local varType = ixAndVar[2].varKind
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algebraicVariables #= Should be alg if more specialised scheme is used.=#, varName)
      #=TODO: Do I need to modify this?=#
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end
  equations = simCode.residualEquations
  equations = createResidualEquationsMDK(stateVariables, algebraicVariables, equations, simCode::SimulationCode.SIM_CODE)
  #= Symbolic names =#
  local stateVariablesSym = [:($(Symbol(v))(t)) for v in stateVariables]
  local algebraicVariablesSym = [Symbol(v) for v in algebraicVariables]
  local parVariablesSym = [Symbol(p) for p in parameters]
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquationsMDK(stateVariables, algebraicVariables, simCode)
  local PARAMETER_EQUATIONS = createParameterEquationsMDK(parameters, simCode)
  #= Formulate the problem as a DAE Problem=#
  program = quote
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations
    function $(Symbol("$(modelName)Model"))(tspan = (0.0, 1.0))
      pars = ModelingToolkit.@parameters begin
        ($(parVariablesSym...), t)
      end
      vars = ModelingToolkit.@variables begin
        ($(stateVariablesSym...), $(algebraicVariablesSym...))
      end
      der = Differential(t)
      eqs = [$(equations...)]
      nonLinearSystem = ModelingToolkit.ODESystem(eqs, t, vars, pars, name=:($(Symbol($modelName))))
      pars = Dict($(PARAMETER_EQUATIONS...), t => tspan[1])
      initialValues = [$(START_CONDTIONS_EQUATIONS...)]
      firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
      problem = ModelingToolkit.ODEProblem(firstOrderSystem, initialValues, tspan, pars)
      return problem
    end
    $(Symbol("$(modelName)Model_problem")) = $(Symbol("$(modelName)Model"))()
    function $(Symbol("$(modelName)Simulate"))(tspan = (0.0, 1.0))     
      solve($(Symbol("$(modelName)Model_problem")))
    end
  end
  return (modelName, program)
end


"
 Creates the residual equations in unsorted order
"
function createResidualEquationsMDK(stateVariables::Array, algebraicVariables::Array, equations::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  local eqs::Array{Expr} = []
  local ht = simCode.crefToSimVarHT
  local lhsVariables = Set([])
  for i in 1:length(stateVariables)
    local variableIdx = ht[stateVariables[i]][1]
    local equationIdx = simCode.matchOrder[variableIdx]
    local equation = equations[equationIdx]
    push!(eqs, residualEqtoJuliaMDK(equation, simCode, equationIdx))
  end
  #= Generate algebraic equations=#
  for var in algebraicVariables
    @info var 
    local variableIdx = ht[var][1]
    local equationIdx = simCode.matchOrder[variableIdx]
    local equation = equations[equationIdx]
    local eqDAEExp = equation.exp
    local eqExp = quote 
        0 ~ $(stripBeginBlocks(expToJuliaExpMDK(eqDAEExp, simCode ;derSymbol=true)))
     end
    push!(eqs, eqExp)
  end
  return eqs
end

"Converts a DAE expression into a Julia expression"
function expToJuliaExpMDK(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE, varSuffix=""; varPrefix="x", derSymbol::Bool=false)::Expr
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
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.PARAMETER(__) => quote
              $(LineNumberNode(@__LINE__, "$varName parameter"))
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.ALG_VARIABLE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName, algebraic"))
              $(Symbol(indexAndVar[2].name))
            end
            _ => begin
              @error "Unsupported varKind: $(varKind)"
              throw()
            end
          end
        else #= Currently only time is a builtin variable. Time is represented as t in the generated code =#
          quote
            t
          end
        end
      end
      DAE.UNARY(operator = op, exp = e1) => begin
        o = DAE_OP_toJuliaOperator(op)
        quote
          $(o)($(expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix)))
        end
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        a = expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix, derSymbol = derSymbol)
        b = expToJuliaExpMDK(e2, simCode, varPrefix=varPrefix, derSymbol = derSymbol)
        o = DAE_OP_toJuliaOperator(op)
        quote
          $o($(a), $(b))
        end
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        quote
          $("(" + BDAE.string(op) + " " + expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix, derSymbol = derSymbol) + ")")
        end
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix, derSymbol = derSymbol) + " " + BDAE.string(op) + " "
            + expToJuliaExpMDK(e2, simCode, varPrefix=varPrefix, derSymbol = derSymbol))
        end
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExpMDK(e1, simCode,
                             varPrefix=varPrefix,derSymbol = derSymbol)
            + " "
            + BDAE.string(op)
            + " "
            + expToJuliaExpMDK(e2, simCode,varPrefix=varPrefix, derSymbol=derSymbol))
        end
      end
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        throw(ErrorException("If expressions not allowed in backend code"))
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst)  => begin
        #Call as symbol is really ugly.. please fix me :(
        DAECallExpressionToMDKCallExpression(tmpStr, explst, simCode, hashTable; varPrefix=varPrefix, derAsSymbol=derSymbol)
      end
      DAE.CAST(ty, exp)  => begin
        quote
          $(generateCastExpressionMDK(ty, exp, simCode, varPrefix))
        end
      end
      _ =>  throw(ErrorException("$exp not yet supported"))
    end
  end
  return expr
end


"""
  Transforms a given equation into MDK Julia code
"""
function residualEqtoJuliaMDK(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, equationIdx::Int64)::Expr
  local ht = SimulationCode.makeIndexVarNameDict(simCode.matchOrder, simCode.crefToSimVarHT)
  @debug "Called residual equation to Julia"
  local result::Expr = @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = daeExp) => begin
      #= TODO: Using Symbolics to rewrite equations... Same limitations apply here as well=#
      local variableIdx = MetaGraphs.get_prop(simCode.equationGraph, equationIdx, :vID)
      local varToSolve = simCode.crefToSimVarHT[ht[variableIdx]][2]
      local derSymbolString = "der_$(varToSolve.name)"
      local solveForState = SimulationCode.isState(varToSolve)
      local varToSolveExpr::Symbol =  solveForState ? Symbol(derSymbolString) : Symbol(varToSolve.name)
      local daeVariables = getVariablesInDAE_Exp(daeExp, simCode, Set([]))
      #= Simple quation rewrite..=#
      eqRewrite = function ()
        eqExpr = quote
          ModelingToolkit.@variables begin
            $(daeVariables...)
            $(varToSolveExpr)
          end
          leq = 0 ~ $(stripBeginBlocks(expToJuliaExpMDK(daeExp, simCode ;derSymbol=true)))
          Symbolics.solve_for([leq], [$varToSolveExpr]; check = false #=Does not check for linearity.. =#)
        end
        local res = @eval  $eqExpr
        return res
      end
      local evEqExpr = eqRewrite()
      if solveForState
        local sol = quote
          der($(Symbol(varToSolve.name))) ~ $(evEqExpr[1])
        end
      else
        local sol = quote
          $(Symbol(varToSolve.name)) ~ $(evEqExpr[1])
        end
      end
      #= Return sol. Assign to results =#
      sol
    end
    _ => begin
      throw("traversalError for $eq")
    end
  end
  return result
end


"""
    Generates the initial value for the equations
    TODO: Currently unable to generate start condition in order 
"""
function createStartConditionsEquationsMDK(states::Array,
                                        algebraics::Array,
                                        simCode::SimulationCode.SimCode)::Array{Expr}
  #  local startEquations = createSortedEquations([algVariables..., stateVariables...], simCode; arrayName = "reals")
  local algInit = getStartConditionsMDK(algebraics, simCode)
  local stateInit = getStartConditionsMDK(states, simCode)
  return vcat(algInit, stateInit)
end



function getStartConditionsMDK(vars::Array, simCode::SimulationCode.SimCode)::Array{Expr}
  local startExprs::Array{Expr} = []
  local residuals = simCode.residualEquations
  local ht::Dict = simCode.crefToSimVarHT
  if length(vars) == 0
    return []
  end
  for var in vars
    (index, simVar) = ht[var]
    varName = simVar.name
    local simVarType = simVar.varKind
    local optAttributes::Option{DAE.VariableAttributes} = simVar.attributes
    #= If no attribute. Let it default to zero =#
    if simVar.attributes == nothing
      push!(startExprs, :($(Symbol(varName)) => 0.0))
      continue
    end
    () = @match optAttributes begin
      SOME(attributes) => begin
        () = @match (attributes.start, attributes.fixed) begin
          (SOME(DAE.CREF(start)), SOME(__)) || (SOME(DAE.CREF(start)), _)  => begin
            push!(startExprs, :($(Symbol("$varName")) => $exp))            
            #= We have a variable of sorts =#
            push!(startExprs,
                  quote
                  $(Symbol("$varName")) => pars[$(Symbol(start.ident))]
                  end)
            ()          
          end            
          (SOME(start), SOME(fixed)) || (SOME(start), _)  => begin
            push!(startExprs,
                  quote
                  $(Symbol("$varName")) => $(expToJuliaExpMDK(start, simCode))
                  end)
            ()
          end
          (NONE(), SOME(fixed)) => begin
            push!(startExprs, :($(varName) => 0.0))
            ()
          end
          (_, _) => ()
        end
        NONE() => ()
      end
    end
  end
  return startExprs
end


function createParameterEquationsMDK(parameters::Array, simCode::SimulationCode.SimCode)::Array{Expr}
  local parameterEquations::Array = []
  local hT = simCode.crefToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    local simVarType::SimulationCode.SimVarType = simVar.varKind
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => exp
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
    #= Solution for https://github.com/SciML/ModelingToolkit.jl/issues/991 =#
    push!(parameterEquations,
          quote
          $(LineNumberNode(@__LINE__, "$param eq"))
          $(Symbol(simVar.name)) => float($((expToJuliaExpMDK(bindExp, simCode))))
          end
          )
  end
  return parameterEquations
end


function generateCastExpressionMDK(ty, exp, simCode, varPrefix)
  return @match ty, exp begin
    (DAE.T_REAL(__), DAE.ICONST(__)) => float(eval(expToJuliaExpMDK(exp, simCode, varPrefix=varPrefix)))
    (DAE.T_REAL(__), DAE.CREF(cref)) where typeof(cref.identType) == DAE.T_INTEGER  => begin
      quote
        float($(expToJuliaExpMDK(exp, simCode, varPrefix=varPrefix)))
      end
    end
    _ => throw("Cast $ty: for exp: $exp not yet supported in codegen!")
  end
end
