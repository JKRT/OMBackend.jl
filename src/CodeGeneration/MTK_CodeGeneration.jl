#=
  Author: John Tinnerholm

TODO: Remember the state derivative scheme. What the heck did I mean  with that?
TODO: Make duplicate code better...
TODO: Make this into it's own module
TODO: Fix the redundant string conversion scheme
TODO: Investigate broken if equations
=#
using ModelingToolkit
import OMBackend

"""
  Generates simulation code targetting modeling toolkit.
  Loop code removed was on old branch.
"""
function generateMTKCode(simCode::SimulationCode.SIM_CODE)
  isCycles = isCycleInSCCs(simCode.stronglyConnectedComponents)
  ODE_MODE_MTK(simCode::SimulationCode.SIM_CODE)
end

function ODE_MODE_MTK(simCode::SimulationCode.SIM_CODE)
  @debug "Runnning: generateMTKCode"
  RESET_CALLBACKS()
  local stringToSimVarHT = simCode.stringToSimVarHT
  local equations::Array = BDAE.RESIDUAL_EQUATION[]
  local exp::DAE.Exp
  local modelName::String = simCode.name
  local parameters::Vector = String[]
  local stateDerivatives::Vector = String[]
  local stateVariables::Vector = String[]
  local algebraicVariables::Vector = String[]
  local discreteVariables::Vector = String[]
  local performIndexReduction = false
  for varName in keys(stringToSimVarHT)
    (idx, var) = stringToSimVarHT[varName]
    local varType = var.varKind
    local varType = var.varKind
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => begin
        #=TODO: Seems to be unable to find variables in some cases...
        should there be an and here? Keep track of variables that are also involved in if/when structures
        =#
        if idx in simCode.matchOrder
          push!(algebraicVariables, varName)
        elseif involvedInEvent(idx, simCode) #= We have a variable that is not contained in continious system =#
          #= Treat discrete variables separate =#
          push!(discreteVariables, varName)
        elseif simCode.isSingular
          #=
            If the variable is not involved in an event and the index is not in match order and
            the system is singular.
            This means that the variable is probably algebraic however, we need to perform index reduction.
          =#
          push!(algebraicVariables, varName)
        end
      end
      #=TODO:johti17 Do I need to modify this?=#
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end
  performIndexReduction = simCode.isSingular
  #= Create equations for variables not in a loop + parameters and stuff=#
  local equations = createResidualEquationsMTK(stateVariables,
                                               algebraicVariables,
                                               simCode.residualEquations,
                                               simCode::SimulationCode.SIM_CODE)
  local parVariablesSym = [Symbol(p) for p in parameters]
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquationsMTK(stateVariables, algebraicVariables, simCode)
  local DISCRETE_START_VALUES = getStartConditionsMTK(discreteVariables, simCode)
  local PARAMETER_EQUATIONS = createParameterEquationsMTK(parameters, simCode)
  local PARAMETER_RAW_ARRAY = createParameterArray(parameters, simCode)
  #= Create callback equations. For MTK we disable the saving function for now. =#
  local CALL_BACK_EQUATIONS = createCallbackCode(modelName, simCode; generateSaveFunction = true)
  #= Symbolic names =#
  local stateVariablesSym = [:($(Symbol(v))(t)) for v in stateVariables]
  local algebraicVariablesSym = [:($(Symbol(v))(t)) for v in algebraicVariables]
  local discreteVariablesSym = [:($(Symbol(v))) for v in discreteVariables]
  #= Reset the callback counter=#
  RESET_CALLBACKS()

  #=
  Formulate the problem as a DAE Problem.
  For this variant we keep t on its own line
  https://github.com/SciML/ModelingToolkit.jl/issues/998
  =#
  #=If our model name is separated by . replace it with __ =#
  local MODEL_NAME = replace(modelName, "." => "__")
  program = quote
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations

    function $(Symbol(MODEL_NAME * "Model"))(tspan = (0.0, 1.0))
      @variables t
      parameters = ModelingToolkit.@parameters begin
        ($(parVariablesSym...), $(discreteVariablesSym...))
      end
      #=
      Only variables that are present in the equation system later should be a part of the variables in the MTK system.
      This means that certain algebraic variables should not be listed among the variables (These are the discrete variables).
      =#
      vars = ModelingToolkit.@variables begin
        ($(stateVariablesSym...), $(algebraicVariablesSym...))
      end
      der = Differential(t)
      eqs = [$(equations...)]
      nonLinearSystem = ModelingToolkit.ODESystem(eqs, t, vars, parameters, name=:($(Symbol($modelName))))
      pars = Dict($(PARAMETER_EQUATIONS...), $(DISCRETE_START_VALUES...))
      #= Initial values for the continious system. =#
      initialValues = [$(START_CONDTIONS_EQUATIONS...)]
      firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
      $(MTK_indexReduction(performIndexReduction))
      #=
      These arrays are introduced to handle the bolted on event handling using callbacks.
      The callback handling for MTK is subject of change should hybrid system be implemented for MTK.
      =#
      local event_p = [$(PARAMETER_RAW_ARRAY...)]
      local discreteVars = collect(values(Dict([$(DISCRETE_START_VALUES...)])))
      local event_vars = vcat(collect(values(Dict([$(START_CONDTIONS_EQUATIONS...)]))),
                              #=Discrete variables=# discreteVars)
      local aux = Array{Array{Float64}}(undef, 2)
      aux[1] = event_p
      #= TODO init them with the initial values of them =#
      aux[2] = event_vars
      problem = ModelingToolkit.ODEProblem(reducedSystem,
                                           initialValues,
                                           tspan,
                                           pars,
                                           callback=$(Symbol("$(MODEL_NAME)CallbackSet"))(aux))
      return (problem, initialValues, reducedSystem, tspan, pars, vars)
    end
    $(CALL_BACK_EQUATIONS)
    ($(Symbol("$(MODEL_NAME)Model_problem")), _, _, _, _,_) = $(Symbol("$(MODEL_NAME)Model"))()
    function $(Symbol("$(MODEL_NAME)Simulate"))(tspan)
      solve($(Symbol("$(MODEL_NAME)Model_problem")), tspan=tspan)
    end
    function $(Symbol("$(MODEL_NAME)Simulate"))(tspan = (0.0, 1.0); solver=Rodas5())
      solve($(Symbol("$(MODEL_NAME)Model_problem")), tspan=tspan, solver)
    end
  end
  #= MODEL_NAME is preprocessed with . replaced with _=#
  return MODEL_NAME, program
end

"""
 This generates code targetting modeling toolkit (but separates algebraic loops from the rest of the system).
 Discrete variables are treated as parameters 
"""
function ODE_MODE_MTK_LOOP(simCode::SimulationCode.SIM_CODE)
  @debug "Runnning: generateMTKCode"
  RESET_CALLBACKS()
  local stringToSimVarHT = simCode.stringToSimVarHT
  local equations::Array = []
  local exp::DAE.Exp
  local modelName::String = simCode.name
  local parameters::Array = []
  local stateDerivatives::Array = []
  local stateVariables::Array = []
  local algebraicVariables::Array = []
  local discreteVariables::Vector = []
  #= Loop arrays=#
  local stateDerivativesLoop::Array = []
  local stateVariablesLoop::Array = []
  local algebraicVariablesLoop::Array = []
  local discreteVariablesLoop::Vector = []
  
  #= 
  We have a cycle! 
  =#
  loop = getCycleInSCCs(simCode.stronglyConnectedComponents)
  
  for varName in keys(stringToSimVarHT)
    (idx, var) = stringToSimVarHT[varName]
    if simCode.matchOrder[idx] in loop
      local varType = var.varKind
      @match varType  begin
        SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
        SimulationCode.STATE(__) => push!(stateVariablesLoop, varName)
        SimulationCode.PARAMETER(__) => push!(parameters, varName)
        SimulationCode.ALG_VARIABLE(__) => begin
          if idx in simCode.matchOrder
            push!(algebraicVariablesLoop, varName)
          else #= We have a variable that is not contained in continious system =#
            #= Treat discrete variables separate =#
            push!(discreteVariablesLoop, varName)
          end
        end
        #=TODO: Do I need to modify this?=#
        SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivativesLoop, varName)
      end
    else #= Someplace else=#
      local varType = var.varKind
      @match varType  begin
        SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
        SimulationCode.STATE(__) => push!(stateVariables, varName)
        SimulationCode.PARAMETER(__) => push!(parameters, varName)
        SimulationCode.ALG_VARIABLE(__) => begin
          if idx in simCode.matchOrder
            push!(algebraicVariables, varName)
          else #= We have a variable that is not contained in continious system =#
            #= Treat discrete variables separate =#
            push!(discreteVariables, varName)
          end
        end
        #=TODO: Do I need to modify this?=#
        SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
      end
    end
  end
  local allVariables = vcat(stateVariables, algebraicVariables, stateVariablesLoop, algebraicVariablesLoop)
  #= Create equations for variables not in a loop + parameters and stuff=#
  local equations = createSortedEquations(stateVariables,
                                          simCode; arrayName = "u")
  local parVariablesSym = [Symbol(p) for p in parameters]  
  local START_CONDTIONS_EQUATIONS = getStartConditions(allVariables, "reals", simCode)
  
  local DISCRETE_START_VALUES = getStartConditionsMTK(discreteVariables, simCode)
  local PARAMETER_EQUATIONS = createParameterEquationsMTK(parameters, simCode)
  local PARAMETER_RAW_ARRAY = createParameterArray(parameters, simCode)
  #= Create callback equations. For MTK we disable the saving function for now. =#
  local CALL_BACK_EQUATIONS = createCallbackCode(modelName, simCode; generateSaveFunction = true)
  #= Symbolic names =#
  
  local stateVariablesSym = [:($(Symbol(v))(t)) for v in stateVariables]
  local algebraicVariablesSym = [:($(Symbol(v))(t)) for v in algebraicVariables]

  local stateVariablesSymNoT = [:($(Symbol(v))) for v in stateVariables]
  local algebraicVariablesSymNoT = [:($(Symbol(v))) for v in algebraicVariables]

  DUMMY_EQUATION_OUTSIDE_LOOP = Expr[]
  for v in vcat(stateVariablesSymNoT, algebraicVariablesSymNoT)
    push!(DUMMY_EQUATION_OUTSIDE_LOOP, :($v ~ $v))
  end
  
  local discreteVariablesSym = [:($(Symbol(v))) for v in discreteVariables]
  RESET_CALLBACKS()
  #= Create loop variables and the loop itself=#
  local equations_loop = createResidualEquationsMTK(stateVariablesLoop,
                                                    algebraicVariablesLoop,
                                                    simCode.residualEquations,
                                                    simCode::SimulationCode.SIM_CODE)
  local algebraicVariablesSymLoop = [:($(Symbol(v))) for v in algebraicVariablesLoop]
  local stateVariablesSymLoop = [:($(Symbol(v))) for v in stateVariablesLoop]
  local START_CONDTIONS_EQUATIONS_LOOP = createStartConditionsEquationsMTK(stateVariablesLoop, algebraicVariablesLoop, simCode)
  #= 
  Formulate the problem as a DAE Problem.ac
  For this variant we keep t on its own line 
  https://github.com/SciML/ModelingToolkit.jl/issues/998
  =#
  program = quote
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations
    using NLsolve
    #= t for time=#
    function $(Symbol(modelName, "AlgebraicLoop"))()
      parameters = ModelingToolkit.@parameters begin
        ($(parVariablesSym...), $(discreteVariablesSym...), t)
      end
      vars = ModelingToolkit.@variables begin
        ($(stateVariablesSym...), $(stateVariablesSymLoop...), $(algebraicVariablesSymLoop...))
      end
      eqs = [$(equations_loop...), $(DUMMY_EQUATION_OUTSIDE_LOOP...)]
      nonLinearSystem = ModelingToolkit.NonlinearSystem(eqs, vars, parameters, name=:($(Symbol($modelName))))
      return nonLinearSystem
    end

    function makeNLProblem()
      loop = $(Symbol(modelName, "AlgebraicLoop"))()
      nlsys_function = generate_function(loop, expression=Val{false})[2]
    end
    $(Symbol(modelName, "NonLinearFunction")) = makeNLProblem()
    function $(Symbol("$(modelName)Model_ODE"))(dx, x, aux, t)
      p = aux[1]
      u = aux[2]
      func!(res, u) = $(Symbol(modelName, "NonLinearFunction"))(res, u, vcat(p,[t]))
      sol = nlsolve(func!, u, ftol=1e-8; method = :newton)
      aux[2] = sol.zero
      #= Equations goes here =#
      $(equations...)
    end
    function $(Symbol("$(modelName)Model_problem"))(tspan)
      aux = [[$(PARAMETER_RAW_ARRAY...)], fill(0.0, $(length(allVariables)))]
      pars = aux[1]
      reals = aux[2]
      $(START_CONDTIONS_EQUATIONS)
      x = fill(0.0, $(length(stateVariables)))
      dx = fill(0.0, $(length(stateVariables)))
      odeF = $(Symbol("$(modelName)Model_ODE"))
      ODEProblem(odeF, x, tspan, aux)
    end
    
    $(CALL_BACK_EQUATIONS)
    $(Symbol("$(modelName)Model_problemInstance")) = $(Symbol("$(modelName)Model_problem"))((0.0,1.))
    function $(Symbol("$(modelName)Simulate"))(tspan)     
      solve($(Symbol("$(modelName)Model_problemInstance")), tspan=tspan)
    end
    function $(Symbol("$(modelName)Simulate"))(tspan = (0.0, 1.0); solver=Tsit5())
      solve($(Symbol("$(modelName)Model_problemInstance")), tspan=tspan, solver)
    end
  end
  
  return (modelName, program)

end

"""
     Creates the residual equations in unsorted order
  """
function createResidualEquationsMTK(stateVariables::Array, algebraicVariables::Array, equations::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  if isempty(equations)
    return Expr[]
  end
  local eqs::Vector{Expr} = Expr[]
  local ht = simCode.stringToSimVarHT
  local lhsVariables = Set([])
  local usedEquations = []
  for i in 1:length(stateVariables)
    local variableIdx = ht[stateVariables[i]][1]
    local equationIdx = simCode.matchOrder[variableIdx]
    local equation = equations[equationIdx]
    push!(usedEquations, equationIdx)
    push!(eqs, residualEqtoJuliaMTK(equation, simCode, equationIdx))
  end
  #=
  Generate algebraic equations.
  It might be the case that the algebraic variables are solved someplace else and not in
  the main set of equations. For instance a particular algebraic variable might be solved in some when equation.
  This sitaution is assumed to have been checked statically prior to this by the frontend.
  If there is no equation for which a algebraic variable is not solved we ignore it.
  =#
  local totalEquations = [i for i in 1:length(equations)]
  local remainingEquations = setdiff(totalEquations, usedEquations)
  #= Additional note. MTK does not require algebraic equations to follow a specific format.
  So we can just generate our remaining equations =#
  for eqIdx in remainingEquations
    local equation = equations[eqIdx]
    local eqDAEExp = equation.exp
    local eqExp = :(0 ~ $(stripBeginBlocks(expToJuliaExpMTK(eqDAEExp, simCode ;derSymbol=false))))
    push!(eqs, eqExp)
  end
  return eqs
end

"""
  Transforms a given equation into MTK Julia code.
  This function uses Symbolics to rewrite the equation into a matching format for MTK.
"""
function residualEqtoJuliaMTK(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, equationIdx::Int64)::Expr
  local ht = SimulationCode.makeIndexVarNameDict(simCode.matchOrder, simCode.stringToSimVarHT)
  @debug "Called residual equation to Julia"
  local result::Expr = @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = daeExp) => begin
      #= TODO: Using Symbolics to rewrite equations... Same limitations apply here as well=#
      local variableIdx = MetaGraphs.get_prop(simCode.equationGraph, equationIdx, :vID)
      local varToSolve = simCode.stringToSimVarHT[ht[variableIdx]][2]
      local derSymbolString = "der_" * varToSolve.name
      local solveForState = SimulationCode.isState(varToSolve)
      #= If we solve for state we solve for der_<var-name>=#
      #= Why did I do this. Ask Lennart! =#
      local varToSolveExpr::Symbol =  solveForState ? Symbol(derSymbolString) : Symbol(varToSolve.name)
#      local varToSolveExpr::Symbol = Symbol(varToSolve.name)
      #= These are the variables present in this specific equation. =#
      local daeVariables = getVariablesInDAE_Exp(daeExp, simCode, Set([]))
      #= Simple quation rewrite..=#
      eqRewrite = function ()
        eqExpr = quote
          ModelingToolkit.@variables begin
            $(daeVariables...)
            $(varToSolveExpr)
          end
          leq = 0 ~ $(stripBeginBlocks(expToJuliaExpMTK(daeExp, simCode ;derSymbol=true)))
          Symbolics.solve_for([leq], [$varToSolveExpr]; check = false #= Does not check for linearity!.. =#)
        end
        local res = @eval $eqExpr
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
function createStartConditionsEquationsMTK(states::Array,
                                        algebraics::Array,
                                        simCode::SimulationCode.SimCode)::Array{Expr}
  #  local startEquations = createSortedEquations([algVariables..., stateVariables...], simCode; arrayName = "reals")
  local algInit = getStartConditionsMTK(algebraics, simCode)
  local stateInit = getStartConditionsMTK(states, simCode)
  return vcat(algInit, stateInit)
end


function getStartConditionsMTK(vars::Array, simCode::SimulationCode.SimCode)::Array{Expr}
  local startExprs::Array{Expr} = []
  local residuals = simCode.residualEquations
  local ht::Dict = simCode.stringToSimVarHT
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
#            push!(startExprs, :($(Symbol("$varName")) => $(exp))) what is this?
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
                  $(Symbol("$varName")) => $(expToJuliaExpMTK(start, simCode))
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

"""
  Creates parameters on a MTK parameters compatible format.
"""
function createParameterEquationsMTK(parameters::Array, simCode::SimulationCode.SimCode)::Array{Expr}
  local parameterEquations::Array = []
  local hT = simCode.stringToSimVarHT
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
            $(Symbol(simVar.name)) => float($((expToJuliaExpMTK(bindExp, simCode))))
          end
          )
  end
  return parameterEquations
end


"""
  Creates a parameter array.
  A parameter array is an array containing the values of the parameters sorted by index.
  The index here is the index assigned by the code generator earlier in the lowering
  of the hybrid DAE.
"""
function createParameterArray(parameters::Vector{T}, simCode::SIM_T) where {T, SIM_T}
  local paramArray = []
  local hT = simCode.stringToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    local simVarType::SimulationCode.SimVarType = simVar.varKind
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => exp
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
    #= Evaluate the parameters. If it is a variable, and can't be evaluated look it up in the parameter dictonary. =#
    local parValue
    try
      #= The boundvalue is known =#
      val = eval(expToJuliaExpMTK(bindExp, simCode))
      parValue = :($(val))
    catch #=If the bound value is a more complex expression. =#
      parValue = :(0) #pars[Num($(param))]) (More complex parameters are yet to be used in the benchmark..)
    end
    push!(paramArray, parValue)
  end
  return paramArray
end

"""
  Generates code for DAE cast expressions for MTK code.
"""
function generateCastExpressionMTK(ty, exp, simCode, varPrefix)
  return @match ty, exp begin
    (DAE.T_REAL(__), DAE.ICONST(__)) => begin
      quote
        float($(expToJuliaExpMTK(exp, simCode, varPrefix=varPrefix)))
      end
    end
    (DAE.T_REAL(__), DAE.CREF(cref)) where typeof(cref.identType) === DAE.T_INTEGER => begin
      quote
        float($(expToJuliaExpMTK(exp, simCode, varPrefix=varPrefix)))
      end
    end
    #= Conversion to a float other alternatives, =#
    (DAE.T_REAL(__), _) => begin
      quote
        float($(expToJuliaExpMTK(exp, simCode, varPrefix=varPrefix)))
      end
    end
    _ => throw("Cast $ty: for exp: $exp not yet supported in codegen!")
  end
end
