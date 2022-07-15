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
=#

#=
  Author: John Tinnerholm
  TODO: Remember the state derivative scheme. What did I mean with that?
  TODO: Make duplicate code better...
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

"""
  The entry point of MTK code generation.
  Either calls ODE_MODE_MTK_PROGRAM_GENERATION
  or do code generation for a model with structural submodels.
"""
function ODE_MODE_MTK(simCode::SimulationCode.SIM_CODE)
    #=If our model name is separated by . replace it with __ =#
  local MODEL_NAME = replace(simCode.name, "." => "__")
  if (isempty(simCode.structuralTransitions) && length(simCode.subModels) < 1)
    #= Generate using the standard name =#
    return ODE_MODE_MTK_PROGRAM_GENERATION(simCode, simCode.name)
  end
  #= Handle structural submodels =#
  #= TODO: Extract the active model =#
  local activeModelSimCode = getActiveModel(simCode)
  local activeModelName = simCode.activeModel
  local structuralModes = Expr[]
  for mode in simCode.subModels
    push!(structuralModes, ODE_MODE_MTK_MODEL_GENERATION(mode, mode.name))
  end
  if isempty(simCode.subModels)
    local modelName = MODEL_NAME * "DEFAULT"
    defaultModel = ODE_MODE_MTK_MODEL_GENERATION(simCode, modelName)
    activeModelName = modelName
    push!(structuralModes, defaultModel)
  end
  local structuralCallbacks = createStructuralCallbacks(simCode, simCode.structuralTransitions)
  local structuralAssignments = createStructuralAssignments(simCode, simCode.structuralTransitions)
  #= Initialize array where the common variables are stored. That is variables all modes have =#
  local commonVariables = createCommonVariables(simCode.sharedVariables)
  #= END =#
  code = quote
    import DAE
    import DataStructures.OrderedCollections
    import SCode
    import OMBackend
    import OMBackend.CodeGeneration
    using ModelingToolkit
    using DifferentialEquations
    $(structuralModes...)
    $(structuralCallbacks...)
    #=
      This function can be used to fetch the top level callbacks that is the collected callbacks of the model.
      Each callback is coupled to each when with recompilation expression.
    =#
    function $(Symbol(MODEL_NAME * "Model"))(tspan = (0.0, 1.0))
      #=  Assign the initial model  =#
      (subModel, initialValues, reducedSystem, _, pars, vars1) = $(Symbol(activeModelName  * "Model"))(tspan)
      #= Assign the structural callbacks =#
      $(structuralAssignments)
      $(commonVariables)
      #= END =#
      callbackConditions = $(if isempty(structuralCallbacks)
                               :(CallbackSet())
                             elseif length(structuralCallbacks) == 1
                               :(CallbackSet(first(callbackSet)))
                             else
                               :(CallbackSet(Tuple(callbackSet), ())) #TODO: I think this only applies to continious callbacks.
                             end)
      #= Create the composite model =#
      compositeProblem = ModelingToolkit.ODEProblem(
        reducedSystem,
        initialValues,
        tspan,
        pars,
        callback = callbackConditions,
      )
      result = $(if simCode.metaModel == nothing
                 :(OMBackend.Runtime.OM_ProblemStructural($(activeModelName),
                                                          compositeProblem,
                                                          structuralCallbacks,
                                                          commonVariables))
                 else
                 :(OMBackend.Runtime.OM_ProblemRecompilation($(activeModelName),
                                                             compositeProblem,
                                                             structuralCallbacks,
                                                             callbackConditions))
                 end)
      return result
    end
    function $(Symbol("$(MODEL_NAME)Simulate"))(tspan = (0.0, 1.0); solver=Rodas5())
      $(Symbol("$(MODEL_NAME)Model_problem")) = $(Symbol("$(MODEL_NAME)Model"))(tspan)
      OMBackend.Runtime.solve($(Symbol("$(MODEL_NAME)Model_problem")), tspan, solver)
    end
  end
  return (MODEL_NAME, code)
end

"""
  Generates a MTK program with a model
"""
function ODE_MODE_MTK_PROGRAM_GENERATION(simCode::SimulationCode.SIM_CODE, modelName)
  local MODEL_NAME = replace(modelName, "." => "__")
  local model = ODE_MODE_MTK_MODEL_GENERATION(simCode, modelName)
  program = quote
    using ModelingToolkit
    using DifferentialEquations
    import ModelingToolkit.IfElse
    indexof(sym,syms) = findfirst(isequal(sym),syms)
    $(model)
    function $(Symbol("$(MODEL_NAME)Simulate"))(tspan = (0.0, 1.0); solver=Rodas5())
      ($(Symbol("$(MODEL_NAME)Model_problem")), ivs, $(Symbol("$(MODEL_NAME)Model_ReducedSystem")), tspan, pars, vars) = $(Symbol("$(MODEL_NAME)Model"))(tspan)
      solve($(Symbol("$(MODEL_NAME)Model_problem")), solver)
    end
  end
  #= MODEL_NAME is preprocessed with . replaced with _=#
  return MODEL_NAME, program
end

"""
  Generates a MTK model
"""
function ODE_MODE_MTK_MODEL_GENERATION(simCode::SimulationCode.SIM_CODE, modelName)
  @debug "Runnning: ODE_MODE_MTK_MODEL"
  RESET_CALLBACKS()
  local stringToSimVarHT = simCode.stringToSimVarHT
  local equations::Vector = BDAE.RESIDUAL_EQUATION[]
  local exp::DAE.Exp
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
      SimulationCode.DISCRETE(__) => push!(discreteVariables, varName)
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
  local EQUATIONS = createResidualEquationsMTK(stateVariables,
                                               algebraicVariables,
                                               simCode.residualEquations,
                                               simCode::SimulationCode.SIM_CODE)
  local parVariablesSym = [Symbol(p) for p in parameters]
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquationsMTK(stateVariables, algebraicVariables, simCode)
  local DISCRETE_START_VALUES = vcat(generateInitialEquations(simCode.initialEquations, simCode; parameterAssignment = true),
                                   getStartConditionsMTK(discreteVariables, simCode))

  local PARAMETER_EQUATIONS = createParameterEquationsMTK(parameters, simCode)
  local PARAMETER_ASSIGNMENTS = createParameterAssignmentsMTK(parameters, simCode)
  local PARAMETER_RAW_ARRAY = createParameterArray(parameters, PARAMETER_ASSIGNMENTS, simCode)
  #= Create callback equations.
    For MTK we disable the saving function for now.
  =#
  local CALL_BACK_EQUATIONS = createCallbackCode(modelName, simCode; generateSaveFunction = false)
  local IF_EQUATION_COMPONENTS::Vector{Tuple{Expr, Expr, Vector{Expr}, Vector{Expr}, Vector{Tuple}}} = createIfEquations(stateVariables, algebraicVariables, simCode)
  #= Symbolic names =#
  local stateVariablesSym = [:($(Symbol(v))(t)) for v in stateVariables]
  local algebraicVariablesSym = [:($(Symbol(v))(t)) for v in algebraicVariables]
  local discreteVariablesSym = [:($(Symbol(v))(t)) for v in discreteVariables]
  #=Preprocess the component of the if equations =#
  local DISCRETE_DUMMY_EQUATIONS = [:(der($(Symbol(dv))) ~ 0) for dv in discreteVariables]
  local IF_EQUATION_EVENTS = [component[1] for component in IF_EQUATION_COMPONENTS]
  local IF_EQUATION_EVENT_DECLARATION = if isempty(IF_EQUATION_EVENTS)
    :(events = [])
  else
    :(events = $(IF_EQUATION_EVENTS...))
  end
  local CONDITIONAL_EQUATIONS = [component[2] for component in IF_EQUATION_COMPONENTS]
  local ifConditionNameAndIV = collect(Iterators.flatten([component[5] for component in IF_EQUATION_COMPONENTS]))
  local ifConditionalVariables = collect(Iterators.flatten([component[4] for component in IF_EQUATION_COMPONENTS]))
  local ZERO_DYNAMICS_COND_EQUATIONS = collect(Iterators.flatten([component[3] for component in IF_EQUATION_COMPONENTS]))
  #= Expand the start conditions with initial equations for the zero dynamic equations for the conditional equations =#  
  local ifConditionalStartEquations = [:($(Symbol(first(v))) => $(last(v))) for v in ifConditionNameAndIV]
  START_CONDTIONS_EQUATIONS = vcat(ifConditionalStartEquations, DISCRETE_START_VALUES, START_CONDTIONS_EQUATIONS)
  #=
    Merge the ifConditional components into the rest of the system and merge the state conditionals with the sates
  =#
  stateVariablesSym = vcat(ifConditionalVariables, discreteVariablesSym, stateVariablesSym)
  EQUATIONS = vcat(EQUATIONS, DISCRETE_DUMMY_EQUATIONS, ZERO_DYNAMICS_COND_EQUATIONS, CONDITIONAL_EQUATIONS)
  EQUATIONS = rewriteEquations(EQUATIONS,
                               t,
                               vcat(stateVariablesSym, algebraicVariablesSym),
                               parVariablesSym)
  #= Reset the callback counter=#
  RESET_CALLBACKS()
  #=
    Formulate the problem as a DAE Problem.
    For this variant we keep it on its own line
    https://github.com/SciML/ModelingToolkit.jl/issues/998
  =#
  #=If our model name is separated by . replace it with __ =#
  local MODEL_NAME = replace(modelName, "." => "__")
  model = quote
    $(CALL_BACK_EQUATIONS)
    function $(Symbol(MODEL_NAME * "Model"))(tspan = (0.0, 1.0))
      ModelingToolkit.@variables t
      D = ModelingToolkit.Differential(t)
      parameters = ModelingToolkit.@parameters begin
        ($(parVariablesSym...))
      end   
      #=
        Only variables that are present in the equation system later should be a part of the variables in the MTK system.
        This means that certain algebraic variables should not be listed among the variables (These are the discrete variables).
      =#
      $(decomposeVariables(stateVariablesSym, algebraicVariablesSym))
      componentVars = []
      for constructor in variableConstructors
        res = eval(ModelingToolkit.Symbolics._parse_vars("CustomCall", Real, constructor()))
        #= t is no longer needed here =#
        push!(componentVars, res[2:end])
      end
      vars =  collect(Iterators.flatten(componentVars))
      #= Initial values for the continious system. =#
      pars = Dict($(PARAMETER_EQUATIONS...))
      startEquationComponents = []
      $(decomposeStartEquations(START_CONDTIONS_EQUATIONS))
      for constructor in startEquationConstructors
        push!(startEquationComponents, constructor())
      end
      initialValues = collect(Iterators.flatten(startEquationComponents))
      #= End construction of initial values =#
      equationComponents = []
      $(decomposeEquations(EQUATIONS, PARAMETER_ASSIGNMENTS))
      for constructor in equationConstructors
        push!(equationComponents, constructor())
      end
      eqs = collect(Iterators.flatten(equationComponents))
      $(IF_EQUATION_EVENT_DECLARATION)
      nonLinearSystem = $(odeSystemWithEvents(!(isempty(ifConditionalStartEquations)), modelName))
      firstOrderSystem = nonLinearSystem
      $(performStructuralSimplify(performIndexReduction)) #TODO. Causes issues for some systems... my own simplify?
      #=
        These arrays are introduced to handle the bolted on event handling using callbacks.
        The callback handling for MTK is subject of change should hybrid system be implemented for MTK.
      =#
      local event_p = [$(PARAMETER_RAW_ARRAY...)]
      local discreteVars = collect(values(ModelingToolkit.OrderedDict($(DISCRETE_START_VALUES...))))
      #= Merge the discrete and event parameters =#
      event_p = vcat(event_p, discreteVars)
      local aux = Vector{Any}(undef, 3)
      #= TODO init them with the initial values of them =#
      aux[1] = event_p
      aux[2] = Float64[]
      aux[3] = reducedSystem
      #=
       This array is used to help indexing state/discrete variables.
       It maps the index given by OMBackend to the actual index in the states. 
       A[i] = <true_index>
      =#
      problem = ModelingToolkit.ODEProblem(reducedSystem,
                                           initialValues,
                                           tspan,
                                           pars,
                                           callback=$(Symbol("$(MODEL_NAME)CallbackSet"))(aux))
      return (problem, initialValues, reducedSystem, tspan, pars, vars)
    end
  end
  return model
end

"""
   Creates equations from the residual equations in unsorted order
"""
function createResidualEquationsMTK(stateVariables::Vector, algebraicVariables::Vector, equations::Vector{BDAE.RESIDUAL_EQUATION}, simCode::SimulationCode.SIM_CODE)::Vector{Expr}
  if isempty(equations)
    return Expr[]
  end
  local eqs::Vector{Expr} = Expr[]
  for eq in equations
    local eqDAEExp = eq.exp
    local eqExp = :(0 ~ $(expToJuliaExpMTK(eqDAEExp, simCode; derSymbol=false)))
    push!(eqs, eqExp)
  end
    return eqs
end

include("mtkExternals.jl")

"""
  Generates the initial value for the equations
  TODO: Currently unable to generate start condition in order
"""
function createStartConditionsEquationsMTK(states::Vector,
                                        algebraics::Vector,
                                        simCode::SimulationCode.SimCode)::Vector{Expr}
  local algInit = getStartConditionsMTK(algebraics, simCode)
  local stateInit = getStartConditionsMTK(states, simCode)
  local initialEquations = simCode.initialEquations
  local ieqInit = generateInitialEquations(initialEquations, simCode)
  #=
    Start with the start conditions above.
    Generate the equations in order afterwards
  =#
  #= Place the initial equations last =#
  return vcat(algInit, stateInit, ieqInit)
end

"""
  Generates initial equations.
  Currently unsorted unless they are sorted before being passed to the simulation code phase.
"""
function generateInitialEquations(initialEqs, simCode::SimulationCode.SimCode; parameterAssignment = true)::Vector{Expr}
  local initialEqsExps = Expr[]
  for ieq in initialEqs
    #= LHS will typically be a variable. Don't have to be though.. =#
    lhs = expToJuliaExpMTK(ieq.lhs, simCode)
    rhs = @match ieq.rhs begin
      DAE.CREF(__) => begin
        #= Evaluate the right hand side at this point =#
        local crefAsStr = string(ieq.rhs)
        local simCodeVar = last(simCode.stringToSimVarHT[crefAsStr])
        local res = if SimulationCode.isStateOrAlgebraic(simCodeVar)
          #= Otherwise get the start attribute =#
          expToJuliaExpMTK(ieq.rhs, simCode)
        else
          evalSimCodeParameter(simCodeVar)
        end
      end
      #= For more complicated expressions, we do local constant folding. =#
      _ => begin
        res = evalDAE_Expression(ieq.rhs, simCode)
        res
      end
    end
    if parameterAssignment
      push!(initialEqsExps,
            quote
              $lhs => $rhs
            end)
    else
      push!(initialEqsExps,
            quote
              $lhs = $rhs
            end)
    end
  end
  return initialEqsExps
end

"""
  Given a vector of variables and the simulation code
  extracts the start attributes to generate initial conditions.
"""
function getStartConditionsMTK(vars::Vector, simCode::SimulationCode.SimCode)::Vector{Expr}
  local startExprs::Vector{Expr} = Expr[]
  local residuals = simCode.residualEquations
  local ht::Dict = simCode.stringToSimVarHT
  local warnings::String = ""
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
            #=
              We have a variable of sorts.
              We look evaluate the value in the pars list.
            =#
            push!(startExprs,
                  quote
                  $(Symbol("$varName")) => pars[$(Symbol(string(start)))]
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
        NONE() => begin
          if ! (simVarType isa SimulationCode.STATE || simVarType isa SimulationCode.DISCRETE)
            warnings *= "\n Assumed starting value of 0.0 for variable: " * varName * "\n"
            push!(startExprs, :($(Symbol(varName)) => 0.0))
          end
        end
      end
    end
  end
  if !isempty(warnings)
    @warn warnings
  end
  return startExprs
end

"""
  Creates the components of the If-Equations
"""
function createIfEquations(stateVariables, algebraicVariables, simCode)
  local ifEquations::Vector{Tuple{Expr, Expr, Vector{Expr}, Vector{Expr}, Vector{Tuple}}} = Tuple{Expr, Expr, Vector{Expr}, Vector{Expr}, Vector{Tuple}}[]
  local i::Int = 0
  for ifEq in simCode.ifEquations
    i += 1
    push!(ifEquations, createIfEquation(stateVariables, algebraicVariables, ifEq, i, simCode))
  end
  return ifEquations
end

"""
This function creates symbolic if equations for use in MTK.
The function returns a tuple, where the first part of the tuple represent the conditions and the affect of the if-equation on the form:
  continuous_events = [
    <Condtion> => <affect> 
    <Condtion> => <affect>
    ....
  ]
Each condition generates one variable with zero dynamics the variable being true or not depending on the branch.
  Example:
  if <condtion> then
    <equations>
  elseif then
    <equations>
  else
    <equations>
  end if;
Would result in:
continuous_events = [
    <condition> => [ifCond1 ~ true, ifCond2 ~ false]
    <condition> => [ifCond1 ~ false, ifCond2 ~ true]
]
An if equation with a single condition would only generate one condtion:
continuous_events = [
    <condition> => [ifCond1 ~ true]
]

The second value in the return tuble represent the if-equations itself:
<lhs> = IfElse.ifelse(<condition>, <value>, IfElse.ifelse(<condition>, <value>, <value>))
  lhs can be one or several variables. (TODO, fix the case for several variables in this kind of branch)

The third part of the tuple contains a set of zero dynamic equations (One for each if equation condition variable)
See the following issue: https://github.com/SciML/ModelingToolkit.jl/issues/1523

The forth part of the tuple contains a vector of symbolic variables.
One for each conditional variable created.
"""
function createIfEquation(stateVariables, algebraicVariables, ifEq::SimulationCode.IF_EQUATION, identifier, simCode)
  local result::Tuple{Expr, Expr, Vector{Expr}, Vector{Expr}, Vector{Tuple}}
  
  generateAffects(n, nConditions, condRes) = begin
    #=
      The cond res is what we are going to evaluate it to.
    =#
    local exprs = Expr[]
    local e
    for i in 1:nConditions
      e = if i == n 
        :($(Symbol(("ifCond$(i)"))) ~ $(condRes))
      else
        :($(Symbol(("ifCond$(i)"))) ~ $(!condRes))
      end
      push!(exprs, e)
    end
    return exprs
  end

  #= TODO. We currently assume residuals that we have made causal and that the original equations are written in a certain form.=#
  deCausalize(eq) = begin
    @match eq.exp begin 
      DAE.BINARY(exp1, _, exp2) => begin
        (:($(expToJuliaExpMTK(exp2, simCode))), :($(expToJuliaExpMTK(exp1, simCode))))
      end
      _ => begin
        throw("Unsupported equation:" * string(eq))
      end
    end
  end

  generateIfExpressions(branches, target) = begin
    if branches[target].targets == -1
      return :($(first(deCausalize(branches[target].residualEquations[1]))))
    end
    #= Otherwise generate code for the other part =#
    local branch = branches[target]
    local cond = :( $(Symbol("ifCond$(target)")) == true )
    #= Assume single equations for now =#
    :(ModelingToolkit.IfElse.ifelse($(cond),
                                    $(first(deCausalize(branch.residualEquations[1]))),
                                    $(generateIfExpressions(branches, branches[target].targets))))
  end
  
  local i::Int = 0
  local nBranches::Int = length(ifEq.branches)
  local conditions = Expr[]
  local ivConditions = Bool[]
  for branch in ifEq.branches
    i += 1
    @match branch begin
      SimulationCode.BRANCH(condition, residuals, -1#=Else=#, targets, _, _, _, _, _) => begin        
      end
      SimulationCode.BRANCH(condition, residuals, _, targets, _, _, _, _, _) => begin
        local mtkCond = transformToMTKConditionEquation(branch.condition, simCode)
        local ivCond = evalInitialCondition(mtkCond)
        local branchesWithConds::Int = nBranches - 1
        local affects::Vector{Expr} = generateAffects(i, branchesWithConds, ivCond)
        local cond = :([$(mtkCond)] => [$(affects...)])
        push!(conditions, cond)
        push!(ivConditions, ivCond)
      end
    end
  end
  local target = 1
  local resEq = ifEq.branches[target].residualEquations[1]
  ifExpressions = :($(last(deCausalize(resEq))) ~ $(generateIfExpressions(ifEq.branches, target)))
  #= Generate zero dynamic equations for the conditions =#
  conditionEquations = Expr[]
  conditionVariables = Expr[]
  conditionVariableNames = Tuple{String, Bool}[]
  for i in 1:length(conditions)
    push!(conditionEquations, :(der($(Symbol("ifCond$(i)"))) ~ 0))
    push!(conditionVariables, :($(Symbol("ifCond$(i)"))(t)))
    push!(conditionVariableNames, ("ifCond$(i)", !(ivConditions[i])))
  end
  conditionExpr = :([$(conditions...)])
  result = (conditionExpr, ifExpressions, conditionEquations, conditionVariables, conditionVariableNames)
  return result
end

"""
  Creates parameters on a MTK parameters compatible format.
"""
function createParameterEquationsMTK(parameters::Vector, simCode::SimulationCode.SimCode)::Vector{Expr}
  local parameterEquations::Vector = []
  local hT = simCode.stringToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    local simVarType::SimulationCode.SimVarType = simVar.varKind
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => begin
        exp
      end
      _ => begin
        throw(ErrorException("Unknown SimulationCode.SimVarType for parameter: " * string(simVarType)))
      end
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
  Creates parameters assignments on a MTK parameters compatible format.
"""
function createParameterAssignmentsMTK(parameters::Vector, simCode::SimulationCode.SimCode)::Vector{Expr}
  local parameterEquations::Vector = []
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
            $(Symbol(simVar.name)) = float($((expToJuliaExpMTK(bindExp, simCode))))
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
function createParameterArray(parameters::Vector{T1}, parameterAssignments::Vector{T2}, simCode::SIM_T) where {T1, T2, SIM_T}
  local paramArray = []
  local hT = simCode.stringToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    @info "Index when creating parameter array for $(string(simVar)) : " index
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
      if val isa Float64
        parValue = :($(val))
      else
        parValue = :($(Symbol(param))) #:(0.0)
      end
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

"""
 This function decomposes the continuous variables.
 This means that if the set of variables are greater than 50 a new inner function is generated as to not
 impact the JIT of the system to much.
"""
function decomposeVariables(stateVariables, algebraicVariables)
  local nStateVars = length(stateVariables)
  local nAlgVars = length(algebraicVariables)
  if  1 < nStateVars < 50 &&  1 < nAlgVars < 50
    expr = quote
      function generateStateVariables()
        $(Tuple([:t, stateVariables...]))
      end
      function generateAlgebraicVariables()
        $(Tuple([:t, algebraicVariables...]))
      end
      variableConstructors = Function[generateStateVariables, generateAlgebraicVariables]
    end
  elseif (1 < nStateVars < 50) && nAlgVars == 0
    expr = quote
      function generateStateVariables()
        $(Tuple([:t, stateVariables...]))
      end
      variableConstructors = Function[generateStateVariables]
    end
  else
    #= Split the array in chunks of 50  for the state and algebraic variables=#
    local stateVectors = collect(Iterators.partition(stateVariables, 50))
    local algVectors = collect(Iterators.partition(algebraicVariables, 50))
    #= For each vector in stateVectors create a constructor for the variables =#
    local i = 1::Int
    local exprs = Expr[]
    constructors = quote
      variableConstructors = Function[]
    end
    push!(exprs, constructors)
    for stateVector in stateVectors
      stateConstructorExpr = quote
        function $(Symbol("generateStateVariables" * string(i)))()
          $(Tuple([:t, stateVector...]))
        end
        push!(variableConstructors, $(Symbol("generateStateVariables" * string(i))))
      end
      push!(exprs, stateConstructorExpr)
      i += 1
    end
    local i = 1
    #= decompose the algebraic variables if needed =#
    for algVector in algVectors
      algConstructorExpr = quote
        function $(Symbol("generateAlgebraicVariables" * string(i)))()
          $(Tuple([:t, algVector...]))
        end
        push!(variableConstructors, $(Symbol("generateAlgebraicVariables" * string(i))))
      end
      push!(exprs, algConstructorExpr)
      i += 1
    end
    #= Generate the composite expression =#
    expr = quote
      $(exprs...)
    end
  end
  return expr
end

"""
  Similar to decompose variables, however, decomposes the set of equations instead.
  This function does so by dividing the total number of equations into separate blocks with 50 equations in each block.
  This function is suppose to be called after decompose variables.
  the parameter assignments contains the set of values assigned to parameters at the beginning of the simulation, furthermore the initialEquation contains the initial equations.
  The reason for including them here is to be able to perform additional simplifications.
"""
function decomposeEquations(equations, parameterAssignments)
  local nStateVars = length(equations)
  local equationVectors = collect(Iterators.partition(equations, 50))
  local exprs = Expr[]
  local constructors = quote
    equationConstructors = Function[]
  end
  push!(exprs, constructors)
  local i = 0
  for equationVector in equationVectors
    eqv = [rewriteEq(i) for i in equationVector]
    equationConstructorExpr = quote
      $(parameterAssignments...) #TODO: It does not seem to work to use parameters as constants, something goes wrong in the substitution.
      function $(Symbol("generateEquations" * string(i)))()
        [$(eqv...)]
      end
      push!(equationConstructors, $(Symbol("generateEquations" * string(i))))
    end
    push!(exprs, equationConstructorExpr)
    i += 1
  end
  expr = quote
    $(exprs...)
  end
  return expr
end

"""
  Similar to decomposeEquations but for start equations.
"""
function decomposeStartEquations(equations)
  local nStateVars = length(equations)
  local equationVectors = collect(Iterators.partition(equations, 50))
  local exprs = Expr[]
  local constructors = quote
    startEquationConstructors = Function[]
  end
  push!(exprs, constructors)
  local i = 0
  for equationVector in equationVectors
    equationConstructorExpr = quote
      function $(Symbol("generateStartEquations" * string(i)))()
        [$(equationVector...)]
      end
      push!(startEquationConstructors, $(Symbol("generateStartEquations" * string(i))))
    end
    push!(exprs, equationConstructorExpr)
    i += 1
  end
  expr = quote
    $(exprs...)
  end
  return expr
end
