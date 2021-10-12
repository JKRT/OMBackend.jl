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

function bDAEVarKindToSimCodeVarKind(backendVar::BDAE.Var)::SimulationCode.SimVarType
  varKind = @match backendVar.varKind begin
    BDAE.STATE(__) => SimulationCode.STATE()
    BDAE.PARAM(__) || BDAE.CONST(__) => SimulationCode.PARAMETER(backendVar.bindExp)
    BDAE.VARIABLE(__) => SimulationCode.ALG_VARIABLE()
    _ => @error("Kind $(typeof(backendVar.varKind)) of backend variable not handled.")
  end
end

"""
`BDAE_identifierToVarString(backendVar::BDAE.Var)`
Converts a `BDAE.Var` to a simcode string, for use in variable names. 
The . separator is replaced with 3 `_`
"""
function BDAE_identifierToVarString(backendVar::BDAE.Var)
  local varName::DAE.ComponentRef = backendVar.varName
  @match varName begin
    DAE.CREF_IDENT(__) => string(varName)
    DAE.CREF_QUAL(__) => string(varName)
    _ => @error("Type $(typeof(varName)) not handled.")
  end
end

"""
  Transform BDAE-Structure to SimulationCode.SIM_CODE
  We only handle one equation system for now.
  mode specifies what mode we should use for code generation.
"""
function transformToSimCode(backendDAE::BDAE.BACKEND_DAE; mode)::SimulationCode.SIM_CODE
  #= Fetch the different components of the model.=#
  local equationSystems::Array = backendDAE.eqs
  local allOrderedVars::Array{BDAE.Var} = [v for es in equationSystems for v in es.orderedVars.varArr]
  local allSharedVars::Array{BDAE.Var} = getSharedVariablesLocalsAndGlobals(backendDAE.shared)
  local allBackendVars = vcat(allOrderedVars, allSharedVars)
  local simVars::Array{SimulationCode.SIMVAR} = allocateAndCollectSimulationVariables(allBackendVars)
  # Assign indices and put all variable into an hash table
  local stringToSimVarHT = createIndices(simVars)
  local equations = [eq for es in equationSystems for eq in es.orderedEqs]
  #= Split equations into three parts. Residuals whenEquations and If-equations =#
  (resEqs::Vector{BDAE.RESIDUAL_EQUATION}, whenEqs::Vector{BDAE.WHEN_EQUATION}, ifEqs::Vector{BDAE.IF_EQUATION}) =
    allocateAndCollectSimulationEquations(equations)
  #= Sorting/Matching for the set of residual equations (This is used for the start condtions) =#
  local eqVariableMapping = createEquationVariableBidirectionGraph(resEqs, allBackendVars, stringToSimVarHT)
  local numberOfVariablesInMapping = length(eqVariableMapping.keys)
  (isSingular, matchOrder, digraph, stronglyConnectedComponents) =
    matchAndCheckStronglyConnectedComponents(equations, eqVariableMapping, numberOfVariablesInMapping, stringToSimVarHT; mode = mode)
  #= 
    The set of if equations needs to be handled in a separate way. 
    Each branch might contain a separate section of variables etc that needs to be sorted and processed. 
    !!It is assumed that the frontend have checked each branch for balance at this point!!
  =#
  simCodeIfEquations::Vector{IF_EQUATION} = constructSimCodeIFEquations(ifEqs,
                                                                        resEqs,
                                                                        allBackendVars,
                                                                        stringToSimVarHT)
  #= Construct SIM_CODE =#
  SimulationCode.SIM_CODE(backendDAE.name,
                          stringToSimVarHT,
                          resEqs,
                          #=TODO fix initial equations here =#
                          BDAE.RESIDUAL_EQUATION[],
                          whenEqs,
                          simCodeIfEquations,
                          isSingular,
                          matchOrder,
                          digraph,
                          stronglyConnectedComponents)
end

"""
  Given a set of BDAE IF_EQUATIONS. 
  Constructs the set of simulation-code if-equations.
  Each if equation can be seen as a small basic block graph.
Currently we merge the other residual equation with the equations of one branch.
TODO:
Possible rework the identifier scheme.
Unique identifier, static variable?
"""
function constructSimCodeIFEquations(BDAE_ifEquations::Vector{BDAE.IF_EQUATION},
                                     resEqs::Vector{BDAE.RESIDUAL_EQUATION},
                                     allBackendVars::Vector,
                                     stringToSimVarHT)::Vector{IF_EQUATION}
  simCodeIfEquations::Vector{IF_EQUATION} = IF_EQUATION[]
  for BDAE_ifEquation in BDAE_ifEquations
    #= Enumerate the branches of the if equation =#
    #= Handle if and else if branches=#
    local conditions = BDAE_ifEquation.conditions
    local condition
    local equations
    local target
    local isSingular
    local matchOrder
    local equationGraph
    local sccs
    local branches::Vector{BRANCH} = BRANCH[]
    local lastConditionIdx = 0
    for conditionIdx in 1:length(conditions)
      condition = listGet(conditions, conditionIdx)
      #= Merge the given residual equations with the equations of this particular branch. =#
      equations::Vector{BDAE.RESIDUAL_EQUATION} = listArray(listAppend(listGet(BDAE_ifEquation.eqnstrue, conditionIdx),
                                                                       arrayList(resEqs)))
      target = conditionIdx + 1
      identifier = conditionIdx
      local eqVariableMapping = createEquationVariableBidirectionGraph(equations, allBackendVars, stringToSimVarHT)
      #= Match and get the strongly connected components =#
      local numberOfVariablesInMapping = length(eqVariableMapping.keys)
      (isSingular, matchOrder, digraph, stronglyConnectedComponents) =
        matchAndCheckStronglyConnectedComponents(eqVariableMapping, numberOfVariablesInMapping, stringToSimVarHT)
      #= Add the branch to the collection. =#
      branch = BRANCH(condition,
                      equations,
                      identifier,
                      target,
                      isSingular,
                      matchOrder,
                      digraph,
                      stronglyConnectedComponents,
                      stringToSimVarHT)
      push!(branches, branch)
      lastConditionIdx = conditionIdx
    end
    #=
     The procedure above added the code for the elseif branches. 
      It is also possible that we have an else branch.
      The else branch is located in the false equations.
      The condition for the else branch to be inactive is active 
      if all preceeding branches failed to evaluate to true.
    =#
    #= Check if we have an else if not we are done.=#
    if listEmpty(BDAE_ifEquation.eqnsfalse)
      break
    end
    condition = DAE.SCONST("ELSE_BRANCH")
    equations = listArray(listAppend(BDAE_ifEquation.eqnsfalse, arrayList(resEqs)))
    lastConditionIdx += 1
    target = lastConditionIdx + 1
    identifier = ELSE_BRANCH #= Indicate else =#
    local eqVariableMapping = createEquationVariableBidirectionGraph(equations, allBackendVars, stringToSimVarHT)
    local numberOfVariablesInMapping = length(eqVariableMapping.keys)
    (isSingular, matchOrder, digraph, stronglyConnectedComponents) =
      matchAndCheckStronglyConnectedComponents(eqVariableMapping, numberOfVariablesInMapping, stringToSimVarHT)
    branch = BRANCH(condition,
                    equations,
                    identifier,
                    target,
                    isSingular,
                    matchOrder,
                    digraph,
                    stronglyConnectedComponents,
                    stringToSimVarHT)
    push!(branches, branch)    
    ifEq = IF_EQUATION(branches)
    push!(simCodeIfEquations, ifEq)
  end
  return simCodeIfEquations
end

"""
Author:johti17
  This function does matching, it also checks for strongly connected components.
If the system is singular we try index reduction before proceeding.
""" 
function matchAndCheckStronglyConnectedComponents(equations, eqVariableMapping, numberOfVariablesInMapping, stringToSimVarHT; mode)::Tuple
  (isSingular::Bool, matchOrder::Vector) = GraphAlgorithms.matching(eqVariableMapping, numberOfVariablesInMapping)
  local digraph::MetaGraphs.MetaDiGraph

  #=
    If the equation is singular this means our system is singular.
    Index reduction might resolve the issues with this system.
  =#
  if isSingular && mode == OMBackend.MTK_MODE
    # @error "The system of equations is Singular"
    # throw("Error: Singular system")
    #= TODO: Try index reduction =#
    digraph = GraphAlgorithms.merge(matchOrder, eqVariableMapping)
    stronglyConnectedComponents::Array = GraphAlgorithms.stronglyConnectedComponents(digraph)
    return (isSingular, matchOrder, digraph, stronglyConnectedComponents)
  end

  if isSingular && mode == DAE_MODE
    throw("TODO: index reduction not implemented for DAE-mode")
    #= Do index reduction here. =#
  end
  #= We will not get here... For now:) =#
  return matchAndCheckStronglyConnectedComponents(eqVariableMapping, numberOfVariablesInMapping, stringToSimVarHT)
end


"""
  This function conducts matching and checks for strongly connected components.
  If the system is singular after matching it throws an error.
"""
function matchAndCheckStronglyConnectedComponents(eqVariableMapping, numberOfVariablesInMapping, stringToSimVarHT)::Tuple
  (isSingular::Bool, matchOrder::Vector) = GraphAlgorithms.matching(eqVariableMapping, numberOfVariablesInMapping)
  if isSingular
    @error "The system of equations is Singular"
    throw("Error: Singular system")
  end
  #= The merge algortihm takes a matching and produces a resulting digraph =#
  local digraph::MetaGraphs.MetaDiGraph = GraphAlgorithms.merge(matchOrder, eqVariableMapping)
  if OMBackend.PLOT_EQUATION_GRAPH
    @info "Dumping the equation graph"
    local labels = makeLabels(digraph, matchOrder, stringToSimVarHT)
    GraphAlgorithms.plotEquationGraphPNG("./digraphOutput.png", digraph, labels)
  end
  stronglyConnectedComponents::Array = GraphAlgorithms.stronglyConnectedComponents(digraph)
  return (isSingular, matchOrder, digraph, stronglyConnectedComponents)
end

"""
  Author: johti17:
  Splits a given set of equations into different types
"""
function allocateAndCollectSimulationEquations(equations::T)::Tuple where {T}
  local isRe(eq) = typeof(eq) == BDAE.RESIDUAL_EQUATION
  local isWhen(eq) = typeof(eq) == BDAE.WHEN_EQUATION
  local isIf(eq) = typeof(eq) == BDAE.IF_EQUATION
  (filter(isRe, equations), filter(isWhen, equations), filter(isIf, equations))
end

"""
Returns the shared global and local variable for the shared data in
an equation system. If no such data is present. Return two empty arrays
"""
function getSharedVariablesLocalsAndGlobals(shared::BDAE.Shared)
  @match shared begin
    BDAE.SHARED(__) => vcat(shared.globalKnownVars, shared.localKnownVars)
    _ => []
  end
end

function allocateAndCollectSimulationVariables(bDAEVariables::Array{BDAE.Var})
  collectVariables(bDAEVariables)
end



"""
  Collect variables from array of BDAE.Var:
  Save the name and it's kind of each variable.
  Index will be set to NONE.
"""
function collectVariables(allBackendVars::Array{BDAE.Var})
  local numberOfVars::Integer = length(allBackendVars)
  local simVars::Array = Array{SimulationCode.SimVar}(undef, numberOfVars)
  for (i, backendVar) in enumerate(allBackendVars)
    local simVarName::String = BDAE_identifierToVarString(backendVar)
    local simVarKind::SimulationCode.SimVarType = bDAEVarKindToSimCodeVarKind(backendVar)
    simVars[i] = SimulationCode.SIMVAR(simVarName, NONE(), simVarKind, backendVar.values)
  end
  return simVars
end
