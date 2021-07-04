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

"""
  Collect variables from array of BDAE.Var:
  Save the name and it's kind of each variable.
  Index will be set to NONE.
"""
function collectVariables(allBackendVars::Array{BDAE.Var})
  local numberOfVars::Integer = length(allBackendVars)
  local simVars::Array = Array{SimulationCode.SimVar}(undef, numberOfVars)
  for (i, backendVar) in enumerate(allBackendVars)
    local simVarName::String = bDAEIdentToSimCodeVarName(backendVar)
    local simVarKind::SimulationCode.SimVarType = bDAEVarKindToSimCodeVarKind(backendVar)
    simVars[i] = SimulationCode.SIMVAR(simVarName, NONE(), simVarKind, backendVar.values)
  end
  return simVars
end

function bDAEVarKindToSimCodeVarKind(backendVar::BDAE.Var)::SimulationCode.SimVarType
  varKind = @match backendVar.varKind begin
    BDAE.STATE(__) => SimulationCode.STATE()
    BDAE.PARAM(__) || BDAE.CONST(__) => SimulationCode.PARAMETER(backendVar.bindExp)
    BDAE.VARIABLE(__) => SimulationCode.ALG_VARIABLE()
    _ => @error("Kind $(typeof(backendVar.varKind)) of backend variable not handled.")
  end
end

function bDAEIdentToSimCodeVarName(backendVar::BDAE.Var)
  local varName::DAE.ComponentRef = backendVar.varName
  @match varName begin
    DAE.CREF_IDENT(__) => string(varName)
    DAE.CREF_QUAL(__) => string(varName)
    _ => @error("Type $(typeof(varName)) not handled.")
  end
end

"""
  Transform BDAE-Structure to SimulationCode.SIM_CODE
  We only handle one equation system for now
"""
function transformToSimCode(backendDAE::BDAE.BACKEND_DAE)::SimulationCode.SIM_CODE
  local equationSystems::Array = backendDAE.eqs
  local allOrderedVars::Array{BDAE.Var} = [v for es in equationSystems for v in es.orderedVars.varArr]
  local allSharedVars::Array{BDAE.Var} = getSharedVariablesLocalsAndGlobals(backendDAE.shared)
  local allBackendVars = vcat(allOrderedVars, allSharedVars)
  local simVars::Array{SimulationCode.SIMVAR} = allocateAndCollectSimulationVariables(allBackendVars)
  # Assign indices and put all variable into an hash table
  local crefToSimVarHT = createIndices(simVars)
  local equations = [eq for es in equationSystems for eq in es.orderedEqs]
  #= Split equations into three parts. Residuals whenEquations and If-equations =#
  (resEqs::Vector{BDAE.RESIDUAL_EQUATION}, whenEqs::Vector{BDAE.WHEN_EQUATION}, ifEqs::Vector{BDAE.IF_EQUATION}) =
    allocateAndCollectSimulationEquations(equations)
  #= Sorting/Matching for the set of residual equations (This is used for the start condtions) =#
  local eqVariableMapping = createEquationVariableBidirectionGraph(resEqs, allBackendVars, crefToSimVarHT)
  local numberOfVariablesInMapping = length(eqVariableMapping.keys)
  (isSingular, matchOrder, digraph, stronglyConnectedComponents) =
    matchAndCheckStronglyConnectedComponents(eqVariableMapping, numberOfVariablesInMapping)
  #= 
    The set of if equations needs to be handled in a separate way. 
    Each branch might contain a separate section of variables etc that needs to be sorted and processed. 
    !!It is assumed that the frontend have checked each branch for balance at this point!!
  =#
  simCodeBranches::Vector{SIM_CODE_IF_EQUATION} = constructSimCodeIFEquations(ifEqs, resEqs, allBackendVars, crefToSimVarHT)
  #= Construct SIM_CODE =#
  SimulationCode.SIM_CODE(backendDAE.name,
                          crefToSimVarHT,
                          resEqs,
                          #=TODO fix initial equations here =#
                          BDAE.RESIDUAL_EQUATION[],
                          whenEqs,
                          simCodeBranches,
                          isSingular,
                          matchOrder,
                          digraph,
                          stronglyConnectedComponents)
end

"""
  Given a set of BDAE IF_EQUATIONS. 
  Constructs the set of simulation code if equations.
  Each if equation can be seen as a small basic block graph.
Currently we merge the other residual equation with the equations of one branch.
"""
function constructSimCodeIFEquations(ifEquations::Vector{BDAE.IF_EQUATION},
                                     resEqs::Vector{BDAE.RESIDUAL_EQUATION},
                                     allBackendVars, crefToSimVarHT)::Vector{SIM_CODE_IF_EQUATION}
  simCodeIfEquations = SIM_CODE_IF_EQUATION[]
  for ifEquation in ifEquations
    #= Enumerate the branches of the if equation =#
    #= Handle if and else if branches=#
    local conditions = ifEquation.conditions
    local condition
    local equations
    local target
    local isSingular
    local matchOrder
    local equationGraph
    local sccs
    local branches::Vector{SIM_CODE_BRANCH} = SIM_CODE_BRANCH[]
    for conditionIdx in 1:length(conditions)
      condition = listGet(conditions, conditionIdx)
      #= Merge the given residual equations with the equations of this particular branch. =#
      equations = listAppend(listGet(ifEquation.eqnstrue, conditionIdx), arrayList(resEqs))
      target = conditionIdx + 1
      local eqVariableMapping = createEquationVariableBidirectionGraph(equations, allBackendVars, crefToSimVarHT)
      #= Match and get the strongly connected components =#
      local numberOfVariablesInMapping = length(eqVariableMapping.keys)
      (isSingular, matchOrder, digraph, stronglyConnectedComponents) =
        matchAndCheckStronglyConnectedComponents(eqVariableMapping, numberOfVariablesInMapping)
      #= Add the branch to the collection. =#
      branch = SIM_CODE_BRANCH(condition, equations, identifier, target, isSingular, matchOrder, digraph, stronglyConnectedComponents)
      push!(branches, branch)
    end
    push!(simCodeIfEquations, branches)
  end
  return simCodeIfEquations
end

"""
  This function does matching, it also checks for strongly connected components
"""
function matchAndCheckStronglyConnectedComponents(eqVariableMapping, numberOfVariablesInMapping)::Tuple
  (isSingular::Bool, matchOrder::Array) = GraphAlgorithms.matching(eqVariableMapping, numberOfVariablesInMapping)
  local digraph::MetaGraphs.MetaDiGraph = GraphAlgorithms.merge(matchOrder, eqVariableMapping)
  if OMBackend.PLOT_EQUATION_GRAPH
    local labels = makeLabels(digraph, matchOrder, crefToSimVarHT)
    GraphAlgorithms.plotEquationGraph("./digraphOutput$(backendDAE.name).pdf", digraph, labels)
  end
  stronglyConnectedComponents::Array = GraphAlgorithms.topological_sort(digraph)
  return (isSingular, matchOrder, digraph, stronglyConnectedComponents)
end

"""
John:
  Splits a given set of equations into different types
"""
function allocateAndCollectSimulationEquations(equations::T)::Tuple where {T}
  isRe(eq) = typeof(eq) == BDAE.RESIDUAL_EQUATION
  isWhen(eq) = typeof(eq) == BDAE.WHEN_EQUATION
  isIf(eq) = typeof(eq) == BDAE.IF_EQUATION
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
