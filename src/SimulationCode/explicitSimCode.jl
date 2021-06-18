import ..OMBackend

"""
  This function transforms the Hybrid-DAE to so called explicit simulation code.
  In explicit simcode the residual equations are sorted vertically and horisontally
 """
function transformToExplicitSimCode(backendDAE::BDAE.BACKEND_DAE)::SimulationCode.EXPLICIT_SIM_CODE
  local equationSystems::Array = backendDAE.eqs
  local allOrderedVars::Array{BDAE.Var} = [v for es in equationSystems for v in es.orderedVars.varArr]
  local allSharedVars::Array{BDAE.Var} = getSharedVariablesLocalsAndGlobals(backendDAE.shared)
  local allBackendVars::Array{BDAE.Var} = vcat(allOrderedVars, allSharedVars)
  local simVars::Array{SimulationCode.SIMVAR} = allocateAndCollectSimulationVariables(allBackendVars)
  local crefToSimVarHT::OrderedDict = createIndices(simVars)
  local equations = [eq for es in equationSystems for eq in es.orderedEqs]
  #= Split equations into three parts. Residuals whenEquations and If-equations =#
  (resEqs,whenEqs,ifEqs) = allocateAndCollectSimulationEquations(equations)
  #=Only assume I have residuals for now=#
  @assert(length(whenEqs) == 0, "IF EQUATION NOT YET SUPPORTED IN EXPLICIT CODE GEN")
  @assert(length(ifEqs) == 0, "WHEN EQUATION NOT YET SUPPORTED IN EXPLICIT CODE GEN")
  local indexToEquation::OrderedDict = createEquationIndicies(resEqs)
  #= Create equation <-> variable mapping =#
  local eqVariableMapping = createEquationVariableBidirectionGraph(equations, allBackendVars, crefToSimVarHT)
  (isSingular::Bool, matchOrder::Array) = GraphAlgorithms.matching(eqVariableMapping, length(eqVariableMapping.keys))
  digraph = GraphAlgorithms.merge(matchOrder, eqVariableMapping)
  if OMBackend.PLOT_EQUATION_GRAPH
    local labels = makeLabels(digraph, matchOrder, crefToSimVarHT)
    GraphAlgorithms.plotEquationGraph("./digraphOutput$(backendDAE.name).pdf", digraph, labels)
  end
  stronglyConnectedComponents::Array = GraphAlgorithms.topological_sort(digraph) 
  #= Reorder the residuals =#
  reOrderedResiduals = []
  reverseTopologicalSort = stronglyConnectedComponents
  loopsExist = LightGraphs.is_cyclic(digraph)
  if ! loopsExist
    for i in reverseTopologicalSort
      push!(reOrderedResiduals, resEqs[i[1]])
    end
  else
    @info "Loop encountered exiting..."
    @error "Unresolved algebraic loop"
    throw()
  end
#= Return explicit simulation code. =#
  return SimulationCode.EXPLICIT_SIM_CODE(backendDAE.name,
                                          crefToSimVarHT,
                                          indexToEquation,
                                          eqVariableMapping,
                                          reOrderedResiduals,
                                          initialEquations,
                                          whenEqs,
                                          ifEqs,
                                          isSingular,
                                          matchOrder,
                                          digraph,
                                          stronglyConnectedComponents)
end

function createEquationIndicies(resEqs)
  local index = 1;
  local ht::OrderedDict = OrderedDict()
  for e in resEqs
    ht[index] = e
    index += 1
  end
  return ht
end

#TODO Code for detecting algebraic loops to be removed.
  # #= Check algebraic loops =#
  # loopsExist = false
  # for i in reverseTopologicalSort
  #   if length(i) > 1
  #     loopsExist = true
  #     break
  #   end
  # end



function dumpInfoOfSort(matchOrder, reverseTopologicalSort, crefToSimVarHT)
  #= Tearing could be added here =#
  @info "Strongly connected components!" reverseTopologicalSort
  @info "Equation solving order:"
  ht = makeIndexVarNameDict(matchOrder, crefToSimVarHT)
  #= Print predefined variables=#
  for i in 1:length(matchOrder)
    if matchOrder[i] == 0
      @info "The following variable is given: $(ht[i])"
    end
  end
  for i in reverseTopologicalSort
    if matchOrder[i] != 0
      variableIdx = MetaGraphs.get_prop(digraph, i, :vID)
      equationIdx = matchOrder[variableIdx]
      @info "$(ht[variableIdx]) is solved in: $(string(indexToEquation[equationIdx]))"
    else
      @info "matchOrder for equation $i => 0. "
    end
  end
end
