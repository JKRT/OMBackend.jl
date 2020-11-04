import ..OMBackend

"""
  TODO:
  This function transforms DAE-mode style simcode into explicit simcode.
  In explicit simcode the residual equations are sorted vertically and horisontally
    #1 Get all variables beloning to a specific residual
 """
function transformToExplicitSimCode(backendDAE::BDAE.BACKEND_DAE)::SimulationCode.EXPLICIT_SIM_CODE
  local equationSystems::Array = backendDAE.eqs
  local allOrderedVars::Array{BDAE.Var} = [v for es in equationSystems for v in es.orderedVars.varArr]
  local allSharedVars::Array{BDAE.Var} = getSharedVariablesLocalsAndGlobals(backendDAE.shared)
  local allBackendVars::Array{BDAE.Var} = vcat(allOrderedVars, allSharedVars)
  local simVars::Array{SimulationCode.SIMVAR} = allocateAndCollectSimulationVariables(allBackendVars)
  local crefToSimVarHT::OrderedDict = createExplicitIndicies(simVars)
  local equations = [eq for es in equationSystems for eq in es.orderedEqs]
  #= Split equations into three parts. Residuals whenEquations and If-equations =#
  (resEqs,whenEqs,ifEqs) = allocateAndCollectSimulationEquations(equations)
  #=Only assume I have residuals for now=#
  @assert(length(whenEqs) == 0, "IF EQUATION NOT YET SUPPORTED IN EXPLICIT CODE GEN")
  @assert(length(ifEqs) == 0, "WHEN EQUATION NOT YET SUPPORTED IN EXPLICIT CODE GEN")
  local indexToEquation::OrderedDict = createEquationIndicies(resEqs)
  #= Create equation <-> variable mapping =#
  local eqVariableMapping = createEquationVariableBidirectionGraph(equations, allBackendVars, crefToSimVarHT)
  @info "eqVariableMapping: $eqVariableMapping"
  (isSingular::Bool, matchOrder::Array) = GraphAlgorithms.matching(eqVariableMapping, length(eqVariableMapping.keys))
  @info "is singular: $isSingular"  
  (g, labels, sortedGraph::OrderedDict) = GraphAlgorithms.merge(matchOrder, eqVariableMapping)
  stronglyConnectedComponents::Array = GraphAlgorithms.tarjan(sortedGraph)
  if OMBackend.PLOT_EQUATION_GRAPH
    @info "Plotting"
    GraphAlgorithms.plotEquationGraph("./digraphOutput$(backendDAE.name).pdf", g, labels)
  end
  #= Reorder the residuals =#
  reOrderedResiduals = []
  reverseTopologicalSort = stronglyConnectedComponents

  #= Check algebraic loops =#
  loopsExist = false
  for i in reverseTopologicalSort
    if length(i) > 1
      loopsExist = true
      break
    end
  end
    #= Tearing could be added here =#
  @info "Strongly connected components!" reverseTopologicalSort
  if ! loopsExist
    for i in reverseTopologicalSort
      push!(reOrderedResiduals, resEqs[i[1]])
    end
  else
    @info "Loop encountered exiting..."
    @error "Unresolved algebraic loop"
  end

  return SimulationCode.EXPLICIT_SIM_CODE(backendDAE.name,
                                         crefToSimVarHT,
                                         indexToEquation,
                                         eqVariableMapping,
                                         reOrderedResiduals,
                                         whenEqs,
                                         ifEqs,
                                         isSingular,
                                         matchOrder,
                                         sortedGraph,
                                         stronglyConnectedComponents)
end

"
 This function
"
function createEquationVariableBidirectionGraph(equations, allBackendVars, crefToSimVarHT)::OrderedDict
  local eqCounter::Int = 1
  local varCounter::Int = 1
  local variableEqMapping = OrderedDict()
  local unknownVariables = filter((x) -> BDAEUtil.isStateOrVariable(x.varKind), allBackendVars)
  nEquations = length(equations)
  nVariables = length(unknownVariables)
  @assert(nEquations == nVariables, "The set of variables != set of equations: #Variables: $nVariables, #Equations $nEquations")

  @info "After getting variables # $nVariables"
  for eq in equations
    #= Fetch all variables beloning to the specific equation =#
    variablesForEq = Backend.BackendEquation.getAllVariables(eq, allBackendVars)
    variableEqMapping["e$(eqCounter)"] = sort(getIndiciesOfVariables(variablesForEq, crefToSimVarHT))
    eqCounter += 1
  end
  @info variableEqMapping
  return variableEqMapping
end

function getIndiciesOfVariables(variables, crefToSimVarHT::OrderedDict{String, Tuple{Integer, SimVar}})::Array
  indicies = []
  @info "crefToSimVar" crefToSimVarHT
  @info "Variables:" variables
  for v in variables
    candidate = crefToSimVarHT[string(v)]
    @info "Candidate:" candidate
    if ! isStateOrAlgebraic(candidate[2])
      continue
    else
      push!(indicies, candidate[1])
    end
  end
  return indicies
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

function createExplicitIndicies(simulationVars::Array{SIMVAR})
  local ht::OrderedDict{String, Tuple{Integer, SimulationCode.SimVar}} = OrderedDict()
  local stateCounter = 0
  local parameterCounter = 0
  local numberOfStates = 0
  for var in simulationVars
    @match var.varKind begin
      SimulationCode.STATE(__) => begin
        stateCounter += 1
        @set! var.index = SOME(stateCounter)
        stVar = SimulationCode.SIMVAR(var.name, var.index, SimulationCode.STATE_DERIVATIVE(var.name), var.attributes)
        push!(ht, var.name => (stateCounter, var))
        #=Adding the state derivative as well =#
        push!(ht, "der($(var.name))" => (stateCounter, stVar))
      end
      SimulationCode.PARAMETER(__) => begin
        parameterCounter += 1
        push!(ht, var.name => (parameterCounter, var))
      end
      _ => continue
    end
  end
  local algIndexCounter::Integer = stateCounter
  for var in simulationVars
    @match var.varKind begin
      SimulationCode.ALG_VARIABLE(__) => begin
        algIndexCounter += 1
        @set! var.index = SOME(algIndexCounter)
        push!(ht, var.name => (var.index.data, var))
      end
      _ => continue
    end
  end
  return ht
end
