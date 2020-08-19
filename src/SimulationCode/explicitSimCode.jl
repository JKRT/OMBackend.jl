
"""
TODO:
This function transforms DAE-mode style simcode into explicit simcode.
In explicit simcode the residual equations are sorted vertically and horisontally
#1 Get all variables beloning to a specific residual
"""
function transformToExplicitSimCode(backendDAE::BDAE.BACKEND_DAE)::SimulationCode.EXPLICIT_SIMCODE
  local equationSystems::Array = backendDAE.eqs
  local allOrderedVars::Array{BDAE.Var} = [v for es in equationSystems for v in es.orderedVars.varArr]
  local allSharedVars::Array{BDAE.Var} = getSharedVariablesLocalsAndGlobals(backendDAE.shared)
  local allBackendVars::Array{BDAE.Var} = vcat(allOrderedVars, allSharedVars)
  local simVars::Array{SimulationCode.SIMVAR} = allocateAndCollectSimulationVariables(allBackendVars)
  local crefToSimVarHT = createExplicitIndicies(simVars)
  local equations = [eq for es in equationSystems for eq in es.orderedEqs]
  #= Split equations into three parts. Residuals whenEquations and If-equations =#
  (resEqs,whenEqs,ifEqs) = allocateAndCollectSimulationEquations(equations)
  #=Only assume I have residuals for now=#
  assert(length(whenEqs) == 0)
  assert(length(ifEqs) == 0)
  local indexToEquation = createEquationIndicies(resEqs)
  #= Create equation <-> variable mapping =#
  local variableEqMapping = createEquationVariableBidirectionGraph(equations, allBackendVars)
  (isSingular::Bool, matchOrder::Array) = GraphAlgorithms.matching(dict, length(dict.keys))
  (_, labels, sortedGraph::OrderedDict) = GraphAlgorithms.merge(machOrder, variableEqMapping)
  stronglyConnectedComponents::Array = GraphAlgorithms.tarjan(sortedGraph)
  return SimulationCode.EXPLICIT_SIMCODE(backendDAE.name,
                                         crefToSimVarHT,
                                         indexToEquation,
                                         variableEqMapping,
                                         resEqs,
                                         whenEqs,
                                         ifEqs,
                                         isSingular,
                                         matchOrder,
                                         sortedGraph,
                                         sortedGraph)
end

function createEquationVariableBidirectionGraph(equations, allBackendVars, crefToSimVarHT)::Dict
  local eqCounter::Int = 1
  local varCounter::Int = 1
  local variableEqMapping = OrderedDict()
  assert(length(equations) == length(allBackendVars), "The set of variables != set of equations")
  for eq in equations
    variablesForEq = getAllVariables(eq, allBackendVars)
    variableEqMapping["e$(eqCounter)"] = getIndiciesOfVariables(variablesForEq,crefToSimVarHT)
    eqCounter += 1
  end
  return variableEqMapping
end

function getIndiciesOfVariables(variables, crefToSimVarHT::OrderedDict{String, Tuple{Integer, SimVar}})::Array
  indicies = []
  for v in variables
    push!(indicies, crefToSimVarHT[string(v)][1])
  end
  return indicies
end

function createEquationIndicies(resEqs)
  local index = 1;
  local ht::OrderedDict()
  for e in resEqs
    ht["e$(index)"] = e
  end
  return ht
end

function createExplicitIndicies(simVar:Array{SimulationCode.SIMVAR})
  local ht::OrderedDict{String, Tuple{Integer, SimulationCode.SimVar}} = Dict()
  local stateCounter = 0
  local parameterCounter = 0
  local numberOfStates = 0
  for var in simulationVars
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
