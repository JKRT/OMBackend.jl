"
Returns true if simvar is either a algebraic or a state variable
"
function isStateOrAlgebraic(simvar::SimVar)::Bool
  return isAlgebraic(simvar) || isState(simvar)
end

"
Returns true if simvar is  an algebraic variable
"
function isAlgebraic(simvar::SimVar)::Bool
  res = @match simvar.varKind begin
    ALG_VARIABLE(__) => true
    _ => false
  end
end

"
Returns true if simvar is  an algebraic variable
"
function isState(simvar::SimVar)::Bool
  res = @match simvar.varKind begin
    STATE(__) => true
    _ => false
  end
end

function dumpVariableEqMapping(mapping::OrderedDict)::String
  local dump = "\n"
  local equations = keys(mapping)
  for e in equations
    variablesAtEq = "{"
    for v in mapping[e]
      variablesAtEq *= "$(v),"
    end
    variablesAtEq *= "}"
    dump *= "Equation $e: involves: $(variablesAtEq)\n"
  end
  return dump
end

"
input digraph
input variablesHT
  cref -> variable information dictonary.
output 
  An array of labels for a directed graph g.
"
function makeLabels(digraph, matchOrder, variablesHT)
  variableIndexToName::OrderedDict = makeIndexVarNameDict(matchOrder, variablesHT)
  labels = []
  for i in 1:length(matchOrder)
    variableIdx = MetaGraphs.get_prop(digraph, i, :vID)
    equationIdx = matchOrder[variableIdx]
    push!(labels, "e$(equationIdx)|$(variableIndexToName[variableIdx])|index_$(i)")
  end
  return labels
end

"
  idx -> var-name.
  Supply matching order and a ht. 
"
function makeIndexVarNameDict(matchOrder, variablesHT)::OrderedDict
  local unknownVariables = filter((x) -> isVariableOrState(x[2].varKind), collect(values(variablesHT)))
  variableIndexToName::OrderedDict = OrderedDict()
  for v in unknownVariables
    variableIndexToName[v[1]] = v[2].name
  end
  return variableIndexToName
end

function isVariableOrState(type::SimVarType)
  return @match type begin
    ALG_VARIABLE(__) => true
    STATE(__) => true
    _ => false
  end
end



"""
   This functions create and assigns indices for variables
   Thus Construct the table that maps variable name to the actual variable.
It executes the following steps:
1. Collect all variables
2. Search all states (e.g. x and y) and give them indices starting at 1 (so x=1, y=2). Then give the corresponding state derivatives (x' and y') the same indices.
3. Remaining algebraic variables will get indices starting with i+1, where i is the number of states.
4. Parameters will get own set of indices, starting at 1.
"""
function createIndices(simulationVars::Array{SimulationCode.SIMVAR})::OrderedDict{String, Tuple{Integer, SimulationCode.SimVar}}
  local ht::OrderedDict{String, Tuple{Integer, SimulationCode.SimVar}} = OrderedDict()
  local stateCounter = 0
  local parameterCounter = 0
  local numberOfStates = 0
  for var in simulationVars
    @match var.varKind begin
      SimulationCode.STATE(__) => begin
        stateCounter += 1
        @set var.index = SOME(stateCounter)
        stVar = SimulationCode.SIMVAR(var.name, var.index, SimulationCode.STATE_DERIVATIVE(var.name), var.attributes)
        push!(ht, var.name => (stateCounter, var))
        #=Adding the state derivative as well=#
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
        var = @set var.index = SOME(algIndexCounter)
        push!(ht, var.name => (var.index.data, var))
      end
      _ => continue
    end
  end
  return ht
end

function createEquationVariableBidirectionGraph(equations, allBackendVars, crefToSimVarHT)::OrderedDict
  local eqCounter::Int = 1
  local varCounter::Int = 1
  local variableEqMapping = OrderedDict()
  local unknownVariables = filter((x) -> BDAEUtil.isVariable(x.varKind), allBackendVars)
  local stateVariables = filter((x) -> BDAEUtil.isState(x.varKind), allBackendVars)
  #= Treat states as solved =#
  nEquations = length(equations)  - length(stateVariables)
  nVariables = length(unknownVariables)
  @assert(nEquations == nVariables, "The set of variables != set of equations: #Variables: $nVariables, #Equations $nEquations")
  for eq in equations
    #= Fetch all variables beloning to the specific equation =#
    variablesForEq = Backend.BDAEUtil.getAllVariables(eq, allBackendVars)
    variableEqMapping["e$(eqCounter)"] = sort(getIndiciesOfVariables(variablesForEq, crefToSimVarHT))
    eqCounter += 1
  end
  return variableEqMapping
end


function getIndiciesOfVariables(variables, crefToSimVarHT::OrderedDict{String, Tuple{Integer, SimVar}})::Array
  indicies = []
  for v in variables
    candidate = crefToSimVarHT[string(v)]
    if isAlgebraic(candidate[2]) #Just add algebraic variables
      push!(indicies, candidate[1])
    elseif isState(candidate[2])
      push!(indicies, candidate[1])
    else
      continue
    end
  end
  return indicies
end
