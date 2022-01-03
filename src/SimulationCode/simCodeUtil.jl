#=
  This file contains various utility functions related to simulation code.
=#
"
Returns true if simvar is either a algebraic or a state variable
"
function isStateOrAlgebraic(simvar::SimVar)::Bool
  return isAlgebraic(simvar) || isState(simvar)
end

"""
Returns true if simvar is an algebraic variable
"""
function isAlgebraic(simvar::SimVar)::Bool
  res = @match simvar.varKind begin
    ALG_VARIABLE(__) => true
    _ => false
  end
end

"""
  Fetches the last identifier of a variable.
That is:
getLastIdentOfVar(Foo.Bar.x) => x
"""
function getLastIdentOfVar(var)::String
  getIdentOfComponentReference(var.varName)
end

"""
  Fetches the last ident of a component reference 
"""
function getIdentOfComponentReference(cr)::String
  return begin
    @match cr begin
      DAE.CREF_QUAL(ident = ident, componentRef = componentRef) => begin
        getIdentOfComponentReference(componentRef)
      end
      DAE.CREF_IDENT(ident) => begin
        ident
      end
      DAE.CREF_ITER(ident = ident) => begin
        throw("Case not handled")
      end
    end
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

"""
  Prints what equation involves which variable.
"""
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

"""
input digraph
input variablesHT
  cref -> variable information dictonary.
output 
  An array of labels for a directed graph g.
"""
function makeLabels(digraph, matchOrder, variablesHT)
  variableIndexToName::OrderedDict = makeIndexVarNameDict(matchOrder, variablesHT)
  labels = []
  for i in 1:length(matchOrder)
    try
      variableIdx = MetaGraphs.get_prop(digraph, i, :vID)
      equationIdx = matchOrder[variableIdx]
      idxToName = variableIndexToName[variableIdx]
      push!(labels, "e$(equationIdx)|$(idxToName)|index_$(i)")
    catch #= For instance the case when a vertex v does not have a prop =#
      idxToName = variableIndexToName[i]
      push!(labels, "e$(NONE)|$(idxToName)|index_$(i)")
    end
  end
  return labels
end


"""
  idx -> var-name.
  Supply matching order and a ht. 
"""
function makeIndexVarNameDict(matchOrder, variablesHT)::DataStructures.OrderedDict
  local unknownVariables = filter((x) -> isVariableOrState(x[2].varKind), collect(values(variablesHT)))
  variableIndexToName::DataStructures.OrderedDict = DataStructures.OrderedDict()
  for v in unknownVariables
    variableIndexToName[v[1]] = v[2].name
  end
  return variableIndexToName
end

"""
  idx -> var-name.
  Supply matching order and a ht. 
"""
function makeIndexVarNameUnorderedDict(matchOrder, variablesHT)::Dict
  local unknownVariables = filter((x) -> isVariableOrState(x[2].varKind), collect(values(variablesHT)))
  variableIndexToName::Dict = DataStructures.OrderedDict()
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
Author: johti17 & Andreas
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

"""
  Given the available residual equations and the set of backend variables 
  create a bidrection graph between these equations and the supplied variables.
Note: If we need to do index reduction there might be empty equations here.

"""
function createEquationVariableBidirectionGraph(equations::RES_T,
                                                allBackendVars::VECTOR_VAR,
                                                stringToSimVarHT)::OrderedDict where{RES_T, VECTOR_VAR}
  local eqCounter::Int = 0
  local variableEqMapping = OrderedDict()
  local unknownVariables = filter((x) -> BDAEUtil.isVariable(x.varKind), allBackendVars)
  local stateVariables = filter((x) -> BDAEUtil.isState(x.varKind), allBackendVars)
  local algebraicAndStateVariables = vcat(unknownVariables, stateVariables)
  #= Treat states as solved =#
  nEquations = length(equations)  - length(stateVariables)
  nVariables = length(unknownVariables)
  #= 
  TODO:
  Assert that the set of known variables has potential equations in which they can be used 
  (Introduction of IF equations made this operation more complicated)
  Assumed to have been checked by the frontend?
  =#
  for eq in equations
    #= Fetch all variables belonging to the specific equation =#
    eqCounter += 1
    variablesForEq = Backend.BDAEUtil.getAllVariables(eq, algebraicAndStateVariables)
    variableEqMapping["e$(eqCounter)"] = sort(getIndiciesOfVariables(variablesForEq, stringToSimVarHT))
  end
  return variableEqMapping
end


"""
  Given a set of variables and a ordered dict that maps the component reference 
  to a simulation code variable returns the indices of these variables.
"""
function getIndiciesOfVariables(variables, stringToSimVarHT::OrderedDict{String, Tuple{Integer, SimVar}})::Array
  local indicies = Int[]
  for v in variables
    candidate = stringToSimVarHT[string(v)]
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

"""
  Returns the residual equation a specfic variable is solved in.
  We search for this equation among the residuals in the context.
  The context should be either the top level simcode or a specific branch of some if equation.
"""
function getEquationSolvedIn(variable::V, context::C) where {V, HT, C}
  local ht = context.stringToSimVarHT
  local variableIdx = ht[variable][1]
  local equationIdx = context.matchOrder[variableIdx]
  #= Return the equation at this specific index =#
  return context.residualEquations[equationIdx]
end
