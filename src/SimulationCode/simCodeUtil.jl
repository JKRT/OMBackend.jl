#=
  This file contains various utility functions related to simulation code.
=#
"""
  Returns true if simvar is either a algebraic or a state variable
"""
function isStateOrAlgebraic(simvar::SimVar)::Bool
  return isAlgebraic(simvar) || isState(simvar)
end

"""
  Returns true if the simulation code variable is discrete.
"""
function isDiscrete(simVar::SimVar)::Bool
  res = @match simVar.varKind begin
    DISCRETE(__) => true
    _ => false
  end
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

function isOCCVar(simVar::SimVar)::Bool
  res = @match simVar.varKind begin
    OCC_VARIABLE(__) => true
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
Author: John & Andreas
   This functions create and assigns indices for variables
   Thus Construct the table that maps variable name to the actual variable.
It executes the following steps:
1. Collect all variables
2. Search all states (e.g. x and y) and give them indices starting at 1 (so x=1, y=2). Then give the corresponding state derivatives (x' and y') the same indices.
3. Remaining algebraic variables will get indices starting with i+1, where i is the number of states.
4. Parameters will get own set of indices, starting at 1.
5. Discrete shares the index with the states and starts at #states + 1
6. OCC Variables also shares the indices with the states and starts at #discretes + 1
"""
function createIndices(simulationVars::Vector{SimulationCode.SIMVAR})::OrderedDict{String, Tuple{Integer, SimulationCode.SimVar}}
  local ht::OrderedDict{String, Tuple{Integer, SimulationCode.SimVar}} = OrderedDict()
  local stateCounter = 0
  local parameterCounter = 0
  local discretes = SimulationCode.SIMVAR[]
  local occVariables = SimulationCode.SIMVAR[]
  local numberOfStates = 0
  for var in simulationVars
    @match var.varKind begin
      SimulationCode.STATE(__) => begin
        stateCounter += 1
        @assign var.index = SOME(stateCounter)
        stVar = SimulationCode.SIMVAR(var.name, var.index, SimulationCode.STATE_DERIVATIVE(var.name), var.attributes)
        push!(ht, var.name => (stateCounter, var))
        #=Adding the state derivative as well=#
        push!(ht, "der($(var.name))" => (stateCounter, stVar))
      end
      SimulationCode.OCC_VARIABLE(__) => begin
        push!(occVariables, var)
      end
      SimulationCode.PARAMETER(__) => begin
        parameterCounter += 1
        push!(ht, var.name => (parameterCounter, var))
      end
      SimulationCode.DISCRETE(__) => begin
        push!(discretes, var)
      end
      _ => continue
    end
  end
  local discreteCounter = stateCounter
  for var in discretes
    discreteCounter += 1
    push!(ht, var.name => (discreteCounter, var))
  end
  local occCounter = discreteCounter
  for var in occVariables
    occCounter += 1
    push!(ht, var.name => (occCounter, var))
  end
  local algIndexCounter::Int = occCounter #Change 2022-09-10
  local algSortingIdx::Int = stateCounter #This idx is used by the backend sorting algorithms
  for var in simulationVars
    @match var.varKind begin
      SimulationCode.ALG_VARIABLE(__) => begin
        algIndexCounter += 1
        algSortingIdx += 1
        @assign var.index = SOME(algIndexCounter)
        @assign var.varKind = ALG_VARIABLE(algSortingIdx)
        push!(ht, var.name => (var.index.data, var))
      end
      _ => continue
    end
  end
  @info "stateCounter" stateCounter
  @info "discreteCounter" discreteCounter
  @info "algIndexCounter" algIndexCounter
  return ht
end

"""
  Given a set of residual equations, a set of if-equations and the set of all backend variables.
  This function creates a bidrectional graph between these equations and the supplied variables.
  (Note: If we need to do index reduction there might be empty equations here).
"""
function createEquationVariableBidirectionGraph(equations::RES_T,
                                                ifEquations::IF_T,
                                                allBackendVars::VECTOR_VAR,
                                                stringToSimVarHT)::OrderedDict where{RES_T, IF_T, VECTOR_VAR}
  local eqCounter::Int = 0
  local variableEqMapping = OrderedDict()
  local unknownVariables = filter((x) -> BDAEUtil.isVariable(x.varKind), allBackendVars)
  #=TODO: The set of discrete variables are currently not in use. =#
  local discreteVariables = filter((x) -> BDAEUtil.isDiscrete(x.varKind), allBackendVars)
  local stateVariables = filter((x) -> BDAEUtil.isState(x.varKind), allBackendVars)
  local algebraicAndStateVariables = vcat(unknownVariables, stateVariables)
  local nDiscretes = length(discreteVariables)
  @info "#stateVariables" length(stateVariables)
  @info "#discretes" nDiscretes
  @info "#algebraic" length(unknownVariables)
  @info "#equations" length(equations)

  for eq in equations
    eqCounter += 1
    variablesForEq = Backend.BDAEUtil.getAllVariables(eq, algebraicAndStateVariables)
    # @debug "Variables in equation:"
    # println("Equation:", string(eq))
    # println("Variables:")
    # for v in variablesForEq
    #   println("\t", string(v))
    # end
    local indices = getIndiciesOfVariables(variablesForEq, stringToSimVarHT)
    # @debug "Indices where:"
    # for idx in indices
    #   println("\t", string(idx))
    # end
    variableEqMapping["e$(eqCounter)"] = sort(indices)
  end
  #=
   There is an additional case to consider.
   If some variables are solved by *some* branch
   (The branches are required to be balanced for ordinary if-equations)
   in an if equation it should be included in the mapping.
  =#
  for ifEq in ifEquations
    #= Select one branch =#
    ifEqBranch = listArray(listGet(ifEq.eqnstrue, 1))
    for eq in ifEqBranch
      eqCounter += 1
      variablesForEq = Backend.BDAEUtil.getAllVariables(eq, algebraicAndStateVariables)
      variableEqMapping["e$(eqCounter)"] = sort(getIndiciesOfVariables(variablesForEq, stringToSimVarHT))
    end
  end
  return variableEqMapping
end

"""
 Same as the other createEquationVariableBidirectionGraph however, here we assume a system that have no if-equations.
"""
function createEquationVariableBidirectionGraph(equations::RES_T,
                                                allBackendVars::VECTOR_VAR,
                                                stringToSimVarHT)::OrderedDict where{RES_T, IF_T, VECTOR_VAR}
  local eqCounter::Int = 0
  local variableEqMapping = OrderedDict()
  local unknownVariables = filter((x) -> BDAEUtil.isVariable(x.varKind), allBackendVars)
  local discreteVariables = filter((x) -> BDAEUtil.isDiscrete(x.varKind), allBackendVars)
  local stateVariables = filter((x) -> BDAEUtil.isState(x.varKind), allBackendVars)
  local algebraicAndStateVariables = vcat(unknownVariables, stateVariables)
  local nDiscretes = length(discreteVariables)
  @debug "#stateVariables" length(stateVariables)
  @debug "#algebraic" length(unknownVariables)
  @debug "#equations" length(equations)
  for eq in equations
    eqCounter += 1
    variablesForEq = Backend.BDAEUtil.getAllVariables(eq, algebraicAndStateVariables)
    variableEqMapping["e$(eqCounter)"] = sort(getIndiciesOfVariables(variablesForEq, stringToSimVarHT))
  end
  return variableEqMapping
end

"""
  Given a set of variables and a dictonary that maps the component reference
  to some simulation code variable.
This function returns the indices of these variables.
*NOTE*:
  That the index of the algebraic variable is treated in a different way here.
  That is, the index of the algebraic variable is offset by the total number of discrete variables
"""
function getIndiciesOfVariables(variables,
                                stringToSimVarHT::OrderedDict{String, Tuple{Integer, SimVar}})
  local indicies = Int[]
  for v in variables
    idx, var  = stringToSimVarHT[string(v)]
    if isAlgebraic(var)
      #= Algebraic variables use a special idx for backend sorting purposes. =#
      push!(indicies, var.varKind.sortIdx)
    elseif isState(var)
      push!(indicies, idx)
    elseif isOCCVar(var)
      push!(indicies, idx)
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

"""
  Creates a OCC graph.
  Returns the graph and the root variables.
(This function also adds info to the model)
"""
function getOCCGraph(flatModel)
  unresolvedFlatModel = OMFrontend.Main.FLAT_MODEL(flatModel.name,
                                                   flatModel.variables,
                                                   flatModel.unresolvedConnectEquations,
                                                   flatModel.initialEquations,
                                                   flatModel.algorithms,
                                                   flatModel.initialAlgorithms,
                                                   ImmutableList.nil,
                                                   NONE(),
                                                   flatModel.DOCC_equations,
                                                   flatModel.unresolvedConnectEquations,
                                                   flatModel.active_DOCC_Equations,
                                                   flatModel.comment)
  local name::String = unresolvedFlatModel.name
  local conns::OMFrontend.Main.Connections
  local conn_eql::List{OMFrontend.Main.Equation}
  local csets::OMFrontend.Main.ConnectionSets.Sets
  local csets_array::Vector{List{OMFrontend.Main.Connector}}
  local ctable::OMFrontend.Main.CardinalityTable.Table
  local broken::OMFrontend.Main.BrokenEdges = ImmutableList.nil
  local rootEquations::List{OMFrontend.Main.Equation} = MetaModelica.nil
  local rootReferenceVariables::Vector{Tuple} = Tuple{OMFrontend.Main.NFComponentRef,
                                                      OMFrontend.Main.NFComponentRef}[]
  (unresolvedFlatModel, conns) = OMFrontend.Main.collect(unresolvedFlatModel)
  (unresolvedFlatModel, conns) = OMFrontend.Main.elaborate(unresolvedFlatModel, conns)
  if OMFrontend.Main.System.getHasOverconstrainedConnectors()
    (_, broken, graph) = OMFrontend.Main.handleOverconstrainedConnections(unresolvedFlatModel, conns, name)
    (roots, _, broken) = OMFrontend.Main.findResultGraph(graph, name)
    rootEquations = OMFrontend.Main.findRootEquations(roots, graph,
                                                      unresolvedFlatModel.equations)
    for re in rootEquations
      push!(rootReferenceVariables,
            (re.lhs, re.rhs))
    end
  end
  OMFrontend.Main.printNFOCConnectionGraph(graph)
  println("Broken equations")
  OMFrontend.Main.printFlatEdges(broken)
  #= Remove the broken edge from the set of edges =#
  @assign graph.connections = arrayList(filter((x)->(!in(x, broken)), listArray(graph.connections)))
  println("Graph with broken removed")
  OMFrontend.Main.printNFOCConnectionGraph(graph)
  #= Convert the branches to regular edges =#
  local uniqueRoots = graph.uniqueRoots
  local definiteRoots = graph.definiteRoots
  local potentialRoots = graph.potentialRoots
  @info "Length uniqueroots" length(uniqueRoots)
  @info "definiteRoots" length(definiteRoots)
  @info "Potential roots" length(potentialRoots)
  println("Print: definiteRoots")
  for r in potentialRoots
    println(OMFrontend.Main.toString(first(r)))
  end
  println("Print: final roots")
  for r in roots
    println(OMFrontend.Main.toString(r))
  end
  #= Get the roots involved in the structural change =#
  rootVariables::List{OMFrontend.Main.ComponentRef} = MetaModelica.list(r for r in roots)
  #  @info "roots" roots
  #= Create a graph that we can search. =#
  local connectionEdges = convertFlatEdgeToEdges(graph.connections)
  @info typeof(connectionEdges)
  local allEdges = listAppend(connectionEdges, graph.branches)
  local searchGraph = createSearchGraph(allEdges)
  return (searchGraph, rootVariables, rootReferenceVariables)
end

"""
 Convert the component references to the backend representation and create an adjecency list representation.
"""
function createSearchGraph(allEdges)
  local edgeSet = Dict()
  local searchGraph = Dict{String, Vector{String}}()
  for edge in allEdges
    @match (e1, e2) = edge
    local s1 = OMFrontend.Main.toString(e1)
    local s2 = OMFrontend.Main.toString(e2)
    edgeSet[s1] = e1
    edgeSet[s2] = e2
  end
  for edge in keys(edgeSet)
    searchGraph[edge] = String[]
  end
  for edge in allEdges
    @match (e1, e2) = edge
    local s1 = OMFrontend.Main.toString(e1)
    local s2 = OMFrontend.Main.toString(e2)
    push!(searchGraph[s1], s2)
    push!(searchGraph[s2], s1)
  end
  return searchGraph
end

"""
  Given a list of flat edges convert them to edges.
"""
function convertFlatEdgeToEdges(connections)
  newEdges = Tuple[]
  for connection in connections
    @match connection begin
      (c1, c2, _)  => begin
        push!(newEdges, (c1, c2))
      end
    end
  end
  return arrayList(newEdges)
end


"""
 This function returns true if a backend variable is in the set of of overconstrained connector variables (occVariables).
"""
function isOverconstrainedConnectorVariable(simVarName::String, occVariables::Vector{String})
  #= Inefficient, can be done better... =#
  local isOCCVar = simVarName in occVariables
  return isOCCVar
end

"""
  Get variables that should be marked as irreductable.
"""
function getIrreductableVars(ifEquations::Vector{BDAE.IF_EQUATION},
                             whenEqs::Vector{BDAE.WHEN_EQUATION},
                             algebraicAndStateVariables::VECTOR_VAR) where {VECTOR_VAR}
  local irreductables = Vector{Any}[]
  for eq in ifEquations
    variablesForEq = Backend.BDAEUtil.getAllVariables(eq, algebraicAndStateVariables)
    push!(irreductables, variablesForEq)
  end
  #TODO: Fix when equations
  # for eq in whenEqs
  #   variablesForEq = Backend.BDAEUtil.getAllVariables(eq, algebraicAndStateVariables)
  #   push!(variablesForEq, irreductables)
  # end
  local irreductablesAsStr = map(x -> string(x), collect(Iterators.flatten(irreductables)))
  irreductablesAsStr = filter(x -> x != "time", irreductablesAsStr)
  return irreductablesAsStr
end
