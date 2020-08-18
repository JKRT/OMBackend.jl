
"
  Regular matching. Does not solve singularities.
  Author: John Tinnerholm
  input:
        dict, adjacency list representation of the equation-variable graph
        n, the number of unknown and equations.
  output:
         assign::Array. assign(j) = i, where j is a variable and i the equation in which it is assigned.
         isSingular::Boolean,          Boolean indicating if the system is singular or not.
"

function matching(dict::OrderedDict,n)
  "Calculates the path for equation i"
  function PF(i)
    eMark[i] = true
    success = false
    for j in dict.vals[i]
      if assign[j] == 0
        assign[j] = i
        success = true
        return success;
      end
    end
    #= Otherwise =#
    success = false
    for j in dict.vals[i]
      if vMark[j] != false
        continue
      end
      vMark[j] = true
      success = PF(assign[j])
      if success
        assign[j] = i
      end
    end
    return success
  end
  #=Entry of algorithm=#
  local assign = [0 for i in 1:n]
  local vMark = []
  local eMark = []
  local isSingular = true
  @info "Starting"
  for i in 1:n
    vMark = [false for j in 1:n]
    eMark = [false for j in 1:n]
    success = PF(i)
    if !success
      isSingular = !success
    end
  end
  return isSingular,assign
end



"
  Author: John Tinnerholm
  Given a order, and a graph represented as an adjacency list creates a new digraph
  representing a causalised system.
"
function merge(matchOrder, graph::OrderedDict)
  "
    Remove function for arrays...
  "
  function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
  end

  local counter = 0
  local values = graph.vals
  local g = SimpleDiGraph()
  local arrRepresentation = []
  for i in matchOrder
    add_vertex!(g)
  end
  for eq in eqOrder
    counter += 1
    c = remove!(values[eq], counter)
    #=Now I need to know. What equation solves for the remaining variables=#
    solvedIn = matchOrder[c]
    for e in solvedIn
      push!(arrRepresentation, [e, eq])
      add_edge!(g, e, eq)
    end
  end
  #=Create the dict repr=#
  dictRepr = OrderedDict()
  for i in 1:length(matchOrder)
    key = i
    dictRepr[key] = []
  end
  for i in arrRepresentation
    key = i[1]
    push!(dictRepr[key], i[2])
  end
  labels = ["f$(matchOrder[i])|z$(i)"  for i in 1:length(matchOrder)]
  return g,labels,dictRepr
end



"
  Author: John Tinnerholm
  Tarjans algorithm
"
function tarjan(g::OrderedDict)::Array
  tarjan(g, length(g.keys))
end

"
 Helper function.
 It is assumed that the dict g is orderd 1->N where 1->N is the indices of the nodes.

 input g::OrderedDict
 input n::Int, the number of vertices
 output sccs: The set of strongly connected components
"
function tarjan(g::OrderedDict, n)::Array
  function strongConnect(v::Int)
    vIndicies[v] = index
    vLowLinks[v] = index
    index += 1
    push!(stack, v)
    vOnStack[v] = true
    v2S = g[v]
    for v2 in v2S
      if vIndicies[v2] == 0
        strongConnect(v2)
        vLowLinks[v] = min(vLowLinks[v], vLowLinks[v2])
      elseif vOnStack[v2]
        vLowLinks[v] = min(vLowLinks[v], vLowLinks[v2])
      end
    end
    if vLowLinks[v] == vIndicies[v]
      stronglyConnectedComponentCandidate = []
      while true
        w = pop!(stack)
        vOnStack[w] = false
        push!(stronglyConnectedComponentCandidate, w)
        if w == v
          break
        end
      end
      push!(sccs, stronglyConnectedComponentCandidate)
    end
  end
  local index = 1
  local sccs = []
  local stack = []
  #=
  Incidices for the vertices 1->n
  0 = undefined
  =#
  local vIndicies = [0 for i in 1:n]
  local vLowLinks = [0 for i in 1:n]
  local vOnStack = [false for i in 1:n]
  for v in g.keys
    #= If v is undefined =#
    if vIndicies[v] == 0
      strongConnect(v)
    end
  end
  return sccs
end
