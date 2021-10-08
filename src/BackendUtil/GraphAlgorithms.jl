
"
  This files contains the various graph algorithms.

  More specficially routines for matching, merging and strongly connected components

  Author: John Tinnerholm

"
module GraphAlgorithms

import LightGraphs
import MetaGraphs
import Cairo
using GraphPlot
using Compose
using DataStructures

"""
  Regular matching. (Does not solve singularities).
  Author: John Tinnerholm
  input:
        dict, adjacency list representation of the equation-variable graph
        n, the number of unknown and equations.
  output:
         assign::Array. assign(j) = i, where j is a variable and i the equation in which it is assigned.
         isSingular::Boolean,          Boolean indicating if the system is singular or not.
"""
function matching(dict::DataStructures.OrderedDict, n::Int)
  #= Global arrays for bookkeeping =#
  local assign = [0 for i in 1:n]
  local vMark = Bool[false for i in 1:n]
  local eMark = Bool[false for i in 1:n]
  """ 
    Calculates the path for equation i. 
    returns true if a path is found. 
  """
  function pathFound(i)
    eMark[i] = true
    local success = false
    local equationsI = dict.vals[i]
    for j in equationsI 
      if assign[j] == 0
        assign[j] = i
        success = true
        return success
      end
    end    
    #= Otherwise =#
    success = false
    for j in equationsI
      if vMark[j] != false
        continue
      end
      vMark[j] = true
      success = pathFound(assign[j])
      if success
        assign[j] = i
        return success
      end
    end
    return success
  end
  #=Entry of algorithm=#
  local isSingular = false
  for i in 1:n
    vMark = [false for j in 1:n]
    eMark = [false for j in 1:n]
    success = pathFound(i)
    if !success
      isSingular = !success
    end
  end
  return isSingular, assign
end

"""
Author: John Tinnerholm
  Given a order, and a graph represented as an adjacency list creates a new digraph
  representing a causalised system.
  input matchOrder, assign array(j) = i The variable j is solved in equation i
  input graph equation -> {Equation -> variables belonging to it}
  output LightGraphs.SimpleDiGraph
"""
function merge(matchOrder::Vector, graph::OrderedDict)::MetaGraphs.MetaDiGraph  "Remove function for arrays.."
  function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
  end 
  #= 
    Convert the given map into an array representation.
    Similar format to the assign matrix but represent the dependencies 
    depends(i) = {Set of variables used in equation i}.
  =#
  local depends = graph.vals
  local g = MetaGraphs.MetaDiGraph()
  #= 
    Create vertices. Each equation in matchorder is one vertex in the graph. 
  =#
  local nMatchOrder = length(matchOrder)
  for eq in 1:nMatchOrder
    LightGraphs.add_vertex!(g)
    MetaGraphs.set_prop!(g, eq, :eID, eq)
  end
  for eq in 1:nMatchOrder
    varIdx = findall(x->x==eq, matchOrder)
    if length(varIdx) == 0
      #=Was zero for eq=#
      varIdx = findall(x->x==eq, matchOrder)
      continue; 
    end
    local depVariables = remove!(depends[eq], first(varIdx))
    #= Solve for the remaining variables =#
    if ! isempty(depVariables)
      for v in depVariables
        MetaGraphs.set_prop!(g, matchOrder[v], :vID, v)
        MetaGraphs.set_prop!(g, eq, :vID, varIdx[1])
        LightGraphs.add_edge!(g, matchOrder[v], eq)
      end
    else
      MetaGraphs.set_prop!(g, eq, :vID, varIdx[1])
    end
  end
  @debug dumpGraphProperties(g)
  return g
end

"
  Dumps the properties of a given MetaDiGraph.
"
function dumpGraphProperties(g::MetaGraphs.MetaDiGraph)
  local nVertices = LightGraphs.vertices(g).stop
  local str = "Meta properties of the graph:\n"
  for i in 1:nVertices
    str *= "Properties: $(MetaGraphs.props(g, i))\n"
  end
  return str
end

"
  Topological sort 
"
function topological_sort(g::LightGraphs.AbstractGraph)::Array
  LightGraphs.topological_sort_by_dfs(g)
end

function stronglyConnectedComponents(g::LightGraphs.AbstractGraph)::Array
  LightGraphs.strongly_connected_components_kosaraju(g)
end


"""
    Plots the given equation graph as a png
"""
function plotEquationGraphPNG(filePath::String, g::LightGraphs.AbstractGraph, labels, dims = (64cm, 64cm)::Tuple)
  plot = gplot(g,
               nodelabel=labels,
               nodefillc="blue",
               nodelabelc="white",
               edgestrokec="black",
               layout=spring_layout)
  draw(Compose.PNG(filePath, dims...), plot)
end

# function plotEquationGraph(filePath::String, g::LightGraphs.AbstractGraph, dims = (64cm, 64cm)::Tuple)
#   plot = gplot(g,
#                nodefillc="blue",
#                nodelabelc="white",
#                edgestrokec="black",
#                layout=spring_layout)
#   draw(Compose.PDF(filePath, dims...), plot)
# end

function connected_components(g::LightGraphs.AbstractGraph)
  LightGraphs.connected_components(g)
end

"""
  Author: John Tinnerholm
  Tarjans algorithm
"""
function tarjan(g::OrderedDict)::Array
  tarjan(g, length(g.keys))
end

"""
 Helper function.
 It is assumed that the dict g is orderd 1->N where 1->N is the indices of the nodes.
 input g::OrderedDict
 input n::Int, the number of vertices
 output sccs: The set of strongly connected components
"""
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




end #= GraphAlgorithms =#
