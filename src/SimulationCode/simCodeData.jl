
"""
Kind of a simulation variable
"""
abstract type SimVarType end

"""
State variable
"""
struct  STATE <: SimVarType
end

"""
 State Derivative
"""
struct  STATE_DERIVATIVE <: SimVarType
  varName::String
end

"""
Algebraic variable
"""
struct ALG_VARIABLE <: SimVarType end

"""
Input variable
"""
struct  INPUT <: SimVarType end

"""
Parameter variable
"""
struct PARAMETER <: SimVarType
  bindExp::Option{DAE.Exp}
end


"""
Abstract type for a simulation variable
"""
abstract type SimVar end

"""
Variable data type used for code generation
"""
struct SIMVAR <: SimVar
  "Readable name of variable"
  name :: String
  "Index of variable, 0 based, type based"
  index::Option{Integer}
  "Kind of variable, one of SimulationCode.SimVarType"
  varKind::SimVarType
  "Variable attributes. Same as in DAE"
  attributes::Option{DAE.VariableAttributes}
end

"Abstract type for simulation code"
abstract type SimCode end

"""
  Root data structure containing information required for code generation to
  generate simulation code for a Modelica model.
"""
struct SIM_CODE <: SimCode
  name::String
  "Mapping of names to the corresponding variable"
  crefToSimVarHT::Dict{String, Tuple{Integer, SimVar}}
  "Different equations stored within simulation code"
  residualEquations::Array{BDAE.RESIDUAL_EQUATION}
  whenEquations::Array{BDAE.WHEN_EQUATION}
  ifEquations::Array{BDAE.IF_EQUATION}
  "True if the system that we are solving is singular"
  isSingular::Bool
    "
    The merged graph. E.g digraph constructed from matching info.
    The indicies are the same as above and they are shared.
    If the system is singular tearing is needed.
   "
  equationGraph::LightGraphs.AbstractGraph
  "
    The reverse topological sort of the equation-graph
  "
  stronglyConnectedComponents::Array
end


"""
  This is the explicit representation of SimCode.
  In this representation all residual equations
  are sorted horisontally and vertically.
  This representation is selected if we do not use DAE-Mode

"""
struct EXPLICIT_SIM_CODE <: SimCode
  "Name of the model"
  name::String
  "
   Mapping of names to the corresponding variable.
   Each variable has a unique index.
   The purpose of this mapping is to have aliases for each unique index.
   We also need to keep track of each simulation variable.
  "
  nameToVar::OrderedDict{String, Tuple{Integer, SimVar}}
  indexToEquation::OrderedDict{Int, BDAE.RESIDUAL_EQUATION}
  "Equation <-> Variable graph (Bidirectional)"
  eqVariableMapping::OrderedDict{String, Array{Int}}

  "Regular equations are encoded as residuals"
  residualEquations::Array{BDAE.RESIDUAL_EQUATION}  
  whenEquations::Array{BDAE.WHEN_EQUATION}
  ifEquations::Array{BDAE.IF_EQUATION}
  "
   If matching resulted in singularity. 
  "
  isSingular::Bool
  "
    The match order:
    Result of assign array, e.g array(j) = equation_i
  "
  matchOrder::Array{Int}

  "
    The merged graph. E.g digraph constructed from matching info.
    The indicies are the same as above and they are shared.
    If the system is singular tearing is needed.
  "
  sortedGraph::LightGraphs.AbstractGraph
  "
    The strongly connected components of the sorted graph.
    This information can be used for tearing.
  "
  stronglyConnectedComponents::Array
end
