
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
struct SIMVAR{T0 <: String, T1 <: Option{Int}, T2 <: SimVarType, T3 <: Option{DAE.VariableAttributes}} <: SimVar
  "Readable name of variable"
  name::T0
  "Index of variable, 0 based, type based"
  index::T1
  "Kind of variable, one of SimulationCode.SimVarType"
  varKind::T2
  "Variable attributes. Same as in DAE"
  attributes::T3
end

"Abstract type for simulation code"
abstract type SimCode end

"""
  Root data structure containing information required for code generation to
  generate simulation code for a Modelica model.
"""
struct SIM_CODE{T0<:String,
                T1<:AbstractDict{String, Tuple{Integer, SimVar}},
                T2<:Vector{BDAE.RESIDUAL_EQUATION},
                T4<:Vector{BDAE.WHEN_EQUATION},
                T5<:Vector{BDAE.IF_EQUATION},
                T6<:Bool,
                T7<:Vector{Int},
                T8<:LightGraphs.AbstractGraph,
                T9<:Vector{Int}} <: SimCode
  name::T0
  "Mapping of names to the corresponding variable"
  crefToSimVarHT::T1
  "Different equations stored within simulation code"
  residualEquations::T2
  "The Initial equations"
  initialEquations::T2
  "When equations"
  whenEquations::T4
  "If Equations"
  ifEquations::T5
  "True if the system that we are solving is singular"
  isSingular::T6
  "
   The match order:
   Result of assign array, e.g array(j) = equation_i
  "
  matchOrder::T7
    "
    The merged graph. E.g digraph constructed from matching info.
    The indicies are the same as above and they are shared.
    If the system is singular tearing is needed.
   "
  equationGraph::T8
  "
    The reverse topological sort of the equation-graph
  "
  stronglyConnectedComponents::T9
end


struct UNSORTED_SIM_CODE <: SimCode
  name::String
  "Mapping of names to the corresponding variable"
  crefToSimVarHT::Dict{String, Tuple{Integer, SimVar}}
  "Different equations stored within simulation code"
  residualEquations::Array{BDAE.RESIDUAL_EQUATION}
  whenEquations::Array{BDAE.WHEN_EQUATION}
  ifEquations::Array{BDAE.IF_EQUATION}
end
