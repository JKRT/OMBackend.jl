
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

const ELSE_BRANCH = -1

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

""" Abstract type for different variants of simulation code """
abstract type SimCode end

"""
  Abstract type for control flow constructs for simulation code
"""
abstract type Construct end

"""
  Represents a branch in simulation code. 
  Contains a single condition, a set of inner equations and a set of possible targets.
  Since each branch potentially contains a set of equations information exist so that the code in each branch can be 
  matched (Similar to the larger system )
"""
struct BRANCH{T1 <: DAE.Exp,
              T2 <: Vector{BDAE.RESIDUAL_EQUATION},
              T3 <: Int, #= Integer code Each branch has one target (next) The ID of one branch is target - 1=#
              T4 <: Bool,
              T5 <: Vector{Int},
              T6 <: Graphs.AbstractGraph,
              T7 <: Vector{Vector{Int}},
              T8 <: AbstractDict{String, Tuple{Integer, SimVar}}} <: Construct
  
  condition::T1
  residualEquations::T2
  identifier::T3 #= A value of -1 indicates that this branch is an else branch =#
  targets::T3
  isSingular::T4
  matchOrder::T5
  equationGraph::T6
  sccs::T7
  stringToSimVarHT::T8
end

"""
  A representation of a simcode IF Equation. 
  Similar to the main simcode module it contains information to make construct a causal representation easier.
"""
struct IF_EQUATION{Branches <: Vector{BRANCH}} <: Construct
  branches::Branches
end


"""
  Root data structure containing information required for code generation to
  generate simulation code for a Modelica model.
"""
struct SIM_CODE{T0<:String,
                T1<:AbstractDict{String, Tuple{Integer, SimVar}},
                T2<:Vector{BDAE.RESIDUAL_EQUATION},
                #= When equations are discret events (In the sense they occur once, the condition is checked per time step. ) =#
                T4<:Vector{BDAE.WHEN_EQUATION},
                #=
                  If equations are represented via a vector of possible branches in which the code can operate. 
                Similar to basic blocks
                =#
                T5<:Vector{IF_EQUATION},
                T6<:Bool,
                T7<:Vector{Int},
                T8<:Graphs.AbstractGraph,
                T9<:Vector,
                T10 <: Vector{BDAE.STRUCTURAL_TRANSISTION},
                T11 <: Vector,
                T12 <: String} <: SimCode
  name::T0
  "Mapping of names to the corresponding variable"
  stringToSimVarHT::T1
  "Different equations stored within simulation code"
  residualEquations::T2
  "The Initial equations"
  initialEquations::T2
  "When equations"
  whenEquations::T4
  "If Equations (Simulation code branches). Each branch contains a condition a set of residual equations and a set of targets"
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
  " The reverse topological sort of the equation-graph "
  stronglyConnectedComponents::T9
  "Contains all structural transistions"
  structuralTransistions::T10  
  "Structural submodels"
  subModels::T11
  "Initial model"
  activeModel::T12
end
