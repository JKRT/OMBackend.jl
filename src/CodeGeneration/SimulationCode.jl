#= /*
* This file is part of OpenModelica.
*
* Copyright (c) 1998-2020, Open Source Modelica Consortium (OSMC),
* c/o Linköpings universitet, Department of Computer and Information Science,
* SE-58183 Linköping, Sweden.
*
* All rights reserved.
*
* THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
* THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
* ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
* RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
* ACCORDING TO RECIPIENTS CHOICE.
*
* The OpenModelica software and the Open Source Modelica
* Consortium (OSMC) Public License (OSMC-PL) are obtained
* from OSMC, either from the above address,
* from the URLs: http:www.ida.liu.se/projects/OpenModelica or
* http:www.openmodelica.org, and in the OpenModelica distribution.
* GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
*
* This program is distributed WITHOUT ANY WARRANTY; without
* even the implied warranty of  MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
* IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
*
* See the full OSMC Public License conditions for more details.
*
=#

"""
  File: SimulationCode.jl
  Data strucutres used for simulation code.
"""
module SimulationCode

import DAE
import BDAE
using MetaModelica
using DataStructures

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
  "Mapping of names to the corresponding variable.
   Each variable has a unique index.
   The purpose of this mapping is to have aliases for each unique index.
   We also need to keep track of each simulation variable.
  "
  nameToVar::OrderedDict{String, Tuple{Integer, SimVar}}
  indexToEquation::OrderedDict{Int, BDAE.RESIDUAL_EQUATION}
  "Equation <-> Variable graph (Bidirectional)"
  variableEqMapping::OrderedDict{String, Int}
  "Regular equations are encoded as residuals"
  residualEquations::Array{BDAE.RESIDUAL_EQUATION}
  whenEquations::Array{BDAE.WHEN_EQUATION}
  ifEquations::Array{BDAE.IF_EQUATION}
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
  sortedGraph::OrderedDict{Int,Array{Int}}
  "
    The strongly connected components of the sorted graph.
    This information can be used for tearing.
  "
  stronglyConnectedComponents::Array
end

#=
We include directly, since Julia does not allow circular imports
=#
include("simCodeDump.jl")

end # module SimulationCode
