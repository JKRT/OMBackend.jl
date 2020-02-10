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
import BackendDAE
using MetaModelica

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
struct PARAMETER <: SimVarType end


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
end


"""
  Root data structure containing information required for code generation to
  generate simulation code for a Modelica model.
"""
struct SIM_CODE
  "Mapping of names to the corresponding variable"
  crefToSimVarHT::Dict{String, Tuple{Integer, SimVarType}} # Index and type
 # "Mapping of index to corresponding variables"
 # indexToVarHT::Dict{Integer, Tuple{String, SimVarType}}
  "Array of Equations"
  equations::Array
end

end # module SimulationCode
