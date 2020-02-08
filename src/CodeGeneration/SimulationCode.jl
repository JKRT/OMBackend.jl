#= /*
* This file is part of OpenModelica.
*
* Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
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
  Some description about this module
"""
module SimulationCode

import DAE
import BackendDAE

using MetaModelica
using ExportAll

"""
Kind of a simulation variable
"""
abstract type SimVarType end

"""
State variable
"""
abstract type STATE <: SimVarType end

"""
Algebraic variable
"""
abstract type  ALG_VARIABLE <: SimVarType end

"""
Input variable
"""
abstract type INPUT <: SimVarType end

"""
Parameter variable
"""
abstract type PARAMETER <: SimVarType end


"""
Some description 2
"""
abstract type SimVar end

"""
Additional information about a variable used in code generation.
"""
struct SIMVAR <: SimVar
  name :: String # TODO: Replace with DAE.ComponentRef
  "Index of variable, 0 based, type based"
  index::Integer
  varKind #= Kind of variable (State, Differentiated State,  Algebraic Variable, Parameter )=# ::SimVarType
end


"""
Root data structure containing information required for code generation to
generate simulation code for a Modelica model.
"""
struct SIM_CODE
  #crefToSimVarHT::Dict{()}
end


@exportAll()

end # module SimulationCode
