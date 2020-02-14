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
  File: SimCodeDump.jl
  Dumping functions for simulation code structures.
"""

using BackendDAE
using SimulationCode

function dumpSimCode(simCode::SimulationCode.SIM_CODE, heading::String)
  print(BackendDAE.DOUBLE_LINE + "\n")
  print("SIM_CODE: " + heading + "\n")
  print(BackendDAE.DOUBLE_LINE + "\n\n")

  print("SimCodeVars" + "\n")
  print(BackendDAE.LINE + "\n")
  BackendDAE.dictPrettyPrint(simCode.crefToSimVarHT)
  print("\n")

  print("SimCodeEquations" + "\n")
  print(BackendDAE.LINE + "\n")
  for eq in simCode.equations
    BackendDAE.printEqTraverse(eq, 0)
    print("\n")
  end
  print("\n")
end

function Base.string(simCode::SimulationCode.SIM_CODE)::String
  str = stringHeading3(simCode.crefToSimVarHT, "SimCodeVars")
  str = str + heading3("SimCodeEquations")
  for eq in simCode.equations
    str = str + string(eq)
  end
  return str + "\n"
end
