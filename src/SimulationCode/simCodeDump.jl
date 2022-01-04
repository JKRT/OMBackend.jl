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
using DataStructures
"""
  File: simCodeDump.jl
  Dumping functions for simulation code structures.
"""
function dumpSimCode(simCode::SimulationCode.SIM_CODE, heading::String)
  print(BDAE.DOUBLE_LINE + "\n")
  print("SIM_CODE: " + heading + "\n")
  print(BDAE.DOUBLE_LINE + "\n\n")

  print("SimCodeVars" + "\n")
  print(BDAE.LINE + "\n")
  dictPrettyPrint(simCode.stringToSimVarHT)
  print("\n")

  print("SimCode residual equations" + "\n")
  print(BDAE.LINE + "\n")
  for eq in simCode.residualEquations
    BDAE.printEqTraverse(eq, 0)
    print("\n")
  end
  print("\n")
end

function string(ht::OrderedDict{T1, Tuple{T2, SimVar}}) where {T1, T2}
  ks = keys(ht)
  res = ""
  for k in ks
    res *= "Variable name:" * k * " Index:" * string(first(ht[k])) * " Attributes:{" * string(last(ht[k])) * "}\n"
  end
  return res
end

function Base.string(simCode::SimulationCode.SIM_CODE)::String
  str = BDAEUtil.stringHeading3("", "SimCode Variables")
  str =  str * string(simCode.stringToSimVarHT)
  str = str + BDAEUtil.heading3("SimCode Equations")
  for eq in simCode.residualEquations
    str = str + string(eq)
  end
  return str + "\n"
end

"""
 Converts a ```backendVar::BDAE.Var``` to the simcode format.
"""
function string(backendVar::BDAE.Var)
  BDAE.string(varName; separator = "___")
end


"""
  This function just forwards the call to BDAE
"""
function string(element)
  BDAE.string(element)
end

function string(v::SIMVAR)
  return v.name  * "," * string(v.index) * ","  * string(v.varKind)
end
