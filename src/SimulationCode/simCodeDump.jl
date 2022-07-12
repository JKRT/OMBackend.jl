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


#=
  TODO: Change to use buffered print instead.
=#

using DataStructures
"""
  File: simCodeDump.jl
  Dumping functions for simulation code structures.
"""
function dumpSimCode(simCode::SimulationCode.SIM_CODE, heading::String = "Simulation-Code")
  local buffer = IOBuffer()
  print(buffer, BDAEUtil.DOUBLE_LINE + "\n")
  print(buffer, "SIM_CODE: " + heading + "\n")
  print(buffer, BDAEUtil.DOUBLE_LINE + "\n\n")

  print(buffer, "Simulation Code Variables:" + "\n")
  print(buffer, BDAEUtil.LINE + "\n")
  local stateVariables = []
  local parameters = []
  local algVariables = []
  local discreteVariables = []
  local stateDerivatives = []
  for varName in keys(simCode.stringToSimVarHT)
    (idx, var) = simCode.stringToSimVarHT[varName]
    local varType = var.varKind
    local varType = var.varKind
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algVariables, varName)
      SimulationCode.DISCRETE(__) => push!(discreteVariables, varName)
    end
  end
  local algAndState = vcat(algVariables, stateVariables)
  println(buffer, "Parameters & Constants:")
  for p in parameters
    println(buffer, p * "| Index:" * string(first(simCode.stringToSimVarHT[p])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  println(buffer, "State Variables:")
  for s in stateVariables
    println(buffer, s * "| Index:" * string(first(simCode.stringToSimVarHT[s])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  println("State Derivatives:")
  for sd in stateDerivatives
    println(buffer, sd * "| Index:" * string(first(simCode.stringToSimVarHT[sd])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  println(buffer, "Algebraic Variables:")
  for a in algVariables
    println(buffer, a * "| Index:" * string(first(simCode.stringToSimVarHT[a])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  println(buffer, "Discrete Variables:")
  for d in discreteVariables
    println(buffer, d * "| Index:" * string(first(simCode.stringToSimVarHT[d])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  print(buffer, "\n")
  print(buffer, "Initial Equations" + "\n")
  println(buffer, BDAEUtil.LINE)
  for ieq in simCode.initialEquations
    print(buffer, string(ieq))
  end
  println(buffer, BDAEUtil.LINE)
  print(buffer, "Residual Equations" + "\n")
  print(buffer, BDAEUtil.LINE + "\n")
  i::Int = 0
  for eq in simCode.residualEquations
    i += 1
    print(buffer, "Index:" * string(i) * "|" * BDAE.string(eq))
  end
  println(buffer, BDAEUtil.LINE)
  println(buffer, "If-equations")
  print(buffer, BDAEUtil.LINE + "\n")
  for ifEq in simCode.ifEquations
    print(buffer, string(ifEq))
  end
  println(buffer, BDAEUtil.LINE)

  println(buffer, "When-Equations")
  println(buffer, BDAEUtil.LINE)
  for wEq in simCode.whenEquations
    print(buffer, string(wEq))
  end
  println(buffer, BDAEUtil.LINE)
  nIfEqs = 0
  for ifEq in simCode.ifEquations
    #Required to be balanced
    nIfEqs += length(first(ifEq.branches).residualEquations)
  end
  nWhenEquations = 0
  for wEq in simCode.whenEquations
    nWhenEquations += length(wEq.whenEquation.whenStmtLst)
    if isSome(wEq.whenEquation.elsewhenPart)
      @match SOME(elseW) = wEq.whenEquation.elsewhenPart
      nWhenEquations += length(elseW.whenEquation.whenStmtLst)
    end
  end
  println(buffer, BDAEUtil.LINE)
  println(buffer, "Simulation Code Statistics:")
  println(buffer, BDAEUtil.LINE)
  println(buffer, "Total Number of Variables:" * string(length(algAndState) + length(discreteVariables)))
  println(buffer, "\tNumber of Algebraic Variables:" * string(length(algVariables)))
  println(buffer, "\tNumber of State Variables:" * string(length(stateVariables)))
  println(buffer, "\tNumber of Discrete Variables:" * string(length(discreteVariables)))

  println(buffer, "Total Number of Equations:" * string(length(simCode.residualEquations) + nIfEqs + nWhenEquations))
  println(buffer, "\tNumber of Residual Equations:" * string(length(simCode.residualEquations)))
  println(buffer, "\tNumber of Equations in If-Equations:" * string(nIfEqs))
  println(buffer, "\tNumber of Equations in When-Equations:" * string(nWhenEquations))
  println(buffer, BDAEUtil.LINE)
  return String(take!(buffer))
end

function string(ht::OrderedDict{T1, Tuple{T2, SimVar}}) where {T1, T2}
  ks = keys(ht)
  res = ""
  for k in ks
    res *= "Name:" * k * "|Index:" * string(first(ht[k])) * "|Attributes:{" * string(last(ht[k])) * "}|\n"
  end
  return res
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

function string(ifEq::IF_EQUATION)
  res = ""
  for branch in ifEq.branches
    res *= "IF (" * string(branch.condition) * ")\n"
    for eq in branch.residualEquations
      res *= string(eq)
    end
    res *= "END\n"
  end
  return res
end
