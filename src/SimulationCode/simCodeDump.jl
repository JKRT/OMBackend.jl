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
  print(BDAEUtil.DOUBLE_LINE + "\n")
  print("SIM_CODE: " + heading + "\n")
  print(BDAEUtil.DOUBLE_LINE + "\n\n")

  print("Simulation Code Variables:" + "\n")
  print(BDAEUtil.LINE + "\n")
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
  println("Parameters & Constants:")
  for p in parameters
    println(p * "| Index:" * string(first(simCode.stringToSimVarHT[p])))
  end
  print(BDAEUtil.LINE + "\n")
  println("State Variables:")
  for s in stateVariables
    println(s * "| Index:" * string(first(simCode.stringToSimVarHT[s])))
  end
  print(BDAEUtil.LINE + "\n")
  println("State Derivatives:")
  for sd in stateDerivatives
    println(sd * "| Index:" * string(first(simCode.stringToSimVarHT[sd])))
  end
  print(BDAEUtil.LINE + "\n")
  println("Algebraic Variables:")
  for a in algVariables
    println(a * "| Index:" * string(first(simCode.stringToSimVarHT[a])))
  end
  print(BDAEUtil.LINE + "\n")
  println("Discrete Variables:")
  for d in discreteVariables
    println(d * "| Index:" * string(first(simCode.stringToSimVarHT[d])))
  end
  print(BDAEUtil.LINE + "\n")
  print("\n")
  print("Initial Equations" + "\n")
  println(BDAEUtil.LINE)
  for ieq in simCode.initialEquations
    print(string(ieq))
  end
  println(BDAEUtil.LINE)
  print("Residual Equations" + "\n")
  print(BDAEUtil.LINE + "\n")
  i::Int = 0
  for eq in simCode.residualEquations
    i += 1
    print("Index:" * string(i) * "|" * BDAE.string(eq))
  end
  println(BDAEUtil.LINE)
  println("If-equations")
  print(BDAEUtil.LINE + "\n")
  for ifEq in simCode.ifEquations
    print(string(ifEq))
  end
  println(BDAEUtil.LINE)

  println("When-Equations")
  println(BDAEUtil.LINE)
  for wEq in simCode.whenEquations
    print(string(wEq))
  end
  println(BDAEUtil.LINE)
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
  println(BDAEUtil.LINE)
  println("Simulation Code Statistics:")
  println(BDAEUtil.LINE)
  println("Total Number of Variables:" * string(length(algAndState) + length(discreteVariables)))
  println("\tNumber of Algebraic Variables:" * string(length(algVariables)))
  println("\tNumber of State Variables:" * string(length(stateVariables)))
  println("\tNumber of Discrete Variables:" * string(length(discreteVariables)))

  println("Total Number of Equations:" * string(length(simCode.residualEquations) + nIfEqs + nWhenEquations))
  println("\tNumber of Residual Equations:" * string(length(simCode.residualEquations)))
  println("\tNumber of Equations in If-Equations:" * string(nIfEqs))
  println("\tNumber of Equations in When-Equations:" * string(nWhenEquations))
  println(BDAEUtil.LINE)
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
