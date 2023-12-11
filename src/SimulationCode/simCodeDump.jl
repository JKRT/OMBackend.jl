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
function dumpSimCode(simCode::SimulationCode.SIM_CODE, heading::String = "Simulation-Code")
  local buffer = IOBuffer()
  print(buffer, BDAEUtil.DOUBLE_LINE + "\n")
  print(buffer, "SIM_CODE: " + heading + "\n")
  print(buffer, BDAEUtil.DOUBLE_LINE + "\n\n")
  local stateVariables = String[]
  local parameters = String[]
  local algVariables = String[]
  local discreteVariables = String[]
  local stateDerivatives = String[]
  local occVariables = String[]
  local dsVariables = String[]
  for f in simCode.functions
    println(buffer, string(f))
  end
  print(buffer, "Simulation Code Variables:" + "\n")
  print(buffer, BDAEUtil.LINE + "\n")
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
      SimulationCode.OCC_VARIABLE(__) => push!(occVariables, varName)
      SimulationCode.DATA_STRUCTURE(__) => push!(dsVariables, varName)
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
  println(buffer, "OCC Variables:")
  for a in occVariables
    println(buffer, a * "| Index:" * string(first(simCode.stringToSimVarHT[a])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  println(buffer, "Discrete Variables:")
  for d in discreteVariables
    println(buffer, d * "| Index:" * string(first(simCode.stringToSimVarHT[d])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  println(buffer, "Datastructure Variables:")
  for d in dsVariables
    println(buffer, d * "| Index:" * string(first(simCode.stringToSimVarHT[d])))
  end
  print(buffer, BDAEUtil.LINE + "\n")
  print(buffer, "\n")
  if !isempty(simCode.initialEquations)
    print(buffer, "Initial Equations" + "\n")
    println(buffer, BDAEUtil.LINE)
    for ieq in simCode.initialEquations
      print(buffer, string(ieq))
    end
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
  if !isempty(simCode.ifEquations)
    println(buffer, "If-equations:")
    print(buffer, BDAEUtil.LINE + "\n")
    for ifEq in simCode.ifEquations
      print(buffer, string(ifEq))
    end
    println(buffer, BDAEUtil.LINE)
  end
  if !isempty(simCode.whenEquations)
    println(buffer, "When-Equations:")
    println(buffer, BDAEUtil.LINE)
    for wEq in simCode.whenEquations
      print(buffer, string(wEq))
    end
    println(buffer, BDAEUtil.LINE)
  end
  if !isempty(simCode.structuralTransitions)
    println(buffer, BDAEUtil.LINE)
    println(buffer, "Structural-Equations:")
    println(buffer, BDAEUtil.LINE)
    for st in simCode.structuralTransitions
      print(buffer, string(st))
    end
    println(buffer, BDAEUtil.LINE)
  end
  if !isempty(simCode.sharedEquations)
    println(buffer, "Shared Variables:")
    println(buffer, BDAEUtil.LINE)
    for sv in simCode.sharedVariables
      print(buffer, string(sv) * "\n")
    end
    println(buffer, BDAEUtil.LINE)
  end
  if !isempty(simCode.sharedEquations)
    println(buffer, "Shared Equations:")
    println(buffer, BDAEUtil.LINE)
    for se in simCode.sharedEquations
      print(buffer, string(se))
    end
    println(buffer, BDAEUtil.LINE)
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
    #if isSome(wEq.whenEquation.elsewhenPart)
    #  @match SOME(elseW) = wEq.whenEquation.elsewhenPart
    #  nWhenEquations += length(elseW.whenEquation.whenStmtLst)
    #end
  end
  println(buffer, BDAEUtil.LINE)
  println(buffer, "Simulation Code Statistics:")
  println(buffer, BDAEUtil.LINE)
  println(buffer, "Total Number of Variables:" * string(length(algAndState) + length(discreteVariables) + length(occVariables)))
  println(buffer, "\tNumber of State Variables:" * string(length(stateVariables)))
  println(buffer, "\tNumber of Algebraic Variables:" * string(length(algVariables)))
  println(buffer, "\tNumber of OCC Variables:" * string(length(occVariables)))
  println(buffer, "\tNumber of Discrete Variables:" * string(length(discreteVariables)))
  println(buffer, "Total Number of Equations:" * string(length(simCode.residualEquations) + nIfEqs + nWhenEquations))
  println(buffer, "\tNumber of Residual Equations:" * string(length(simCode.residualEquations)))
  println(buffer, "\tNumber of Equations in If-Equations:" * string(nIfEqs))
  println(buffer, "\tNumber of Equations in When-Equations:" * string(nWhenEquations))
  println(buffer, "Total Number of Functions:" * string(length(simCode.functions)))
  println(buffer, BDAEUtil.LINE)
  #local varsInEvent = 0
  #println(buffer, "Variables involved in events:")
  print(buffer, BDAEUtil.DOUBLE_LINE + "\n")
  println(buffer, "END SIM_CODE")
  print(buffer, BDAEUtil.DOUBLE_LINE + "\n")
  if !isempty(simCode.subModels)
    for (i, sm) in enumerate(simCode.subModels)
      println(buffer, "\n\n")
      println(buffer, dumpSimCode(sm, "Structural-Sub-model #" * string(i)))
      println(buffer, "\n\n")
    end
  end

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

function Base.string(ieq::SimulationCode.DYNAMIC_OVERCONSTRAINED_CONNECTOR_EQUATION)
  BDAE.string(ieq.structuralDOCC_equation)
end

function string(st::IMPLICIT_STRUCTURAL_TRANSISTION)
  string(st.structuralWhenEquation)
end

function Base.string(simStructChange::SimulationCode.EXPLICIT_STRUCTURAL_TRANSISTION)
  local structuralChange = simStructChange.structuralTransition
  local str = "STRUCTURAL_TRANSITION: "
  str = str * "FROM: <" * structuralChange.fromState * "> TO: <" * structuralChange.toState * "> WHEN: " * string(structuralChange.transistionCondition) * "\n"
  return str
end

function string(f::EXTERNAL_MODELICA_FUNCTION)
  local buffer = IOBuffer()
  println(buffer, "function EXTERNAL " * f.name)
  for arg in f.inputs
    println(buffer, " " * string(arg))
  end
  for arg in f.outputs
    println(buffer, " " * string(arg))
  end
  println(buffer, "calling externally defined function: " * f.libInfo)
  println(buffer, "end " * f.name)
  return String(take!(buffer))
end

function string(f::MODELICA_FUNCTION)
  local buffer = IOBuffer()
  println(buffer, "function " * f.name)
  for arg in f.inputs
    println(buffer, "input " * string(arg))
  end
  for arg in f.outputs
    println(buffer, "output " * string(arg))
  end
  for l in f.locals
    println(buffer, " local:" * string(l))
  end
  for s in f.statements
    println(buffer, " " * string(s))
  end
  println(buffer, "end " * f.name)
  return String(take!(buffer))
end
