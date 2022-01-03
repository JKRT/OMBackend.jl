#=
# This file is part of OpenModelica.
#
# Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
# c/o Linköpings universitet, Department of Computer and Information Science,
# SE-58183 Linköping, Sweden.
#
# All rights reserved.
#
# THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
# THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
# ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
# RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
# ACCORDING TO RECIPIENTS CHOICE.
#
# The OpenModelica software and the Open Source Modelica
# Consortium (OSMC) Public License (OSMC-PL) are obtained
# from OSMC, either from the above address,
# from the URLs: http:www.ida.liu.se/projects/OpenModelica or
# http:www.openmodelica.org, and in the OpenModelica distribution.
# GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
#
# This program is distributed WITHOUT ANY WARRANTY; without
# even the implied warranty of  MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
# IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
#
# See the full OSMC Public License conditions for more details.
#
  Author: John Tinnerholm, john.tinnerholm@liu.se
=#

"""
  Creates the structural callbacks 
"""
function createStructuralCallbacks(simCode, structuralTransistions::Vector{ST}) where {ST}
  local structuralCallbacks = Expr[]
  for structuralTransisiton in structuralTransistions
    push!(structuralCallbacks, createStructuralCallback(simCode, structuralTransisiton))
  end
  return structuralCallbacks
end

"""
  Creates a single structural callback.
"""
function createStructuralCallback(simCode, structuralTransistion::ST) where {ST}
  local cond = transformToZeroCrossingCondition(structuralTransistion.transistionCondition)
  local callbackName = createCallbackName(structuralTransistion)
  quote
    function $(Symbol(callbackName))(destinationSystem)
      #= Represent structural change. =#
      local structuralChange = OMBackend.Runtime.StructuralChange($(structuralTransistion.toState), false, destinationSystem)
      #= The affect simply activates the structural callback informing us to generate code for a new system =#
      function affect!(integrator)
        structuralChange.structureChanged = true
      end
      function condition(u, t, integrator)
        return $(expToJuliaExpMTK(cond, simCode))
      end    
      local cb = ContinuousCallback(condition, affect!)
      return (cb, structuralChange)
    end
  end
end

function createStructuralAssignments(simCode, structuralTransistions::Vector{ST}) where {ST}
  local structuralAssignments = Expr[]
    for structuralTransisiton in structuralTransistions
      push!(structuralAssignments, createStructuralAssignment(simCode, structuralTransisiton))
    end
  result = quote
    structuralCallbacks = OMBackend.Runtime.StructuralChange[]
    callbackSet = []
    $(structuralAssignments...)
  end
  return result
end

"""
  This function creates a structural assignment.
  That is the constructor for a structural callback guiding structural change.
"""
function createStructuralAssignment(simCode, structuralTransistion::ST) where {ST}
  local callbackName = createCallbackName(structuralTransistion)
  local toState = structuralTransistion.toState
  local toStateProblem = Symbol(toState * "Problem")
  local toStateModel = Symbol(toState * "Model")
  local integratorCallbackName = structuralTransistion.fromState * structuralTransistion.toState * "_CALLBACK"
  local structuralChangeStructure = structuralTransistion.fromState * structuralTransistion.toState * "_STRUCTURAL_CHANGE"
  quote
    ($(toStateProblem), _, _, _, _, _) = ($(toStateModel))(tspan)
    ($(Symbol(integratorCallbackName)), $(Symbol(structuralChangeStructure))) = $(Symbol(callbackName))($(toStateProblem))
    push!(structuralCallbacks, $(Symbol(structuralChangeStructure)))
    push!(callbackSet, ($(Symbol(integratorCallbackName))))
  end
end

function createCallbackName(structuralTransisiton::ST) where {ST}
  return "structuralCallback" * structuralTransisiton.fromState * structuralTransisiton.toState
end                                                             


"""
  Creates the variables that are shared between the structural modes of a model.
"""
function createCommonVariables(commonVariables)
  commonVariablesExpr = Expr[]
  for variable in commonVariables
    res = :(
      push!(commonVariables, $variable)
    )
    push!(commonVariablesExpr, res)
  end
  return quote
    commonVariables = String[]
    $(commonVariablesExpr...)
  end
end
