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
  local idx = 1
  for structuralTransisiton in structuralTransistions
    push!(structuralCallbacks, createStructuralCallback(simCode, structuralTransisiton, idx))
    idx += 1
  end
  return structuralCallbacks
end

"""
  Creates a single structural callback for an explicit transition.
"""
function createStructuralCallback(simCode, simCodeStructuralTransistion::SimulationCode.EXPLICIT_STRUCTURAL_TRANSISTION, idx)
  local structuralTransistion = simCodeStructuralTransistion.structuralTransistion
  local cond = transformToZeroCrossingCondition(structuralTransistion.transistionCondition)
  local callbackName = createCallbackName(structuralTransistion, 0)
  quote
    function $(Symbol(callbackName))(destinationSystem)
      #= Represent structural change. =#
      local structuralChange = OMBackend.Runtime.StructuralChange($(structuralTransistion.toState), false, destinationSystem)
      #= The affect simply activates the structural callback informing us to generate code for a new system =#
      function affect!(integrator)
        structuralChange.structureChanged = true
      end
      function condition(x, t, integrator)
        return $(expToJuliaExp(cond, simCode))
      end    
      local cb = ContinuousCallback(condition, affect!)
      return (cb, structuralChange)
    end
  end
end

"""
  Creates an implicit structural callback where the final state is unknown.
  These structural callback can only occur as a part of a when equation.
  The when equation might also make other changes to the variables before recompilation.
TODO:
Also make sure to create possible other elements in the structural when equation
"""
function createStructuralCallback(simCode, simCodeStructuralTransistion::SimulationCode.IMPLICIT_STRUCTURAL_TRANSISTION, idx)
  local structuralTransistion = simCodeStructuralTransistion.structuralWhenEquation
  local callbackName = createCallbackName(structuralTransistion, idx)
  local whenCondition = structuralTransistion.whenEquation.condition
  local zeroCrossingCond = transformToZeroCrossingCondition(whenCondition)
  #= 
    Construct the specified structural change as a pair. 
    For now assume that we only change some parameter.
  =#
  local stmtLst = structuralTransistion.whenEquation.whenStmtLst
  (whenOperators, recompilationDirective) = createStructuralWhenStatements(stmtLst, simCode)
  local componentToModify = string(recompilationDirective.componentToChange)
  local newValue = expToJuliaExpMTK(recompilationDirective.newValue, simCode)
  modification = quote
    ($(componentToModify), $(newValue))
  end
  @match SOME(metaModel) = simCode.metaModel
  quote
    function $(Symbol(callbackName))()
      #= The recompilation directive. =#
      local modification = $(modification)
      #= Represent structural change. =#
      #=NOTE: For the implicit callbacks the model is assumed to be the same=#
      local structuralChange = OMBackend.Runtime.StructuralChangeRecompilation($(simCode.name),
                                                                               false,
                                                                               $(metaModel),
                                                                               modification,
                                                                               $(simCode.stringToSimVarHT))
      #= The affect simply activates the structural callback informing us to generate code for a new system =#
      function affect!(integrator)
        #= Expand the when operators =#
        $(whenOperators...)
        #= Change the struct to mark that a structural change have occured=#
        structuralChange.structureChanged = true
      end
      function condition(x, t, integrator)
        str = $callbackName
        return $(expToJuliaExp(zeroCrossingCond, simCode))
      end    
      local cb = ContinuousCallback(condition, affect!)
      return (cb, structuralChange)
    end
  end
end

"""
"""
function createStructuralAssignments(simCode, structuralTransistions::Vector{ST}) where {ST}
  local structuralAssignments = Expr[]
  local idx = 1
  for structuralTransisiton in structuralTransistions
    @match structuralTransisiton begin
      SimulationCode.EXPLICIT_STRUCTURAL_TRANSISTION(__) => begin
        push!(structuralAssignments, createStructuralAssignment(simCode, structuralTransisiton))
      end
      SimulationCode.IMPLICIT_STRUCTURAL_TRANSISTION(__) => begin
        push!(structuralAssignments, createStructuralAssignment(simCode, structuralTransisiton, idx))
      end
    end
    idx += 1
  end
  result = quote
    structuralCallbacks = OMBackend.Runtime.AbstractStructuralChange[]
    callbackSet = []
    $(structuralAssignments...)
  end
  return result
end

"""
  This function creates a structural assignment.
  That is the constructor for a structural callback guiding structural change.
"""
function createStructuralAssignment(simCode, simCodeStructuralTransistion::SimulationCode.EXPLICIT_STRUCTURAL_TRANSISTION)
  local structuralTransistion = simCodeStructuralTransistion.structuralTransistion
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

"""
Creates a structural assignment for an implicit structural transisiton.
These are numbered from 1->N
"""
function createStructuralAssignment(simCode, simCodeStructuralTransistion::SimulationCode.IMPLICIT_STRUCTURAL_TRANSISTION, idx::Int)
  local structuralTransistion = simCodeStructuralTransistion.structuralWhenEquation
  local callbackName = createCallbackName(structuralTransistion, idx)
  local integratorCallbackName = callbackName * "_CALLBACK"
  local structuralChangeStructure = callbackName * "_STRUCTURAL_CHANGE"
  quote
    ($(Symbol(integratorCallbackName)), $(Symbol(structuralChangeStructure))) = $(Symbol(callbackName))()
    push!(structuralCallbacks, $(Symbol(structuralChangeStructure)))
    push!(callbackSet, ($(Symbol(integratorCallbackName))))
  end
end

function createCallbackName(structuralTransisiton::BDAE.STRUCTURAL_TRANSISTION, idx = 0)
  return "structuralCallback" * structuralTransisiton.fromState * structuralTransisiton.toState
end

"""
  Creates a structural callback for the when equation.
  The name is up to change.
"""
function createCallbackName(structuralTransisiton::BDAE.STRUCTURAL_WHEN_EQUATION, idx::Int)
  return string("structuralCallbackWhenEquation", idx)
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

"""
  Generates statements for the structural when equation construct.
  This function returns a tuple where the first part is a vector of statements occuring in the when equation.
  The second part is the recompilation statement itself that specifies what structural changes are to occur.
  This last part is then used by the runtime to just-in-time compile the model when the event occurs.
"""
function createStructuralWhenStatements(@nospecialize(whenStatements::List{BDAE.WhenOperator}),
                                        simCode::SimulationCode.SIM_CODE)
  local res::Vector{Expr} = []
  local recompilationOperator
  for wStmt in  whenStatements
    @match wStmt begin
      BDAE.ASSIGN(__) => begin
        exp1 = expToJuliaExp(wStmt.left, simCode, varPrefix="integrator.u")
        exp2 = expToJuliaExp(wStmt.right, simCode)
        push!(res, :($(exp1) = $(exp2)))
      end
      BDAE.RECOMPILATION(__) => begin
        recompilationOperator = wStmt
      end
      BDAE.REINIT(__) => begin
        throw("Reinit is not allowed in a structural when equation")
      end
      _ => begin
        throw("Unsupported statement in the when equation:" * string(wStmt))
      end
    end
  end
  return (res, recompilationOperator)
end
