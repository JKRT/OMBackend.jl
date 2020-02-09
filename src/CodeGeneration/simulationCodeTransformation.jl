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

import DAE

using MetaModelica
using BackendDAE
using BackendEquation
using SimulationCode


"""
  Collect variables from array of BackendDAE.Var:
  Save the name and it's kind of each variable.
  Index will be set to NONE.
"""
function collectVariables(allBackendVars::Array{BackendDAE.Var})
  local numberOfVars::Integer = length(allBackendVars)
  local simVars::Array = Array{SimulationCode.SimVar}(undef, numberOfVars)
  for (i, backendVar) in enumerate(allBackendVars)
    local simVarName::String = bDAEIdentToSimCodeVarName(backendVar)
    local simVarKind::SimulationCode.SimVarType = bDAEVarKindToSimCodeVarKind(backendVar)
    simVars[i] = SimulationCode.SIMVAR(simVarName, NONE(), simVarKind)
  end
  return simVars
end

function collectEquations(equations::Array)
end
using BackendDAE

function bDAEVarKindToSimCodeVarKind(backendVar::BackendDAE.Var)::SimulationCode.SimVarType
  varKind = @match backendVar.varKind begin
    BackendDAE.STATE(__) => SimulationCode.STATE()
    BackendDAE.PARAM(__) => SimulationCode.ALG_VARIABLE()
    BackendDAE.PARAM(__) || BackendDAE.CONST(__) => SimulationCode.PARAMETER()
    _ => @error("Kind $(typeof(backendVar.varKind)) of backend variable not handled.")
  end
end

function bDAEIdentToSimCodeVarName(backendVar::BackendDAE.Var)
  local varName = backendVar.varName
  @match varName begin
    DAE.CREF_IDENT(__) => string(backendVar.varName.ident)
    _ => @error("Type $(typeof(varName)) not handled.")
  end
end

"""
  Transform BackendDAEStructure to SimulationCode.SIM_CODE
"""
function transformToSimCode(backendDAE::BackendDAE.BACKEND_DAE)::SimulationCode.SIM_CODE
  local equationSystems::Array = backendDAE.eqs
  # Idea:
  # Collect variables and generate index
  # Have seperate indeces (starting with 1) for state derivatives, parameters and everything else
  # If x is a state we need two variables der(x) and x and both need the same index
  # dx[1] aka der(x) and x[1] aka x
  # Parameter a would then be in p[1]
  # Collect variables and sort for states
  local allOrderedVars::Array{BackendDAE.Var} =
    let
      [v for es in equationSystems for v in es.orderedVars.varArr]
    end

  local allSharedVars::Array{BackendDAE.Var} = getSharedVariablesLocalsAndGlobals(backendDAE.shared) #=TODO: One equation sys for now=#

  local allBackendVars = vcat(allOrderedVars, allSharedVars)
  #==#
  local simVars::Array{SimulationCode.SIMVAR} =
    allocateAndCollectSimulationVariables(allBackendVars)
  # Assign indices and put all variable into an hash table
  local crefToSimVarHT = createIndices(simVars)

  #= kabdelhak: is this efficient or inefficient? i have no idea but it works for now =#
  local equations = BackendEquation.concenateEquations(backendDAE.eqs)

  local simulationEquations = allocateAndCollectSimulationEquations(equations)
  #= Construct SIM_CODE =#
  simCode = SimulationCode.SIM_CODE(crefToSimVarHT, equations)
  return simCode
end

"""
  Assign the indices.
  Construct the HashTable.
"""
function createIndices(simulationVars::Array{SimulationCode.SIMVAR})::Dict{String, SimulationCode.SimVar}
  local crefToSimVarHT = Dict{String, SimulationCode.SimVar}()
  for var in simulationVars
    crefToSimVarHT[var.name] = var
  end
  (crefToSimVarHT)
end

"""
  Does nothing for now
"""
function allocateAndCollectSimulationEquations(equations)
 equations
end

"""
Returns the shared global and local variable for the shared data in
an equation system. If no such data is present. Return two empty arrays
"""
function getSharedVariablesLocalsAndGlobals(shared::Shared)
  @match shared begin
    BackendDAE.SHARED(__) => vcat(shared.globalKnownVars, shared.localKnownVars)
    _ => []
  end
end

function allocateAndCollectSimulationVariables(bDAEVariables::Array{BackendDAE.Var})
  collectVariables(bDAEVariables)
end
