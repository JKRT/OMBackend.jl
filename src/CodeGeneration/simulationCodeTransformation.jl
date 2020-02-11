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
using Setfield

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

function bDAEVarKindToSimCodeVarKind(backendVar::BackendDAE.Var)::SimulationCode.SimVarType
  varKind = @match backendVar.varKind begin
    BackendDAE.STATE(__) => SimulationCode.STATE()
    BackendDAE.PARAM(__) || BackendDAE.CONST(__) => SimulationCode.PARAMETER()
    BackendDAE.VARIABLE(__) => SimulationCode.ALG_VARIABLE()
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
  local equations = let
      [eq for es in equationSystems for eq in es.orderedEqs]
  end
  @info equations
  local simulationEquations = allocateAndCollectSimulationEquations(equations)
  #= Construct SIM_CODE =#
  simCode = SimulationCode.SIM_CODE(backendDAE.name, crefToSimVarHT, simulationEquations)
  return simCode
end

"""
-Andreas:
  Assign the indices.
  Construct the HashTable.
  1. Collect all variables
2. Search all states (e.g. x and y) and give them indices starting at 1 (so x=1, y=2). Then give the corresponding state derivatives (x' and y') the same indices.
3. Remaining algebraic variables will get indices starting with i+1, where i is the number of states.
4. Parameters will get own set of indices, starting at 1.

"""
function createIndices(simulationVars::Array{SimulationCode.SIMVAR})::Dict{String, Tuple{Integer, SimulationCode.SimVarType}}
  @info simulationVars
  local ht::Dict{String, Tuple{Integer, SimulationCode.SimVarType}} = Dict()
  local stateCounter = 0
  local parameterCounter = 0
  local numberOfStates = 0
  for var in simulationVars
    @match var.varKind begin
      SimulationCode.STATE(__) => begin
        stateCounter += 1
        @set var.index = SOME(stateCounter)
        push!(ht, var.name => (stateCounter, var.varKind))
        #=Adding the state derivative as well=#
        push!(ht, "der($(var.name))" => (stateCounter, SimulationCode.STATE_DERIVATIVE(var.name)))
      end
      SimulationCode.PARAMETER(__) => begin
        parameterCounter += 1
        push!(ht, var.name => (parameterCounter, var.varKind))
      end
      _ => continue
    end
  end
  for var in simulationVars
    @match var.varKind begin
      SimulationCode.ALG_VARIABLE(__) => begin
        var = @set var.index = SOME(stateCounter + 1)
        push!(ht, var.name => (var.index.data, var.varKind))
      end
      _ => continue
    end
  end
  @info ht
  return ht
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
