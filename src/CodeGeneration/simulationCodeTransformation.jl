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

using MetaModelica
using SimulationCode
using BackendDAE

"""
Collect variables from array of BackendDAE.Var:
Save the name and it's kind of each variable.
Index will be set to -2.
"""
function collectVariables(allBackendVars::Array{BackendDAE.Var})

  numberOfVars = length(allBackendVars)
  simVars = Array{SimulationCode.SimVars}(undef, numberOfVars)

  for (i,backendVar) in enumerate(allBackendVars)
    simVarName = ""
    if typeof(backendVar.varName) == DAE.CREF_IDENT
      simVarName = string(backendVar.varName.ident)
    else
      @error("Type $(typeof(backendvar)) not handled.")
    end

    if typeof(backendVar.varKind) == BackendDAE.STATE
      simVarKind = SimulationCode.STATE
    elseif typeof(backendVar.varKind) == BackendDAE.PARAM ||
           typeof(backendVar.varKind) == BackendDAE.CONST
      simVarKind = SimulationCode.PARAMETER
    elseif typeof(backendVar.varKind) == BackendDAE.PARAM
      simVarKind = SimulationCode.ALG_VARIABLE
    else
      @error("Kind $(typeof(backendVar.varKind)) of backend variable not handled.")
    end

    simVar = SimulationCode.SIMVAR(name=simVarName, index=-2, varKind=simVarKind)
    simVars[i] = simVar
  end

  return simVars
end


"""
Transform BackendDAEStructure to SimulationCode.SIM_CODE
"""
function transformToSimCode(backendDAE::BackendDAE.BackendDAEStructure)::SimulationCode.SIM_CODE

  # Idea:
  # Collect variables and generate index
  # Have seperate indeces (starting with 1) for state derivatives, parameters and everything else
  # If x is a state we need two variables der(x) and x and both need the same index
  # dx[1] aka der(x) and x[1] aka x
  # Parameter a would then be in p[1]

  # Collect variables and sort for states
  allBackendVars = Array{BackendDAE.Var}(backendDAE.eqs.orderdVars, backendDAE.shared.globalKnownVars #=, backendDAE.shared.locaKnownVars=#)
  simVars = Array{SimulationCode.SIMVAR}(undef, numberOfVars)
  simVars = collectVariables(allBackendVars)

  # Assign indices and put all variable into an hash table
  crefToSimVarHT = createIndices(simVars)

  # Collect equations
  equations = collectEquations(backendDAE)

  simCode = SimulationCode.SIM_CODE(equations, crefToSimVarHT)

  return simCode
end
