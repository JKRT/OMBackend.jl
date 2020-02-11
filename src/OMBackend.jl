#=
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
*/ =#

#= Setup to support multiple modules by adding them to the load path=#
const CURRENT_DIRECTORY = @__DIR__
const BACKEND_DIRECTORY = CURRENT_DIRECTORY * "/Backend"
const CODE_GENERATION_DIRECTORY = CURRENT_DIRECTORY * "/CodeGeneration"
if ! (CURRENT_DIRECTORY in LOAD_PATH)
  @info("Setting up loadpath..")
  push!(LOAD_PATH, CURRENT_DIRECTORY, BACKEND_DIRECTORY, CODE_GENERATION_DIRECTORY)
  @info("Done setting up loadpath: $LOAD_PATH")
end

@info("initialize the backend API")

module OMBackend

using MetaModelica
using ExportAll
using Absyn

import BackendDAE
import BackendDAECreate
import BackendDump
import Causalize
import DAE
import Prefix
import SCode
import SimulationCode
import SimCodeDump
import CodeGeneration

function translate(frontendDAE::DAE.DAElist)
  local bDAE = lower(frontendDAE)
  local simCode = generateSimulationCode(bDAE)
  generateTargetCodeAndSimulate(simCode)

  return true
end

"""
 Transforms given Frontend DAE IR to causalized backend DAE IR (BDAE IR)
"""
function lower(frontendDAE::DAE.DAElist)
  local bDAE::BackendDAE.BackendDAEStructure
  local simCode::SIM_CODE
  #= Create Backend structure from Frontend structure =#
  bDAE = BackendDAECreate.lower(frontendDAE)
  BackendDump.dumpBackendDAEStructure(bDAE, "translated");
  #= detect state variables =#
  bDAE = Causalize.detectStates(bDAE)
  BackendDump.dumpBackendDAEStructure(bDAE, "states marked");
  #= causalize system, for now DAEMode =#
  bDAE = Causalize.daeMode(bDAE)
  BackendDump.dumpBackendDAEStructure(bDAE, "residuals");
  return bDAE
end

"""
  Transforms causalized BDAE IR to simulation code
"""
function generateSimulationCode(bDAE::BackendDAE.BackendDAEStructure)::SimulationCode.SIM_CODE
  simCode = CodeGeneration.transformToSimCode(bDAE)
  SimCodeDump.dumpSimCode(simCode, "transformed simcode")
  (simCode)
end


function generateTargetCodeAndSimulate(simCode::SimulationCode.SIM_CODE)
  #= Target code =#
  (fileName, functions)= CodeGeneration.generateCode(simCode)
  @info functions
  @info fileName
  CodeGeneration.writeDAE_equationsToFile(fileName, functions)
  #include(fileName)
end

end #=OMBackend=#
