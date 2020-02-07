#=
* This file is part of OpenModelica.
*
* Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
* c/o Link�pings universitet, Department of Computer and Information Science,
* SE-58183 Link�ping, Sweden.
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
const EXAMPLE_DAE = CURRENT_DIRECTORY * "/ExampleDAE"

if ! (CURRENT_DIRECTORY in LOAD_PATH)
  @info("Setting up loadpath..")
  push!(LOAD_PATH, CURRENT_DIRECTORY, BACKEND_DIRECTORY, EXAMPLE_DAE)
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
import ExampleDAEs

function translate()
  local lst::DAE.DAElist = ExampleDAEs.HelloWorld_DAE
  execute_translation_steps(lst)
end

function translate(lst_DAE_IR::DAE.DAElist)
  execute_translation_steps(lst_DAE_IR)
end

function execute_translation_steps(lst::DAE.DAElist)
  local bDAE::BackendDAE.BackendDAEStructure
  #= Create Backend structure from Frontend structure =#
  dae = BackendDAECreate.lower(lst)
  BackendDump.dumpBackendDAEStructure(dae, "translated");

  #= detect state variables =#
  dae = Causalize.detectStates(dae)
  BackendDump.dumpBackendDAEStructure(dae, "states marked");

  #= causalize system, for now DAEMode =#
  dae = Causalize.daeMode(dae)
  BackendDump.dumpBackendDAEStructure(dae, "residuals");

  #= create simCode -> target code =#
end

end #=OMBackend=#
