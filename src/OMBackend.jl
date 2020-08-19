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

#= Setup to support multiple modules by adding them to the load path =#

# const BACKEND_DIRECTORY = realpath(CURRENT_DIRECTORY * "/Backend")
# const CODE_GENERATION_DIRECTORY = realpath(CURRENT_DIRECTORY * "/CodeGeneration")
# const EXAMPLE_DAE_DIRECTORY = realpath(CURRENT_DIRECTORY * "./../test/ExampleDAE")

# @info "Starting.."
# @info LOAD_PATH
# if ! (CURRENT_DIRECTORY in LOAD_PATH)
#   @debug("Setting up loadpath..")
#   push!(LOAD_PATH, CURRENT_DIRECTORY, BACKEND_DIRECTORY, CODE_GENERATION_DIRECTORY, EXAMPLE_DAE_DIRECTORY)
#   @debug("Done setting up loadpath: $LOAD_PATH")
# end


# @info("initialize backend API")
# @info "Our current loadpath: $LOAD_PATH"

module OMBackend

const CURRENT_DIRECTORY = @__DIR__
include("$CURRENT_DIRECTORY/FrontendUtil/FrontendUtil.jl")
include("$CURRENT_DIRECTORY/Backend/Backend.jl")
include("$CURRENT_DIRECTORY/SimulationCode/SimulationCode.jl")
include("$CURRENT_DIRECTORY/CodeGeneration/CodeGeneration.jl")
include("$CURRENT_DIRECTORY/../test/ExampleDAE/ExampleDAEs.jl")
include("backendAPI.jl")

end #=OMBackend=#
