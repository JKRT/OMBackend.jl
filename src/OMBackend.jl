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
const CURRENT_DIRECTORY = @__DIR__
const BACKEND_DIRECTORY = realpath(CURRENT_DIRECTORY * "/Backend")
const CODE_GENERATION_DIRECTORY = realpath(CURRENT_DIRECTORY * "/CodeGeneration")

const EXAMPLE_DAE_DIRECTORY = realpath(CURRENT_DIRECTORY * "./../test/ExampleDAE")

@info "Starting.."
@info LOAD_PATH
if ! (CURRENT_DIRECTORY in LOAD_PATH)
  @debug("Setting up loadpath..")
  push!(LOAD_PATH, CURRENT_DIRECTORY, BACKEND_DIRECTORY, CODE_GENERATION_DIRECTORY, EXAMPLE_DAE_DIRECTORY)
  @debug("Done setting up loadpath: $LOAD_PATH")
end


@info("initialize backend API")
@info "Our current loadpath: $LOAD_PATH"

module OMBackend

using MetaModelica
using ExportAll
using Absyn

import BDAE
import BDAECreate
import Causalize
import DAE
import Prefix
import SCode
import Base.Meta
@info "Test"
import SimulationCode
import CodeGeneration
import ExampleDAEs

global EXAMPLE_MODELS = Dict("HelloWorld" => ExampleDAEs.helloWorld_DAE)

function info()
  println("OMBackend.jl")
  println("A Julia backend for the Equation Oriented Language Modelica!")
  println("Run any test module by executing runExample(<model-name>)")
  println("Available Example models include:")
  for (k,v) in EXAMPLE_MODELS
    println(k)
  end
  println("Cheers // The Devs")
end

function runExample(modelName::String)
  println("Translating Hybrid DAE of: $modelName")
  translate(EXAMPLE_MODELS[modelName])
  println("Translation done:")
  println("Simulating with default settings:")
end

"""
Contains expressions of models in memory.
Do NOT mutate in other modules!
//John
"""
global COMPILED_MODELS = Dict()

function translate(frontendDAE::DAE.DAE_LIST)::Tuple{String, String}
  local bDAE = lower(frontendDAE)
  local simCode = generateSimulationCode(bDAE)
  generateTargetCode(simCode)
end

"""
 Transforms given Frontend DAE IR to causalized backend DAE IR (BDAE IR)
"""
function lower(frontendDAE::DAE.DAE_LIST)::BDAE.BDAEStructure
  local bDAE::BDAE.BDAEStructure
  local simCode::SIM_CODE
  @assert typeof(listHead(frontendDAE.elementLst)) == DAE.COMP
  #= Create Backend structure from Frontend structure =#
  bDAE = BDAECreate.lower(frontendDAE)
  @debug(BDAE.stringHeading1(bDAE, "translated"));
  bDAE = Causalize.detectIfEquations(bDAE)
  @debug(BDAE.stringHeading1(bDAE, "if equations transformed"));
  bDAE = Causalize.detectStates(bDAE)
  @debug(BDAE.stringHeading1(bDAE, "states marked"));
  bDAE = Causalize.residualizeEveryEquation(bDAE)
  @debug(BDAE.stringHeading1(bDAE, "residuals"));
  return bDAE
end

"""
  Transforms causalized BDAE IR to simulation code for DAE-mode
"""
function generateSimulationCode(bDAE::BDAE.BDAEStructure)::SimulationCode.SimCode
  simCode = CodeGeneration.transformToSimCode(bDAE)
  @debug BDAE.stringHeading1(simCode, "SIM_CODE: transformed simcode")
  return simCode
end

"""
  Generates code interfacing DifferentialEquations.jl
  The resulting code is saved in a dictonary which contains functions that where simulated
  this session. Returns the generated modelName and corresponding generated code
"""
function generateTargetCode(simCode::SimulationCode.SIM_CODE)
  #= Target code =#
  (modelName, modelCode) = CodeGeneration.generateCode(simCode)
  @debug "Functions:" modelCode
  @debug "Model:" modelName
  COMPILED_MODELS[modelName] = modelCode
  return (modelName, modelCode)
end

function writeModelToFile(modelName::String)
  CodeGeneration.writeDAE_equationsToFile(fileName, functions)
end

"""
  Evaluates the in memory representation of a named model
"""
function simulateModel(modelName::String, tspan=(0.0, 1.0))
  local modelCode = COMPILED_MODELS[modelName]
  @debug "Generated modelCode : $modelCode"
  local res = Meta.parse("begin $modelCode end") #Hack
  @debug res
  eval(res)
  eval(Meta.parse("$(modelName)Simulate($(tspan))"))
end

end #=OMBackend=#
