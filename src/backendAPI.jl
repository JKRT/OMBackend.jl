using MetaModelica
using ExportAll
using Absyn

import .Backend.BDAE
import .Backend.BDAECreate
import .Backend.Causalize
import ..CodeGeneration
import OMBackend.ExampleDAEs
import .FrontendUtil.DAE
import .FrontendUtil.Prefix
import .SimulationCode
import Base.Meta
import SCode
@info "Test"
global EXAMPLE_MODELS = Dict("HelloWorld" => OMBackend.ExampleDAEs.helloWorld_DAE)

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
  simulateModel(modelName)
  println("Model now in memory")
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
  simCode = SimulationCode.transformToSimCode(bDAE)
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
