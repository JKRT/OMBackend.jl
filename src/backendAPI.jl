using MetaModelica
using ExportAll
using Absyn
#= For interactive evaluation. =#
using ModelingToolkit
#import Symbolics
using SymbolicUtils


import .Backend.BDAE
import .Backend.BDAECreate
import .Backend.BDAEUtil
import .Backend.Causalize
import ..CodeGeneration
import OMBackend.ExampleDAEs
import .SimulationCode
import Base.Meta
import SCode
import JuliaFormatter
import Plots

global EXAMPLE_MODELS = Dict("HelloWorld" => OMBackend.ExampleDAEs.helloWorld_DAE
                             , "LotkaVolterra" => OMBackend.ExampleDAEs.lotkaVolterra_DAE
                             , "VanDerPol" => OMBackend.ExampleDAEs.vanDerPol_DAE
                             , "BouncingBall" => OMBackend.ExampleDAEs.bouncingBall_DAE
                             , "Influenza" => OMBackend.ExampleDAEs.influenza_DAE)
#= Settings =#
@enum BackendMode begin
  DAE_MODE = 1
  ODE_MODE = 2
  MODELING_TOOLKIT_MODE = 3
end

function info()
  println("OMBackend.jl")
  println("A Julia backend for the Equation Oriented Language Modelica!")
  println("Run any test module by executing runExample(<model-name>)")
  println("Available Example models include:")
  for (k,v) in EXAMPLE_MODELS
    println(k)
  end
end

function plotGraph(shouldPlot::Bool)
  global PLOT_EQUATION_GRAPH = shouldPlot
end

function runExample(modelName::String)
  println("Translating Hybrid DAE of: $modelName")
  translate(EXAMPLE_MODELS[modelName])
  println("Translation done:")
  simulateModel(modelName)
  println("Simulating with default settings:")
  println("Model now in memory")
end


function runExampleExplicit(modelName::String)
  println("Translating Hybrid DAE of: $modelName")
  translate(EXAMPLE_MODELS[modelName], ODE_MODE)
  println("Translation done:")
  println("Simulating with default settings:")
  simulateModel(modelName)
  println("Model now in memory")
end


"""
Contains expressions of models currently in memory.
Do NOT mutate in other modules!
//John
"""
global COMPILED_MODELS = Dict()

""" 
  Active simulation functions 
  Do NOT mutate in other modules!
"""
global COMPILED_SIMULATION_CODE = Dict()

function translate(frontendDAE::DAE.DAE_LIST; BackendMode = DAE_MODE)::Tuple{String, Expr}
  local bDAE = lower(frontendDAE)
  local simCode
  if BackendMode == DAE_MODE
    simCode = generateSimulationCode(bDAE)
    return generateTargetCode(simCode)
  elseif BackendMode == ODE_MODE
    simCode = generateExplicitSimulationCode(bDAE)
    return ("Error", :())
  elseif BackendMode == MODELING_TOOLKIT_MODE
    @debug "Experimental: Generates and runs code using modelling toolkit"
    simCode = generateSimulationCode(bDAE)
    return generateMDKTargetCode(simCode)
  else
    @error "No mode specificed: valid modes are:"
    println("ODE_MODE")
    println("DAE_MODE")
  end
end

"""
 Transforms given DAE-IR/Hybrid DAE to backend DAE-IR (BDAE-IR)
"""
function lower(frontendDAE::DAE.DAE_LIST)::BDAE.BDAEStructure
  local bDAE::BDAE.BDAEStructure
  local simCode::SIM_CODE
  @debug "Length of frontend DAE:" length(frontendDAE.elementLst)
  @assert typeof(listHead(frontendDAE.elementLst)) == DAE.COMP
  #= Create Backend structure from Frontend structure =#
  bDAE = BDAECreate.lower(frontendDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "translated"));
  bDAE = Causalize.detectIfExpressions(bDAE)  
  @debug(BDAEUtil.stringHeading1(bDAE, "if equations transformed"));
  bDAE = Causalize.detectStates(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "states marked"));
  bDAE = Causalize.residualizeEveryEquation(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "residuals"));
  return bDAE
end

"""
  Transforms  BDAE-IR to simulation code for DAE-mode
"""
function generateSimulationCode(bDAE::BDAE.BDAEStructure)::SimulationCode.SimCode
  simCode = SimulationCode.transformToSimCode(bDAE)
  @debug BDAE.stringHeading1(simCode, "SIM_CODE: transformed simcode")
  return simCode
end

"""
  Transforms BDAE-IR to simulation code for MDK mode
"""
function generateUnsortedSimulationCode(bDAE::BDAE.BDAEStructure)::SimulationCode.UNSORTED_SIM_CODE
  simCode = SimulationCode.transformToSimCodeNoSort(bDAE)
  @debug BDAE.stringHeading1(simCode, "SIM_CODE: transformed simcode")
  return simCode
end

"""
  Transforms BDAE-IR to simulation code for DAE-mode
"""
function generateExplicitSimulationCode(bDAE::BDAE.BDAEStructure)
  simCode = SimulationCode.transformToExplicitSimCode(bDAE)
  @info "Code generation for ODE-mode not yet supported! Exiting.."
#  @debug BDAE.stringHeading1(simCode, "SIM_CODE: transformed simcode")
#  return simCode
end

"""
  Generates code interfacing DifferentialEquations.jl
  The resulting code is saved in a dictonary which contains functions that where simulated
  this session. Returns the generated modelName and corresponding generated code
"""
function generateTargetCode(simCode::SimulationCode.SIM_CODE)
  #= Target code =#
  (modelName::String, modelCode::Expr) = CodeGeneration.generateCode(simCode)
  @debug "Functions:" modelCode
  @debug "Model:" modelName
  COMPILED_MODELS[modelName] = modelCode
  return (modelName, modelCode)
end


"""
  Generates code interfacing ModelingToolkit.jl
  The resulting code is saved in a table which contains functions that where simulated
  this session. Returns the generated modelName and corresponding generated code
"""
function generateMDKTargetCode(simCode::SimulationCode.SIM_CODE)
  #= Target code =#
  (modelName::String, modelCode::Expr) = CodeGeneration.generateMDKCode(simCode)
  @debug "Functions:" modelCode
  @debug "Model:" modelName
  COMPILED_MODELS[modelName] = modelCode
  return (modelName, modelCode)
end

"""
  Generates code interfacing DifferentialEquations.jl
  The resulting code is saved in a dictonary which contains functions that where simulated
  this session. Returns the generated modelName and corresponding generated code
"""
function generateTargetCode(simCode::SimulationCode.EXPLICIT_SIM_CODE)
  #= Target code =#
  @info "Code generation for explicit simcode is not yet supported"
  try
    (modelName, modelCode) = CodeGeneration.generateCode(simCode)
  catch
    @info "ODE mode failed"
  end
end

function writeModelToFile(modelName::String)
  model = COMPILED_MODELS[modelName]
  fileName = "$modelName.jl"
  try
    formattedModel = JuliaFormatter.format_text(model)
    CodeGeneration.writeDAE_equationsToFile(modelName, model)
  catch
    @info "Failed writing $model to file: $fileName"
  end
end

"
  Prints a model. 
  If the specified model exists. Print it to stdout.
"
function printModel(modelName::String; keepComments = true)
    try
      local model::Expr = COMPILED_MODELS[modelName]
      #= Remove all the redudant blocks from the model =#
      if keepComments == false
        strippedModel = CodeGeneration.stripComments(model)
      end
      local strippedModel = CodeGeneration.stripBeginBlocks(model)
      local modelStr::String = "$strippedModel"
      formattedResults = JuliaFormatter.format_text(modelStr;
                                                    remove_extra_newlines = true,
                                                    always_use_return = true)
      println(formattedResults)
    catch e 
      @error "Model: $(modelName) is not compiled. Available models are: $(availableModels())"
      throw("Error printing model")
    end
end


"
    Prints compiled models to stdout
"
function availableModels()::String
  str = "Compiled models:\n"
    for m in keys(COMPILED_MODELS)
      str *= "  $m\n"
    end
  return str
end

"""
  Simulates model interactivly. 

TODO:
    (Currently does a redudant string conversion)
"""
function simulateModel(modelName::String; tspan=(0.0, 1.0))
  local modelCode::Expr = COMPILED_MODELS[modelName]
  local modelCodeStr = ""
  try
    modelCodeStr::String = "$modelCode"
    local parsedModel = Meta.parse.(modelCodeStr)
    @eval $parsedModel
    local modelRunnable = Meta.parse("OMBackend.$(modelName)Simulate($(tspan))")
    #= Run the model with the supplied tspan. =#
    @eval Main $modelRunnable
  catch err
    @info "Interactive evaluation failed: $err"
    println(modelCodeStr)
    @info "Dump of model-code"
    Base.dump(parsedModel)
    throw(err)
  end
end

"
  The default plot function of OMBackend.
  All labels of the variables and the name is given by default
"
function plot(sol::CodeGeneration.OMSolution)
  local nsolution = sol.diffEqSol
  local t = nsolution.t
  local rescols = collect(eachcol(transpose(hcat(nsolution.u...))))
  labels = permutedims(sol.idxToName.vals)
  Plots.plot(t, rescols; labels=labels)
end


"
  An alternative plot function in OMBackend.
  All labels of the variables and the name is given by default
"
function plot(sol)
  Plots.plot(sol)
end


"""
  Turns on logging
TODO:
  Make more fine grained
"""
function turnOnLogging()
  ENV["JULIA_DEBUG"] = "OMBackend"
end

function turnOnLoggingCodeGeneration()
  ENV["JULIA_DEBUG"] = "CodeGeneration"
end

