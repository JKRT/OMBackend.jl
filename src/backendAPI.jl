using MetaModelica
using ExportAll
using Absyn
#= For interactive evaluation. =#
using ModelingToolkit
using SymbolicUtils


import .Backend.BDAE
import .Backend.BDAECreate
import .Backend.BDAEUtil
import .Backend.Causalize
import ..CodeGeneration
import OMBackend.ExampleDAEs
import .SimulationCode
import ..Runtime
import Base.Meta
import SCode
import JuliaFormatter
import Plots
import REPL
import OMBackend

const latexSymbols = REPL.REPLCompletions.latex_symbols
#= Settings =#
@enum BackendMode begin
  DAE_MODE = 1
  ODE_MODE = 2 #Currently not in operation
  MTK_MODE = 3
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

"""
  If this function is called a DAG representing the equations is plotted to your working directory.
"""
function plotGraph(shouldPlot::Bool)
  global PLOT_EQUATION_GRAPH = shouldPlot
end


"""
Contains expressions of models currently in memory.
Do NOT mutate in other modules!
//John
"""
const COMPILED_MODELS = Dict()

""" 
  MTK models
"""
const COMPILED_MODELS_MTK = Dict()

function translate(frontendDAE::DAE.DAE_LIST; BackendMode = DAE_MODE)::Tuple{String, Expr}
  local bDAE = lower(frontendDAE)
  local simCode
  if BackendMode == DAE_MODE
    simCode = generateSimulationCode(bDAE)
    return generateTargetCode(simCode)
  elseif BackendMode == MTK_MODE
    @debug "Experimental: Generates and runs code using modelling toolkit"
    simCode = generateSimulationCode(bDAE)
    return generateMTKTargetCode(simCode)
  else
    @error "No mode specificed: valid modes are:"
    println("DAE_MODE")
    println("ModelingToolkit_MODE")
  end
end

"""
 Transforms given DAE-IR/Hybrid DAE to backend DAE-IR (BDAE-IR)
"""
function lower(frontendDAE::DAE.DAE_LIST)::BDAE.BACKEND_DAE
  local bDAE::BDAE.BACKEND_DAE
  local simCode::SIM_CODE
  @debug "Length of frontend DAE:" length(frontendDAE.elementLst)
  @assert typeof(listHead(frontendDAE.elementLst)) == DAE.COMP
  #= Create Backend structure from Frontend structure =#
  bDAE = BDAECreate.lower(frontendDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "translated"));
  #= Expand arrays =#
  (bDAE, expandedVars) = Causalize.expandArrayVariables(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "Array variables expanded"));
  #= Expand Array variables in equation system=#
  bDAE = Causalize.detectAndReplaceArrayVariables(bDAE, expandedVars)
  @debug(BDAEUtil.stringHeading1(bDAE, "Equation system variables expanded"));
  #= Transform if expressions to if equations =#
  @debug(BDAEUtil.stringHeading1(bDAE, "if equations transformed"));
  bDAE = Causalize.detectIfExpressions(bDAE)
  #= Mark state variables =#
  bDAE = Causalize.detectStates(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "states marked"));
  bDAE = Causalize.residualizeEveryEquation(bDAE)
  #= =#
  @debug(BDAEUtil.stringHeading1(bDAE, "residuals"));
  return bDAE
end

"""
  Transforms  BDAE-IR to simulation code for DAE-mode
"""
function generateSimulationCode(bDAE::BDAE.BACKEND_DAE)::SimulationCode.SimCode
  simCode = SimulationCode.transformToSimCode(bDAE)
  @debug BDAEUtil.stringHeading1(simCode, "SIM_CODE: transformed simcode")
  return simCode
end

"""
  Transforms BDAE-IR to simulation code for DAE-mode
"""
function generateExplicitSimulationCode(bDAE::BDAE.BACKEND_DAE)
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
function generateMTKTargetCode(simCode::SimulationCode.SIM_CODE)
  #= Target code =#
  (modelName::String, modelCode::Expr) = CodeGeneration.generateMTKCode(simCode)
  @debug "Functions:" modelCode
  @debug "Model:" modelName
  COMPILED_MODELS_MTK[modelName] = modelCode
  return (modelName, modelCode)
end

function writeModelToFile(modelName::String, filePath::String; keepComments = true, formatFile = true)
  model = COMPILED_MODELS[modelName]
  fileName = "$modelName.jl"
  try
    if keepComments == false
      strippedModel = CodeGeneration.stripComments(model)
    end
    local strippedModel = CodeGeneration.stripBeginBlocks(model)
    local modelStr::String = "$strippedModel"
    formattedModel = if formatFile
      JuliaFormatter.format_text(modelStr,
                                 remove_extra_newlines = true,
                                 always_use_return = true)
    else
      modelStr
    end
    CodeGeneration.writeDAE_equationsToFile(filePath, formattedModel,)
  catch e
    @info "Failed writing $model to file: $fileName"
    @error e
  end
end

"""
  Prints a model. 
  If the specified model exists. Print it to stdout.
"""
function printModel(modelName::String; MTK = false, keepComments = true, keepBeginBlocks = true)
  try
    local model::Expr
    model =  if !MTK
      COMPILED_MODELS[modelName]
    else
      COMPILED_MODELS_MTK[modelName]
    end
      strippedModel = "$model"
      #= Remove all the redudant blocks from the model =#
      if keepComments == false
        strippedModel = CodeGeneration.stripComments(model)
      end
      if keepBeginBlocks == false
        strippedModel = CodeGeneration.stripBeginBlocks(model)
      end
      
      local modelStr::String = "$strippedModel"
      formattedResults = JuliaFormatter.format_text(modelStr;
                                                    remove_extra_newlines = true,
                                                    always_use_return = false)
      println(formattedResults)
    catch e 
      @error "Model: $(modelName) is not compiled. Available models are: $(availableModels())"
      throw("Error printing model")
    end
end

"""
    Prints compiled models to stdout
"""
function availableModels()::String
  str = "Compiled models (DAE-MODE):\n"
    for m in keys(COMPILED_MODELS)
      str *= "  $m\n"
    end

  str = "Compiled models (MTK-MODE):\n"
    for m in keys(COMPILED_MODELS_MTK)
      str *= "  $m\n"
    end
  
  return str
end

"""
  Simulates model interactivly. 
  
"""
function simulateModel(modelName::String; MODE = DAE_MODE ,tspan=(0.0, 1.0))
  local modelCode::Expr
  if MODE === DAE_MODE
    modelCode = COMPILED_MODELS[modelName]
    try
      @eval $(:(import OMBackend))
      @eval $modelCode
      local modelRunnable = Meta.parse("OMBackend.$(modelName)Simulate($(tspan))")
      #= Run the model with the supplied tspan. =#
      @eval Main $modelRunnable
    catch err
      @info "Interactive evaluation failed: $err"
      println(modelCodeStr)
      @info "Dump of model-code"
      #    Base.dump(parsedModel) TODO
      throw(err)
    end
  elseif MODE == MTK_MODE
    #= This does a redudant string conversion for now due to modeling toolkit being as is...=#
      modelCode = COMPILED_MODELS_MTK[modelName]
      local modelCodeStr = ""
    try
      @eval $(:(import OMBackend))
      modelCodeStr::String = "$modelCode"
      local parsedModel = Meta.parse.(modelCodeStr)
      @eval $parsedModel
      local modelRunnable = Meta.parse("OMBackend.$(modelName)Simulate($(tspan))")
      #= Run the model with the supplied tspan. =#
      @eval Main $modelRunnable
    catch err
      @info "Interactive evaluation failed: $err with mode: $(MODE)"
      @info err
      #= Base.dump(parsedModel) =#
      throw(err)
    end
  else
    throw("Unsupported mode")
  end
end

"
  The default plot function of OMBackend.
  All labels of the variables and the name is given by default
"
function plot(sol::Runtime.OMSolution)
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
"""
function turnOnLogging(mod = "OMBackend"::String)
  ENV["JULIA_DEBUG"] = mod
end
