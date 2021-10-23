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
import OMFrontend

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
  MTK models
"""
const COMPILED_MODELS_MTK = Dict()

function translate(frontendDAE::DAE.DAE_LIST; BackendMode = MTK_MODE)::Tuple{String, Expr}
  local bDAE = lower(frontendDAE)
  local simCode
  if BackendMode == DAE_MODE
<<<<<<< HEAD
    simCode = generateSimulationCode(bDAE; mode = DAE_MODE)
    return generateTargetCode(simCode)
=======
    throw("DAE mode is removed.")
>>>>>>> Major changes. Remove DAE-mode. Now only supports MTK mode
  elseif BackendMode == MTK_MODE
    @debug "Experimental: Generates and runs code using modelling toolkit"
    simCode = generateSimulationCode(bDAE; mode = MTK_MODE)
    return generateMTKTargetCode(simCode)
  else
    @error "No mode specificed: valid modes are:"
    println("DAE_MODE")
    println("ModelingToolkit_MODE")
  end
end

"""
`function dumpInitialSystem(frontendDAE::DAE.DAE_LIST)`
 Dumps a textual representation of the initial system.
"""
function dumpInitialSystem(frontendDAE::DAE.DAE_LIST)::String
  str =  "Length of frontend DAE:" * string(length(frontendDAE.elementLst)) * "\n"
  bDAE = BDAECreate.lower(frontendDAE)
  str *= BDAEUtil.stringHeading1(bDAE, "translated")
  return str
end

"""
`function printInitialSystem(frontendDAE::DAE.DAE_LIST)`
 Dumps a textual representation of the initial system.
"""
function printInitialSystem(frontendDAE::DAE.DAE_LIST)
  print(dumpInitialSystem(frontendDAE::DAE.DAE_LIST))
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
  #= We always residualize since residuals are easier to work with =#
  bDAE = Causalize.residualizeEveryEquation(bDAE)
  @debug(BDAEUtil.stringHeading1(bDAE, "residuals"));
  return bDAE
end


"""
     Transforms given FlatModelica to backend DAE-IR (BDAE-IR).
  """
function lower(frontendDAE::OMFrontend.Main.FlatModel)
  local bDAE::BDAE.BACKEND_DAE
  local simCode::SIM_CODE
  @debug "Length of frontend DAE:" length(frontendDAE.elementLst)
  bDAE = BDAECreate.lower(frontendDAE)
end

"""
  Transforms  BDAE-IR to simulation code for DAE-mode
"""
function generateSimulationCode(bDAE::BDAE.BACKEND_DAE; mode)::SimulationCode.SimCode
  simCode = SimulationCode.transformToSimCode(bDAE; mode = mode)
  @debug BDAEUtil.stringHeading1(simCode, "SIM_CODE: transformed simcode")
  return simCode
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
  modelName = replace(modelName, "." => "__")
  COMPILED_MODELS[modelName] = modelCode
  return (modelName, modelCode)
end


"""
`generateMTKTargetCode(simCode::SimulationCode.SIM_CODE)`
  Generates code interfacing ModelingToolkit.jl
  The resulting code is saved in a table which contains functions that where simulated
  this session. Returns the generated modelName and corresponding generated code
"""
function generateMTKTargetCode(simCode::SimulationCode.SIM_CODE)
  #= Target code =#
  (modelName::String, modelCode::Expr) = CodeGeneration.generateMTKCode(simCode)
  @debug "Functions:" modelCode
  @debug "Model:" modelName
  modelName = replace(modelName, "." => "__")
  COMPILED_MODELS_MTK[modelName] = modelCode
  return (modelName, modelCode)
end

"""
  Writes a model to file by default the file is formatted and comments are kept.
"""
function writeModelToFile(modelName::String, filePath::String; keepComments = true, formatFile = true, mode = MTK_MODE)
  model = COMPILED_MODELS_MTK[modelName]
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
function printModel(modelName::String; MTK = true, keepComments = true, keepBeginBlocks = true)
  try
    local model::Expr
    model = COMPILED_MODELS_MTK[modelName]
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
    Prints available compiled models to stdout
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
`simulateModel(modelName::String; MODE = DAE_MODE ,tspan=(0.0, 1.0))`
  Simulates model interactivly.
"""
function simulateModel(modelName::String; MODE = MTK_MODE ,tspan=(0.0, 1.0))
  #= Strings containing . need to be in a format suitable for Julia =#
  modelName = replace(modelName, "." => "__")
  local modelCode::Expr  
  if MODE == MTK_MODE
    #= This does a redudant string conversion for now due to modeling toolkit being as is...=#
      modelCode = COMPILED_MODELS_MTK[modelName]
      local modelCodeStr = ""
    try
      @eval $(:(import OMBackend))
      strippedModel = CodeGeneration.stripBeginBlocks(modelCode)
      modelCodeStr::String = "$strippedModel"
      local parsedModel = Meta.parse(modelCodeStr)
      @eval $parsedModel
      local modelRunnable = Meta.parse("OMBackend.$(modelName)Simulate($(tspan))")
      #= Run the model with the supplied tspan. =#
      @eval Main $modelRunnable
    catch err
      @info "Interactive evaluation failed: $err with mode: $(MODE)"
      @info err
      println(modelCodeStr)
      @info modelCodeStr
      throw(err)
    end
  else
    throw("Unsupported mode")
  end
end

"
`plot(sol::Runtime.OMSolution)`
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
`function plot(sol)`
  An alternative plot function in OMBackend.
  All labels of the variables and the name is given by default
"
function plot(sol)
  Plots.plot(sol)
end


"""
`turnOnLogging(mod = "OMBackend"::String)`\n
Turns on logging. An optional parameter `mod` can be used to specify which model should be logged for more granuality.
"""
function turnOnLogging(mod = "OMBackend"::String)
  ENV["JULIA_DEBUG"] = mod
end
