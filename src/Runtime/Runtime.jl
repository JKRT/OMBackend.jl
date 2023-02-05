#=
  Simulation runtime implemented based on the integrator interface
=#
module Runtime
include("RuntimeUtil.jl")
import .RuntimeUtil
import Absyn
import SCode
import OMBackend
import OMBackend.SimulationCode
import OMBackend.CodeGeneration
import OMFrontend
import DAE

import ModelingToolkit
import ModelingToolkit.IfElse

using DataStructures
using ModelingToolkit
using DifferentialEquations
using MetaModelica

abstract type AbstractOMSolution end

"""
  Wrapper
"""
struct OMSolution{T1, T2} <: AbstractOMSolution
  "Solution given by DifferentialEquations.jl"
  diffEqSol::T1
  "Various metadata for the specific model"
  idxToName::T2
end

"""
  Wrapper object for equation based models that contain several solutions
"""
struct OMSolutions{T1, T2} <: AbstractOMSolution
  "Set of solutions given by DifferentialEquations.jl"
  diffEqSol::Vector{T1}
  "Various metadata for the specific model"
  idxToName::T2
end


abstract type AbstractStructuralChange end

"""
  Wrapper callback for a static structural change.
  That is the model
  we are simulating changes during the simulation but the future model can be predicted statically
"""
mutable struct StructuralChange{SYS} <: AbstractStructuralChange
  "The name of the next mode"
  name::String
  "Indicates if the structure has changed"
  structureChanged::Bool
  "The system we are switching to."
  system::SYS
end

"""
  Wrapper callback for structural change that triggers a recompilation.
"""
mutable struct StructuralChangeRecompilation{MOD <: Tuple} <: AbstractStructuralChange
  "The name of the next mode"
  name::String
  "Indicates if the structure has changed"
  structureChanged::Bool
  "The meta model. That is a SCode representation of the model itself"
  metaModel::SCode.CLASS
  "The modification to be applied during recompilation"
  modification::MOD
  """
    The symbol table for the old model.
    This is used to map indices of variables when the structure of the model changes
  """
  stringToSimVarHT
end

"""
 Wrapper callback for a dynamic connection reconfiguration
"""
mutable struct StructuralChangeDynamicConnection <: AbstractStructuralChange
  "The name of the next mode"
  name::String
  "Indicates if the structure has changed"
  structureChanged::Bool
  "The meta model. A flat representation of the model itself."
  #= Would it be better to modify the SCode instead? Less code to change?=#
  flatModel::OMFrontend.Main.FLAT_MODEL
  "The index of the specific dynamic connection equation."
  index::Int
  """
    The symbol table for the old model.
    This is used to map indices of variables when the structure of the model changes
  """
  stringToSimVarHT
  """
    If equations are to be added or removed.
  """
  activeEquations::Bool
end

mutable struct OM_ProblemStructural{T0 <: String, T1, T2, T3}
  "The name of the active mode"
  activeModeName::T0
  "The problem we are currently solving"
  problem::T1
  "The set of structural callbacks"
  structuralCallbacks::T2
  """ Variables that all modes have in common """
  commonVariables::T3
end

mutable struct OM_ProblemRecompilation{T0 <: String, T1, T2, T3}
  "The name of the active mode"
  activeModeName::T0
  "The problem we are currently solving"
  problem::T1
  "The set of structural callbacks"
  structuralCallbacks::T2
  "The set of callback conditons"
  callbackConditions::T3
end

mutable struct OM_Problem{T0, T1}
  problem::T0
end

#=
The current scheme for structural change.
Callbacks are created in the model that encompasses the two submodels.

These callbacks contains a boolean field that indicate if the structure has changed.
It also contains a field indicating what system we should switch to.

During solving, this field is set by the callback.
In the solver loop we iterate through these callbacks.

If a structural change was detected we act on it and change to the system pointed to by the callback.
Depending on the encompassing system we either just in time recompile the new system or we switch to the new.
Saving our current time step and reinitialize our new changed system
.
We should also statically detect if VSS simulation is needed since it is more resource heavy than regular simulation
=#

"""
  Custom solver function for Modelica code with structuralCallbacks to monitor the solving process
  (Using the integrator interface) from DifferentialEquations.jl
"""
function solve(omProblem::OM_ProblemStructural, tspan, alg; kwargs...)
  @info "Calling omProblem::OM_ProblemStructural"
  local problem = omProblem.problem
  local structuralCallbacks = omProblem.structuralCallbacks
  local commonVariableSet = omProblem.commonVariables
  local symsOfInitialMode = getSyms(problem)
  local activeModeName = omProblem.activeModeName
  #= The issue here is that the symbols might differ due to a change in the cref scheme =#
  local indicesOfCommonVariablesForStartingMode = getIndicesOfCommonVariables(symsOfInitialMode, commonVariableSet; destinationPrefix = activeModeName)
  #= Create integrator =#
  integrator = init(problem, alg, kwargs...)
  #@info "Value of tspan[2]" tspan[2]
  add_tstop!(integrator, tspan[2])
  local oldSols = []
  #= Run the integrator=#
  @label START_OF_INTEGRATION
  for i in integrator
    #= Check structural callbacks in order =#
    retCode = check_error(integrator)
    for cb in structuralCallbacks
      if cb.structureChanged
        #= Find the correct variables and map them between the two models =#
        local newSystem = cb.system
        indicesOfCommonVariables = getIndicesOfCommonVariables(getSyms(newSystem),
                                                               commonVariableSet;
                                                               destinationPrefix = cb.name)
        newU0 = Float64[i.u[idx] for idx in indicesOfCommonVariablesForStartingMode]
        #= Save the old solution together with the name and the mode that was active =#
        push!(oldSols, (integrator.sol, getSyms(problem), activeModeName))
        #= Now we have the start values for the next part of the system=#
        integrator = init(cb.system,
                          alg;
                          #                         t0 = i.t,
                          #                          u0 = newU0,
                          dt = 0.0000000001,
                          #                          tstart = tspan[1],
                          #                          tstop = tspan[2],
                          kwargs...)
        #=
          Reset with the new values of u0
          and set the active moed to the mode we are currently using.
        =#
        activeModeName = cb.name
        reinit!(integrator, newU0; t0 = i.t, reset_dt = true)
        cb.structureChanged = false
        #= goto to save preformance =#
        @goto START_OF_INTEGRATION
      end
    end
  end
  #= The solution of the integration procedure =#
  local solution = integrator.sol
  push!(oldSols, solution)
  return oldSols
end

global SHOULD_DO_REINITIALIZATION = true

"""
  Custom solver function for Modelica code with structuralCallbacks to monitor the solving process
  (Using the integrator interface) from DifferentialEquations.jl
"""
function solve(omProblem::OM_ProblemRecompilation, tspan, alg; kwargs...)
  local problem = omProblem.problem
  local structuralCallbacks = omProblem.structuralCallbacks
  local callbackConditions = omProblem.callbackConditions
  local activeModeName = omProblem.activeModeName
  #local integrator = init(problem, alg, dtmax = 0.001, kwargs...)
  local integrator = init(problem, alg,  kwargs...)
  local oldSols = []
  local solutions = []
  #= Run the integrator=#
  @label START_OF_INTEGRATION
  while true
    local i = integrator
    if i.t + i.dt >= tspan[2]
      @info "Ending simulation"
      Base.invokelatest(solve!, i)
      break
    end
    @info "Integration step was" i.t
    @info "dt was:" i.dt
    #= Check structural callbacks in order =#
    retCode = check_error(integrator)
    @info "Checking structural callbacks"
    for j in 1:length(structuralCallbacks)
      local cb = structuralCallbacks[j]
      if cb.structureChanged && i.t < tspan[2]
        @info "Callback triggered at:" integrator.t
        #= Recompile the system =#
        local oldHT = cb.stringToSimVarHT
        local newU0
        if ! SHOULD_DO_REINITIALIZATION
          @time (newProblem, newSymbolTable, initialValues, sc) = recompilation(cb.name, cb, integrator.u, tspan, callbackConditions)
          #= Assuming the indices are the same (Which is not always true) =#
          local symsOfOldProblem = getSyms(problem)
          local symsOfNewProlem = getSyms(newProblem)
          newU0 = RuntimeUtil.createNewU0(symsOfOldProblem,
                                          symsOfNewProlem,
                                          newSymbolTable,
                                          initialValues,
                                          integrator,
                                          sc)
          println("U0 before variables are set:")
          local discrete_events = REDUCED_SYSTEM.discrete_events
          @info "discrete_events:" discrete_events
          #= If there are discrete events these should be evaluated before proceeding =#
          events = if length(discrete_events) > 0
            RuntimeUtil.evalDiscreteEvents(discrete_events, newU0, i.t)
          else
            []
          end
          for e in events
            newU0[e[1]] = e[2]
          end
          println("Generate new u0:")
          println(newU0)
          @info "New u0 generated"
          # # TMP for System 10 With optimization
          global activeU0 = newU0

          #= Save the old solution together with the name and the mode that was active =#
          push!(solutions, integrator.sol)
          push!(oldSols, (integrator.sol, getSyms(problem), activeModeName))
          #= Now we have the start values for the next part of the system=#
          integrator = init(newProblem,
                            alg,
                            kwargs...)
          reinit!(integrator,
                  newU0;
                  t0 = i.t,
                  reset_dt = true)
        else #= Francesco proposed solution. =#
          #=
            Rerun OCC algorithms.
            Find the root variables of the OCC graph + the variables they assign and the root sources.
            That is the reference variables for the roots.
          =#
          @time (rootIndices, variablesToSetIdx, rootSources, variablestoReset) = returnRootIndices(cb.name,
                                                                                  cb,
                                                                                  integrator.u,
                                                                                  tspan,
                                                                                  problem)
          @assert length(rootIndices) == length(keys(rootSources)) "Root sources and indices must have the same length. Length was $(length(rootIndices)) == $(length(keys(rootSources)))"
          @info "Root indices" rootIndices
          @info "variablesToSetIdx" variablesToSetIdx
          @info "rootsources" rootSources
          local newU0::Vector{Float64} = Float64[v for v in integrator.u]
          local stateVars = states(OMBackend.LATEST_REDUCED_SYSTEM)
          #= This is bad, do not use strings thise way. =#
          local stateVarsAsStr = [replace(string(s), "(t)" => "") for s in stateVars]
          local OM_NameToMTKIdx = Dict()
          local rootKeys = [k for k in keys(rootSources)]
          local rootValues = [v for v in values(rootSources)]
          local variablesToSet = collect(Iterators.flatten([v for v in values(variablestoReset)]))
          @info "variablesToSet" variablesToSet
          local rootKeysToMTKIdx = indexin(rootKeys, stateVarsAsStr)
          @info "rootKeysToMTKIdx" rootKeysToMTKIdx
          local rootValsToMTKIdx = indexin(rootValues, stateVarsAsStr)
          @info "rootValsToMTKIdx" rootValsToMTKIdx
          local variablesToResetMTKIdx = indexin(variablesToSet, stateVarsAsStr)
          local rootToEquationMap::Dict{String, Symbolics.Equation} = Dict()
          @info "variablesToResetMTKIdx" variablesToResetMTKIdx
          for (i, rk) in enumerate(rootKeys)
            OM_NameToMTKIdx[rk] = rootKeysToMTKIdx[i]
          end
          for (i, rv) in enumerate(rootValues)
            OM_NameToMTKIdx[rv] = rootValsToMTKIdx[i]
          end
          for (i, vr) in enumerate(variablesToSet)
            OM_NameToMTKIdx[vr] = variablesToResetMTKIdx[i]
          end
          @info "OM_NameToMTKIdx" OM_NameToMTKIdx
          #=
           Start by setting the root variables we got from returnRootIndices.
           Each of these variables have a reference variable.
          =#
          for k in rootKeys
            rootStart = OM_NameToMTKIdx[k]
            rootSource = OM_NameToMTKIdx[rootSources[k]]
            @info "integrator.u[$(rootStart)] = integrator.u[$(rootSource)]"
            integrator.u[rootStart] = integrator.u[rootSource]
            rootToEquationMap[k] = ~(0, stateVars[rootStart] - stateVars[rootSource])
          end
          #=
            Get the variables of the system
          =#
          local equationDeps = ModelingToolkit.equation_dependencies(OMBackend.LATEST_REDUCED_SYSTEM)
          #= TODO: Not good should strive to fix the internal representations. =#
          local allEquationsAsStr = map((x)-> begin
                                          if !isempty(x)
                                              [replace(string(y), "(t)" => "") for y in x]
                                          end
                                        end, equationDeps)
          local rootKeys = [string(k) for k in keys(rootSources)]
          @info "rootKeys" rootKeys
          local equationToAddMap = Dict()
          #=
            Find the equations to replace
          =#
          global ALLEQUATIONSASSTR = allEquationsAsStr
          local oldEquations = equations(OMBackend.LATEST_REDUCED_SYSTEM)
          local assignedRoots = String[]
          for (i, eqVariables) in enumerate(allEquationsAsStr)
            if eqVariables !== nothing && length(eqVariables) == 2
              firstV = first(eqVariables)
              secondV = last(eqVariables)
              check1 = firstV in rootKeys && secondV in rootValues
              check2 = firstV in rootValues && secondV in rootKeys
              check3 = length(ModelingToolkit.difference_vars(oldEquations[i])) == 2
              if (check1 || check2) && check3
                @info "Root already assigned. Root was at $(i)-->" eqVariables
                push!(assignedRoots, last(eqVariables))
              end
            end
          end
          @info "assignedRoots" assignedRoots
          for (i, eqVariables) in enumerate(allEquationsAsStr)
            if eqVariables !== nothing && length(eqVariables) == 2
              local cand = first(eqVariables)
              if cand in rootKeys && ! (cand in assignedRoots)
                @info "Eq" cand
                equationToAddMap[cand] = i
              end
            end
          end
          #= The equation indices are the locations in which we are to insert our new equations=#
          @info "equationToAddMap" equationToAddMap
          #= One assignment need to be changed. =#
          for k in rootKeys
            @info "Value of rootStart" k
            varsToSet = variablestoReset[k]
            val = OM_NameToMTKIdx[k]
            @info "varsToSet" varsToSet
            for v in varsToSet
              vIdx = OM_NameToMTKIdx[v]
              @info "integrator[$(vIdx)] = integrator[$(val)]"
              integrator.u[vIdx] = integrator.u[val]
            end
          end
          push!(solutions, integrator.sol)
          push!(oldSols, (integrator.sol, getSyms(problem), activeModeName))
          newEquations = Symbolics.Equation[]
          @info "equationToAddMap" equationToAddMap
          @info "rootToEquationMap" rootToEquationMap
          for eq in oldEquations
            println(eq)
            push!(newEquations, eq)
          end
          for eq in keys(rootToEquationMap)
            if eq in assignedRoots
              @info "Eq that was in root was" eq
              continue
            end
            local replacementIdx = equationToAddMap[eq] #Hack
            local newEquation = rootToEquationMap[eq]
            @info "Setting newEquation[$(replacementIdx)] = $(eq)"
            newEquations[replacementIdx] = newEquation
            @info  "New equations at $(replacementIdx)" newEquation
          end
          #newEquations[4] = rootToEquationMap["G2"]
          println("New equations:")
          for eq in newEquations
            println(eq)
          end
          @time newSystem = ODESystem(newEquations,
                                independent_variable(OMBackend.LATEST_REDUCED_SYSTEM),
                                states(OMBackend.LATEST_REDUCED_SYSTEM),
                                parameters(OMBackend.LATEST_REDUCED_SYSTEM);
                                      name = Symbol(cb.name))
          newSystem = OMBackend.CodeGeneration.structural_simplify(newSystem)
          global NEW_SYSTEM = newSystem
          @time newProblem = ModelingToolkit.ODEProblem(
            newSystem,
            integrator.u,
            tspan,
            problem.p,
            #=
              TODO currently only handles a single structural callback.
            =#
            callback = callbackConditions #Should be changed look at method below
          )
          @time integrator = init(newProblem,
                                  alg;
                                  kwargs...)
          @time reinit!(integrator,
                        integrator.u;
                        t0 = i.t - i.dt,
                        reset_dt = true)
        end
        #= End recompilation =#
        #=
          Reset with the new values of u0
          and set the active mode to the mode we are currently using.
        =#
        #= Point the problem to the new problem =#
        problem = newProblem
        #= Note that the structural event might be triggered again here (We kill it later) =#
        @info "Stepping... at $(integrator.t)"
        Base.invokelatest(step!, integrator, true)
        #step!(integrator, true)
        @info "After step"
        #=
         We reset the structural change pointer again here just to make sure
         that we do not trigger the callback again.
        =#
        cb.structureChanged = false
        #= goto to save preformance =#
        @debug "Reset!"
        @goto START_OF_INTEGRATION
      end
    end
    @info "taking a step"
    Base.invokelatest(step!, i, true)
    #step!(integrator, true)
  end
  push!(solutions, integrator.sol)
  return solutions
end

"""
  Recompile the metamodel with some component changed.
  Returns a tuple of a new problem together with a new symbol table.
inputs:
  Name of the active model
  The structural callback causing the recompilation
  The current u values of the integrator
  The provided timespan
  The callback conditions (This to make sure that the new model have the same callbacks)
"""
function recompilation(activeModeName,
                       structuralCallback::StructuralChangeRecompilation,
                       integrator_u,
                       tspan,
                       callbackConditions)::Tuple
  #=  Recompilation =#
  #= Have the SCode =#
  #= - 1) Fetch the parameter from the structural callback =#
  local metaModel = structuralCallback.metaModel
  local modification = structuralCallback.modification
  local inProgram = MetaModelica.list(metaModel)
  local elementToChange = first(modification)
  #= Get the symbol table, using this we evaluate the new value based on the value of some parameter. =#
  local newValueExpr = Meta.parse(last(modification))
  newValueExpr = quote
    stringToSimVarHT = $(structuralCallback.stringToSimVarHT)
    $(newValueExpr)
  end
  local newValue = eval(newValueExpr)
  #= - 2) Change the parameters in the SCode via API (As specified by the modification)=#
  #=  2.1 Change the parameter so that it is the same as the modifcation. =#
  newProgram = MetaModelica.list(RuntimeUtil.setElementInSCodeProgram!(elementToChange, newValue, metaModel))
  local classToInstantiate = activeModeName
  #=- 3) Call the frontend + the backend + JIT compile Julia code in memory =#
  local flatModelica = first(OMFrontend.instantiateSCodeToFM(classToInstantiate, newProgram))
  #= 3.1 Run the backend=#
  (resultingModel, simulationCode) = runBackend(flatModelica, classToInstantiate)
  #= 4.0 Revaulate the model=#
  local modelName = replace(activeModeName, "." => "__") * "Model"
  @eval $(resultingModel)
  modelCall = quote
    $(Symbol(modelName))($(tspan))
  end
  (problem, callbacks, initialValues, reducedSystem, tspan, pars, vars) = @eval $(modelCall)
  #= Reconstruct the composite problem to keep the callbacks =#
  compositeProblem = ModelingToolkit.ODEProblem(
    reducedSystem,
    initialValues,
    tspan,
    pars,
    #=
      TODO currently only handles a single structural callback.
    =#
      callback = callbackConditions #Should be changed look at method below
  )
  #=Changed System=#
  #= 4.1 Update the structural callback with the new situation =#
  @match SOME(newMetaModel) = simulationCode.metaModel
  structuralCallback.metaModel = newMetaModel
  structuralCallback.stringToSimVarHT = simulationCode.stringToSimVarHT
  #= 4.x) Assign this system to newSystem. =#
  return (compositeProblem, simulationCode.stringToSimVarHT, initialValues, false)
end

"""
  Structural callback for dynamic connection handling.
  Returns (problem, symbol table, initial values).
"""
function recompilation(activeModeName,
                       structuralCallback::StructuralChangeDynamicConnection,
                       integrator_u,
                       tspan,
                       callbackConditions)
  @info "recompilation"
  #= Fetch the model that we were generating from memory. =#
  local flatModel = OMBackend.CodeGeneration.FLAT_MODEL
  local unresolvedConnectEquations = flatModel.unresolvedConnectEquations
  #= Get the relevant equation =#
  local indexOfEquation = structuralCallback.index
  local equationIf = MetaModelica.listGet(flatModel.DOCC_equations, indexOfEquation)
  @assert length(equationIf.branches) == 1
  println("Assertions fullfilled")
  if ! structuralCallback.activeEquations
    println("First branch")
    equationsToAdd = first(equationIf.branches).body
    newFlatModel = RuntimeUtil.createNewFlatModel(flatModel, unresolvedConnectEquations, equationsToAdd)
  else
    newFlatModel = RuntimeUtil.createNewFlatModel(flatModel, indexOfEquation, unresolvedConnectEquations)
  end
  (resultingModel, simulationCode) = runBackend(newFlatModel, activeModeName)
  println("New model generated")
  local model = replace(activeModeName, "." => "__")
  local modelName = string(model, "Model")
  #local result = OMBackend.modelToString(model; MTK = true, keepComments = false, keepBeginBlocks = false)
  println("We have a new model!\n");
  resultingModel = OMBackend.CodeGeneration.stripComments(resultingModel)
  resultingModel = OMBackend.CodeGeneration.stripBeginBlocks(resultingModel)
  result = "$resultingModel"
  @eval $(resultingModel)
  OMBackend.writeStringToFile(string("modfied", modelName * ".jl"), result)
  local newParsedModel = Meta.parse(result)
  @eval $(newParsedModel)
  local modelCall = quote
    $(Symbol(modelName))($(tspan))
  end
  (problem, callbacks, initialValues, reducedSystem, tspan, pars, vars) = @eval $(modelCall)
  @info "We have constructed a new problem..."
  global PROBLEM = problem
  global REDUCED_SYSTEM = reducedSystem
  #= Update meta model somehow=#
  return (problem,
          simulationCode.stringToSimVarHT,
          problem.u0,
          true)
end

"""
  This function returns the root indices of the OCC graph.
"""
function returnRootIndices(activeModeName,
                structuralCallback::StructuralChangeDynamicConnection,
                integrator_u,
                tspan,
                callbackConditions)
  @info "Calling returnRootIndices"
  local flatModel = OMBackend.CodeGeneration.FLAT_MODEL
  (variablestoReset, rootSources) = RuntimeUtil.resolveDOOCConnections(flatModel, flatModel.name)
  @info "variablestoReset" variablestoReset
  @info rootSources
  local rootVariables = keys(variablestoReset)
  @info "Root variables" rootVariables
  local ht = structuralCallback.stringToSimVarHT
  rootIndices = Int[]
  variablesToSet = []
  variablesToSetIdx = Vector{Int}[]
  for v in rootVariables
    indexOfRoot = first(ht[v])
    @info "Index of $(v) was:" indexOfRoot
    push!(rootIndices, indexOfRoot)
    push!(variablesToSet, values(variablestoReset[v]))
  end
  for variables in variablesToSet
    tmp = Int[]
    for v in variables
      idx = first(ht[v])
      push!(tmp, idx)
    end
    push!(variablesToSetIdx, tmp)
  end
  return (rootIndices, variablesToSetIdx, rootSources, variablestoReset)
end

"""
 Runs the backend.
  Translates the flat model to Flat Modelica.
  Generates the simulation code.
  Creates the new model.
Returns, a tuple of the new model and the simulation code of this model.
"""
function runBackend(flatModelica, classToInstantiate)
  local bdae = OMBackend.lower(flatModelica)
  local simulationCode = OMBackend.generateSimulationCode(bdae; mode = OMBackend.MTK_MODE)
  local newModel = OMBackend.CodeGeneration.ODE_MODE_MTK_MODEL_GENERATION(simulationCode, classToInstantiate)
  return (newModel, simulationCode)
end

"""
  Solving procedure without structural callbacks.
"""
function solve(omProblem::OM_Problem, tspan, alg; kwargs...)
  local problem = omProblem.problem
  #= Create integrator =#
  integrator = init(problem, alg, stop_at_next_tstop = true, kwargs...)
  add_tstop!(integrator, tspan[2])
  for i in integrator
  end
  #= Return the final solution =#
  return integrator.sol
end

"""
    Fetches the symbolic variables from a problem.
"""
function getSyms(problem::ODEProblem)
  return problem.f.syms
end

"""
  Fetches  the symbolic variables from a solution
"""
function getSymsFromSolution(sol)
  return sol.prob.f.syms
end


"""
  Get a vector of indices of the variables between syms1 and syms2.
  The destination mode might have a prefix for it's symbols.
  The prefix string is used to give the common variables the correct prefix
  between transistions.
"""
function getIndicesOfCommonVariables(syms1::Vector{Symbol}, names::Vector{String}; destinationPrefix::String = "")
  #= The common variables have the name without the prefix of the destination system =#
  local syms2 = [Symbol(destinationPrefix * "_" * name * "(t)") for name in names]
  local indicesOfCommonVariables = Int[]
  local idxDict1 = DataStructures.OrderedDict()
  local idxDict2 = DataStructures.OrderedDict()
  local counter = 1
  for sym in syms1
    idxDict1[sym] = counter
    counter += 1
  end
  counter = 1
  for sym in syms2
    idxDict2[sym] = counter
    counter += 1
  end
  local commonVariables = ∩(keys(idxDict1), keys(idxDict1))
  (smallestKeyset, dict) = if length(keys(idxDict1)) < length(keys(idxDict2))
    keys(idxDict1), idxDict2
  else
    keys(idxDict2), idxDict1
  end
  for key in smallestKeyset
    if haskey(dict, key)
      push!(indicesOfCommonVariables, dict[key])
    end
  end
  return indicesOfCommonVariables
end

end
