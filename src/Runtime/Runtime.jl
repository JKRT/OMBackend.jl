#=
  Simulation runtime implemented based on the integrator interface
=#
module Runtime
include("RuntimeUtil.jl")
import .RuntimeUtil
import Absyn
import SCode
import OMBackend
import OMFrontend
import ..CodeGeneration

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
  Wrapper callback for structural change that triggers a recompilation
"""
mutable struct StructuralChangeRecompilation{MOD} <: AbstractStructuralChange
  "The name of the next mode"
  name::String
  "Indicates if the structure has changed"
  structureChanged::Bool
  "The meta model"
  metaModel::SCode.CLASS
  modification::MOD
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

mutable struct OM_ProblemRecompilation{T0 <: String, T1, T2}
  "The name of the active mode"
  activeModeName::T0
  "The problem we are currently solving"
  problem::T1
  "The set of structural callbacks"
  structuralCallbacks::T2
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
Saving our current time step and reinitialize our new changed system.

We should also statically detect if VSS simulation is needed since it is more resource heavy than regular simulation
=#

"""
  Custom solver function for Modelica code with structuralCallbacks to monitor the solving process
  (Using the integrator interface) from DifferentialEquations.jl
"""
function solve(omProblem::OM_ProblemStructural, tspan, alg; kwargs...)
  local problem = omProblem.problem
  local structuralCallbacks = omProblem.structuralCallbacks
  local commonVariableSet = omProblem.commonVariables
  local symsOfInitialMode = getSyms(problem)
  local activeModeName = omProblem.activeModeName
  #= The issue here is that the symbols might differ due to a change in the cref scheme =#
  local indicesOfCommonVariablesForStartingMode = getIndicesOfCommonVariables(symsOfInitialMode, commonVariableSet; destinationPrefix = activeModeName)
  #= Create integrator =#
  integrator = init(problem, alg, dtmax = 0.01, kwargs...)
  #@info "Value of tspan[2]" tspan[2]
  add_tstop!(integrator, tspan[2])
  local oldSols = Tuple[]
  #= Run the integrator=#
  @label START_OF_INTEGRATION
  #@info "START OF INTEGRATION"
  for i in integrator
    #= Check structural callbacks in order =#
    retCode = check_error(integrator)
    #@info "Value of retCide" retCode
    #@info "SOLVER STEP:"
    for cb in structuralCallbacks
      if cb.structureChanged
        #@info("STRUCTURE CHANGED!")
        #@info("Status of i: $(i)")
        #= Find the correct variables and map them between the two models  =#
        local newSystem = cb.system
        indicesOfCommonVariables = getIndicesOfCommonVariables(getSyms(newSystem),
                                                               commonVariableSet;
                                                               destinationPrefix = cb.name)
        #@info "DONE COMPUTING INDICES"
        newU0 = Float64[i.u[idx] for idx in indicesOfCommonVariables]
        #= Save the old solution together with the name and the mode that was active =#
        push!(oldSols, (integrator.sol, getSyms(problem), activeModeName))
        #= Now we have the start values for the next part of the system=#
        integrator = init(cb.system,
                          alg;
                          t0 = i.t,
                          u0 = newU0,
                          dt = 0.01,
                          dtmax = 0.01,
                          tstart = tspan[1],
                          tstop = tspan[2],
                          kwargs...)
        #=
        Reset with the new values of u0
        and set the active moed to the mode we are currently using.
        =#
        activeModeName = cb.name
        reinit!(integrator, newU0; t0 = i.t, reset_dt = true)
        cb.structureChanged = false
        #@info "!!DONE CHANGING THE STRUCTURE!! Restarting"
        #= goto to save preformance =#
        @goto START_OF_INTEGRATION
      end
    end
  end
  #= The solution of the integration procedure =#
  local solution = integrator.sol
  #@info "Solution:" solution  
  #= The final solution =#
  #= in oldSols we have the old solution. =#
  local startingSol = if ! isempty(oldSols)
    first(first(oldSols))
  else #= If we have no old solution we are done =#
    return solution
  end
#  @info "Starting solution" startingSol
  local newTimePoints = vcat(startingSol.t, solution.t) #TODO should be a loop here.
  #= We are creating a new Vector{Vector{Float64}}: =#
  #= Each solution in oldSol is solution to some subsolution of the VSS, before structural change =#
  #= If the dimensions between the solutions are not equivivalent we should only keep the columns we need =#
#  indicesOfCommonVariables = getIndicesOfCommonVariables(getSyms(problem), getSyms(cb.system))
  #= The starting point is the initial variables of our starting mode =#
  local newUs = [startingSol[i,:]  for i in indicesOfCommonVariablesForStartingMode]
  #= Now we need to merge these with the latest solution =#
  #@info "Common variable set:" commonVariableSet
  #@info "DestinationPrefix:" activeModeName
  #@info "Syms from the solution:" getSymsFromSolution(solution)
  local tmp = getIndicesOfCommonVariables(getSymsFromSolution(solution), commonVariableSet; destinationPrefix = activeModeName)
  #@info "Value of tmp:" tmp
  for i in 1:length(commonVariableSet)
    newUs[i] = vcat(newUs[i], solution[tmp[i],:])
  end
  #=Convert into a matrix and then into a vector of vector again to get the right dimensions. =#
  newUs = transpose(hcat(newUs...))
  newUs = [newUs[:,i] for i in 1:size(newUs,2)]
  #= For anyone reading this.. this the transformations above could have been done better! =#
  #= Should be the common varibles =#  
  local sol = SciMLBase.build_solution(solution.prob, solution.alg, newTimePoints, newUs)
  #= Return the final solution =#
  return sol
end


"""
  Custom solver function for Modelica code with structuralCallbacks to monitor the solving process
  (Using the integrator interface) from DifferentialEquations.jl
"""
function solve(omProblem::OM_ProblemRecompilation, tspan, alg; kwargs...)
  local problem = omProblem.problem
  local structuralCallbacks = omProblem.structuralCallbacks
  local activeModeName = omProblem.activeModeName
  integrator = init(problem, alg, dtmax = 0.01, kwargs...)
  @info "Value of tspan[2]" tspan[2]
  local oldSols = []
  #= Run the integrator=#
  @label START_OF_INTEGRATION
  @info "START OF INTEGRATION"
  for i in integrator
    #= Check structural callbacks in order =#
    retCode = check_error(integrator)
    for cb in structuralCallbacks
      if cb.structureChanged
        @info "STRUCTURAL CHANGE DETECTED"
        #= Recompile the system =#
        @info "RECOMPILING"
        local newSystem = recompilation(cb.name, cb, integrator.u, tspan)
        @info "RECOMPILATION DONE"
        #= End recompilation =#
        @info "DONE COMPUTING INDICES"
        #= Assuming the indices are the same (Which is not neccesary the same) =#
        local newU0 = Float64[integrator.u[idx] for idx in 1:10]        
        #= Save the old solution together with the name and the mode that was active =#
        push!(oldSols, (integrator.sol, getSyms(problem), activeModeName))
        #= Now we have the start values for the next part of the system=#
        @info "Making a new integrator"
        integrator = init(newSystem,
                          alg;
                          t0 = i.t,
                          u0 = newU0,
                          tstop = tspan[2],
                          dt = 0.01,
                          dtmax = 0.01,
                          kwargs...)
        @info "We have a new integrator"
        #=
        Reset with the new values of u0
        and set the active moed to the mode we are currently using.
       =#
        add_tstop!(integrator, tspan[2])
        activeModeName = cb.name
        reinit!(integrator, newU0; t0 = i.t, reset_dt = true)
        cb.structureChanged = false
 #       @info "!!DONE CHANGING THE STRUCTURE!! Restarting"
        #= goto to save preformance =#
        @goto START_OF_INTEGRATION
      end
    end
  end
  sol = integrator.sol
  push!(oldSols, sol)
  return oldSols
end

"""
  Recompiles the metamodel with some component changed
"""
function recompilation(activeModeName, structuralCallback, u, tspan)
  #=  Recompilation =#
  #= Have the SCode =#        
  #= - 1) Fetch the parameter from the structural callback =#
  local metaModel = structuralCallback.metaModel
  local modification = structuralCallback.modification
  local inProgram = MetaModelica.list(metaModel)
  local elementToChange = first(modification)
  local newValue = last(modification)
  #= - 2) Change the parameters in the SCode via API =#
  #= !!!TODO!!! =#
  #=  2.1 Change the parameter so that it is the same as the modifcation. =#
  newProgram = MetaModelica.list(RuntimeUtil.setElementInSCodeProgram!(elementToChange, newValue, metaModel))
  local classToInstantiate = activeModeName
  #=- 3) Call the frontend + the backend + JIT compile Julia code in memory =#
  local flatModelica = first(OMFrontend.instantiateSCodeToFM(classToInstantiate, newProgram))
  @info "Modified model!"
  println(OMFrontend.toString(flatModelica))
  @info "DONE"
  local bdae = OMBackend.lower(flatModelica)
  local simulationCode = OMBackend.generateSimulationCode(bdae; mode = OMBackend.MTK_MODE)
  #= Here I also have a new index mapping. That is I know what the new var -> u is. =#
  local resultingModel = OMBackend.CodeGeneration.ODE_MODE_MTK_MODEL_GENERATION(simulationCode, classToInstantiate)
  local modelName = replace(activeModeName, "." => "__") * "Model"
  @eval $(resultingModel)
  modelCall = quote
    $(Symbol(modelName))($(tspan))
  end
  (problem, initialValues, reducedSystem, tspan, pars, vars) = @eval $(modelCall)
  #= - 4) Assign this system to newSystem. =#
  return problem
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
