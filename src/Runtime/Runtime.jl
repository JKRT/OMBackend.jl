#=
  Simulation runtime implemented based on the integrator interface
=#
module Runtime

using DataStructures
using DiffEqBase
using DifferentialEquations
using ModelingToolkit

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
  Custom solver function to monitor the solving process
  (Using the integrator interface) from DifferentialEquations.jl
"""

function solve(problem, tspan, alg, structuralCallbacks, commonVariableDict; kwargs...)
  @info "Calling solve!"
  #= Create integrator =#
  integrator = init(problem, alg, dtmax = 0.1, kwargs...)
  add_tstop!(integrator, tspan[2])
  oldSols = []
  #= Run the integrator=#
  @label START_OF_INTEGRATION
  println(integrator.t)
  @show integrator
  for i in integrator
    #= Check structural callbacks in order =#
    println("Stepping:")
#    println("Status of i: $(i)")
    retCode = check_error(integrator)
    for cb in structuralCallbacks
      if cb.structureChanged
        println("STRUCTURE CHANGED!")
        println("Status of i: $(i)")
        #= Find the correct variables and map them between the two models  =#
        indicesOfCommonVariables = getIndicesOfCommonVariables(getSyms(problem), getSyms(cb.system))
        newU0 = Float64[i.u[idx] for idx in indicesOfCommonVariables]
        #= Now we have the start values for the next part of the system=#
        push!(oldSols, integrator.sol)
        integrator = init(cb.system,
                          alg;
                          t0 = i.t,
                          u0 = newU0,
                          tstop = tspan[2],
                          dtmax = 0.1,
                          kwargs...)
        #= Reset with the new values of u0 =#
        reinit!(integrator, newU0; t0 = i.t, reset_dt = true)
        cb.structureChanged = false
        println("!!DONE CHANGING THE STRUCTURE!! Restarting")
        #= goto to save preformance =#
        @goto START_OF_INTEGRATION
      end
    end
  end
  #= The solution of the integration procedure =#
  local solution = integrator.sol
  #= The final solution =#
  #= in oldSols we have the old solution. =#
  local newTimePoints = vcat(first(oldSols).t, solution.t) #TODO should be a loop here.

  #=TODO: We have a set of common variables. Find the indices for this set. =#
  #= We are creating a new Vector{Vector{Float64}}: =#
  #= Each solution in oldSol is solution to some subsolution of the VSS, before structural change =#
  #= If the dimensions between the solutions are not equivivalent we should only keep the columns we need =#
#  indicesOfCommonVariables = getIndicesOfCommonVariables(getSyms(problem), getSyms(cb.system))
  local startingSol = first(oldSols)
  local newUs = [startingSol[i,:]  for i in [1,2,5,6]] #TODO be clever here. 
  #= Now we need to merge these with the new solution =#
  @info "New US before vcat:" newUs[1]
  a = vcat(newUs[1], solution[1,:])
  b = vcat(newUs[2], solution[2,:])
  c = vcat(newUs[3], solution[3,:])
  d = vcat(newUs[4], solution[4,:])

  #=Convert into a matrix and then into a vector of vector again to get the right dimensions. =#
  newUs = transpose(hcat([a,b,c,d]...))
  newUs = [newUs[:,i] for i in 1:size(newUs,2)]
  #= For anyone reading this.. this could have been done better! =#
  
  #= Should be the common varibles =#  
  local sol = SciMLBase.build_solution(solution.prob, solution.alg, newTimePoints, newUs)
  #= Return the final solution =#
  return sol
end

"""
  Solving procedure without structural callbacks.
"""
function solve(problem, tspan, alg; kwargs...)
  @info "Calling solve!"
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
function getSyms(problem)
  return problem.f.syms
end


"""
  Get the indices of the variables common between syms1 and syms2
"""
function getIndicesOfCommonVariables(syms1, syms2)
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
  # @info "idxDict1" idxDict1
  # @info "idxDict2" idxDict2
  # @info commonVariables
  # @info smallestKeyset
  for key in smallestKeyset
    if haskey(dict, key)
      push!(indicesOfCommonVariables, dict[key])
    end
  end
  return indicesOfCommonVariables
end

end
