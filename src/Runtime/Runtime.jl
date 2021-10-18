#=
  Simulation runtime implemented based on the integrator interface
=#
module Runtime

using ModelingToolkit
using DiffEqBase
using DifferentialEquations

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
function solve(problem, tspan, alg, structuralCallbacks; kwargs...)
  @info "Calling solve!"
  #= Create integrator =#
  integrator = init(problem, alg, stop_at_next_tstop = true, kwargs...)
  add_tstop!(integrator, tspan[2])
  oldSols = []
  #= Run the integrator=#
  @label START_OF_INTEGRATION
  println(integrator.t)
  @show integrator
  for i in integrator
    #= Check structural callbacks in order =#
    println("Stepping:")
    println("Status of i: $(i)")
    retCode = check_error(integrator)
    for cb in structuralCallbacks
      if cb.structureChanged
        println("STRUCTURE CHANGED!")
        println("Status of i: $(i)")
        #= TODO: Here we need to find the correct variables and map them between the two models. These indices need to be saved. =#
        newU0 = [i.u[1], i.u[2], i.u[5], i.u[6]]
        push!(oldSols, integrator.sol)
        #= =#
        #reinit!(integrator; t0 = i.t, reset_dt = true)
        integrator = init(cb.system,
                          alg;
                          t0 = i.t,
                          u0 = newU0,
                          tstop = tspan[2])
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

  local startingSol = first(oldSols)
  local newUs = [startingSol.u[1], startingSol.u[2], startingSol.u[5], startingSol.u[6]]

  #= Now we need to merge these with the new solution =#
  newUs[1] = vcat(newUs[1], solution.u[1])
  newUs[2] = vcat(newUs[2], solution.u[2])
  newUs[3] = vcat(newUs[3], solution.u[3])
  newUs[4] = vcat(newUs[4], solution.u[4])
  @info "New time points" newTimePoints
  @info "Old u values:" first(oldSols).u
  @info "New u values:" solution.u
  @info "Merged u values:" newUs
  #= Should be the common varibles =#
  local vars = newUs #[solution[i][2] for i in 1:length(solution)]
  local N = length(vars)
  println("Vars>", vars)
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


end
