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

#= TODO: add StructuralChange callbacks here somehow so we can monitor the execution =#

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
    println(i)
    retCode = check_error(integrator)
    for cb in structuralCallbacks
      if cb.structureChanged
        println("STRUCTURE CHANGED!")
        println("Status of i: $(i)")
        newU0 = [i.u[1], i.u[2], i.u[5], i.u[6]]
        reinit!(integrator; t0 = i.t, reset_dt = true)
        push!(oldSols, integrator.sol)
        @info "Value of i.t" i.t
        integrator = init(cb.system,
                          alg;
                          t0 = i.t,
                          u0 = newU0,                          
                          tstop = tspan[2])
        #= Reset with the new values of u0 =#
        reinit!(integrator, newU0; t0 = i.t, reset_dt = true)
        cb.structureChanged = false
        println("!!DONE CHANGING THE STRUCTURE!! Restarting")
        @goto START_OF_INTEGRATION
      end
    end
  end
  #= Return the final solution =#
  return [oldSols..., integrator.sol]
end

"""
Solve without structural callbacks
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
