using Revise
using ModelingToolkit
import OMBackend
include("Pendulum.jl")
include("FreeFall.jl")

#=
This model contains variables but no equations of its own.
The first set of equations are those of pendulum
=#

mutable struct StructuralChange
  structureChanged::Bool
  system
end

"""
  Structural callback. f is the structure we are changing into.
"""
function structuralCallback(f)
  structuralChange = StructuralChange(false, f)
  function affect!(integrator)
    println("Detecting structure change!")
    structuralChange.structureChanged = true
  end
  function condition(u, t, integrator)
    return t - 5
    end
  cb = ContinuousCallback(condition, affect!)
  return (cb, structuralChange)
end


function BreakingPendulum(tspan)
  (_, initialValues, reducedSystem, tspan, pars) = PendulumModel(tspan)
  (freeFall,_,_,_,_) = FreeFallModel(tspan)
  (sCB, changeStructure1) = structuralCallback(freeFall)
  #= Our initial model is the pendulum model (but with a structural change specified.) =#
   breakingPendulumModel = ModelingToolkit.ODEProblem(
    reducedSystem,
    initialValues,
    tspan,
    pars,
    callback = CallbackSet(sCB),
   )
  local commonVariableSet = [Symbol("x(t)"), Symbol("y(t)"), Symbol("vx(t)"), Symbol("vy(t)")]
  return (breakingPendulumModel, [changeStructure1], commonVariableSet)
end

(problem, structuralCallbacks) = BreakingPendulum((0.0, 7.))
finalSolution = OMBackend.Runtime.solve(problem, (0.0, 7.), Rodas5(), structuralCallbacks, commonVariableSet)
