using Revise
using ModelingToolkit
import OMBackend
include("Pendulum.jl")
include("BouncingBall.jl")

#=
This model contains variables but no equations of its own.
The first set of equations are those of pendulum
=#

mutable struct StructuralChange #= We keep a pointer to this structure. It is changed outside the callback) =#
  structureChanged::Bool
  system
end

"""
  Structural callback. f is the structure we are changing into.
"""
function structuralCallback(f)
  #= Represent structural change. =#
  local structuralChange = StructuralChange(false, f)
  function affect!(integrator)
    println("Detecting structure change!")
    structuralChange.structureChanged = true
  end
  function condition(u, t, integrator)
    return t - 3
    end    
  local cb = ContinuousCallback(condition, affect!)
  return (cb, structuralChange)
end

function BreakingPendulum(tspan)
  (pendulum, initialValues, reducedSystem, tspan, pars, vars1) = PendulumModel(tspan)
  (bouncingBall, _ ,_ ,_ ,_, vars2) = BouncingBallModel(tspan)
  (sCB, changeStructure1) = structuralCallback(bouncingBall)
   breakingPendulumModel = ModelingToolkit.ODEProblem(
    reducedSystem,
    initialValues,
    tspan,
    pars,
    callback = CallbackSet(sCB),
   )

  #=Some inline testing=#
  indices = OMBackend.Runtime.getIndicesOfCommonVariables(OMBackend.Runtime.getSyms(bouncingBall), OMBackend.Runtime.getSyms(pendulum))
  #= To keep track of common variables. We will know the common variables a priori =#
  local commonVariableSet = [Symbol("x(t)"), Symbol("y(t)"), Symbol("vx(t)"), Symbol("vy(t)")]
  return (breakingPendulumModel, [changeStructure1], commonVariableSet)
end


(problem, structuralCallbacks, commonVariableSet) = BreakingPendulum((0.0, 7.))
finalSolution = OMBackend.Runtime.solve(problem, (0.0, 7.), Rodas5(), structuralCallbacks, commonVariableSet)
