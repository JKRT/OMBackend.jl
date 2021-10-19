using Revise
using ModelingToolkit
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
  @info "Indices values:" indices
  #= To keep track of common variables. =#
  local commonVariableDict = ["PendulumModel" => vars1, "BouncingBall" => vars2]
  return (breakingPendulumModel, [changeStructure1], commonVariableDict)
end


(problem, structuralCallbacks, commonVariableDict) = BreakingPendulum((0.0, 7.))
import OMBackend
finalSolution = OMBackend.Runtime.solve(problem, (0.0, 7.), Rodas5(), structuralCallbacks, commonVariableDict)
