using DiffEqBase
using DifferentialEquations
using Plots
using Sundials
using Revise

using OM

#= Just assuming that it is indeed compiled:) =#
flyingBall = OM.flatten("FlyingBall", "./FlyingGodBall.mo")

function BouncingBallRealsStartConditions(p, t0)
  local x0 = Array{Float64}(undef, 2)
  local dx0 = Array{Float64}(undef, 2)
  x0[2] = (1.0) #= h =#
  return x0, dx0
end

function BouncingBallRealsDifferentialVars()
    return Bool[1, 1]
end

function BouncingBallRealsDAE_equations(res, dx, x, p, t) #=time=#
    res[1] = ((dx[1]) - ((-(p[1])))) #= g =#
    res[2] = ((dx[2]) - (x[1])) #= v =#
end


#= Assuming it has been generated =#
function BouncingBallRealsDAE_equations2(res, dx, x, p, t) #=time=#
    res[1] = ((dx[1]) - (((p[1])))) #= g =#
    res[2] = ((dx[2]) - (x[1])) #= v =#
end

function BouncingBallRealsParameterVars()
    p = Array{Float64}(undef, 2)
    p[1] = (9.81) #= g =#
    p[2] = (0.7) #= e =#
    return p
end

global RES
global EVENT_TIME_VALUES

function BouncingBallRealsCallbackSet(p)#=parameters=#
    #= WHEN EQUATIONS. Discrete callbacks =#
    condition1(x, t, integrator) = begin
        return x[2]# - 0
    end
    affect1!(integrator) =
        begin
          integrator.u[1] = -((p[2]) * ((integrator.u[1]))) #= e =#
          terminate!(integrator)
          global EVENT_TIME_VALUES = integrator.u
          global RES = true
        end
    cb1 = ContinuousCallback(condition1,
                             affect1!,
                             rootfind=true,
                             save_positions=(false,true),
                             affect_neg! = affect1!)
    #= IF EQUATIONS =#
    return CallbackSet(cb1)
end

function BouncingBallRealsSimulate(tspan = (0.0, 1.0))
    # Define problem
    p = BouncingBallRealsParameterVars()
    (x0, dx0) = BouncingBallRealsStartConditions(p, tspan[1])
    differential_vars = BouncingBallRealsDifferentialVars()
    #= Pass the residual equations =#
    problem = DAEProblem(
        BouncingBallRealsDAE_equations,
        dx0,
        x0,
        tspan,
        p,
        differential_vars = differential_vars,
        callback = BouncingBallRealsCallbackSet(p),
    )
    # Solve with IDA:)
  solution = solve(problem, IDA())
  if RES
    @info "Result after change"
    @info "Event time values:" EVENT_TIME_VALUES
    @info "Solution:" solution
    #= Compile other models in memory. Assume I have done so now =#
    old_solution = solution
    problem = DAEProblem(
      BouncingBallRealsDAE_equations2,
      dx0,
      EVENT_TIME_VALUES,
      (old_solution.t[length(old_solution.t)], 1.0),
      p,
      differential_vars
    )
    # Solve with IDA:)
    new_solution = solve(problem, IDA())
  end
  return solution
end
