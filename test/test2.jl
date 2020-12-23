using DiffEqBase
using DifferentialEquations
using Plots
using Sundials
using Revise


function BouncingBallRealsDAE_equations(res, dx, x, p, t) #=time=#
    res[1] = dx[1] - p[1] #= g =#
    res[2] = dx[2] - x[1] #= v =#
end

function condition(u,t,integrator) # Event when event_f(u,t) == 0
  u[1]
end

function affect!(integrator)
  integrator.u[2] = -integrator.u[2]
end

function BouncingBallRealsParameterVars()
    p = Array{Float64}(undef, 1)
    p[1] = 9.81 #= g =#
    return p
end

cb = ContinuousCallback(condition, affect!)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
p = BouncingBallRealsParameterVars()
prob = DAEProblem(BouncingBallRealsDAE_equations,Array{Float64}(undef, 2),Array{Float64}(undef, 2), u0,tspan,callback=cb)
sol = solve(prob,IDA())
plot(sol)
