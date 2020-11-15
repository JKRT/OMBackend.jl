
using DifferentialEquations
using DiffEqBase
using Sundials
using Plots

function BouncingBallStartConditions(p, t0)
  local x0 = Array{Float64}(undef, 5)
  local dx0 = Array{Float64}(undef, 5)
  x0[5] #= flying =# = (true)
  x0[2] #= h =# = (1.0)
  return x0, dx0
end

function BouncingBallDifferentialVars()
  return Bool[1, 1, 0, 0, 0]
end

function BouncingBallDAE_equations(res, dx, x, p, t #=time=#)
  res[1] = ((dx[2]  #= der(h) =#) - (x[1] #= v =#))
  res[2] = (dx[1]  #= der(v) =#)
            if 0 < x[5]  #= flying =#
              - p[1] #= g =#
            else
              0.0
            end
  res[3] = ((x[4] #= impact =#) - ((x[2] #= h =#) <= (0.0)))
end

function BouncingBallParameterVars()
  p = Array{Float64}(undef, 2)
  p[1] #= g =# = 9.81
  p[2] #= e =# = 0.7
  return p
end

function BouncingBallCallbackSet()
  condition(x,t,integrator) = begin
    (#= Array exp=# reduce(|, [(((x[2] #= h =#) <= (0.0)) && ((x[1] #= v =#) <= (0.0))),(x[4] #= impact =#),]))
  end
  affect!(integrator) = begin
    x[3] =
      if 0 < x[6]
        (- ((p[2] #= e =#) * x[3]))
      else
        0.0
      end
  end
  integrator.x[1] = integrator.x[3] #= v_new =#
  integrator.x[5]  = integrator.x[3]  > 0.0
  #=The Discrete Callback corresponds to when. Discrete Hybrid Systems =#
  #= The continious callbacks would be applied to if-eqs=#
  cb1 = DiffEqBase.DiscreteCallback(condition, affect!)
  return DiffEqBase.CallbackSet(cb1)
end

function BouncingBallSimulate(tspan = (0.0, 3.0))
  # Define problem
  p_is = BouncingBallParameterVars()
  (x0, dx0) = BouncingBallStartConditions(p_is, tspan[1])
  @info x0
  differential_vars = BouncingBallDifferentialVars()
  #= Pass the residual equations =#
  problem = DAEProblem(BouncingBallDAE_equations, dx0, x0, tspan, p_is, differential_vars=differential_vars, callback=CallbackSet())
  # Solve with IDA:)
  solution = solve(problem, IDA())
  return solution
end
