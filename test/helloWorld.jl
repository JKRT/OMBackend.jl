
using DifferentialEquations
using Sundials
using Plots

"""
DAE equation for HelloWord.mo
    der(x) = a*x

"""
function helloWord_DAE_equations(res, dx #=The state derivatives =#, x #= State & alg variables =#, p, t #=time=#)
  res[1] = p[1]*x[1] - dx[1]
end


"""
true if a certain index has a state or not in the system of equations.
"""
function helloWorld_differential_vars()
  return [true]
end

"""
For all states give them initial values:
For instance given 10 state variables,
"""
function helloWorld_startConditions(p, t0)
  x0 = [1.0]
  dx0 = [p[1]*x0[1]]
  return x0, dx0
end

"""
Assign the set values to the parameters
"""
function helloWorld_parameter_vars()
  return p=[0.5]
end


"""
Simulate helloWord

# To test run
```julia-repl
julia> tspan = (0.0, 1.0)
julia> solution = helloWorld_simulate(tspan)
julia> plotSolution(solution, tspan)
```
"""
function helloWorld_simulate(tspan = (0.0, 1.0))
  # Define problem
  p_is = helloWorld_parameter_vars()
  println(p_is)
  (x0, dx0) = helloWorld_startConditions(p_is, tspan[1])
  println(x0, dx0)
  differential_vars = helloWorld_differential_vars()
  println(differential_vars)
  #= Pass the residual equations =#
  problem = DAEProblem(helloWord_DAE_equations, dx0, x0, tspan, p_is, differential_vars=differential_vars)
  # Solve with IDA:)
  solution = solve(problem, IDA())
  return solution
end

plot(helloWorld_simulate())
