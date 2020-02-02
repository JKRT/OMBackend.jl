
using DifferentialEquations
using Sundials

"""
DAE equation for HelloWord.mo
    der(x) = a*x

"""
function helloWord_DAE_equations(res, dx, x, p, t)
  res[1] = p[1]*x[1] - dx[1]
end


"""
Boolean array telling if variable of array x is a state
"""
function helloWorld_differential_vars()
  return [true]
end


function helloWorld_startConditions()
  x0 = [1.0, 0, 0]
  dx0 = [-0.04, 0.04, 0.0]

  return x0, dx0
end

function helloWorld_simulate(tspan = (0.0, 1.0))

  # Define problem
  (x0, dx0) = helloWorld_startConditions()
  differential_vars = helloWorld_differential_vars()
  probblem = DAEProblem(f, dx0, x0, tspan, differential_vars=differential_vars)

  # Solve with IDA
  solution = solve(problem, IDA())

  return solution
end
