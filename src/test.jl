
  using DifferentialEquations
  using Sundials
  using Plots

    function testStartConditions(p, t0)
      x0 = [1.0]
      dx0 = [p[1]*x0[1]]
      return x0, dx0
    end
  function testDifferentialVars()
      return Any[true]
     end
    function testDAE_equations(res, dx #=The state derivatives =#, x #= State & alg variables =#, p, t #=time=#)
      res[1] = dx[1] - (- p[1]) * x[1]

     end

function testParameterVars
  return p=[0.5]
end

 function testSimulate(tspan = (0.0, 1.0))
  # Define problem
  p_is = testparameterVars()
  (x0, dx0) =testStartConditions(p_is, tspan[1])
  differential_vars = testDifferential_vars()
  #= Pass the residual equations =#
  problem = DAEProblem(testDAE_equations, dx0, x0, tspan, p_is, differential_vars=differential_vars)
  # Solve with IDA:)
  solution = solve(problem, IDA())
  return solution
end
