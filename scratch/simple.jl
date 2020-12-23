using DiffEqBase
using DifferentialEquations
using Plots
using Sundials
using Revise

function SimpleCirquitStartConditions(p, t0)  
  local x0 = zeros(2)
  local dx0 = zeros(2)
  A = 12.0;
  C = 1.0;
  L = 0.01;
  R1 = 10.0;
  R2 = 100.0;
  f = 1.0;
   
  i2 = x0[1]
  u2 = x0[2]

  u3 = R2 * i2
  u = A * sin(2 * 3.1415 * f * t0);
  u1 = u - u2 
  i1 = u1/R1 
  i = i2 + i1
  u4 = u3 - u
  
  dx0[1] = u4/L
  dx0[2] = i1/C
  return (x0, dx0)
end
function SimpleCirquitDifferentialVars()
  return Bool[1, 1]
end

global counter = 0 

function SimpleCirquitDAE_equations(res, dx, x, p, t)
  global counter += 1
  A = 12.0;
  C = 1.0;
  L = 0.01;
  R1 = 10.0;
  R2 = 100.0;
  f = 1.0;
   
  i2 = x[1]
  u2 = x[2]

  u3 = R2 * i2
  u = A * sin(2 * 3.1415 * f * t);
  u1 = u - u2 
  i1 = u1/R1 
  i = i2 + i1
  u4 = u3 - u
  @info "***************"
  @info "res-1: $(res[1])"
  @info "res-2: $(res[2])"
  @info "x1: $(x[1])"
  @info "x2: $(x[2])"
  @info "dx: $dx"
  @info "x: $x"
  # @info "time: $t"
  #= Modifications =#
 # dx[1] = i1/C
  #= Modifications =#
  res[2] = C * dx[2] - i1 #dx[2] = der(u2)
  res[1] = L * dx[1] - u4 #dx[1] #der(i2)
  if counter == 4
    throw("STOP!")
  end
end

function SimpleCirquitParameterVars()  
  p = Array{Float64}(undef, 6)
  p[1] =  1.0
  p[4] = 1.0
  p[3] = 0.01
  p[2] = 12.0
  p[6] = 10.0  
  p[5] = 100.0
  return p
end

function SimpleCirquitCallbackSet(p)        
  return CallbackSet()
end    

function SimpleCirquitSimulate(tspan = (0.0, 1.0))
  p = SimpleCirquitParameterVars()
  (x0, dx0) = SimpleCirquitStartConditions(p, tspan[1])
  differential_vars = SimpleCirquitDifferentialVars()
  problem = DAEProblem(
    SimpleCirquitDAE_equations,
    dx0,
    x0,
    tspan,
    p;
    differential_vars = differential_vars,
    callback = SimpleCirquitCallbackSet(p),
  )
  
  solution = solve(problem, IDA())
  
  return solution
end
