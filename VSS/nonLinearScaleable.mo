model nonlinearScalable
  parameter Integer N=10;
  parameter Real a=0.5;
  Real x[N](each start=2.5);
  Real y(start=0, fixed=true);
equation
  for i in 1:N loop
    N+1 = exp(time*i*a+x[i]) + sum(x);
  end for;
  der(y) = sum(x)*time;
end nonlinearScalable;
