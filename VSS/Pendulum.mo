model Pendulum
  parameter Real x0 = 10;
  parameter Real y0 = 10;
  parameter Real g = 9.81;
  parameter Real length = sqrt(x0^2 + y0^2);
  /* Common variables */
  Real x(start = x0);
  Real y(start = y0);
  Real vx;
  Real vy;
 /* Model specific variables */
  Real phi(start = 1., fixed = true);
  Real phid;
equation
  der(phi) = phid;
  der(x) = vx;
  der(y) = vy;
  x = length * sin(phi);
  y = -length * cos(phi);
  der(phid) =  -g / length * sin(phi);
end Pendulum;