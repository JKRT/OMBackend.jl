model FreeFall
  Real x;
  Real y;
  Real vx;
  Real vy;
  parameter Real g = 9.81;
equation
  der(x) = vx;
  der(y) = vy;
  der(vx) = 0.0;
  der(vy) = -g;
end FreeFall;
