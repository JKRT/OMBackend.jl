model FreeFall
  Real x;
  Real y;
  Real vx;
  Real vy;
equation
  der(x) = vx;
  der(y) = vy;
  der(vx) = 0.0;
  der(vy) = -9.81;
end FreeFall;
