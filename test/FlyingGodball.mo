model FlyingBall
  parameter Real e=0.7;
  parameter Real g=9.81;
  Real h(start = 1);
  Real v;
equation 
  der(h) = v;
  der(v) = g;
end FlyingBall;