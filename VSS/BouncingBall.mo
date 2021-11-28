model BouncingBall
  parameter Real e=0.7;
  parameter Real g=9.81;
  Real x;   
  Real y(start = 1.0);
  Real vx;
  Real vy;
equation
  der(x) = vx;
  der(y) = vy;
  der(vy) = -g;
  der(vx) = 0.0;  
  when y <= 0 then
    reinit(vy, -e*pre(vy));
  end when;
end BouncingBall;