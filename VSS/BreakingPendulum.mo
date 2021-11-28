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

model Pendulum
  parameter Real x0 = 10;
  parameter Real y0 = 10;
  parameter Real g = 9.81;
  parameter Real L = sqrt(x0^2 + y0^2);
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
  x = L * sin(phi);
  y = -L * cos(phi);
  der(phid) =  -g / L * sin(phi);
end Pendulum;

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
end BouncingBall

model BreakingPendulum
/* Common variables */
Real x;
Real y;
Real vx;
Real vy;
/*
  The  structure keyword
(structure requires that variables declared in a model
that contain them must also be existant in the submodules).
This means that the variables present in the breaking pendulum
need also to be present in the different modes the system can be in.
*/
structural FreeFall freeFall;
structural BouncingBall bouncingBall;
structural Pendulum pendulum;

equation
/* New Modelica construct the switch equation */
   switch-when t == 5
     case Pendulum then FreeFall;
     case _  FreeFall; //Or bouncing ball.
   end switch-when;
end BreakingPendulum;
