#= Breaking Pendulum mockup =#
include("FreeFall.jl")
include("Pendulum.jl")

function BreakingPendulium()
    #= 
      1. Start by simulating the Pendulum model 
      Add a callback to the Pendulum model. 
      This callback specifices that after 5 seconds the Pendulum breaks.
      
      2. When this stuructural change happens terminate integration. 
         reflatten equations 

       3.   Rstart the simulation
             If the same number of equations {
                 reinit (Like when equations)
             }
             else {
                We need to save the values of the relevant variables and start solving from this point.
             } 
    =#

    #= Since we know the Breaking Pendulum has 2 stuructural components we to create them =#
    FreeFall = FreeFallModel()
    Pendulum = PendulumModel()
    #= It is assumed that we have statically examined that the rules for structural holds =#
    
end
