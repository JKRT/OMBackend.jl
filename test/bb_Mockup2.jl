using DifferentialEquations
using DiffEqBase
using Plots
using Sundials
#using Revise

function BouncingBallStartConditions(p, t0)
    local x0 = Array{Float64}(undef, 5)
    local dx0 = Array{Float64}(undef, 5)
    x0[5] = (true) #= flying =#
    x0[2] = (1.0) #= h =#
    return x0, dx0
end

function BouncingBallDifferentialVars()
    return Bool[1, 1, 0, 0, 0]
end

function BouncingBallDAE_equations(res, dx, x, p, t) #=time=#
    res[1] = ((dx[2]) - (x[1])) #= v =#
    res[2] = ((dx[1]) - (p[1])) #= ifEq_tmp0 =#
    res[3] = ((x[4]) - ((x[2]) <= (0.0))) #= h =#
end

function BouncingBallParameterVars()
    p = Array{Float64}(undef, 2)
    p[1] = (9.81) #= g =#
    p[2] = (0.7) #= e =#
    return p
end

function BouncingBallCallbackSet(p)#=parameters=#
    #= WHEN EQUATIONS. Discrete callbacks =#
    condition1(x, t, integrator) =
        Bool(floor((reduce(|, [(((x[2]) <= (0.0)) && ((x[1]) <= (0.0))), Bool(floor(abs((x[4]))))])))) #= impact =#
    affect1!(integrator) =
        let
            #v = (integerator.u[3]) #= v_new =#
            (integrator.u[5]) = ((integrator.u[3]) > (0.0)) #= v_new =#
            (integrator.u[3]) = (
                if ((integerator.u[4]) && !integrator.uprev[4]) #= impact =#
                    ((-((p[2]) * ((integrator.uprev[1])))))
                else #= e =#
                    (0.0)
                end
            )
        end
    cb1 = DiscreteCallback(condition1, affect1!)

    #= IF EQUATIONS =#
    condition2(x, t, integrator) = Bool(floor((x[5]))) #= flying =#
    affect2!(integrator) =
        let
            integrator.u[6] = - ((-(p[1]))) #= g =#

        end
    cb2 = ContinuousCallback(condition2, affect2!)
    condition3(x, t, integrator) = !Bool(floor((x[5]))) #= flying =#
    affect3!(integrator) =
        let
            integrator.u[6] = (0.0) #= ifEq_tmp0 =#

        end
    cb3 = ContinuousCallback(condition3, affect3!)

    return CallbackSet(cb1, cb2, cb3)
end

function BouncingBallSimulate(tspan = (0.0, 3.0))
    # Define problem
    p_is = BouncingBallParameterVars()
    (x0, dx0) = BouncingBallStartConditions(p_is, tspan[1])
    differential_vars = BouncingBallDifferentialVars()
    #= Pass the residual equations =#
    problem = DAEProblem(
        BouncingBallDAE_equations,
        dx0,
        x0,
        tspan,
        p_is,
        differential_vars = differential_vars,
        callback = BouncingBallCallbackSet(p_is),
    )
    # Solve with IDA:)
    solution = solve(problem, IDA())
    return solution
end
