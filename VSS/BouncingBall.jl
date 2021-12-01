begin
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations
    function BouncingBallModel(tspan = (0.0, 1.0))
        @variables t        #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:98 =#
        parameters = ModelingToolkit.@parameters((e, g)) #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:99 =#
        vars = ModelingToolkit.@variables((x(t), y(t), vx(t), vy(t))) #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:106 =#
        der = Differential(t)
        eqs = [der(x) ~ vx, der(y) ~ vy, der(vx) ~ 0.0, der(vy) ~ -g]
        nonLinearSystem = ModelingToolkit.ODESystem(
            eqs,
            t,
            vars,
            parameters,
            name = :($(Symbol("BouncingBall"))),
        )
        pars = Dict(e => float(0.7), g => float(9.81))
        initialValues = [x => 0.0, y => 1.0, vx => 0.0, vy => 0.0]
        firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
        reducedSystem = firstOrderSystem
        local event_p = [0.7, 9.81]
        local discreteVars = collect(values(Dict([])))
        local event_vars = vcat(
            collect(values(Dict([x => 0.0, y => 1.0, vx => 0.0, vy => 0.0]))),
            discreteVars,
        )
        local aux = Array{Array{Float64}}(undef, 2)
        aux[1] = event_p
        aux[2] = event_vars
        problem = ModelingToolkit.ODEProblem(
            reducedSystem,
            initialValues,
            tspan,
            pars,
            callback = BouncingBallCallbackSet(aux),
        )
        return (problem, initialValues, reducedSystem, tspan, pars, vars)
    end
    begin
        saved_values_BouncingBall = SavedValues(Float64, Tuple{Float64,Array})
        function BouncingBallCallbackSet(aux)
            local p = aux[1]
            local reals = aux[2]
            begin
                condition1 = ((x, t, integrator) -> x[2] - 0.0)
                affect1! = (integrator -> integrator.u[4] = -(p[1] * integrator.u[4]))
                cb1 = ContinuousCallback(
                    condition1,
                    affect1!,
                    rootfind = true,
                    save_positions = (true, true),
                    affect_neg! = affect1!,
                )
            end
            begin
                savingFunction(u, t, integrator) =
                    let
                        (t, deepcopy(integrator.p))
                    end
                cb2 = SavingCallback(savingFunction, saved_values_BouncingBall)
            end
            return CallbackSet(cb1, cb2)
        end
    end
    (BouncingBallModel_problem, _, _, _, _, _) = BouncingBallModel()
    function BouncingBallSimulate(tspan)
        return solve(BouncingBallModel_problem, tspan = tspan)
    end
    function BouncingBallSimulate(tspan = (0.0, 1.0); solver = Rodas5())
        return solve(BouncingBallModel_problem, tspan = tspan, solver)
    end
end