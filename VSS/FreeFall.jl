begin
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations
    function FreeFallModel(tspan = (0.0, 1.0))
        @variables t        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:101 =#
        parameters = ModelingToolkit.@parameters((g,)) #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:102 =#
        vars = ModelingToolkit.@variables((x(t), y(t), vx(t), vy(t))) #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:109 =#
        der = Differential(t)
        eqs = [der(x) ~ vx, der(y) ~ vy, der(vx) ~ 0.0, der(vy) ~ -g]
        nonLinearSystem = ModelingToolkit.ODESystem(
            eqs,
            t,
            vars,
            parameters,
            name = :($(Symbol("FreeFall"))),
        )
        pars = Dict(g => float(9.81))
        initialValues = [x => 0.0, y => 0.0, vx => 0.0, vy => 0.0]
        firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
        reducedSystem = firstOrderSystem
        local event_p = [9.81]
        local discreteVars = collect(values(Dict([])))
        local event_vars = vcat(
            collect(values(Dict([x => 0.0, y => 0.0, vx => 0.0, vy => 0.0]))),
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
            callback = FreeFallCallbackSet(aux),
        )
        return (problem, initialValues, reducedSystem, tspan, pars, vars)
    end
    begin
        saved_values_FreeFall = SavedValues(Float64, Tuple{Float64,Array})
        function FreeFallCallbackSet(aux)
            local p = aux[1]
            local reals = aux[2]
            begin
                savingFunction(u, t, integrator) =
                    let
                        (t, deepcopy(integrator.p))
                    end
                cb1 = SavingCallback(savingFunction, saved_values_FreeFall)
            end
            return CallbackSet(cb1)
        end
    end
    (FreeFallModel_problem, _, _, _, _, _) = FreeFallModel()
    function FreeFallSimulate(tspan)
        return solve(FreeFallModel_problem, tspan = tspan)
    end
    function FreeFallSimulate(tspan = (0.0, 1.0); solver = Rodas5())
        return solve(FreeFallModel_problem, tspan = tspan, solver)
    end
end