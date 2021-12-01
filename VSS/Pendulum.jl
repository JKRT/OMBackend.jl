begin
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations
    function PendulumModel(tspan = (0.0, 1.0))
        @variables t        #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:98 =#
        parameters = ModelingToolkit.@parameters((x0, y0, g, L)) #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:99 =#
        vars = ModelingToolkit.@variables((x(t), y(t), phi(t), phid(t), vx(t), vy(t))) #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:106 =#
        der = Differential(t)
        eqs = [
            der(x) ~ vx,
            der(y) ~ vy,
            der(phi) ~ phid,
            der(phid) ~ -((g * sin(phi)) / L),
            0 ~ x - L * sin(phi),
            0 ~ y - -(L * cos(phi)),
        ]
        nonLinearSystem = ModelingToolkit.ODESystem(
            eqs,
            t,
            vars,
            parameters,
            name = :($(Symbol("Pendulum"))),
        )
        pars = Dict(
            x0 => float(10.0),
            y0 => float(10.0),
            g => float(9.81),
            L => float(sqrt(x0^2.0 + y0^2.0)),
        )
        initialValues =
            [vx => 0.0, vy => 0.0, x => pars[x0], y => pars[y0], phi => 1.0, phid => 0.0]
        firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
        reducedSystem = ModelingToolkit.dae_index_lowering(firstOrderSystem)
        local event_p = [10.0, 10.0, 9.81, 0]
        local discreteVars = collect(values(Dict([])))
        local event_vars = vcat(
            collect(
                values(
                    Dict([
                        vx => 0.0,
                        vy => 0.0,
                        x => pars[x0],
                        y => pars[y0],
                        phi => 1.0,
                        phid => 0.0,
                    ]),
                ),
            ),
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
            callback = PendulumCallbackSet(aux),
        )
        return (problem, initialValues, reducedSystem, tspan, pars, vars)
    end
    begin
        saved_values_Pendulum = SavedValues(Float64, Tuple{Float64,Array})
        function PendulumCallbackSet(aux)
            local p = aux[1]
            local reals = aux[2]
            begin
                savingFunction(u, t, integrator) =
                    let
                        (t, deepcopy(integrator.p))
                    end
                cb1 = SavingCallback(savingFunction, saved_values_Pendulum)
            end
            return CallbackSet(cb1)
        end
    end
    (PendulumModel_problem, _, _, _, _, _) = PendulumModel()
    function PendulumSimulate(tspan)
        return solve(PendulumModel_problem, tspan = tspan)
    end
    function PendulumSimulate(tspan = (0.0, 1.0); solver = Rodas5())
        return solve(PendulumModel_problem, tspan = tspan, solver)
    end
end