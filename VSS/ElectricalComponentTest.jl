begin
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations
    function ElectricalComponentTest__SimpleCircuitModel(tspan = (0.0, 1.0))
        @variables t        #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:98 =#
        parameters = ModelingToolkit.@parameters((R1_R, C_C, R2_R, L_L, AC_A, AC_w)) #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:99 =#
        vars = ModelingToolkit.@variables((
            C_v(t),
            L_i(t),
            R1_v(t),
            R1_i(t),
            R1_p_v(t),
            R1_p_i(t),
            R1_n_v(t),
            R1_n_i(t),
            C_i(t),
            C_p_v(t),
            C_p_i(t),
            C_n_v(t),
            C_n_i(t),
            R2_v(t),
            R2_i(t),
            R2_p_v(t),
            R2_p_i(t),
            R2_n_v(t),
            R2_n_i(t),
            L_v(t),
            L_p_v(t),
            L_p_i(t),
            L_n_v(t),
            L_n_i(t),
            AC_v(t),
            AC_i(t),
            AC_p_v(t),
            AC_p_i(t),
            AC_n_v(t),
            AC_n_i(t),
            G_p_v(t),
            G_p_i(t),
        )) #= C:\Users\johti17\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\MTK_CodeGeneration.jl:106 =#
        der = Differential(t)
        eqs = [
            der(C_v) ~ C_i / C_C,
            der(L_i) ~ L_v / L_L,
            0 ~ (R2_p_v + -AC_p_v) - 0.0,
            0 ~ (R2_p_v + -R1_p_v) - 0.0,
            0 ~ (C_p_v + -R1_n_v) - 0.0,
            0 ~ (R2_n_v + -L_p_v) - 0.0,
            0 ~ (AC_n_v + -C_n_v) - 0.0,
            0 ~ (AC_n_v + -L_n_v) - 0.0,
            0 ~ (AC_n_v + -G_p_v) - 0.0,
            0 ~ (R1_n_i + C_p_i) - 0.0,
            0 ~ (L_p_i + R2_n_i) - 0.0,
            0 ~ (((G_p_i + C_n_i) + L_n_i) + AC_n_i) - 0.0,
            0 ~ ((R1_p_i + R2_p_i) + AC_p_i) - 0.0,
            0 ~ R1_R * R1_i - R1_v,
            0 ~ R1_v - (R1_p_v - R1_n_v),
            0 ~ 0.0 - (R1_p_i + R1_n_i),
            0 ~ R1_i - R1_p_i,
            0 ~ C_v - (C_p_v - C_n_v),
            0 ~ 0.0 - (C_p_i + C_n_i),
            0 ~ C_i - C_p_i,
            0 ~ R2_R * R2_i - R2_v,
            0 ~ R2_v - (R2_p_v - R2_n_v),
            0 ~ 0.0 - (R2_p_i + R2_n_i),
            0 ~ R2_i - R2_p_i,
            0 ~ L_v - (L_p_v - L_n_v),
            0 ~ 0.0 - (L_p_i + L_n_i),
            0 ~ L_i - L_p_i,
            0 ~ AC_v - AC_A * sin(AC_w * t),
            0 ~ AC_v - (AC_p_v - AC_n_v),
            0 ~ 0.0 - (AC_p_i + AC_n_i),
            0 ~ AC_i - AC_p_i,
            0 ~ G_p_v - 0.0,
        ]
        nonLinearSystem = ModelingToolkit.ODESystem(
            eqs,
            t,
            vars,
            parameters,
            name = :($(Symbol("ElectricalComponentTest.SimpleCircuit"))),
        )
        pars = Dict(
            R1_R => float(10.0),
            C_C => float(0.01),
            R2_R => float(100.0),
            L_L => float(0.1),
            AC_A => float(1.0),
            AC_w => float(1.0),
        )
        initialValues = [
            R1_v => 0.0,
            R1_i => 0.0,
            R1_p_v => 0.0,
            R1_p_i => 0.0,
            R1_n_v => 0.0,
            R1_n_i => 0.0,
            C_i => 0.0,
            C_p_v => 0.0,
            C_p_i => 0.0,
            C_n_v => 0.0,
            C_n_i => 0.0,
            R2_v => 0.0,
            R2_i => 0.0,
            R2_p_v => 0.0,
            R2_p_i => 0.0,
            R2_n_v => 0.0,
            R2_n_i => 0.0,
            L_v => 0.0,
            L_p_v => 0.0,
            L_p_i => 0.0,
            L_n_v => 0.0,
            L_n_i => 0.0,
            AC_v => 0.0,
            AC_i => 0.0,
            AC_p_v => 0.0,
            AC_p_i => 0.0,
            AC_n_v => 0.0,
            AC_n_i => 0.0,
            G_p_v => 0.0,
            G_p_i => 0.0,
            C_v => 0.0,
            L_i => 0.0,
        ]
        firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
        reducedSystem = firstOrderSystem
        local event_p = [10.0, 0.01, 100.0, 0.1, 1.0, 1.0]
        local discreteVars = collect(values(Dict([])))
        local event_vars = vcat(
            collect(
                values(
                    Dict([
                        R1_v => 0.0,
                        R1_i => 0.0,
                        R1_p_v => 0.0,
                        R1_p_i => 0.0,
                        R1_n_v => 0.0,
                        R1_n_i => 0.0,
                        C_i => 0.0,
                        C_p_v => 0.0,
                        C_p_i => 0.0,
                        C_n_v => 0.0,
                        C_n_i => 0.0,
                        R2_v => 0.0,
                        R2_i => 0.0,
                        R2_p_v => 0.0,
                        R2_p_i => 0.0,
                        R2_n_v => 0.0,
                        R2_n_i => 0.0,
                        L_v => 0.0,
                        L_p_v => 0.0,
                        L_p_i => 0.0,
                        L_n_v => 0.0,
                        L_n_i => 0.0,
                        AC_v => 0.0,
                        AC_i => 0.0,
                        AC_p_v => 0.0,
                        AC_p_i => 0.0,
                        AC_n_v => 0.0,
                        AC_n_i => 0.0,
                        G_p_v => 0.0,
                        G_p_i => 0.0,
                        C_v => 0.0,
                        L_i => 0.0,
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
            callback = ElectricalComponentTest__SimpleCircuitCallbackSet(aux),
        )
        return (problem, initialValues, reducedSystem, tspan, pars, vars)
    end
    begin
        var"saved_values_ElectricalComponentTest.SimpleCircuit" =
            SavedValues(Float64, Tuple{Float64,Array})
        function ElectricalComponentTest__SimpleCircuitCallbackSet(aux)
            local p = aux[1]
            local reals = aux[2]
            begin
                savingFunction(u, t, integrator) =
                    let
                        (t, deepcopy(integrator.p))
                    end
                cb1 = SavingCallback(
                    savingFunction,
                    var"saved_values_ElectricalComponentTest.SimpleCircuit",
                )
            end
            return CallbackSet(cb1)
        end
    end
    (ElectricalComponentTest__SimpleCircuitModel_problem, _, _, _, _, _) =
        ElectricalComponentTest__SimpleCircuitModel()
    function ElectricalComponentTest__SimpleCircuitSimulate(tspan)
        return solve(ElectricalComponentTest__SimpleCircuitModel_problem, tspan = tspan)
    end
    function ElectricalComponentTest__SimpleCircuitSimulate(
        tspan = (0.0, 1.0);
        solver = Rodas5(),
    )
        return solve(
            ElectricalComponentTest__SimpleCircuitModel_problem,
            tspan = tspan,
            solver,
        )
    end
end