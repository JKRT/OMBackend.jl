begin
    

    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:166 =#
    using DiffEqBase
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:167 =#
    using DifferentialEquations
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:168 =#
    using Plots
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:169 =#
    using Sundials
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:170 =#
    begin
        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:112 =#
        function BouncingBallRealsStartConditions(p, t0)
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:112 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:113 =#
            local x0 = zeros(2)
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:114 =#
            local dx0 = zeros(2)
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:115 =#
            begin
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:207 =#
                begin
                    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:463 =#
                end
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:208 =#
                begin
                    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:495 =#
                    begin
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:479 =#
                        #= h:479 =#
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:480 =#
                        x0[2] = begin
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:374 =#
                                1.0
                            end
                    end
                end
            end
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:116 =#
            return (x0, dx0)
        end
    end
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:171 =#
    begin
        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:125 =#
        begin
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:126 =#
            function BouncingBallRealsDifferentialVars()
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:126 =#
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:127 =#
                return Bool[1, 1]
            end
        end
    end
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:172 =#
    begin
        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:120 =#
        function BouncingBallRealsDAE_equations(res, dx, x, p, t)
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:120 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:121 =#
            begin
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:256 =#
                res[1] = begin
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:420 =#
                        begin
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\codeGenerationUtil.jl:255 =#
                                dx[1]
                            end - begin
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:412 =#
                                -(begin
                                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:394 =#
                                        #= g parameter:394 =#
                                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:395 =#
                                        p[1]
                                    end)
                            end
                    end
            end
            begin
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:256 =#
                res[2] = begin
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:420 =#
                        begin
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\codeGenerationUtil.jl:255 =#
                                dx[2]
                            end - begin
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:390 =#
                                #= v state:390 =#
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:391 =#
                                x[1]
                            end
                    end
            end
        end
    end
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:173 =#
    begin
        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:132 =#
        function BouncingBallRealsParameterVars()
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:132 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:133 =#
            p = Array{Float64}(undef, 2)
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:134 =#
            begin
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:240 =#
                #= g:240 =#
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:241 =#
                p[1] = begin
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:374 =#
                        9.81
                    end
            end
            begin
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:240 =#
                #= e:240 =#
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:241 =#
                p[2] = begin
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:374 =#
                        0.7
                    end
            end
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:135 =#
            return p
        end
    end
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:174 =#
    begin
        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:139 =#
        function BouncingBallRealsCallbackSet(p)
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:139 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:140 =#
            #= WHEN EQUATIONS:140 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:141 =#
            begin
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:283 =#
                function condition1(x, t, integrator)
                    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:283 =#
                    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:284 =#
                    begin
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:420 =#
                        begin
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:390 =#
                                #= h state:390 =#
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:391 =#
                                x[2]
                            end - begin
                                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:374 =#
                                0.0
                            end
                    end
                end
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:286 =#
                function affect1!(integrator)
                    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:286 =#
                    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:287 =#
                    Expr[quote
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:347 =#
    integrator.u[1] = begin
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:412 =#
            -(begin
                    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:420 =#
                    begin
                            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:394 =#
                            #= e parameter:394 =#
                            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:395 =#
                            p[2]
                        end * begin
                            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\codeGenerationUtil.jl:263 =#
                            integrator.u[1]
                        end
                end)
        end
end]
                end
                #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:289 =#
                cb1 = ContinuousCallback(condition1, affect1!, rootfind = true, save_positions = (false, true), affect_neg! = affect1!)
            end
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:142 =#
            #= IF EQUATIONS:142 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:143 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:144 =#
            return CallbackSet(begin
                        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:193 =#
                        begin
                            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:186 =#
                            cb1
                        end
                    end)
        end
    end
    #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:175 =#
    begin
        #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:148 =#
        function BouncingBallRealsSimulate(tspan = (0.0, 1.0))
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:148 =#
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:150 =#
            p = BouncingBallRealsParameterVars()
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:151 =#
            (x0, dx0) = BouncingBallRealsStartConditions(p, tspan[1])
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:152 =#
            differential_vars = BouncingBallRealsDifferentialVars()
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:154 =#
            problem = DAEProblem(BouncingBallRealsDAE_equations, dx0, x0, tspan, p, differential_vars = differential_vars, callback = BouncingBallRealsCallbackSet(p))
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:158 =#
            solution = solve(problem, IDA())
            #= C:\Users\John\Projects\Programming\JuliaPackages\OM.jl\OMBackend.jl\src\CodeGeneration\CodeGeneration.jl:159 =#
            return solution
        end
    end
end
