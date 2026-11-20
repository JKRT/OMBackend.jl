"""
The omc route: omc's frontend (OMCompiler.jl, the Rust omc loaded into Julia) -> omc's flat model -> OMJLIR
(OMBackendOMCompilerExt) -> OMBackend. Not part of `runtests.jl`: it needs an OpenModelica installation with
`libOpenModelicaCompiler` (a build with `-DOM_OMC_ENABLE_RUST=ON`, or the Rust nightly of build.openmodelica.org).
CI's job `omc-route` runs it.

    OPENMODELICAHOME=<installation> julia --project=<env with OMBackend, OMCompiler, Test> test/omcRouteTests.jl

Small models with analytic values: a state, an event, an algebraic loop, an array, a function.
"""

using Test
import OMBackend
import OMCompiler

const OMHOME::String = get(ENV, "OPENMODELICAHOME", "")
isempty(OMHOME) && error("omcRouteTests: set OPENMODELICAHOME to an OpenModelica installation with libOpenModelicaCompiler")

"""
    RouteCase(source, tspan, checks)

A model's source, its time span, and `checks`: (variable, value at the end of the run).
"""
struct RouteCase
  source::String
  tspan::Tuple{Float64, Float64}
  checks::Vector{Tuple{String, Float64}}
end

const CASES::Vector{Pair{String, RouteCase}} = [
  "Decay" => RouteCase("model Decay Real x(start = 1, fixed = true); equation der(x) = -x; end Decay;",
                       (0.0, 1.0), [("x", exp(-1.0))]),
  #= h(0) = 1, falls freely: the first impact at sqrt(2 / 9.81) = 0.4515 s, then bounces with e = 0.7. =#
  "Ball" => RouteCase("""model Ball
                           parameter Real e = 0.7;
                           Real h(start = 1, fixed = true);
                           Real v(start = 0, fixed = true);
                         equation
                           der(h) = v;
                           der(v) = -9.81;
                           when h <= 0 then
                             reinit(v, -e * pre(v));
                           end when;
                         end Ball;""",
                      (0.0, 0.6), [("v", -0.7 * -9.81 * sqrt(2 / 9.81) - 9.81 * (0.6 - sqrt(2 / 9.81)))]),
  "Loop" => RouteCase("model Loop Real x; Real y; equation x + y = 2 + time; x - y = sin(time); end Loop;",
                      (0.0, 1.0), [("x", (3 + sin(1.0)) / 2), ("y", (3 - sin(1.0)) / 2)]),
  "Rates" => RouteCase("""model Rates
                            parameter Real k[3] = {1, 2, 3};
                            Real x[3](each start = 1, each fixed = true);
                          equation
                            der(x) = -k .* x;
                          end Rates;""",
                       (0.0, 1.0), [("x[1]", exp(-1.0)), ("x[2]", exp(-2.0)), ("x[3]", exp(-3.0))]),
  "Square" => RouteCase("""package Square
                             function sq
                               input Real u;
                               output Real y;
                             algorithm
                               y := u * u;
                             end sq;
                             model M
                               Real x(start = 0, fixed = true);
                             equation
                               der(x) = sq(time);
                             end M;
                           end Square;""",
                        (0.0, 1.0), [("x", 1 / 3)]),
]

"""
The class to instantiate for a case: the model, or `M` inside a package.
"""
function className(name::String, source::String)::String
  return startswith(source, "package") ? name * ".M" : name
end

@testset "omc route" begin
  OMCompiler.init(omhome = OMHOME)
  for (name::String, case::RouteCase) in CASES
    @testset "$name" begin
      local loaded::String = OMCompiler.evaluate("loadString(\"" * replace(case.source, "\"" => "\\\"") * "\")")
      @test strip(loaded) == "true"
      local cls::String = className(name, case.source)
      local (model::OMCompiler.Node, functions::Vector) = OMCompiler.instantiate(cls)
      local ir::OMBackend.OMJLIR.Model = OMBackend.lowerOMCFlatModel(model, functions)
      OMBackend.translate(ir)
      local sol = OMBackend.simulate(cls; tspan = case.tspan, reltol = 1e-8, abstol = 1e-10)
      @test string(sol.retcode) == "Success"
      for (variable::String, expected::Float64) in case.checks
        local values::Vector{Float64} = OMBackend.variableValues(sol, variable)
        @test isapprox(values[end], expected; rtol = 1e-5, atol = 1e-8)
      end
    end
  end
end
