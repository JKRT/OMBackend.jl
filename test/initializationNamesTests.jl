#= The initialization matches parameters and pins to the variables by their names (`variableName`), not their printing:
   an array element's symbol prints as `var"k[1]"`, SimCode names it `k[1]` (Models/InitializationNames.mo). =#
using Test
import OMBackend
import OMFrontend

@testset "The initialization's variables by name" begin
  local file::String = joinpath(@__DIR__, "Models", "InitializationNames.mo")
  local (fm, cache) = OMFrontend.instantiateSCodeToFM("InitializationNames.FreeArrayParameter",
                                                      OMFrontend.translateToSCode(OMFrontend.parseFile(file)))
  OMBackend.translate(fm; functionList = OMFrontend.cacheToFunctionList(cache))
  local sol = OMBackend.simulateModel(OMBackend.canonicalName("InitializationNames.FreeArrayParameter");
                                      tspan = (0.0, 1.0), reltol = 1e-8, abstol = 1e-8)
  @test sol.retcode == OMBackend.CodeGeneration.ModelingToolkit.SciMLBase.ReturnCode.Success
  #= k = {3, 4} from y = k .* x = {3, 8} at x = {1, 2}: x = {exp(-3 t), 2 exp(-4 t)} =#
  @test OMBackend.getVariableValues(sol, "x[1]")[end] ≈ exp(-3) rtol = 1e-6
  @test OMBackend.getVariableValues(sol, "x[2]")[end] ≈ 2 * exp(-4) rtol = 1e-6
end
