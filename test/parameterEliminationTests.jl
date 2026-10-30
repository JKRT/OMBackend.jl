#= A read by a literal subscript names one element: an enumeration's or a Boolean's literal (`p[E.b]`, `q[true]`) as an
   Integer's, so the parameters it reads are folded (SimCode eliminateConstantParameters; only an array read by a
   subscript that reads a variable keeps its elements, dynamicallyIndexedElements!). A small parameter array read by a
   variable subscript (`r[i]`) is lowered to a selection over `i` before: its value checked. Models/EnumIndexedParameters.mo. =#
using Test
import OMBackend
import OMFrontend

const EIP_SOLUTION = OMBackend.CodeGeneration.ModelingToolkit.SciMLBase.AbstractODESolution

"The simulation code and the final value of `v` of `model` of Models/EnumIndexedParameters.mo, simulated to 1."
function _eipRun(model::String, v::String)::Tuple{OMBackend.SimCode.SIM_CODE, Float64}
  local file::String = joinpath(@__DIR__, "Models", "EnumIndexedParameters.mo")
  local (fm::OMFrontend.Frontend.FlatModel, cache::OMFrontend.Frontend.FunctionTreeImpl.Tree) =
    OMFrontend.instantiateSCodeToFM(model, OMFrontend.translateToSCode(OMFrontend.parseFile(file)))
  OMBackend.translate(fm; functionList = OMFrontend.cacheToFunctionList(cache))
  local simCode::OMBackend.SimCode.SIM_CODE = OMBackend.LAST_SIM_CODE[]
  local sol::EIP_SOLUTION = OMBackend.simulateModel(OMBackend.canonicalName(model); tspan = (0.0, 1.0))
  return (simCode, OMBackend.getVariableValues(sol, v)[end])
end

"The keys of the parameters of the simulation code `sc`, sorted."
_eipParameters(sc::OMBackend.SimCode.SIM_CODE)::Vector{String} =
  sort(String[k for (k, (_, v)) in sc.stringToSimVarHT if v.varKind isa OMBackend.SimCode.PARAMETER])

@testset "Parameters read by literal and by variable subscripts" begin
  local (literal::OMBackend.SimCode.SIM_CODE, x::Float64) = _eipRun("EnumIndexedParameters.Read", "x")
  @test isempty(_eipParameters(literal))
  @test x ≈ 12.0
  local (_, y::Float64) = _eipRun("EnumIndexedParameters.Dynamic", "y")
  @test isapprox(y, 2.0; atol = 1e-6)
end
