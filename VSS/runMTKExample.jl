#= Demo file for the EDF presentation =#
using Revise
import OM
import OMBackend
"""
Given a set of models and a file  run the models in the file.
"""
function writeSetOfModelsToFile(models, file)
  results = []
  for model in models
    @info "Exporting : $model to file"
    res = OM.exportModelToFile(model, "$(file).mo"; targetFile="$file.jl", mode = OMBackend.MTK_MODE)
    push!(results, res)
  end
  @info "Done exporting models to file"
end
tst = ["ElectricalComponentTest.SimpleCircuit"]
F = "ElectricalComponentTest"
@info "Write model to file"
@time writeSetOfModelsToFile(tst, F)
@info "Including the Simple circuit"
@time include("ElectricalComponentTest.jl")

@info "Simulating the simple circuit"
@time res = solve(ElectricalComponentTest__SimpleCircuitModel_problem; tspan=(0.0, 1.0));
@info "Simulating the simple circuit. Note max(Δt) = 0.1"
@time res = solve(ElectricalComponentTest__SimpleCircuitModel_problem; tspan=(0.0, 10.0), dtmax=0.1);
@info "Plotting the circuit"
plot(res; legend=false)
