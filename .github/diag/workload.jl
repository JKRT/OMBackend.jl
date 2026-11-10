#=
  The extension's translate workload in a normal process, translate only, under the watchdog (branch
  diag/windows-ext-hang only; run when the precompile stalled inside the workload). Run under driver.jl with a
  second thread for the watchdog:

    OMBACKEND_DIAG_WATCHDOG=300,420 julia -t1,1 --project=<env> workload.jl
=#
include(joinpath(@__DIR__, "..", "..", "ext", "OMBackendOMFrontendExt", "diag.jl"))
diagProcess("script: start")
import OMBackend
diagLibraries("script: imported OMBackend")
import OMFrontend
diagLibraries("script: imported OMFrontend")

"Run the workload of the loaded extension with its own watchdog (its stage lines name the model and the phase)."
function runWorkload()::Nothing
  local ext::Union{Nothing, Module} = Base.get_extension(OMBackend, :OMBackendOMFrontendExt)
  ext === nothing && error("OMBackendOMFrontendExt is not loaded")
  local watchdog = ext.diagWatchdogStart()
  try
    ext._precompileWarmupTranslate(; simulate = false)
  finally
    ext.diagWatchdogStop(watchdog)
  end
  return nothing
end

runWorkload()
diagStage("script: end")
