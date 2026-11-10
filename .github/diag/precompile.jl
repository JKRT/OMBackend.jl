#=
  Pkg.precompile() of the active environment, as one step of the Windows hang diagnostic (branch
  diag/windows-ext-hang only). Run under driver.jl.

    DIAG_PRECOMPILE_TWO_THREADS   the workers run with -t1,1 (a second thread for the workload's watchdog)
=#
import Pkg

"Precompile the environment; a failure is printed and decides the exit code."
function precompileEnvironment()::Int
  try
    if haskey(ENV, "DIAG_PRECOMPILE_TWO_THREADS")
      Pkg.precompile(; configs = `-t1,1` => Base.CacheFlags())
    else
      Pkg.precompile()
    end
  catch err
    println("[precompile] failed:")
    showerror(stdout, err)
    println()
    return 1
  end
  return 0
end

exit(precompileEnvironment())
