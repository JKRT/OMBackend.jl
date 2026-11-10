#=
  Diagnostics of the Windows hang while OMBackendOMFrontendExt precompiles (branch diag/windows-ext-hang only,
  not for merging). Everything is off unless OMBACKEND_DIAG is set. Base only: the file is included at the top of
  the extension, before its imports, and by the scripts of .github/diag.

    OMBACKEND_DIAG                     stage lines on (stderr, and the file OMBACKEND_DIAG_LOG when set)
    OMBACKEND_DIAG_BACKEND_FIRST       the extension imports OMBackend before OMFrontend
    OMBACKEND_DIAG_SKIP_EXT_WORKLOAD   the extension's precompile workload is skipped (the core's stays)
    OMBACKEND_DIAG_ERROR_MODE          Windows: SetErrorMode(0x8003) before the imports (no hard-error popups,
                                       as test/runtests.jl does): a library that fails to load is an error
    OMBACKEND_DIAG_WATCHDOG=a,b        the workload's watchdog: backtraces of all tasks after a seconds, again
                                       after b seconds and then the process ends (needs a second thread)
=#

"The last stage line of this process (the watchdog prints it)."
const DIAG_LAST_STAGE = Ref{String}("")

"True when the diagnostics are on (read at each call: a constant would be kept in the package image)."
diagOn()::Bool = haskey(ENV, "OMBACKEND_DIAG")

"""
    diagLine(msg::String)::String

Write one diagnostic line to stderr, without Julia's IO (usable from the watchdog's thread), and append it to the
file OMBACKEND_DIAG_LOG: Pkg keeps a precompile worker's output until the worker ends, the file is there at once.
"""
function diagLine(msg::String)::String
  local line::String = string("[diag ", Libc.strftime("%H:%M:%S", time()), " pid ", getpid(), "] ",
                              replace(msg, r"\r?\n" => " / "))
  ccall(:jl_safe_printf, Cvoid, (Cstring, Cstring...), "%s\n", line)
  local path::String = get(ENV, "OMBACKEND_DIAG_LOG", "")
  if !isempty(path)
    try
      open(path, "a") do io::IO
        println(io, line)
      end
    catch err
      ccall(:jl_safe_printf, Cvoid, (Cstring, Cstring...), "%s\n", string("[diag] no stage file: ", err))
    end
  end
  return line
end

"A stage of the load or of the workload: the line that was reached last tells where a hang is."
function diagStage(msg::String)::Nothing
  diagOn() || return nothing
  DIAG_LAST_STAGE[] = diagLine(msg)
  return nothing
end

"""
    diagLibraries(at::String)::Nothing

A stage line with the packages of the DLL clash that are loaded (Plots and Glib against OMRuntimeExternalC) and
the gettext, iconv, Glib and Modelica runtime libraries in the process, with their paths: which libintl won.
"""
function diagLibraries(at::String)::Nothing
  diagOn() || return nothing
  local watched::Tuple{Vararg{Symbol}} = (:OMBackend, :OMFrontend, :OMParser, :Plots, :GR, :GR_jll, :Glib_jll,
                                          :Gettext_jll, :GettextRuntime_jll, :Libiconv_jll, :OMRuntimeExternalC)
  local loaded::Vector{Symbol} = Symbol[nameof(m) for m::Module in Base.loaded_modules_array() if nameof(m) in watched]
  local libraries::Vector{String} = String[]
  try
    for path::String in Base.Libc.Libdl.dllist()
      if occursin(r"intl|iconv|gio-|glib-|gobject-|Modelica|libGR\b|gks"i, basename(path)) ||
         occursin("OMRuntimeExternalC", path)
        push!(libraries, path)
      end
    end
  catch err
    push!(libraries, string("dllist failed: ", err))
  end
  diagStage(string(at, "; loaded: ", isempty(loaded) ? "none" : join(loaded, " "),
                   "; libraries: ", isempty(libraries) ? "none" : join(libraries, " | ")))
  return nothing
end

"Windows' error mode of the process (SEM_* bits; 0x8003 = no hard-error, GP fault or open-file popups), else 0."
function diagErrorMode()::UInt32
  @static if Sys.iswindows()
    return ccall((:GetErrorMode, "kernel32"), stdcall, UInt32, ())
  else
    return UInt32(0)
  end
end

"The stage line of the process: threads, pools, Windows' error mode, and whether Julia is generating an image."
function diagProcess(at::String)::Nothing
  diagOn() || return nothing
  @static if Sys.iswindows()
    if haskey(ENV, "OMBACKEND_DIAG_ERROR_MODE")
      local before::UInt32 = ccall((:SetErrorMode, "kernel32"), stdcall, UInt32, (UInt32,), 0x00008003)
      diagStage(string(at, "; SetErrorMode(0x8003), the mode was 0x", string(before; base = 16)))
    end
  end
  diagStage(string(at, "; threads default=", Threads.nthreads(:default), " interactive=", Threads.nthreads(:interactive),
                   " pool=", Threads.threadpool(), " generating_output=", ccall(:jl_generating_output, Cint, ()),
                   " error mode=0x", string(diagErrorMode(); base = 16), " julia=", VERSION))
  return nothing
end

"The watchdog of a workload: the flag that ends it and its task (nothing: not started)."
struct DiagWatchdog
  done::Threads.Atomic{Bool}
  task::Union{Nothing, Task}
end

"End the process without running exit hooks or unloading libraries (on Windows these wait for the loader lock)."
function diagTerminate(code::Int)::Nothing
  @static if Sys.iswindows()
    ccall((:TerminateProcess, "kernel32"), stdcall, Cint, (Ptr{Cvoid}, Cuint),
          ccall((:GetCurrentProcess, "kernel32"), stdcall, Ptr{Cvoid}, ()), code % Cuint)
  else
    ccall(:_exit, Cvoid, (Cint,), code)
  end
  return nothing
end

"One second off the Julia runtime (GC safe: a collection on the workload's thread does not wait for the sleep)."
function diagSleepSecond()::Nothing
  @static if Sys.iswindows()
    @ccall gc_safe = true Sleep(1000::UInt32)::Cvoid
  else
    @ccall gc_safe = true usleep(1000000::UInt32)::Int32
  end
  return nothing
end

"The last stage and the backtrace of every task that has not finished, straight to stderr."
function diagDump(why::String)::Nothing
  diagLine(string("watchdog: ", why, "; last stage: ", DIAG_LAST_STAGE[]))
  ccall(:jl_print_task_backtraces, Cvoid, (Cint,), 0)
  diagLine("watchdog: backtraces written")
  return nothing
end

"""
    diagWatchdogLoop(done::Threads.Atomic{Bool}, dumpAfter::Float64, exitAfter::Float64)::Nothing

The watchdog's loop, on another thread than the workload: a dump after `dumpAfter` seconds, a second one after
`exitAfter` seconds, then the process ends with code 3. The loop does not allocate between the dumps and does not
yield, so it runs while the workload's thread is stuck in compiled code.
"""
function diagWatchdogLoop(done::Threads.Atomic{Bool}, dumpAfter::Float64, exitAfter::Float64)::Nothing
  local t0::Float64 = time()
  local dumped::Bool = false
  while !done[]
    diagSleepSecond()
    local elapsed::Float64 = time() - t0
    if elapsed > exitAfter
      diagDump("the workload has not ended after $(round(Int, elapsed)) s, second dump, ending the process")
      diagTerminate(3)
    elseif !dumped && elapsed > dumpAfter
      dumped = true
      diagDump("the workload has not ended after $(round(Int, elapsed)) s")
    end
  end
  return nothing
end

"""
    diagWatchdogStart()::DiagWatchdog

Start the watchdog when OMBACKEND_DIAG_WATCHDOG gives its two limits and the other thread pool has a thread
(`julia -t1,1`: the workload runs in the interactive pool, the watchdog in the default one).
"""
function diagWatchdogStart()::DiagWatchdog
  local done::Threads.Atomic{Bool} = Threads.Atomic{Bool}(false)
  local spec::String = get(ENV, "OMBACKEND_DIAG_WATCHDOG", "")
  (diagOn() && !isempty(spec)) || return DiagWatchdog(done, nothing)
  local limits::Vector{Float64} = Float64[parse(Float64, s) for s::SubString{String} in split(spec, ',')]
  local dumpAfter::Float64 = limits[1]
  local exitAfter::Float64 = limits[2]
  local task::Union{Nothing, Task} = nothing
  if Threads.threadpool() === :interactive && Threads.nthreads(:default) > 0
    task = Threads.@spawn :default diagWatchdogLoop(done, dumpAfter, exitAfter)
  elseif Threads.threadpool() === :default && Threads.nthreads(:interactive) > 0
    task = Threads.@spawn :interactive diagWatchdogLoop(done, dumpAfter, exitAfter)
  end
  diagLine(task === nothing ? "watchdog: no second thread, not started" :
           "watchdog: started, dump after $(dumpAfter) s, end after $(exitAfter) s")
  return DiagWatchdog(done, task)
end

"End the watchdog and wait for its task (a task left running keeps a precompile worker from finishing)."
function diagWatchdogStop(watchdog::DiagWatchdog)::Nothing
  watchdog.done[] = true
  local task::Union{Nothing, Task} = watchdog.task
  task === nothing || wait(task)
  return nothing
end
