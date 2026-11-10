#=
  Run one step of the Windows hang diagnostic (branch diag/windows-ext-hang only) and stop it when it stalls.

    julia driver.jl <label> <stallSeconds> <capSeconds> -- <command...>

  The command runs with OMBACKEND_DIAG=1 and OMBACKEND_DIAG_LOG=$DIAG_DIR/<label>.stages.log (the stage lines of
  ext/OMBackendOMFrontendExt/diag.jl). A stall: the stage file exists, its last line is not a final stage, and no
  line was added for <stallSeconds>. On a stall, or after <capSeconds> in all: a snapshot of the process that wrote
  the last stage line (snapshot.ps1 on Windows), that process is killed so its parent reports what it printed, and
  the command gets three minutes to end before its process tree is killed.

    DIAG_DIR          the directory of the stage files and snapshots (default: the working directory)
    DIAG_FORGET_EXT   delete the extension's precompile cache first, so the command compiles it again
    GITHUB_OUTPUT     receives where=<none|import|workload|end> (the kind of the last stage) and stalled=<true|false>

  Exit code: the command's, or 124 when it was stopped.
=#

"The last line of the stage file that is not empty (\"\": no file, or no line yet)."
function lastStage(path::String)::String
  isfile(path) || return ""
  local lines::Vector{String} = filter!(!isempty, readlines(path))
  return isempty(lines) ? "" : last(lines)
end

"True for the stage line that ends a step: nothing can stall after it but the command itself."
isFinalStage(line::String)::Bool = occursin(r"(ext: module end|load: end|script: end)$", line)

"The kind of the last stage: where a step stopped (none: the extension's module was never entered)."
function stageKind(line::String)::String
  isempty(line) && return "none"
  isFinalStage(line) && return "end"
  occursin("] workload:", line) && return "workload"
  return "import"
end

"The pid in a stage line (`[diag HH:MM:SS pid N] ...`), nothing when there is none."
function stagePid(line::String)::Union{Nothing, Int}
  local m::Union{Nothing, RegexMatch} = match(r"pid (\d+)\]", line)
  return m === nothing ? nothing : parse(Int, m.captures[1])
end

"Run a command for its side effect; a failure is reported, not thrown (the diagnostic goes on)."
function attempt(cmd::Cmd)::Nothing
  try
    run(pipeline(ignorestatus(cmd); stdout = stdout, stderr = stderr))
  catch err
    println("[driver] ", cmd, " failed: ", err)
  end
  return nothing
end

"The state of the stalled process, from outside: files `<prefix>.*` (processes, threads, windows, screenshot)."
function snapshot(prefix::String, pid::Int)::Nothing
  println("[driver] snapshot of pid ", pid, " -> ", prefix, ".*")
  if Sys.iswindows()
    attempt(`powershell -NoProfile -NonInteractive -ExecutionPolicy Bypass -File $(joinpath(@__DIR__, "snapshot.ps1")) $prefix $pid`)
  else
    open(prefix * ".processes.txt", "w") do io::IO
      try
        run(pipeline(ignorestatus(`ps -o pid,ppid,stat,time,etime,command -p $pid`); stdout = io, stderr = io))
      catch err
        println(io, "ps failed: ", err)
      end
    end
  end
  return nothing
end

"Kill one process (not its children)."
function killProcess(pid::Int)::Nothing
  println("[driver] killing pid ", pid)
  attempt(Sys.iswindows() ? `taskkill /F /PID $pid` : `kill -9 $pid`)
  return nothing
end

"Kill the command and, on Windows, the processes it started."
function killTree(proc::Base.Process)::Nothing
  println("[driver] killing the command's process tree")
  if Sys.iswindows()
    attempt(`taskkill /F /T /PID $(getpid(proc))`)
  else
    kill(proc, Base.SIGKILL)
  end
  return nothing
end

"Wait until the process has ended, at most `seconds`; true when it has."
function waitFor(proc::Base.Process, seconds::Float64)::Bool
  local t0::Float64 = time()
  while process_running(proc) && time() - t0 < seconds
    sleep(2)
  end
  return !process_running(proc)
end

"Append `key=value` to the file GITHUB_OUTPUT names (the outputs of a workflow step)."
function stepOutput(key::String, value::String)::Nothing
  println("[driver] ", key, "=", value)
  local path::String = get(ENV, "GITHUB_OUTPUT", "")
  isempty(path) && return nothing
  open(path, "a") do io::IO
    println(io, key, "=", value)
  end
  return nothing
end

"The extension's precompile cache in the first depot is deleted: the next load or Pkg.precompile compiles it."
function forgetExtension()::Nothing
  local dir::String = joinpath(DEPOT_PATH[1], "compiled", "v$(VERSION.major).$(VERSION.minor)", "OMBackendOMFrontendExt")
  println("[driver] forgetting ", dir, isdir(dir) ? "" : " (not there)")
  rm(dir; force = true, recursive = true)
  return nothing
end

"""
    drive(label::String, stall::Float64, cap::Float64, command::Vector{String})::Int

Run `command` as the step `label` (see the head of the file); the exit code of the step.
"""
function drive(label::String, stall::Float64, cap::Float64, command::Vector{String})::Int
  local dir::String = get(ENV, "DIAG_DIR", pwd())
  mkpath(dir)
  local stages::String = joinpath(dir, label * ".stages.log")
  rm(stages; force = true)
  haskey(ENV, "DIAG_FORGET_EXT") && forgetExtension()
  local env::Dict{String, String} = Dict{String, String}(ENV)
  env["OMBACKEND_DIAG"] = "1"
  env["OMBACKEND_DIAG_LOG"] = stages
  println("[driver] ", label, ": ", join(command, " "), " (stall ", stall, " s, cap ", cap, " s)")
  local proc::Base.Process = run(pipeline(setenv(Cmd(command), env); stdout = stdout, stderr = stderr); wait = false)
  local t0::Float64 = time()
  local changed::Float64 = t0
  local size::Int = -1
  local stopped::Bool = false
  while process_running(proc)
    sleep(5)
    local now::Float64 = time()
    local line::String = lastStage(stages)
    if isfile(stages) && filesize(stages) != size
      size = filesize(stages)
      changed = now
      println("[driver] ", round(Int, now - t0), " s, stage: ", line)
    end
    local stalled::Bool = !isempty(line) && !isFinalStage(line) && now - changed > stall
    if (stalled || now - t0 > cap) && process_running(proc)
      stopped = true
      println("[driver] ", stalled ? "STALL: no stage line for $(round(Int, now - changed)) s" :
              "CAP: $(round(Int, now - t0)) s", "; last stage: ", isempty(line) ? "none" : line)
      local pid::Union{Nothing, Int} = stagePid(line)
      snapshot(joinpath(dir, label), pid === nothing ? getpid(proc) : pid)
      pid === nothing || killProcess(pid)
      waitFor(proc, 180.0) || killTree(proc)
      waitFor(proc, 60.0)
      break
    end
  end
  local final::String = lastStage(stages)
  stepOutput("where", stageKind(final))
  stepOutput("stalled", string(stopped))
  println("[driver] ", label, ": ", stopped ? "stopped" : "exit code $(proc.exitcode)", " after ", round(Int, time() - t0),
          " s; last stage: ", isempty(final) ? "none" : final)
  return stopped ? 124 : Int(proc.exitcode)
end

"The arguments of the script: label, the two limits, `--`, the command."
function main(args::Vector{String})::Int
  local sep::Union{Nothing, Int} = findfirst(==("--"), args)
  (sep == 4 && length(args) > 4) || error("usage: driver.jl <label> <stallSeconds> <capSeconds> -- <command...>")
  return drive(args[1], parse(Float64, args[2]), parse(Float64, args[3]), args[5:end])
end

exit(main(ARGS))
