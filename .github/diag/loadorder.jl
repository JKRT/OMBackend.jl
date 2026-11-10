#=
  Load packages in the order given, in a normal process, with a stage line and the gettext/Glib libraries of the
  process after each one (branch diag/windows-ext-hang only). Run under driver.jl:

    julia --project=<env> loadorder.jl OMBackend OMFrontend
    julia --project=<env> loadorder.jl OMFrontend OMBackend

    DIAG_SELFTEST_STALL   after the first stage line, stall: the self-test of the driver's stall path and of the
                          snapshot (on Windows behind a message box: the snapshot has to show the dialog)
=#
include(joinpath(@__DIR__, "..", "..", "ext", "OMBackendOMFrontendExt", "diag.jl"))

"Import `names` one at a time; the stage after each tells which import a hang or an error is in."
function loadInOrder(names::Vector{String})::Nothing
  diagProcess("load: start, order " * join(names, ", "))
  if haskey(ENV, "DIAG_SELFTEST_STALL")
    @static if Sys.iswindows()
      diagStage("load: self-test, a message box")
      ccall((:MessageBoxW, "user32"), stdcall, Cint, (Ptr{Cvoid}, Cwstring, Cwstring, Cuint),
            C_NULL, "diag self-test text", "diag self-test title", 0)
    end
    diagStage("load: self-test, sleeping")
    sleep(600)
  end
  diagLibraries("load: before the imports")
  for name::String in names
    try
      Core.eval(Main, Expr(:import, Expr(:., Symbol(name))))
      diagLibraries("load: imported " * name)
    catch err
      diagStage("load: import " * name * " threw: " * first(sprint(showerror, err), 2000))
      diagLibraries("load: after the failed import of " * name)
    end
  end
  local loaded::Vector{Symbol} = Symbol[nameof(m) for m::Module in Base.loaded_modules_array()]
  diagStage("load: extension OMBackendOMFrontendExt " * (:OMBackendOMFrontendExt in loaded ? "loaded" : "not loaded"))
  diagStage("load: end")
  return nothing
end

loadInOrder(String[String(a) for a in ARGS])
