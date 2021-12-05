@info("OMBackend: Starting build script")

push!(LOAD_PATH, "@v#.#", "@stdlib")
@info("Current loadpath: $LOAD_PATH")

using Pkg

function buildDeps()
  Pkg.add("ExportAll")
  Pkg.add("Sundials")
  Pkg.add("CSV")
  Pkg.add("DifferentialEquations")
  Pkg.add("Setfield")
  Pkg.add("DataStructures")
  Pkg.add("Plots")
  Pkg.add("GraphPlot")
  Pkg.add("Cairo")
  Pkg.add("Compose")
  Pkg.add("LightGraphs")
  Pkg.add("MetaGraphs")
  Pkg.add("JuliaFormatter")
  Pkg.add("Reduce")
  Pkg.add("MacroTools")
  #= This packages are available using the OpenModelica Julia registry=#
  Pkg.add("ImmutableList")
  Pkg.add("MetaModelica")
  Pkg.add("Absyn")
  Pkg.add("SCode")
  Pkg.add("DAE")
  Pkg.add("DoubleEnded")
  #Pkg.add("OpenModelicaParser") #= Add if we do not have it =#
  Pkg.add("OMFrontend")
  @info("All dependencies succesfull")
  Pkg.build("OMFrontend") #= This will call build on the parser implicitly. =#
  @info("OMFrontend built successfully")
end

buildDeps()
@info("OMBackend: Finished build script")
