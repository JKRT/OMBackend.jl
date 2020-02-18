@info("OMBackend: Starting build script")

push!(LOAD_PATH, "@v#.#", "@stdlib")
@info("Current loadpath: $LOAD_PATH")

using Pkg

# Add dependencies
function buildDeps()
  Pkg.add("ExportAll")
  Pkg.add("Sundials")
  Pkg.add("CSV")
  Pkg.add("DifferentialEquations")
  Pkg.add("Setfield")
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/ImmutableList.jl"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/MetaModelica.jl"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/Absyn.jl"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/SCode.jl"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/DoubleEnded.jl"))
  @info("Build all dependencies succesfull")
end

buildDeps()

@info("OMBackend: Finished build script")
