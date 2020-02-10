@info("OMBackend: Starting build script")

push!(LOAD_PATH, "@v#.#", "@stdlib")
@info("Current loadpath: $LOAD_PATH")

using Pkg

# Add dependencies
function buildDeps()
  Pkg.add("ExportAll")
  Pkg.add("Sundials")
  Pkg.add("DifferentialEquations")
  Pkg.add("Setfield")
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/MetaModelica.jl.git"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/Absyn.jl.git"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/ImmutableList.jl.git"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/SCode.jl.git"))
  Pkg.add(PackageSpec(url="https://github.com/OpenModelica/DoubleEnded.jl.git"))
  @info("Build all dependencies succesfull")
end

buildDeps()

@info("OMBackend: Finished build script")
