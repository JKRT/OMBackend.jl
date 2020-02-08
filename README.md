


# Dependencies
* Julia 1.1, 1.2 or 1.3
* ExportAll.jl
* MetaModelica.jl
* Setfield.jl
* OMCompiler.jl
* DifferentialEquations.jl

# Installation
Install dependencies with
```julia
julia> import Pkg
julia> Pkg.build("OMBackend")
```
or
```julia
julia> include("deps/build.jl")
```
Then precompile with
```julia
julia> using OMBackend
```

# Test
Assuming that you are developing the OMBackend package
```julia
julia> dev .
julia> using OMBackend
julia> BackendDAE.translate()
```

# Package structure
* The main interface should be governed via the file all.jl in ./src/Main
* Utility functions for DAE traversal can be found in the Util package and used by executing import Util
