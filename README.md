


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
(v1.1) pkg> activate .
julia> using OMBackend
```

# Run tests
Assuming you have activated OMBackend
```julia
julia> include("test/runtests.jl")
```
