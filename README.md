# About OMBackend.jl
OMBackend.jl is one component of the new OpenModelica Compiler infrastructure for Julia.
It is able to transform given DAE IR and simulate using DifferentialEquations.jl as it's backend.

# Dependencies
* Julia 1.1, 1.2 or 1.3
* ExportAll.jl
* MetaModelica.jl
* Setfield.jl
* OMCompiler.jl
* DifferentialEquations.jl

Experimental support for Julia 1.5

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
