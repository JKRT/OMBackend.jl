[![Github Action CI](https://github.com/JKRT/OMBackend.jl/workflows/CI/badge.svg)](https://github.com/JKRT/OMBackend.jl/actions) [![License: OSMC-PL](https://img.shields.io/badge/license-OSMC--PL-lightgrey.svg)](OSMC-License.txt)
# About OMBackend.jl
OMBackend.jl is one component of the new OpenModelica Compiler infrastructure for Julia.
It is able to transform a given Hybrid DAE and simulate it using DifferentialEquations.jl.

# Dependencies
* Julia 1.7
* ExportAll.jl
* MetaModelica.jl
* Setfield.jl
* OMCompiler.jl
* DifferentialEquations.jl

And more..

# Installation
Install dependencies with:
```julia
julia> import Pkg
julia> Pkg.build("OMBackend")
```
or:
```julia
julia> include("deps/build.jl")
```
Then precompile with:
```julia
(v1.7) pkg> activate .
julia> using OMBackend
```

# Running tests

In the folder containing OMBackend:
```julia
julia> activate .
```
Then run:

```julia
(OMBackend) pkg> test
```

As an alternative  assuming you have activated OMBackend:
```julia
julia> include("test/runtests.jl")
```

# Example use 
Assuming you use DAE.jl and a suitable frontend you can use OMBackend to simulate your Modelica models.

```
model BouncingBallReals
  parameter Real e=0.7;
  parameter Real g=9.81;
  Real h(start=1);
  Real v;
equation 
  der(h) = v;
  der(v) = -g;  
  when h <= 0 then
    reinit(v, -e*pre(v));
  end when;
end BouncingBallReals;
```

Simply pass the given DAE to the function translate. 

```
julia> OMBackend.translate(BouncingBallReals)
julia> OMBackend.simulate("BouncingBallReals", tspan = (0.0, 2.5))
```

![image](https://user-images.githubusercontent.com/8775827/99516636-b6914280-298e-11eb-85cf-c9041314e9b4.png)
