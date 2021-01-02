![Build Status](https://travis-ci.com/JKRT/OMBackend.jl.svg?token=spA4q6rnAXNpiqHMwpnL&branch=master)
# About OMBackend.jl
OMBackend.jl is one component of the new OpenModelica Compiler infrastructure for Julia.
It is able to transform a given Hybrid DAE and simulate it using DifferentialEquations.jl.

Currently only so called DAE-mode is supported, using the Sundials package with the IDA solver. 

# Dependencies
* Julia 1.5
* ExportAll.jl
* MetaModelica.jl
* Setfield.jl
* OMCompiler.jl
* DifferentialEquations.jl

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
(v1.6) pkg> activate .
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
