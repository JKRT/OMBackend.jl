


# Dependencies
* Julia 1.1, 1.2 or 1.3
* ExportAll.jl
* MetaModelica.jl
* Setfield.jl
* OMCompiler.jl
* DifferentialEquations.jl

# Installation
Install dependencies:
```julia
(v1.3) pkg> add ExportAll
(v1.3) pkg> add Sundials
(v1.3) pkg> add DifferentialEquations
(v1.3) pkg> add Setfield
(v1.3) pkg> add https://github.com/OpenModelica/ImmutableList.jl.git
(v1.3) pkg> add https://github.com/OpenModelica/MetaModelica.jl.git
(v1.3) pkg> add https://github.com/OpenModelica/Absyn.jl.git
(v1.3) pkg> add https://github.com/OpenModelica/SCode.jl.git
(v1.3) pkg> add https://github.com/OpenModelica/DoubleEnded.jl.git
```

# Test
Assuming that you are developing the OMBackend package
```julia
(v1.3) pkg> test OMBackend
```

# Package structure
* The main interface should be governed via the file all.jl in ./src/Main
* Utility functions for DAE traversal can be found in the Util package and used by executing import Util
