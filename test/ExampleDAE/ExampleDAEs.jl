module ExampleDAEs

using Absyn
using MetaModelica
using ImmutableList
import SCode

using ..FrontendUtil

import .DAE
import .Prefix
import .ClassInf

include("helloWorld.jl")
include("bouncingBall.jl")
include("vanDerPol.jl")
include("influenca.jl")
include("lotkaVolterra.jl")

export helloWorld_DAE
export bouncingBall_DAE
export vanDerPol_DAE
export influenca_DAE
export lotkaVolterra_DAE

end #= ExampleDAEs =#
