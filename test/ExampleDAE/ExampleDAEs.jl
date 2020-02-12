module ExampleDAEs

using Absyn
using MetaModelica
using BackendDAE
using ImmutableList

import DAE
import Prefix
import SCode
import ClassInf

include("helloWorld.jl")
include("bouncingBall.jl")
include("vanDerPol.jl")
include("influenca.jl")

export helloWorld_DAE
export bouncingBall_DAE
export vanDerPol_DAE
export influenca_DAE

end #= ExampleDAEs =#
