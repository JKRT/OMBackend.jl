module ExampleDAEs

using Absyn
using MetaModelica
using BackendDAE
using ImmutableList

import DAE
import Prefix
import SCode

include("helloWorld.jl")
include("bouncingBall.jl")

export HelloWorld_DAE
export BouncingBall_DAE

end #= ExampleDAEs =#
