#=
  The backend of OMBackend..
=#

module Backend

using ..FrontendUtil

include("BDAE.jl")
include("BackendEquation.jl")
include("BDAEUtil.jl")
include("BDAECreate.jl")

include("Causalize.jl")

#include("backendDump.jl")

export BDAE, BDAECreate, BackendEquation, BDAEUtil, Causalize, BackendDAEExp

end
