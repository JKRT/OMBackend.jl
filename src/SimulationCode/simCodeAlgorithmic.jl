abstract type ModelicaFunction end

struct MODELICA_FUNCTION <: ModelicaFunction
  name::String
  #= Vector of integer, real string etc=#
  inputs::Vector
  outputs::Vector
  locals::Vector
  statements::Vector{DAE.Statement}
end

struct EXTERNAL_MODELICA_FUNCTION <: ModelicaFunction
  name::String
  inputs::Vector
  outputs::Vector
  libInfo::String
end
