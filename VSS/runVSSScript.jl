#= Script to compile the two models above =#

using Revise
import OM
import OMBackend
import OMFrontend

"""
Translates models. After this we have the model in memory.
"""
function translateModelsMTK(models, file)
  for model in models
    @info "Translating : $model"
    @time OM.translateModel(model, "$(file).mo", mode = OMBackend.MTK_MODE)
  end
end

function dumpModelsMTK(models, file)
  local res
  #= Get the simulation code =#
  local scode = OM.translateToSCode("$(file).mo")
  for model in models
    res = OMFrontend.instantiateSCodeToDAE(model, scode)
    OMBackend.printInitialSystem(res[1])
  end
end


function flattenPendulum()
  local tst = ["Pendulum"]
  local F = "Pendulum"
#  @info "Dumping the models"
  dumpModelsMTK(tst, F)
#  #=lets try to run=#
  translateModelsMTK(tst, F)
end

function flattenFreeFall()
  local tst = ["FreeFall"]
  local F = "FreeFall"
  dumpModelsMTK(tst, F)
#  #=lets try to run=#
  translateModelsMTK(tst, F)
end

function flattenBouncingBall()
  local tst = ["BouncingBall"]
  local F = "BouncingBall"
  dumpModelsMTK(tst, F)
#  #=lets try to run=#
  translateModelsMTK(tst, F)
end

flattenFreeFall()
flattenPendulum()
flattenBouncingBall()

OMBackend.writeModelToFile("FreeFall", "FreeFall.jl"; keepComments = false, formatFile = true, mode = OMBackend.MTK_MODE)
OMBackend.writeModelToFile("Pendulum", "Pendulum.jl"; keepComments = false, formatFile = true, mode = OMBackend.MTK_MODE)
OMBackend.writeModelToFile("BouncingBall", "BouncingBall.jl"; keepComments = false, formatFile = true, mode = OMBackend.MTK_MODE)
