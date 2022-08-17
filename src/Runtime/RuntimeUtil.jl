module RuntimeUtil
import Absyn
import SCode
import OMFrontend
import OMFrontend.Main.SCodeUtil

using MetaModelica

"""
  Wrapper to a function in SCode util.
  inIdent is a string since frontend automatically remove certain parameters.
  TODO: Also search for components in the innermost class.
"""
function getElementFromSCodeProgram(inIdent::String, inClass::SCode.Element)
  result = SCodeUtil.getElementNamed(inIdent, inClass)
  return result
end

"""
  Given a name sets that element to a new value.
  It then returns the modified SCodeProgram.
Currently it is assumed to be at the top level of the class.
TODO: Fix for sublevels as well
"""
function setElementInSCodeProgram!(inIdent::String, newValue::T, inClass::SCode.Element) where {T}
  #= Get all elements from the class together with the corresponding names =#
  local elementToReplace = getElementFromSCodeProgram(inIdent, inClass)::SCode.Element
  local elements = listArray(SCodeUtil.getClassElements(inClass))::Vector{SCode.Element}
  local i = 1
  local indexOfElementToReplace = 0
  local outClass = inClass
  for element in elements
    if SCodeUtil.elementNameEqual(element, elementToReplace)
      indexOfElementToReplace = i
      break
    end
    i += 1
  end
  local modification = SCodeUtil.getComponentMod(elementToReplace)
  @assign modification.binding = makeCondition(newValue)
  @assign elementToReplace.modifications = modification
  elements[indexOfElementToReplace] = elementToReplace
  @assign outClass.classDef.elementLst = arrayList(elements)
  return outClass
end

makeCondition(val::Bool) = begin
  SOME(Absyn.BOOL(val))
end
makeCondition(val::Int) = begin
  SOME(Absyn.INTEGER(val))
end
makeCondition(val::Real) = begin
  SOME(Absyn.REAL(string(val)))
end
makeCondition(val) = begin
  throw("Only primitive values {Integer, Boolean, Real} are currently supported in a Recompilation call")
end

"""
  Converts a symbol (of a MTK variable) to a string.
"""
function convertSymbolToString(symbol::Symbol)
  res = replace(String(symbol), "(t)" => "")
  #= Remove prefixes in front of variables =#
  return res
end

"""
  Converts a list of symbols to a list of strings
"""
function convertSymbolsToStrings(symbols::Vector{Symbol})
  map(convertSymbolToString, symbols)
end

"""
  This function maps variables between two models during a structural change with recompilation.
  It returns a new vector of u₀ variables to initialize the new model.
  We do so by assigning the old values when the structural change occured for all variables
  that occured in the model before the structural change.
TODO:
Remove the special case.
The reason for it is in some example the containing model changes name wheras in others it does not.
Furthermore, the way indices are handled need to be fixed since the index after reduction in MTK is not the same as
the index in the simulation code stage of the backend.
For now we return the array for the special case with dynamic overconstrained connectors.
"""
function createNewU0(symsOfOldProblem::Vector{Symbol},
                     symsOfNewProblem::Vector{Symbol},
                     newHT,
                     initialValues,
                     integrator,
                     specialCase)
#  @debug "Length of old and new" length(symsOfOldProblem) length(symsOfNewProblem)
  local newU0 = Float64[last(initialValues[idx]) for idx in 1:length(symsOfNewProblem)]
  local variableNamesOldProblem = RuntimeUtil.convertSymbolsToStrings(symsOfOldProblem)
  local variableNamesNewProblem = RuntimeUtil.convertSymbolsToStrings(symsOfNewProblem)
  @debug variableNamesOldProblem
  @debug variableNamesNewProblem
  #= It was assumed to only be real variable not discretes, which might have other indices? =#
  local variableNamesWithoutPrefixesOP
  local variableNamesWithoutPrefixesNP
  if ! specialCase
    variableNamesWithoutPrefixesOP = [replace(k, r".*_" => "") for k in variableNamesOldProblem]
    variableNamesWithoutPrefixesNP = [replace(k, r".*_" => "") for k in variableNamesNewProblem]
  else
    return integrator.u
  end
  local largestProblem = if length(variableNamesOldProblem) > length(variableNamesOldProblem)
    variableNamesWithoutPrefixesOP
  else
    variableNamesWithoutPrefixesNP
  end
  for v in largestProblem
    local varNameWithoutPrefix = v
    if varNameWithoutPrefix in variableNamesWithoutPrefixesOP && varNameWithoutPrefix in variableNamesWithoutPrefixesNP
      local oldIndices = findall((x)-> x == varNameWithoutPrefix, variableNamesWithoutPrefixesOP)
      @assert(length(oldIndices) == 1, "Zero or more than one variable with that name. Size of oldIndicies was $(length(oldIndices)). Name was $(v)")
      idxOldVar = first(oldIndices)
      #= Locate the index of a variable with that name in the set of new variables=#
      local indices = findall((x)-> x == varNameWithoutPrefix, variableNamesWithoutPrefixesNP)
      #= I assume here that there are no duplicate variables =#
      @assert(length(indices) == 1, "Zero or more than one variable with that name. Size of indices was $(length(indices)). Name was $(v)")
      local varNameNew = variableNamesNewProblem[first(indices)]
      local idxNewVar = getIdxFromEntry(newHT[varNameNew])
      newU0[idxNewVar] = integrator.u[idxOldVar]
    end
  end
  return newU0
end

"""
  Gets the index from an entry to the symbol table
"""
function getIdxFromEntry(entry::Tuple)::Int
  first(entry)
end


"""
  Creates a new flat model.
  This model either has the equation specified in the if or these equations are to be added to the model
  when this condition is true.
  The name of this model
"""
function createNewFlatModel(flatModel,
                            unresolvedEquations,
                            newEquations)
  local newFlatModel =
    OMFrontend.Main.FLAT_MODEL(flatModel.name,
                               flatModel.variables,
                               flatModel.equations,
                               flatModel.initialEquations,
                               flatModel.algorithms,
                               flatModel.initialAlgorithms,
                               nil,
                               NONE(),
                               nil,
                               nil,
                               Bool[], #= TODO: This equation might need to be changed. =#
                               flatModel.comment)
  # println("Equations to Add:")
  # println("********************************************************************")
  # @debug length(newEquations)
  # for e in newEquations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")

  # println("Existing equations:")
  # println("********************************************************************")
  # @debug "Length of existing system:" length(newFlatModel.equations)
  # for e in newFlatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  @assign newFlatModel.equations = listAppend(unresolvedEquations, newEquations)
  # println("Unresolved System:")
  # println("********************************************************************")
  # for e in newFlatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  #=
    1. Reresolve the connect equations.
    2. Perform constant evaluation.
    3. Run the simplification pass.
  =#
  newFlatModel = OMFrontend.Main.resolveConnections(newFlatModel, newFlatModel.name)
  newFlatModel = OMFrontend.Main.evaluate(newFlatModel)
  newFlatModel = OMFrontend.Main.simplifyFlatModel(newFlatModel)
  # println("Final System")
  # @debug "Length of final system:" length(newFlatModel.equations)
  # println("********************************************************************")
  # for e in newFlatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  return newFlatModel
end

"""
  Creates a new flat model with a set of connection equation removed.
  Note that the flat model passed to this function does not have any active equations.
  This means that this flat model is a new flat model with the unbreakable branches removed.
"""
function createNewFlatModel(flatModel,
                            idx,
                            unresolvedEquations)
  local aDoccs = flatModel.active_DOCC_Equations
  aDoccs[idx] = false
  local newFlatModel =
    OMFrontend.Main.FLAT_MODEL(flatModel.name,
                               flatModel.variables,
                               flatModel.unresolvedConnectEquations,
                               flatModel.initialEquations,
                               flatModel.algorithms,
                               flatModel.initialAlgorithms,
                               nil,
                               NONE(),
                               flatModel.DOCC_equations,
                               flatModel.unresolvedConnectEquations,
                               aDoccs,
                               flatModel.comment)
  # println("********************************************************************")  
  # println("Existing equations:")
  # println("********************************************************************")
  # @debug "Length of existing system:" length(flatModel.equations)
  # for e in flatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  
  # println("********************************************************************")
  # println("New System:")
  # println("********************************************************************")
  # @debug "Length of new system:" length(newFlatModel.unresolvedConnectEquations)
  # for e in newFlatModel.unresolvedConnectEquations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  #=
    1. Reresolve the connect equations.
    2. Perform constant evaluation.
    3. Run the simplification pass.
  =#
  newFlatModel = OMFrontend.Main.resolveConnections(newFlatModel, newFlatModel.name)
  newFlatModel = OMFrontend.Main.evaluate(newFlatModel)
  newFlatModel = OMFrontend.Main.simplifyFlatModel(newFlatModel)
  # println("Final System")
  # @debug "Length of final system:" length(newFlatModel.equations)
  # println("********************************************************************")
  # for e in newFlatModel.equations
  #   println(OMFrontend.Main.toString(e))
  # end
  # println("********************************************************************")
  return newFlatModel
end

end #= RuntimeUtil =#
