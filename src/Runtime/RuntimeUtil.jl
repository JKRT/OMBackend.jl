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
"""
function createNewU0(symsOfOldProblem::Vector{Symbol},
                     symsOfNewProblem::Vector{Symbol},
                     oldHT,
                     newHT,
                     initialValues,
                     integrator)
  local newU0 = Float64[last(initialValues[idx]) for idx in 1:length(symsOfNewProblem)]
  local variableNamesOldProblem  = RuntimeUtil.convertSymbolsToStrings(symsOfOldProblem)
  local variableNamesNewProblem  = RuntimeUtil.convertSymbolsToStrings(symsOfNewProblem)
  local variableNamesWithoutPrefixesOP = [replace(k, r".*_" => "") for k in variableNamesOldProblem]
  local variableNamesWithoutPrefixesNP = [replace(k, r".*_" => "") for k in variableNamesNewProblem]
  local largestProblem = if length(variableNamesOldProblem) > length(variableNamesOldProblem)
    variableNamesWithoutPrefixesOP
  else
    variableNamesWithoutPrefixesNP
  end
  for v in largestProblem
    local varNameWithoutPrefix = v
    if varNameWithoutPrefix in variableNamesWithoutPrefixesOP && varNameWithoutPrefix in variableNamesWithoutPrefixesNP
      local oldIndices = findall((x)-> x == varNameWithoutPrefix, variableNamesWithoutPrefixesOP)
      @assert(length(oldIndices) == 1, "Zero or more than one variable with that name!")
      idxOldVar = first(oldIndices)
      #= Locate the index of a variable with that name in the set of new variables=#
      local indices = findall((x)-> x == varNameWithoutPrefix, variableNamesWithoutPrefixesNP)
      #= I assume here that there are no duplicate variables =#
      @assert(length(indices) == 1, "Zero or more than one variable with that name!")
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

end #= RuntimeUtil =#
