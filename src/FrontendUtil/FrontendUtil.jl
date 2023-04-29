module FrontendUtil

import Absyn
import OMFrontend
using MetaModelica

include("AbsynUtil.jl")
include("Util.jl")
include("Prefix.jl")

export Util
export AbsynUtil
export Prefix

"""
 This function handles certain builtin functions.
 For now it only removes the smooth function
"""
function handleBuiltin(fm::OMFrontend.Main.FLAT_MODEL)
  fm = removeSmoothOperator(fm)
  return fm
end

"""
  Removes the smooth operator from the set of equations.
  (From the specification the application of this function seems to be optional)
See:
https://build.openmodelica.org/Documentation/ModelicaReference.Operators.%27smooth()%27.html
"""
function removeSmoothOperator(fm::OMFrontend.Main.FLAT_MODEL)
  local equations = fm.equations
  equations = OMFrontend.Main.mapList(equations, removeSmooth)
  @assign fm.equations = equations
  return fm
end

"""
  Removes the smooth operator from an equation, returns the argument to smooth.
"""
function removeSmooth(eq::OMFrontend.Main.Equation)
  @match eq begin
    OMFrontend.Main.EQUATION_EQUALITY(lhs = e1, rhs = OMFrontend.Main.CALL_EXPRESSION(call)) where string(OMFrontend.Main.functionName(call)) == "smooth"  => begin
      local arguments = OMFrontend.Main.arguments(call)
      @match x <| y <| nil = arguments
      local newEq = eq
      println("Changed eq:" * OMFrontend.Main.toString(newEq))
      @assign newEq.rhs = y
      println("Changed eq:" * OMFrontend.Main.toString(newEq))
      newEq
    end
    _ => eq
  end
end

end
