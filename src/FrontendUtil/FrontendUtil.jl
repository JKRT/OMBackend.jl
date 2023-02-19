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
  Removes the smooth operator from a equation, returns the argument to smooth
"""
function removeSmooth(eq::OMFrontend.Main.Equation)
  println("Calling remove smooth with " * OMFrontend.Main.toString(eq) * " and type " * string(typeof(eq)))
  if eq isa OMFrontend.Main.EQUATION_EQUALITY
    println(typeof(eq.rhs))
  end
  @match eq begin
    OMFrontend.Main.EQUATION_EQUALITY(lhs = e1, rhs = OMFrontend.Main.CALL_EXPRESSION(OMFrontend.Main.TYPED_CALL(__))) => begin
      @match OMFrontend.Main.CALL_EXPRESSION(OMFrontend.Main.TYPED_CALL(fn, ty, var, arguments, attributes))  = eq.rhs
      @match Absyn.IDENT(n) = OMFrontend.Main.name(fn)
      @match x <| y <| nil = arguments
      println(n)
      if n == "smooth"
        local newEq = eq
        @assign newEq.rhs = y
#        println("Changed eq:" * OMFrontend.Main.toString(newEq))
        newEq
      else
        eq
      end
    end
    _ => eq
  end
end

end
