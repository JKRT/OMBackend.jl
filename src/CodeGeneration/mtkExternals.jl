#= The ModelingToolkit.jl package is licensed under the MIT "Expat" License:
# Copyright (c) 2018-20: Christopher Rackauckas, Julia Computing.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  SOFTWARE.
=#  

#= So we know about t an der in the global scope. This is needed for the rules below to match correctly. =#
@variables t
const D = Differential(t)

"""
  Temporary rewrite function. Not very pretty...
  Original code by Chris R. Modifed it to fix terms of type X * D(Y).
  The solution to solve it is not pretty and is probably quite flaky.
"""
function move_diffs(eq::Equation; rewrite)
  # Do not modify `D(x) ~ ...`, already correct
  # Ignore `δ(x) ~ ...` for now
  res = if !(eq.lhs isa Term && operation(eq.lhs) isa Differential) &&
    !(eq.lhs isa Term && operation(eq.lhs) isa Difference)
    _eq = eq.rhs-eq.lhs
    rhs = rewrite(_eq)
    if rhs === nothing
      eq
    end
    lhs = _eq - rhs
    if !(lhs isa Number) && (operation(lhs) isa Differential)
      lhs ~ -rhs
    elseif !(lhs isa Number) && (operation(lhs) == *)
      #= This code is probably quite flaky, however, it should not be needed after similar things are introduced in MTK. =#
      newRhs = substitute(rhs, rhs => rhs / arguments(lhs)[1])
      newLhs = substitute(lhs, lhs => arguments(lhs)[2])
     #lhs = lhs_
      ~(newLhs, -newRhs)
    else
      -lhs ~ rhs
    end
  else
    eq
  end
  return res
end

function makeODESystem(deqs, iv, vars, pars; name)
  # Equation is not in the canonical form
  # Use a rule to find all g(u)*D(u) terms
  # and move those to the lhs
  D = Differential(iv)
  r1 = SymbolicUtils.@rule ~~a * D(~~b) * ~~c => 0
  r2 = SymbolicUtils.@rule D(~~b) => 0
  debugRewrite(deqs, iv, vars, pars)
  remove_diffs = SymbolicUtils.Postwalk(SymbolicUtils.Chain([r1,r2]))
  deqs = [move_diffs(eq, rewrite = remove_diffs) for eq in deqs]
  debugRewrite(deqs, iv, vars, pars)
  res = ModelingToolkit.ODESystem(deqs, iv, vars, pars; name = name)
  return res
end

function debugRewrite(deqs, t, vars, parameters)
  print("@variables")
  print("(")
  for v in vars
    print(v, ",")
  end
  println(")")
  print("eqs =")
  println("[")
  for eq in deqs
    println(eq, ",")
  end
  println("]")
  print("@parameters")
  print("(")
  for p in parameters
    print(p, ",")
  end
  println(")")
end
