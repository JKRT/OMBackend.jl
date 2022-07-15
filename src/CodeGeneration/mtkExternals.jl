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

#=
  So we know about t an der in the global scope.
  This is needed for the rules below to match correctly.
=#
@variables t
const D = Differential(t)

"""
  Temporary rewrite function. Not very pretty...
  Original code by Chris R. Expanded to fix terms of type X * D(Y).
  The solution to solve it is not pretty and is probably flaky.
"""
function move_diffs(eq::Equation; rewrite)
  # Do not modify `D(x) ~ ...`, already correct
  res =
    if !(eq.lhs isa Term && operation(eq.lhs) isa Differential) && !(eq.lhs isa Term && operation(eq.lhs) isa Difference)
      _eq = eq.rhs-eq.lhs
      rhs = rewrite(_eq)
      if rhs === nothing
        eq
      end
      lhs = _eq - rhs
      if !(lhs isa Number) && (operation(lhs) isa Differential)
        lhs ~ -rhs
      elseif !(lhs isa Number) && (operation(lhs) == *)
        #=
          This code is probably quite flaky, however, it should not be needed after similar things are introduced in MTK.
          TODO: Handle this more elegantly using rewrite rules instead.
        =#
        local newRhs
        local newLhs
        for arg in arguments(lhs)
#          @info "Arg was:" arg typeof(arg)
          if (arg isa Number || arg isa Sym) || (arg isa Term && !(operation(arg) isa Differential))
#            @info "Rewrite case:" lhs
            newRhs = substitute(rhs, rhs => rhs / arg)
            rhs = newRhs
            newLhs = substitute(lhs, lhs => lhs / arg)
            lhs = newLhs
          end
        end
        tmp = ~(newLhs, -newRhs)
#        @info "After Rewrite" tmp
        tmp
      else
        -lhs ~ rhs
      end
    else
      eq
    end
  return res
end

"""
  If index reduction is needed we currently do a hack.
  That is we do an eval of the system, this is computationally expensive, however, it seems to work for now.
  Related to issue https://github.com/SciML/ModelingToolkit.jl/issues/1493
"""
function makeODESystem(deqs, iv, vars, pars, idxReduction; name, continuous_events = [])
  # Equation is not in the canonical form
  # Use a rule to find all g(u)*D(u) terms
  # and move those to the lhs
  D = Differential(iv)
  r1 = SymbolicUtils.@rule ~~a * D(~~b) * ~~c => 0
  r2 = SymbolicUtils.@rule D(~~b) => 0
  #@info "BEFORE REWRITING"
  #println(debugRewrite(deqs, iv, vars, pars; separator = "\n"))
  remove_diffs = SymbolicUtils.Postwalk(SymbolicUtils.Chain([r1,r2]))
  usedStates = Set()
  local rewrittenDeqs = Symbolics.Equation[]
  local req
  for eq in deqs
#    if !(operation(eq.lhs) isa Differential)
      req = move_diffs(eq, rewrite = remove_diffs)
#    end
    #    @info "Left hand side of the equation" req.lhs
    if req.lhs isa Real
      push!(rewrittenDeqs, req)
    elseif !(req.lhs in usedStates)
      #      @info "Not a duplicate" req.lhs
      #      @info "Used equations are" usedStates
      push!(rewrittenDeqs, req)
      push!(usedStates, req.lhs)
    else
      #      @info "Duplicate:" req.lhs
      push!(rewrittenDeqs, eq)
    end
  end
  #= Special computationally heavy routine for systems that need index reduction.. =#
  #@info "AFTER REWRITING"
  #println(debugRewrite(rewrittenDeqs, iv, vars, pars; separator = ",\n"))

  if idxReduction
    @info "We need index reduction"
    fail()
    #@info "We are doing index reduction"
    #= Convert the system to a string=#
    res = debugRewrite(rewrittenDeqs, iv, vars, pars)
    #    println(res)
    #= Parse and evaluate it =#
    expr = Meta.parse(res)
    eval(expr)
    #= Generate a new system with metadata...=#
    res = ModelingToolkit.ODESystem(eqs2, iv, vars2, pars; name = name)
    res2 = ModelingToolkit.dae_index_lowering(res)
    #structural_simplify(res2; simplify = true, allow_symbolic = true, allow_parameter = true)
    return res2
  end
  #println(debugRewrite(deqs, iv, vars, pars))
  res = if ! isempty(continuous_events)
    ModelingToolkit.ODESystem(rewrittenDeqs, iv, vars, pars; name = name, continuous_events = continuous_events)
  else
    ModelingToolkit.ODESystem(rewrittenDeqs, iv, vars, pars; name = name)
  end
  return res
end


"""
  Rewrite equations that do not conform to the requirements of MTK
"""
function rewriteEquations(edeqs, iv, eVars, ePars)
  local der = ModelingToolkit.Differential(t)
  #= Remove the t's =#
  eVars = [Symbol(replace(string(i), "(t)" => "")) for i in eVars]
  preEval = quote 
    vars = ModelingToolkit.@variables begin
      $(eVars...)
    end
    pars = ModelingToolkit.@parameters begin
      $(ePars...)
    end
  end
  #= Make the derivative symbol known =#
  eval(preEval)
  eval(:(der = ModelingToolkit.Differential(t)))
  eval(:(import ModelingToolkit.IfElse))
  local deqs = [eval(i) for i in edeqs]
  #= Rewrite equations =#
  D = Differential(iv)
  local r1 = SymbolicUtils.@rule ~~a * D(~~b) * ~~c => 0
  local r2 = SymbolicUtils.@rule D(~~b) => 0
  local remove_diffs = SymbolicUtils.Postwalk(SymbolicUtils.Chain([r1,r2]))
  local usedStates = Set()
  local rewrittenDeqs = Symbolics.Equation[]
  local req
  for eq in deqs
#    if !(operation(eq.lhs) isa Differential)
      req = move_diffs(eq, rewrite = remove_diffs)
#    end
    #    @info "Left hand side of the equation" req.lhs
    if req.lhs isa Real
      push!(rewrittenDeqs, req)
    elseif !(req.lhs in usedStates)
      #      @info "Not a duplicate" req.lhs
      #      @info "Used equations are" usedStates
      push!(rewrittenDeqs, req)
      push!(usedStates, req.lhs)
    else
      #      @info "Duplicate:" req.lhs
      push!(rewrittenDeqs, eq)
    end
  end
  return rewrittenDeqs
end

function debugRewrite(deqs, t, vars, parameters; separator = ",")
  local buffer = IOBuffer()
  print(buffer, "@variables t;")
  print(buffer, "vars2 = @variables")
  print(buffer, "(")
  for v in vars
    print(buffer, v, separator)
  end
  print(buffer, ");")
  print(buffer, "eqs2 =")
  print(buffer, "[")
  for eq in deqs
    print(buffer, replace(string(eq), "(t)" => "", "Differential" => "D"), separator)
  end
  print(buffer, "];")
  return String(take!(buffer))
end

rewriteEq(eq) = begin
  local eqStr = string(eq)
  res = Meta.parse(replace(eqStr, "Differential(t)" => "D"))
  res
end
