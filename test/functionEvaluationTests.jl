#= The partial evaluation of Modelica function calls (SimCode functionEvaluation): an environment of abstract values
   by the scalars' references. A for loop's iterator is one of them (it was a String-keyed environment's Symbol key:
   a MethodError the interpretation took for "not followed"); an assignment to an array forgets its elements, one
   to a record its fields. =#
using Test
import OMBackend

@testset "Function evaluation by references" begin
  local SC = OMBackend.SimCode
  local IR = OMBackend.OMJLIR
  local src = SC.NO_SOURCE
  local var = (n::Symbol, ty::SC.SType = SC.TYPE_REAL()) -> SC.EXP_CREF(SC.SimCref(n, Int[], IR.ComponentRef(n)), ty)
  local element = (n::Symbol, k::Int) -> SC.EXP_CREF(SC.SimCref(n, Int[k], IR.elementCref(IR.ComponentRef(n), Int[k])), SC.TYPE_REAL())
  local fv = (n::Symbol, ty::SC.SType = SC.TYPE_REAL()) -> SC.FunctionVariable(String(n), ty, SC.Exp[], nothing, Pair{String, SC.Exp}[], IR.ComponentRef(n))
  local add = (a::SC.Exp, b::SC.Exp) -> SC.BINARY(a, SC.OP_ADD, b)
  local values = (f::SC.MODELICA_FUNCTION, args::Vector{Any}) ->
    SC.functionValues(f, args, SC.Interpretation(Dict{String, SC.MODELICA_FUNCTION}(f.name => f), Dict{Any, Any}(), 0, 0))

  #= y := 0; for i in 1:3 loop y := y + i; end for; a := {x, 2}; a[1] := 5; z := a[1] + a[2] (x unknown) =#
  local loop = SC.FOR(SC.Iterator(:i, SC.RANGE(SC.TYPE_ARRAY(SC.TYPE_INTEGER(), [3]), SC.ICONST(1), nothing, SC.ICONST(3))),
                      SC.Statement[SC.ASSIGN(var(:y), add(var(:y), var(:i, SC.TYPE_INTEGER())), src)], src)
  local arrayTy = SC.TYPE_ARRAY(SC.TYPE_REAL(), [2])
  local f = SC.MODELICA_FUNCTION("f", SC.FunctionVariable[fv(:x)], SC.FunctionVariable[fv(:y), fv(:z)],
                                 SC.FunctionVariable[SC.FunctionVariable("a", arrayTy, SC.Exp[SC.ICONST(2)], nothing, Pair{String, SC.Exp}[], IR.ComponentRef(:a))],
                                 SC.Statement[SC.ASSIGN(var(:y), SC.RCONST(0.0), src), loop,
                                              SC.ASSIGN(var(:a, arrayTy), SC.ARRAY_EXP(arrayTy, true, SC.Exp[var(:x), SC.RCONST(2.0)]), src),
                                              SC.ASSIGN(element(:a, 1), SC.RCONST(5.0), src),
                                              SC.ASSIGN(var(:z), add(element(:a, 1), element(:a, 2)), src)], nothing)
  local outs = values(f, Any[SC.UNKNOWN])
  @test outs[1] === 6.0
  @test outs[2] === 7.0
  #= an assignment to the whole array forgets the elements assigned before =#
  local g = SC.MODELICA_FUNCTION("g", SC.FunctionVariable[fv(:x)], SC.FunctionVariable[fv(:z)],
                                 SC.FunctionVariable[SC.FunctionVariable("a", arrayTy, SC.Exp[SC.ICONST(2)], nothing, Pair{String, SC.Exp}[], IR.ComponentRef(:a))],
                                 SC.Statement[SC.ASSIGN(element(:a, 1), SC.RCONST(5.0), src),
                                              SC.ASSIGN(var(:a, arrayTy), SC.ARRAY_EXP(arrayTy, true, SC.Exp[var(:x), SC.RCONST(2.0)]), src),
                                              SC.ASSIGN(var(:z), element(:a, 1), src)], nothing)
  @test values(g, Any[SC.UNKNOWN])[1] === SC.UNKNOWN
  @test values(g, Any[3.0])[1] === 3.0
end
