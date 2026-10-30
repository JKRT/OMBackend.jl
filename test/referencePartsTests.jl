#= OMJLIR's references by their parts: a name is derived from them (`flatName`, `legacyName`), never taken apart.
   A slice of an array split along its inner dimensions (`T_SLICE`), and the algorithm lifter's questions by parts:
   which references overlap, the wholes an assignment's target is a part of, the arrays equations constrain. The
   reduced models' references (state sets, derivatives, array form) are checked in indexReductionTests
   (`referenceProblems`). =#
using Test
import OMBackend

@testset "References by their parts" begin
  local IRp = OMBackend.OMJLIR
  local P = OMBackend.Passes
  local part = (n::Symbol, subs::Vector{Int} = Int[], ty::IRp.ModelicaType = IRp.T_REAL()) ->
    IRp.CrefPart(n, IRp.Subscript[IRp.INDEX(IRp.INTEGER(k)) for k in subs], ty)
  local ref = (ps::IRp.CrefPart...) -> IRp.ComponentRef(IRp.CrefPart[ps...])

  @testset "a slice's inner subscripts in its type" begin
    local T = IRp.T_ARRAY(IRp.T_REAL(), IRp.Dimension[IRp.DIM_SIZE(3), IRp.DIM_SIZE(3)])
    local G = IRp.CrefPart(Symbol("a|b"), IRp.Subscript[], IRp.T_COMPONENTS(:C, [:a, :b]))
    local slice = P.IndexReduction.sliceCref(ref(G, part(:f), part(:T, Int[], T)), [2, 1])
    @test slice.name === Symbol("a|b.f.T[2,1]") && slice.parts[end].ty == IRp.T_SLICE([2, 1], T)
    @test P.legacyName(slice) == "a|b_f_T[2][1]"
    #= an element of the slice: the instance's subscript on the group's part (named as the member), a derivative's
       on its last part =#
    local element = IRp.ComponentRef(IRp.CrefPart[IRp.CrefPart(G.name, IRp.Subscript[IRp.INDEX(IRp.INTEGER(2))], G.ty); slice.parts[2:end]])
    @test element.name === Symbol("b.f.T[2,1]") && P.unsubscriptedName(element) === slice.name
    @test IRp.withLastSubscripts(IRp.derivativeCref(slice), IRp.Subscript[IRp.INDEX(IRp.INTEGER(2))]).name === Symbol("\$DER.a|b.f.T[2,1][2]")
  end

  @testset "references overlap by their parts" begin
    @test IRp.crefsOverlap(ref(part(:buf)), ref(part(:buf, [2])))
    @test IRp.crefsOverlap(ref(part(:r)), ref(part(:r), part(:x)))
    @test IRp.crefsOverlap(ref(part(:r, [1])), ref(part(:r, [1]), part(:x)))
    @test !IRp.crefsOverlap(ref(part(:buf, [1])), ref(part(:buf, [12])))
    @test !IRp.crefsOverlap(ref(part(:x)), ref(part(:xy)))
    @test !IRp.crefsOverlap(ref(part(:m, [1])), ref(part(:m, [1, 2])))
  end

  @testset "an algorithm's wholes; the arrays equations constrain" begin
    @test P.wholesOf(ref(part(:r), part(:a, [2]), part(:b))) == Symbol[:r, Symbol("r.a"), Symbol("r.a[2]")]
    @test P.wholesOf(ref(part(:x, [2]))) == Symbol[:x]
    local eq = IRp.SCALAR_EQUATION(IRp.CREF(ref(part(:c, [1]), part(:T)), IRp.T_REAL()), IRp.REAL(1.0), IRp.T_REAL(),
                                   IRp.NO_SOURCE, IRp.DEFAULT_ATTRIBUTES)
    local constrained = P.constrainedBy(IRp.Equation[eq])
    @test Symbol("c[1].T") in constrained.names && :c in constrained.arrays && !(Symbol("c[1].T") in constrained.arrays)
  end
end
