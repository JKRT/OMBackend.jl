"""
  Generic BDAE expression
"""
abstract type Exp
struct ICONST <: Exp
  integer #= Integer constants =#::Integer
end
struct RCONST <: Exp
  real #= Real constants =#::Real
end
struct SCONST <: Exp
  string #= String constants =#::String
end
struct BCONST <: Exp
  bool #= Bool constants =#::Bool
end
struct CLKCONST <: Exp
  clk #= Clock kinds =#::ClockKind
end
struct ENUM_LITERAL <: Exp
  name::Absyn.Path
  index::Integer
end
struct CREF <: Exp
  componentRef::ComponentRef
  ty::Type
end
struct BINARY <: Exp
  exp1::Exp
  operator::Operator
  exp2::Exp
end
struct UNARY <: Exp
  operator::Operator
  exp::Exp
end
struct LBINARY <: Exp
  exp1::Exp
  operator::Operator
  exp2::Exp
end
struct LUNARY <: Exp
  operator::Operator
  exp::Exp
end
struct RELATION <: Exp
  exp1::Exp
  operator::Operator
  exp2::Exp
  index #= Use -1 as a default; other indexes are used in the backend for some silly reasons =#::Integer
  optionExpisASUB::Option{Tuple{Exp, Integer, Integer}}
end

struct CALL <: Exp
  path::Absyn.Path
  expLst::List{Exp}
  attr::CallAttributes
end
struct RECORD <: Exp
  path::Absyn.Path
  exps #= component values =#::List{Exp}
  comp #= component name =#::List{String}
  ty::Type
end
struct PARTEVALFUNCTION <: Exp
  path::Absyn.Path
  expList::List{Exp}
  ty::Type
  origType::Type
end
struct ARRAY <: Exp
  ty::Type
  scalar #= scalar for codegen =#::Bool
  array #= Array constructor, e.g. {1,3,4} =#::List{Exp} <: Exp
end
struct MATRIX <: Exp
  ty::Type
  integer #= Size of the first dimension =#::Integer
  matrix::List{List{Exp}}
end
struct RANGE <: Exp
  ty #= the (array) type of the expression =#::Type
  start #= start value =#::Exp
  step #= step value =#::Option{Exp}
  stop #= stop value =#::Exp
end
struct TUPLE <: Exp
  PR #= PR. Tuples, used in func calls returning several
  arguments =#::List{Exp}
end
struct CAST <: Exp
  ty #= This is the full type of this expression, i.e. ET_ARRAY(...) for arrays and matrices =#::Type
  exp::Exp
end
struct ASUB <: Exp
  exp::Exp
  sub::List{Exp}
end
struct TSUB <: Exp
  exp::Exp
  ix::Integer
  ty::Type
end
struct RSUB <: Exp
  exp::Exp
  ix::Integer
  #=  Used when generating code for MetaModelica records
  =#
  fieldName::String
  ty::Type
end
struct SIZE <: Exp
  exp::Exp
  sz::Option{Exp}
end
