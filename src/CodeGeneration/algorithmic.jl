#=
Code generation for algorithmic Modelica.
author:johti17
=#

"""
      Generates algorithmic Modelica Code.
      Returns the generated Julia code + the names of the functions that has been generated.
TODO:
  - Solve function arguments in a more elegant way.
  - Add support for for-loops.
  - Solve the lowering of external functions without string splitting.
"""
function generateFunctions(functions::Vector{SimulationCode.ModelicaFunction})::Tuple{Vector{Expr}, Vector{String}}
  local jFuncs = Expr[]
  local names = String[]
  for func in functions
    local inputs = generateIOL(func.inputs)
    local outputs = generateIOL(func.outputs)
    local f
    inputsJL = if length(inputs) > 1
      tuple(inputs...)
    elseif length(inputs) == 1
      inputs[1]
    else
    end
    @match func begin
      SimulationCode.MODELICA_FUNCTION(__) => begin
        local locals = generateLocals(func.locals)
        local statements = generateStatements(func.statements)
        if inputsJL isa Tuple
          f = quote
            function $(Symbol(func.name))($(inputsJL...))
              $(locals...)
              $(statements...)
              return $(if length(outputs) > 1
                         tuple(outputs...)
                       elseif length(outputs) == 1
                         outputs[1]
                       else
                         nothing
                       end)
            end
          end
        else
          f = quote
            function $(Symbol(func.name))($(inputsJL))
              $(locals...)
              $(statements...)
              return $(if length(outputs) > 1
                         tuple(outputs...)
                       elseif length(outputs) == 1
                         outputs[1]
                       else
                         nothing
                       end)
            end
          end
        end
      end
      SimulationCode.EXTERNAL_MODELICA_FUNCTION(__) => begin
        if inputsJL isa Tuple
          local extCall = Meta.parse(func.libInfo)
          extCall = namespaceifyExternalFunction(extCall)
          f = quote
            function $(Symbol(func.name))($(inputsJL...))
              #= Here a call to the external function should be placed =#
              $(extCall)
              $(if length(outputs) > 1
                  tuple(outputs...)
                elseif length(outputs) == 1
                  outputs[1]
                else
                  nothing
                end)
            end
          end
        else
          local extCall = Meta.parse(func.libInfo)
          extCall = namespaceifyExternalFunction(extCall)
          f = quote
            function $(Symbol(func.name))($(inputsJL))
              #= Here a call to the external function should be placed =#
              $(extCall)
              $(if length(outputs) > 1
                  tuple(outputs...)
                elseif length(outputs) == 1
                  outputs[1]
                else
                  nothing
                end)
            end
          end
        end
      end
    end
    push!(jFuncs, f)
    push!(names, func.name)
  end
  return jFuncs, names
end

function generateIOL(inputs::Vector)
  local jInputs = Symbol[]
  for i in inputs
    local s = DAE_VAR_ToJulia(i)
    #= Complex type, prefixed with void* =#
    push!(jInputs, s)
  end
  return jInputs
end

"""
`generateSignatureForRegistration(inputs::Vector{DAE.VAR})`
This function generates the input signature for calls to Symbolics.register
"""
function generateSignatureForRegistration(inputs::Vector{DAE.VAR})
  local jInputs = Expr[]
  for i in inputs
    @match i.ty begin
      DAE.T_COMPLEX(__) => begin
        local s = DAE_VAR_ToJulia(i)
        #= Complex type, prefixed with void* =#
        push!(jInputs, Expr(:(::), s, :(Ptr{Nothing})))
      end
      _ => begin
        local s = DAE_VAR_ToJulia(i)
        push!(jInputs, Expr(:(::), s, :(Any)))
      end
    end
  end
  return jInputs
end

function generateLocals(inputs::Vector)
  local jInputs = Expr[]
  for i in inputs
    local s = DAE_VAR_ToJulia(i)
    push!(jInputs, Expr(:local, s))
  end
  return jInputs
end

function generateStatements(statements::Union{List{DAE.Statement}, Vector{DAE.Statement}})
  local jStmts = Expr[]
  for s in statements
    stmt = generateStatement(s)
    push!(jStmts, stmt)
  end
  return jStmts
end

function generateStatement(s::DAE.Statement)
  throw("Unsupported stmt:" * string(s))
end

function generateStatement(stmt::DAE.STMT_ASSIGN)
  local lhs = string(stmt.exp1)
  local rhs = expToJuliaExpAlg(stmt.exp)
  return :($(Symbol(lhs)) = $(rhs))
end

function generateStatement(stmt::DAE.STMT_TUPLE_ASSIGN)
  local lhs = string(stmt.expExpLst)
  local rhs = expToJuliaExpAlg(stmt.exp)
  return :($(Symbol(lhs)) = $(rhs))
end

function generateStatement(stmt::DAE.STMT_ASSIGN_ARR)
  #return string(stmt.lhs) * ":=" * string(stmt.exp) * "|" * string(stmt.type_)
  throw("generateStatement(stmt::DAE.STMT_ASSIGN_ARR) not impl")
end

function generateStatement(stmt::DAE.STMT_WHILE, simCode)
  local buffer = IOBuffer()
  local cond = expToJuliaExpAlg(stmt.exp)
  local stmts = generateStatements(stmt.statementLst)
  quote
    while ($(cond))
      $(stmts...)
    end
  end
end

"""
  Generates for statements.
  Currently not implemented.
"""
function generateStatement(stmt::DAE.STMT_FOR)
  throw("TODO: For not implemented yet")
end

"""
  Generates If statements
"""
function generateStatement(stmt::DAE.STMT_IF)::Expr
  local cond = expToJuliaExpAlg(stmt.exp)
  local stmts = generateStatements(stmt.statementLst)
  local res = @match stmt.else_ begin
    DAE.NOELSE(__) => begin
      local expr = Expr(:if, cond)
      local blck = Expr(:block)
      for stmt in stmts
        push!(blck.args, stmt)
      end
      push!(expr.args. blck)
      expr
    end
    DAE.ELSE(__) => begin
      local expr = Expr(:if, cond)
      local blck = Expr(:block)
      for stmt in stmts
        push!(blck.args, stmt)
      end
      push!(expr.args, blck)
      local elseStmts = generateStatement(stmt.else_)
      push!(expr.args, elseStmts)
      expr
    end
    DAE.ELSEIF(__) => begin
      local expr = Expr(:if, cond)
      local blck = Expr(:block)
      for stmt in stmts
        push!(blck.args, stmt)
      end
      push!(expr.args. blck)
      local elseIfs = generateStatement(stmt.else_)
      Expr(:if, cond, stmts, elseIfs)
    end
  end
  return res
end

"""
For the else branch we generate a block and add the statements of the ELSE to this block.
Should never be called at the top level.
"""
function generateStatement(stmt::DAE.ELSE)::Expr
  local block = Expr(:block)
  stmts = generateStatements(stmt.statementLst)
  for stmt in stmts
    push!(block.args, stmt)
  end
  return block
end

"""
For the elseif branch we create an elseif expression.
Similar to the else this should never be called from the top level.
"""
function generateStatement(stmt::DAE.ELSEIF)::Expr
  local cond = expToJuliaExpAlg(stmt.exp)
  local stmts = generateStatments(stmt.statementLst)
  Expr(:elseif, cond, stmts)
end

"""
  Maps a DAE expression to a Julia expression for algorithmic code in Modelica Functions(!).
  Since functions do not use the model HT the original name is preserved for algorithmic generation.
For algorithmic code outside Modelica functions do not call this function.
"""
function expToJuliaExpAlg(@nospecialize(exp::DAE.Exp))::Expr
  local expr::Expr = begin
    local int::Int64
    local real::Float64
    local bool::Bool
    local tmpStr::String
    local cr::DAE.ComponentRef
    local e1::DAE.Exp
    local e2::DAE.Exp
    local e3::DAE.Exp
    local expl::List{DAE.Exp}
    local lstexpl::List{List{DAE.Exp}}
    @match exp begin
      DAE.BCONST(bool) => quote $bool end
      DAE.ICONST(int) => quote $int end
      DAE.RCONST(real) => quote $real end
      DAE.SCONST(tmpStr) => quote $tmpStr end
      DAE.CREF(Absyn.IDENT("time"), _) => begin
        quote t end
      end
      #= Array accesses =#
      DAE.CREF(DAE.CREF_IDENT(ident, identType, subscriptLst), _) where !isempty(subscriptLst)  => begin
        local arrAccess = SimulationCode.string(exp)
        expr = Meta.parse(arrAccess)
      end
      DAE.CREF(cr, _)  => begin
        local varName = SimulationCode.string(cr)
        quote
          $(Symbol(varName))
        end
      end
      DAE.UNARY(operator = op, exp = e1) => begin
        o = CodeGeneration.DAE_OP_toJuliaOperator(op)
        :($(o)($(expToJuliaExpAlg(e1))))
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        local lhs = expToJuliaExpAlg(e1)
        local rhs = expToJuliaExpAlg(e2)
        local op = CodeGeneration.DAE_OP_toJuliaOperator(op)
        :($op($(lhs), $(rhs)))
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        local operand = expToJuliaExpAlg(e1)
        local op = CodeGeneration.DAE_OP_toJuliaOperator(op)
        :($op($(op)))
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        local lhs = expToJuliaExpAlg(e1)
        local rhs = expToJuliaExpAlg(e2)
        local op = CodeGeneration.DAE_OP_toJuliaOperator(op)
        :($op($(lhs), $(rhs)))
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        local lhs = expToJuliaExpAlg(e1)
        local rhs = expToJuliaExpAlg(e2)
        local op = CodeGeneration.DAE_OP_toJuliaOperator(op)
        :($op($(lhs), $(rhs)))
      end
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        local cond = expToJuliaExpAlg(e1)
        local thenExp = expToJuliaExpAlg(e2)
        local elseExp = expToJuliaExpAlg(e2)
        quote
          if $(cond)
            $(thenExp)
          else
            $(elseExp)
          end
        end
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst, attr = attr)  => begin
        local expr = Expr(:call, if !(attr.builtin) #Note extend this when needed
                            Symbol(tmpStr)
                          else
                            Symbol(tmpStr)
                          end)
        local args = map(explst) do arg
          expToJuliaExpAlg(arg)
        end
        expr.args = vcat(expr.args, args)
        quote
          $(expr)
        end
      end
      DAE.CALL(path, expLst, attr) => begin
        local expr = Expr(:call, if !(attr.builtin)
                            Symbol(string(path))
                          else
                            Symbol(string(path))
                          end)
        local args = map(expLst) do arg
          expToJuliaExpAlg(arg)
        end
        expr.args = vcat(expr.args, args)
        expr
      end
      DAE.CAST(ty, exp)  => begin
        quote
          $(generateCastExpressionMTK(ty, exp, simCode, varPrefix))
        end
      end
      _ =>  throw(ErrorException("$exp not yet supported"))
    end
  end
  return expr
end

"""
  Converts a DAE var to an equvivalent Julia repr.
  Simple for now.
"""
function DAE_VAR_ToJulia(v::DAE.VAR)
  local vName = string(v.componentRef)
  Symbol(vName)
end

"""
  Adds OMRuntimeExternalC as a prefix to external calls
"""
function namespaceifyExternalFunction(expr::Expr)
  res = if expr.head == :(=)
    local callExpr = last(expr.args)
    @match Expr(:call, [funcName, y...,z]) = callExpr
    exp = Expr(:call, Expr(:(.), Symbol("OMRuntimeExternalC"), QuoteNode(funcName)), y...,z)
    #= Add the prefixes call to the right hand side of the expression. =#
    expr.args[2] = exp
    expr
  else #Otherwise a side effect call or a call that returns directly.
    @assert expr.head === :call "Invalid call passed to namespaceifyExternalFunction"
    @match Expr(:call, [funcName, y...,z]) = expr
    Expr(:call, Expr(:(.), Symbol("OMRuntimeExternalC"), QuoteNode(funcName)))
  end
  return res
end
