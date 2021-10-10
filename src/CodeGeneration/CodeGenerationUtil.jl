#= /*
* This file is part of OpenModelica.
*
* Copyright (c) 1998-2020, Open Source Modelica Consortium (OSMC),
* c/o Linköpings universitet, Department of Computer and Information Science,
* SE-58183 Linköping, Sweden.
*
* All rights reserved.
*
* THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
* THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
* ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
* RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
* ACCORDING TO RECIPIENTS CHOICE.
*
* The OpenModelica software and the Open Source Modelica
* Consortium (OSMC) Public License (OSMC-PL) are obtained
* from OSMC, either from the above address,
* from the URLs: http:www.ida.liu.se/projects/OpenModelica or
* http:www.openmodelica.org, and in the OpenModelica distribution.
* GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
*
* This program is distributed WITHOUT ANY WARRANTY; without
* even the implied warranty of  MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
* IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
*
* See the full OSMC Public License conditions for more details.
*
=#

#=
# Author: John Tinnerholm (johti17)
=#
import MacroTools


"""
Return string containing the OSMC copyright stuff.
"""
function copyRightString()

  strOut = string("#= /*
* This file is part of OpenModelica.
*
* Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
* c/o Linköpings universitet, Department of Computer and Information Science,
* SE-58183 Linköping, Sweden.
*
* All rights reserved.
*
* THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
* THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
* ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
* RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
* ACCORDING TO RECIPIENTS CHOICE.
*
* The OpenModelica software and the Open Source Modelica
* Consortium (OSMC) Public License (OSMC-PL) are obtained
* from OSMC, either from the above address,
* from the URLs: http:www.ida.liu.se/projects/OpenModelica or
* http:www.openmodelica.org, and in the OpenModelica distribution.
* GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
*
* This program is distributed WITHOUT ANY WARRANTY; without
* even the implied warranty of  MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
* IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
*
* See the full OSMC Public License conditions for more details.
*
=#")
  return strOut
end


"""
johti17
    Transform a condition into a zero crossing function.
    For instance y < 10 -> y - 10
TODO:
    Assumes a Real-Expression.
    Also assume that relation expressions are written as <= That is we go from positive 
    to negative.. 
    Fix this.
"""
function transformToZeroCrossingCondition(conditonalExpression::DAE.Exp)
  res = @match conditonalExpression begin
    DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
      DAE.BINARY(lhs, DAE.SUB(DAE.T_REAL_DEFAULT), rhs)
    end
    DAE.LUNARY(operator = op, exp = e1)  => begin
      DAE.UNARY(DAE.SUB(DAE.T_REAL_DEFAULT), e1)
    end
    DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
      DAE.BINARY(e1, DAE.SUB(DAE.T_REAL_DEFAULT), e2)
    end
    DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
      DAE.BINARY(e1, DAE.SUB(DAE.T_REAL_DEFAULT), e2)
    end
    _ => begin
      conditonalExpression
    end
  end
  return res
end

"""
  Flattens a vector of expressions.
"""
function flattenExprs(eqs::Vector{Expr})
  quote
    $(eqs...)
  end
end


"""
 Convert DAE.Exp into a Julia string. 
"""
function expToJL(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE; varPrefix="x")::String
  hashTable = simCode.stringToSimVarHT
  str = begin
    local int::ModelicaInteger
    local real::ModelicaReal
    local bool::Bool
    local tmpStr::String
    local cr::DAE.ComponentRef
    local e1::DAE.Exp
    local e2::DAE.Exp
    local e3::DAE.Exp
    local expl::List{DAE.Exp}
    local lstexpl::List{List{DAE.Exp}}
    @match exp begin
      DAE.ICONST(int) => string(int)
      DAE.RCONST(real)  => string(real)
      DAE.SCONST(tmpStr)  => (tmpStr)
      DAE.BCONST(bool)  => string(bool)
      DAE.ENUM_LITERAL((Absyn.IDENT(str), int))  => str + "()" + string(int) + ")"
      DAE.CREF(cr, _)  => begin
        varName = SimulationCode.string(cr)
        builtin = if varName == "time"
          true
        else
          false
        end
        if ! builtin
          #= If we refeer to time, we simply return t instead of a concrete variable =#
          indexAndVar = hashTable[varName]
          varKind::SimulationCode.SimVarType = indexAndVar[2].varKind
          @match varKind begin
            SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
            SimulationCode.STATE(__) => "$varPrefix[$(indexAndVar[1])] #= $varName =#"
            SimulationCode.PARAMETER(__) => "p[$(indexAndVar[1])] #= $varName =#"
            SimulationCode.ALG_VARIABLE(__) => "$varPrefix[$(indexAndVar[1])] #= $varName =#"
            SimulationCode.STATE_DERIVATIVE(__) => "dx[$(indexAndVar[1])] #= der($varName) =#"
          end
        else #= Currently only time is a builtin variabe. Time is represented as t in the generated code =#
          "t"
        end
      end
      DAE.UNARY(operator = op, exp = e1) => begin
        ("(" + SimulationCode.string(op) + " " + expToJL(e1, simCode) + ")")
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expToJL(e1, simCode, varPrefix=varPrefix) + " " + SimulationCode.string(op) + " " + expToJL(e2, simCode, varPrefix=varPrefix))
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        ("(" + SimulationCode.string(op) + " " + expToJL(e1, simCode, varPrefix=varPrefix) + ")")
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        (expToJL(e1, simCode, varPrefix=varPrefix) + " " + SimulationCode.string(op) + " " + expToJL(e2, simCode, varPrefix=varPrefix))
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        (expToJL(e1, simCode, varPrefix=varPrefix) + " " + SimulationCode.string(op) + " " + expToJL(e2, simCode,varPrefix=varPrefix))
      end
      #=TODO?=#
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        "if" + expToJL(e1, simCode, varPrefix=varPrefix) + "\n" + expToJL(e2, simCode,varPrefix=varPrefix) + "else\n" + expToJL(e3, simCode,varPrefix=varPrefix) + "\nend"
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = expl)  => begin
        #=
          TODO: Keeping it simple for now, we assume we only have one argument in the call
          We handle derivitives seperatly
        =#
        varName = SimulationCode.string(listHead(expl))
        (index, type) = hashTable[varName]
        @match tmpStr begin
        "der" => "dx[$index]  #= der($varName) =#"
          "pre" => begin
            indexForVar = hashTable[varName][1]
            "(integrator.u[$(indexForVar)])"
          end
          "edge" =>  begin
             indexForVar = hashTable[varName][1]
             string(tuple(map((x) -> expToJL(x, simCode, varPrefix=varPrefix), expl)...)...) + " && ! integrator.uprev[$(indexForVar)]"
          end
          _  =>  begin
            tmpStr *= string(tuple(map((x) -> expToJL(x, simCode, varPrefix=varPrefix), expl)...)...)
          end
        end
      end
      DAE.CAST(exp = e1)  => begin
         expToJL(e1, simCode)
      end
      DAE.ARRAY(DAE.T_ARRAY(DAE.T_BOOL(__)), scalar, array) => begin
        local arrayExp = "#= Array exp=# reduce(|, ["
        for e in array
          arrayExp *= expToJL(e, simCode) + ","
        end
        arrayExp *= "])"
      end
      _ =>  throw(ErrorException("$exp not yet supported"))
    end
  end
  return "(" + str + ")"
end


function DAE_OP_toJuliaOperator(op::DAE.Operator)
    return @match op begin
      DAE.ADD() => :+
      DAE.SUB() => :-
      DAE.MUL() => :*
      DAE.DIV() => :/
      DAE.POW() => :^
      DAE.UMINUS() =>  :-
      DAE.UMINUS_ARR() => :-
      DAE.ADD_ARR() => :+
      DAE.SUB_ARR() => :-
      DAE.MUL_ARR() => :*
      DAE.DIV_ARR() => :/
      DAE.MUL_ARRAY_SCALAR() => :*
      DAE.ADD_ARRAY_SCALAR() => :+
      DAE.SUB_SCALAR_ARRAY() =>  :-
      DAE.MUL_SCALAR_PRODUCT() => :*
      DAE.MUL_MATRIX_PRODUCT() => :*
      DAE.DIV_ARRAY_SCALAR() => :/
      DAE.DIV_SCALAR_ARRAY() => :/
      DAE.POW_ARRAY_SCALAR() => :^
      DAE.POW_SCALAR_ARRAY() => :^
      DAE.POW_ARR() => :^
      DAE.POW_ARR2() => :^
      DAE.AND() => :&&
      DAE.OR() => :(||)
      DAE.NOT() => :(!)
      DAE.LESS() => :<
      DAE.LESSEQ() => :<=
      DAE.GREATER() => :>
      DAE.GREATEREQ() => :>=
      DAE.EQUAL() => :(=)
      DAE.NEQUAL() => :(!=)
      DAE.USERDEFINED() => throw("Unknown operator: Userdefined")
      _ => throw("Unknown operator")
    end
end


"
  TODO: Keeping it simple for now, we assume we only have one argument in the call
  We handle derivitives seperatly
"
function DAECallExpressionToJuliaCallExpression(pathStr::String, expLst::List, simCode, ht; varPrefix=varPrefix)::Expr
  @match pathStr begin
    "der" => begin
      varName = SimulationCode.string(listHead(expLst))
      (index, _) = ht[varName]
      quote
        dx[$(index)] #= der($varName) =#
      end
    end
    "pre" => begin
      varName = SimulationCode.string(listHead(expLst))
      (index, _) = ht[varName]
      indexForVar = ht[varName][1]
      quote 
        (integrator.u[$(indexForVar)])
      end
    end
    _  =>  begin
      funcName = Symbol(pathStr)
      argPart = tuple(map((x) -> expToJuliaExp(x, simCode, varPrefix=varPrefix), expLst)...)
      quote 
        $(funcName)($(argPart...))
      end
    end
  end
end

"
  TODO: Keeping it simple for now, we assume we only have one argument in the call..
  Also the der as symbol is really ugly..
"
function DAECallExpressionToMTKCallExpression(pathStr::String, expLst::List,
                                              simCode::SimulationCode.SimCode, ht; varPrefix=varPrefix, derAsSymbol=false)::Expr
  @match pathStr begin
    "der" => begin
      varName = SimulationCode.string(listHead(expLst))
      (index, _) = ht[varName]
      if derAsSymbol
        quote
          $(Symbol("der_$(varName)"))
        end
      else
        quote
          der($(Symbol(varName)))
        end
      end
    end
    _  =>  begin
      funcName = Symbol(pathStr)
      argPart = tuple(map((x) -> expToJuliaExpMTK(x, simCode), expLst)...)
      quote 
        $(funcName)($(argPart...))
      end
    end
  end
end

"
  Removes all comments from a given exp
"
function stripComments(ex::Expr)::Expr
  return Base.remove_linenums!(ex)
end

"""
Transforms:
  <name>[<index>] -> <name>_index
"""
function arrayToSymbolicVariable(arrayRepr::Expr)::Expr
  MacroTools.postwalk(arrayRepr) do x
    MacroTools.@capture(x, T_[index_]) || return let
      x
    end
    return let
      local newVarStr::String = "$(T)_$(index)"
      local newVar = Symbol(newVarStr)
      return newVar
    end
  end
end

"
 Removes all the redudant blocks from a generated expression
"
function stripBeginBlocks(e)::Expr
  MacroTools.postwalk(e) do x
    return MacroTools.unblock(x)
  end
end

"""
Transforms:s
  <name>_index -> <name>[index]
"""
const pattern = r".*_[0-9]+"
function symbolicVariableToArrayRef(e::Expr)::Expr
  MacroTools.postwalk(e) do x
      x isa Symbol || return let
        x
    end
    return let
      if match(pattern, "$x") == nothing
        return x
      end
      local splitted = split("$x", "_")
      local res = Meta.parse("$(splitted[1])[$(splitted[2])]")
      return res
    end
  end
end

"""
  Utility function, traverses a DAE exp. Variables are saved in the supplied variables array
  (Note that variables here refers to parameters as well)
"""
function getVariablesInDAE_Exp(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE, variables::Set)
  local  hashTable = simCode.stringToSimVarHT
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
    #= These are not varaiables, so we simply return what we have collected thus far. =#
    DAE.BCONST(bool) => variables
    DAE.ICONST(int) => variables
    DAE.RCONST(real) => variables
    DAE.SCONST(tmpStr) => variables
    DAE.CREF(cr, _) where SimulationCode.string(cr)=="time" => begin
      push!(variables, Symbol("t"))
    end
    DAE.CREF(cr, _)  => begin
      varName = SimulationCode.string(cr)
      indexAndVar = hashTable[varName]
      push!(variables, Symbol(varName))
      varKind::SimulationCode.SimVarType = indexAndVar[2].varKind
      @match varKind begin
        SimulationCode.STATE(__) || SimulationCode.PARAMETER(__) || SimulationCode.ALG_VARIABLE(__) => begin
          push!(variables, Symbol(varName))
        end           
        _ => begin
          @error "Unsupported varKind: $(varKind)"
          throw()
        end
      end
    end
    DAE.UNARY(operator = op, exp = e1) => begin
      getVariablesInDAE_Exp(e1, simCode, variables)
    end
    DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) || DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) || DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
      getVariablesInDAE_Exp(e1, simCode, variables)
      getVariablesInDAE_Exp(e2, simCode, variables)        
    end
    DAE.LUNARY(operator = op, exp = e1)  => begin
      getVariablesInDAE_Exp(e1, simCode, variables)
    end
    DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
      throw(ErrorException("If expressions not allowed in backend code"))
    end
    #= Should not introduce anything new..  I am a idiot - John 2021=#
    DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst)  => begin
      #TODO only assumes one argument
      getVariablesInDAE_Exp(listHead(explst), simCode, variables)
    end
    DAE.CAST(ty, exp)  => begin
      getVariablesInDAE_Exp(exp, simCode, variables)
    end
    _ =>  throw(ErrorException("$exp not yet supported"))
  end
end

function isCycleInSCCs(sccs)
  for sc in sccs
    if length(sc) > 1
      return true
    end
  end
  return false
end

function getCycleInSCCs(sccs)
  for sc in sccs
    if length(sc) > 1
      return sc
    end
  end
  return []
end
