#=
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

function transformToMTKConditionEquation(cond::DAE.Exp, simCode)
  res = @match cond begin
    DAE.RELATION(e1, DAE.LESS(__), e2) => begin
      :($(expToJuliaExpMTK(e1, simCode)) - $(expToJuliaExpMTK(e2, simCode)) ~ 0)
    end
    #= TODO: Is this correct? =#
    DAE.RELATION(e1, DAE.GREATER(__), e2) => begin
      :(- $(expToJuliaExpMTK(e1, simCode)) - $(expToJuliaExpMTK(e2, simCode)) ~ 0)
    end
    _ => begin
      throw("Operator: " * "'" * string(cond.operator) * "' in: " * string(cond) * " is not supported")
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
Returns:
stateVariables, algebraicVariables, stateVariablesLoop, algebraicVariablesLoop
"""
function separateVariables(simCode)::Tuple
  local stringToSimVarHT = simCode.stringToSimVarHT
  local parameters::Vector = []
  local stateDerivatives::Vector = []
  local stateVariables::Vector = []
  local algebraicVariables::Vector = []
  local discreteVariables::Vector = []
  #= Loop arrays=#
  local stateDerivativesLoop::Vector = []
  local stateVariablesLoop::Vector = []
  local algebraicVariablesLoop::Vector = []
  local discreteVariablesLoop::Vector = []
  #= Separate the variables =#
  for varName in keys(stringToSimVarHT)
    (idx, var) = stringToSimVarHT[varName]
    if simCode.matchOrder[idx] in loop
      local varType = var.varKind
      @match varType  begin
        SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
        SimulationCode.STATE(__) => push!(stateVariablesLoop, varName)
        SimulationCode.PARAMETER(__) => push!(parameters, varName)
        SimulationCode.ALG_VARIABLE(__) => begin
          if idx in simCode.matchOrder
            push!(algebraicVariablesLoop, varName)
          else #= We have a variable that is not contained in continious system =#
            #= Treat discrete variables separate =#
            push!(discreteVariablesLoop, varName)
          end
        end
        #=TODO: Do I need to modify this?=#
        SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivativesLoop, varName)
      end
    else #= Someplace else=#
      local varType = var.varKind
      @match varType  begin
        SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
        SimulationCode.STATE(__) => push!(stateVariables, varName)
        SimulationCode.PARAMETER(__) => push!(parameters, varName)
        SimulationCode.ALG_VARIABLE(__) => begin
          if idx in simCode.matchOrder
            push!(algebraicVariables, varName)
          else #= We have a variable that is not contained in continious system =#
            #= Treat discrete variables separate =#
            push!(discreteVariables, varName)
          end
        end
        #=TODO: Do I need to modify this?=#
        SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
      end
    end
  end
  return (stateVariables, algebraicVariables, stateVariablesLoop, algebraicVariablesLoop)
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

"""
  Get the name of the active model.
  The default is the original model.
"""
function getActiveModel(simCode)
  local activeModelName = simCode.activeModel
  for sm in simCode.subModels
    if sm.name == activeModelName
      return sm
    end
  end
  return activeModelName
end


"""
  TODO: Keeping it simple for now, we assume we only have one argument in the call..
  Also the der as symbol is really ugly..
"""
function DAECallExpressionToMTKCallExpression(pathStr::String, expLst::List,
                                              simCode::SimulationCode.SimCode, ht; varPrefix=varPrefix, derAsSymbol=false)::Expr
  @match pathStr begin
    "der" => begin
      varName = SimulationCode.string(listHead(expLst))
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


"""
  Removes all comments from a given exp
"""
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

"""
  Converts a DAE expression into a MTK expression.
"""
function expToJuliaExpMTK(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE, varSuffix=""; varPrefix="x", derSymbol::Bool=false)::Expr
  hashTable = simCode.stringToSimVarHT
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
      DAE.CREF(cr, _)  => begin
        varName = SimulationCode.string(cr)
        builtin = if varName == "time"
          true
        else
          false
        end
        if ! builtin
          #= If we refer to time, we simply return t instead of a concrete variable =#
          indexAndVar = hashTable[varName]
          varKind::SimulationCode.SimVarType = indexAndVar[2].varKind
          @match varKind begin
            SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
            SimulationCode.STATE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName state"))
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.PARAMETER(__) => quote
              $(LineNumberNode(@__LINE__, "$varName parameter"))
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.ALG_VARIABLE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName, algebraic"))
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.DISCRETE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName, discrete"))
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.OCC_VARIABLE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName, occ variable"))
              $(Symbol(indexAndVar[2].name))
            end
            _ => begin
              @error "Unsupported varKind: $(varKind)"
              fail()
            end
          end
        else #= Currently only time is a builtin variable. Time is represented as t in the generated code =#
          quote
            t
          end
        end
      end
      DAE.UNARY(operator = op, exp = e1) => begin
        o = DAE_OP_toJuliaOperator(op)
        quote
          $(o)($(expToJuliaExpMTK(e1, simCode, varPrefix=varPrefix)))
        end
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        local lhs = expToJuliaExpMTK(e1, simCode, varPrefix=varPrefix, derSymbol = derSymbol)
        local rhs = expToJuliaExpMTK(e2, simCode, varPrefix=varPrefix, derSymbol = derSymbol)
        local op = DAE_OP_toJuliaOperator(op)
        :($op($(lhs), $(rhs)))
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        local operand = expToJuliaExpMTK(e1, simCode, varPrefix=varPrefix, derSymbol = derSymbol)
        local op = DAE_OP_toJuliaOperator(op)
        :($op($(op)))
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        local lhs = expToJuliaExpMTK(e1, simCode, varPrefix=varPrefix, derSymbol = derSymbol)
        local rhs = expToJuliaExpMTK(e2, simCode, varPrefix=varPrefix, derSymbol = derSymbol)
        local op = DAE_OP_toJuliaOperator(op)
        :($op($(lhs), $(rhs)))
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        local lhs = expToJuliaExpMTK(e1, simCode, varPrefix=varPrefix,derSymbol = derSymbol)
        local rhs = expToJuliaExpMTK(e2, simCode,varPrefix=varPrefix, derSymbol = derSymbol)
        local op = DAE_OP_toJuliaOperator(op)
        :($op($(lhs), $(rhs)))
      end
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        throw(ErrorException("If expressions not allowed in backend code"))
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst)  => begin
        #Call as symbol is really ugly.. please fix me :(
        DAECallExpressionToMTKCallExpression(tmpStr, explst, simCode, hashTable; varPrefix=varPrefix, derAsSymbol=derSymbol)
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
  If the system needs to conduct index reduction make sure to inform MTK.
(We avoid structurally simplify for now since that might interfer with some other algorithms)
"""
function performStructuralSimplify(simplify)::Expr
  if (simplify)
    :(reducedSystem = ModelingToolkit.structural_simplify(firstOrderSystem; simplify = false, allow_parameter=true))#ModelingToolkit.dae_index_lowering(firstOrderSystem))
  else
    :(reducedSystem = firstOrderSystem)#ModelingToolkit.structural_simplify(firstOrderSystem; simplify = false, allow_parameter=true))
  end
end

"""
  Generates different constructors for the ODESystem depending on given parameters.
"""
function odeSystemWithEvents(hasEvents, modelName)
  if hasEvents
    :(ODESystem(eqs, t, vars, parameters;
              name=:($(Symbol($modelName))),
              continuous_events = events))
  else
    :(ODESystem(eqs, t, vars, parameters;
              name=:($(Symbol($modelName)))))
  end
end


"""
  Given the variable idx and simCode statically decide if this variable is involved in some event.
"""
function involvedInEvent(idx, simCode)
  return false
end


"""
  This functions evaluate a single DAE-constant:{Bool, Integer, Real, String}.
  If the argument  to this function is not a constant it throws an error.
"""
function evalDAEConstant(daeConstant::DAE.Exp, simCode)
  @match daeConstant begin
    DAE.BCONST(bool) => bool
    DAE.ICONST(int) => int
    DAE.RCONST(real) => real
    DAE.SCONST(tmpStr) => tmpStr
    #= Try to evaluate the expression =#
    DAE.BINARY(__) => begin
      evalDAE_Expression(daeConstant, simCode)
    end
    _ => begin
      local str = string(daeConstant)
      throw("$(str) is not a constant")
    end
  end
end

function evalDAEConstant(daeConstant::DAE.Exp)
  @match daeConstant begin
    DAE.BCONST(bool) => bool
    DAE.ICONST(int) => int
    DAE.RCONST(real) => real
    DAE.SCONST(tmpStr) => tmpStr
    #= Try to evaluate the expression =#
    _ => begin
      local str = string(daeConstant)
      throw("$(str) is not a constant")
    end
  end
end


"""
  Evaluates a simulation code parameter.
  Fails if the function is not a parameter.
"""
function evalSimCodeParameter(v::V, simCode) where V
  @match SimulationCode.SIMVAR(name, _, SimulationCode.PARAMETER(SOME(bindExp)), _) = v
  local val = evalDAEConstant(bindExp, simCode)
  return val
end

"""
 Evalutates the components in a DAE expression (Currently if the components are parameters)
"""
function evalDAE_Expression(expr, simCode)::Expr
  function replaceParameterVariable(exp, ht)
    if Util.isCref(exp)
      local simVar = last(simCode.stringToSimVarHT[string(exp)])
      if ! SimulationCode.isStateOrAlgebraic(simVar)
        @match SimulationCode.SIMVAR(name, _, SimulationCode.PARAMETER(SOME(bindExp)), _) = simVar
        return (bindExp, true, ht)
      else
        @error "Not supported states and algebraics in initial equations is not supported"
        fail()
      end
      println("End one level")
    end
    (exp, true, ht)
  end
  local a = 0
  #= Replaces all known variables in the daeExp =#
  local daeExp = first(Util.traverseExpBottomUp(expr, replaceParameterVariable, 0))
  local jlExpr = expToJuliaExpMTK(daeExp, simCode)
  local evaluatedJLExpr = eval(jlExpr)
  return quote $(evaluatedJLExpr) end
end

"""
  Decide the iv of the condition.
  This currently assumes that the simulation starts at 0.0 (which might not be the case).
This function needs to be improved so that it also evaluates static parameters.
"""
function evalInitialCondition(mtkCond)
  try
    local mtkCondE = eval(mtkCond)
    local v = substitute(mtkCondE.lhs, t => 0.0)
    local ivCond = (v == 0)
    if ivCond == false
      return ivCond
    end      
    if length(v.val.arguments) > 1
      return true
    end
    return ivCond
  catch #= Unable to evaluate at this point in time. =#
    return true
  end
end

function generateIfExpressions(branches, target, resEqIdx::Int, identifier::Int, simCode)
  if branches[target].targets == -1
    return :($(first(deCausalize(branches[target].residualEquations[resEqIdx], simCode))))
  end
  #= Otherwise generate code for the other part =#
  local branch = branches[target]
  local cond = :( $(Symbol("ifCond$(identifier)")) == true )
  #= Assume single equations for now =#
  quote
    ModelingToolkit.IfElse.ifelse($(cond),
                                  $(first(deCausalize(branch.residualEquations[resEqIdx], simCode))),
                                  $(generateIfExpressions(branches, branches[target].targets, resEqIdx, identifier, simCode)))
  end
end

#= TODO. We currently assume residuals that we have made causal and that the original equations are written in a certain form.=#
function deCausalize(eq, simCode)
  @match eq.exp begin
    DAE.BINARY(DAE.RCONST(0.0), _, exp2) => begin
      (:($(expToJuliaExpMTK(exp2, simCode))), :($(expToJuliaExpMTK(eq.exp.exp1, simCode))))
    end
    DAE.BINARY(exp1, _, DAE.RCONST(0.0)) => begin
      (:($(expToJuliaExpMTK(exp1, simCode))), :($(expToJuliaExpMTK(eq.exp.exp2, simCode))))
    end
    DAE.BINARY(exp1, _, exp2) => begin
      (:($(expToJuliaExpMTK(exp2, simCode))), :($(expToJuliaExpMTK(exp1, simCode))))
    end
    _ => begin
      throw("Unsupported equation:" * string(eq))
    end
  end
end


#= TODO. REMOVE =#
function deCausalize2(eq, simCode)
  @match eq.exp begin
    DAE.BINARY(DAE.RCONST(0.0), _, exp2) => begin
      expToJuliaExpMTK(exp2, simCode)
    end
    DAE.BINARY(exp1, _, DAE.RCONST(0.0)) => begin
      expToJuliaExpMTK(exp1, simCode)
    end
    DAE.BINARY(exp1, _, exp2) => begin
      expToJuliaExpMTK(eq.exp, simCode)
    end
    _ => begin
      throw("Unsupported equation:" * string(eq))
    end
  end
end
