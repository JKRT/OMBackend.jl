#= /*
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
=#

module BDAEUtil

using ExportAll
using MetaModelica
using Setfield

import ..BDAE
import ..BackendEquation
import OMBackend
import ..Util
import Absyn
import DAE

"""
This function converts an array of variables to the BDAE variable structure
"""
function convertVarArrayToBDAE_Variables(vars::Array{BDAE.Var})::BDAE.Variables
  local variables::BDAE.Variables = begin
    BDAE.VARIABLES([i for i in vars])
  end
  return variables
end

function createEqSystem(vars::BDAE.Variables, eqs::Array)::BDAE.EQSYSTEM
  (BDAE.EQSYSTEM(vars,
                 eqs,
                 NONE(),
                 NONE(),
                 NONE(),
                 nil,
                 BDAE.UNKNOWN_PARTITION(),
                 BackendEquation.emptyEqns()))
end

"""
  Traverse and update a given structure BDAE.BDAEStructure given a traversalOperation and optional arguments
"""
function mapEqSystems(dae::BDAE.BACKEND_DAE, traversalOperation::Function, args...)
  dae = begin
    local eqs::Array{BDAE.EqSystem, 1}
    @match dae begin
      BDAE.BACKEND_DAE(eqs = eqs) => begin
        for i in 1:arrayLength(eqs)
          eqs[i] = traversalOperation(eqs[i], args...)
        end
        @assign dae.eqs = eqs
        (dae)
      end
      _ => begin
        (dae)
      end
    end
  end
end

function mapEqSystems(dae::BDAE.BACKEND_DAE, traversalOperation::Function)
  dae = begin
    local eqs::Array{BDAE.EqSystem, 1}
    @match dae begin
      BDAE.BACKEND_DAE(eqs = eqs) => begin
        for i in 1:arrayLength(eqs)
          eqs[i] = traversalOperation(eqs[i])
        end
        @assign dae.eqs = eqs
        (dae)
      end
      _ => begin
        (dae)
      end
    end
  end
end

function mapEqSystemEquations(syst::BDAE.EQSYSTEM, traversalOperation::Function)
  syst = begin
    local eqs::Array{BDAE.Equation,1}
    @match syst begin
      BDAE.EQSYSTEM(orderedEqs = eqs) => begin
        for i in 1:arrayLength(eqs)
          eqs[i] = traversalOperation(eqs[i])
        end
        @assign syst.orderedEqs = eqs
        (syst)
      end
    end
  end
end


function mapEqSystemEquationsNoUpdate(syst::BDAE.EQSYSTEM, traversalOperation::Function, extArg)
  extArg = begin
    local eqs::Array{BDAE.Equation,1}
    @match syst begin
      BDAE.EQSYSTEM(orderedEqs = eqs) => begin
        for i in 1:arrayLength(eqs)
          extArg = traversalOperation(eqs[i], extArg)
        end
        (extArg)
      end
    end
  end
end

function mapEqSystemVariablesNoUpdate(syst::BDAE.EQSYSTEM, traversalOperation::Function, extArg)
  extArg = begin
    local varArr::Array{BDAE.Var,1}
    @match syst begin
      BDAE.EQSYSTEM(orderedVars = BDAE.VARIABLES(varArr = varArr)) => begin
        for i in 1:arrayLength(varArr)
          extArg = traversalOperation(varArr[i], extArg)
        end
        (extArg)
      end
    end
  end
  return extArg
end

function traverseEquationExpressions(eq::BDAE.Equation,
                                     traversalOperation::Function,
                                     extArg::T)::Tuple{BDAE.Equation,T} where{T}
   (eq, extArg) = begin
     local lhs::DAE.Exp
     local rhs::DAE.Exp
     local cref::DAE.ComponentRef
     @match eq begin
       BDAE.EQUATION(lhs, rhs) => begin
         (lhs, extArg) = Util.traverseExpTopDown(lhs, traversalOperation, extArg)
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @assign eq.lhs = lhs
         @assign eq.rhs = rhs
         (eq, extArg)
       end
       BDAE.SOLVED_EQUATION(componentRef = cref, exp = rhs) => begin
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @assign eq.rhs = rhs;
         (eq, extArg)
       end
       BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @assign eq.exp = rhs;
         (eq, extArg)
       end
       BDAE.IF_EQUATION(__) => begin
         for eqLst in eq.eqnstrue
           for equation in eqLst
             traverseEquationExpressions(equation, traversalOperation, extArg)
           end
         end
         for equation in eq.eqnsfalse
           traverseEquationExpressions(equation, traversalOperation, extArg)
         end        
         #         throw("Failed to traverse if equation!")
         (eq, extArg)
       end
       _ => begin
         (eq, extArg)
       end
     end
   end
end

"""
Directly maps the DAE type to the BDAE type.
Before casualisation we do not know if variables are state or not.
"""
function DAE_VarKind_to_BDAE_VarKind(kind::DAE.VarKind)::BDAE.VarKind
  @match kind begin
    DAE.VARIABLE(__) => BDAE.VARIABLE()
    DAE.DISCRETE(__) => BDAE.DISCRETE()
    DAE.PARAM(__) => BDAE.PARAM()
    DAE.CONST(__) => BDAE.CONST()
  end
end

function isStateOrVariable(kind::BDAE.VarKind)
  res = @match kind begin
    BDAE.VARIABLE(__) => true
    BDAE.STATE(__) => true
    _ => false
  end
  return res
end

function isVariable(kind::BDAE.VarKind)
  res = @match kind begin
    BDAE.VARIABLE(__) => true
    _ => false
  end
  return res
end

function isState(kind::BDAE.VarKind)
  res = @match kind begin
    BDAE.STATE(__) => true
    _ => false
  end
  return res
end


function isWhenEquation(eq::BDAE.Equation)
  @match eq begin
    BDAE.WHEN_EQUATION(__) => true
    _ => false
  end
end

"""
    kabdelhak:
    Detects if a given expression is a der() call and adds the corresponding
    cref to a hashmap
"""
function detectStateExpression(exp::DAE.Exp, stateCrefs::Dict{DAE.ComponentRef, Bool})
  local cont::Bool
  local outCrefs = stateCrefs
  (outCrefs, cont) = begin
    local state::DAE.ComponentRef
    @match exp begin
      DAE.CALL(Absyn.IDENT("der"), DAE.CREF(state) <| _ ) => begin
        #= Adds a state with boolean value that does not matter, 
           it is later  BDAE.BACKEND_DAE(eqs = eqs) checked if it exists at all =#
        outCrefs[state] = true
        (outCrefs, true)
      end
      _ => begin
        (outCrefs, true)
      end
    end
  end
  return (exp, cont, outCrefs)
end

#=TODO. Did I do something stupid down below here.. ?=#

function countAllUniqueVariablesInSetOfEquations(eqs::Vector{RES_EQ}, vars::Vector{VAR}) where {RES_EQ, VAR}
  vars = Set()
  for eq in eqs
    varsForEq = getAllVariables(eq, vars)
    for v in varsForEq
      push!(vars, v)
    end
  end
  return length(vars)
end

"""
  Author:johti17
  input: Backend Equation, eq
  input: All existing variables
  output All variable in that specific equation
"""
function getAllVariables(eq::BDAE.RESIDUAL_EQUATION, vars::Array{BDAE.Var})::Array{DAE.ComponentRef}
  local componentReferences::List = Util.getAllCrefs(eq.exp)
  local stateCrefs = Dict{DAE.ComponentRef, Bool}()
  (_, stateElements)  = traverseEquationExpressions(eq, detectStateExpression, stateCrefs)
  local stateElementArray = collect(keys(stateElements))
  local componentReferencesNotStates::Array = [componentReferences...]
  local componentReferencesArr::Array = [componentReferences..., stateElementArray...]
  variablesInEq::Array = []
  for var in vars
    local vn = var.varName
    if vn in componentReferencesNotStates && isVariable(var.varKind)
      push!(variablesInEq, vn)
    elseif vn in stateElementArray
      push!(variablesInEq, vn)
    else
      continue
    end
  end
  return variablesInEq
end

"""
  Author:johti17
  input: Backend Equation, eq
  input: All existing variables
  output All variable in that specific equation except the state variables
"""
function getAllVariablesExceptStates(eq::BDAE.RESIDUAL_EQUATION, vars::Array{BDAE.Var})::Array{DAE.ComponentRef}
  local componentReferences::List = Util.getAllCrefs(eq.exp)
  local componentReferencesArr::Array = [componentReferences...]
  local varNames = [v.varName for v in vars]
  variablesInEq::Array = []
  for vn in varNames
    if vn in componentReferencesArr
      push!(variablesInEq, vn)
    end
  end
  return variablesInEq
end

function isArray(cref::DAE.ComponentRef)::Bool
  @match cref begin
    DAE.OPTIMICA_ATTR_INST_CREF(__) || DAE.WILD(__) => false
    _ => begin
      typeof(cref.identType) == DAE.T_ARRAY
    end    
  end
end

function getSubscriptAsIntArray(dims)::Array
  local dimIndices = []
  for d in dims
    if ! (typeof(d) == DAE.DIM_INTEGER)
      throw("Non integers dimensions for arrays are not supported by OMBackend. Variable was $(string(v))")
    else
      push!(dimIndices, d.integer)
    end
  end
  return dimIndices
end

function getSubscriptAsUnicodeString(subscriptLst)::String
  local subscripts = subscriptLst
  local subscriptStr = ""
  for s in subscripts
    @assert(typeof(s) == DAE.INDEX, "DAE.INDEX: is expected for $(s)")
    #= Here I assume integer index!=#
    local indexAsInt::Integer = s.exp.integer
    local result = 0
    local tmp = 0            
    subscriptStr *= getIndexAsAUnicodeString(s)
  end
  return subscriptStr
end

function getIndexAsAUnicodeString(idx::DAE.INDEX)
  local indexAsInt::Integer = idx.exp.integer
  local subscriptStr = ""
  return getIntAsUnicodeSubscript(indexAsInt)
end

"""
input: 100 
output \"₁₀₀\"
"""
function getIntAsUnicodeSubscript(i::Integer)
  local subscriptStr = ""
  while i > 0
    tmp = i % 10
    subscriptStr  *= OMBackend.latexSymbols["\\_" + string(tmp)]
    i =  i ÷ 10
  end
  #= Since the code above generates the string in the reverse order it needs to be re reversed=#
  return reverse(subscriptStr)
end

include("backendDump.jl")
@exportAll()
end
