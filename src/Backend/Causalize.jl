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

module Causalize

using MetaModelica
using Setfield
using ExportAll

import Absyn
import BDAE
import BDAEUtil
import BackendEquation
import DAE

"""
    Variable can be: Variable, Discrete, Constant and Parameters
    From this create Algebraic and State Variables (Known variable)
    Traverse all equations and locate the variables that are derived.
    These we mark as states
"""
function detectStates(dae::BDAE.BDAEStructure)
  BDAEUtil.mapEqSystems(dae, detectStatesEqSystem)
end

"""
  Replaces all if expressions with a temporary variable.
  Adds an equation assigning this variable to the set of equations.
"""
function detectIfEquations(dae::BDAE.BDAEStructure)
  BDAEUtil.mapEqSystems(dae, detectIfEquationsEqSystem)
end

"""
    kabdelhak:
    Detects all states in the system by looking for component references in
    der() calls.
    Updates all variables with those component references to
    varKind BDAE.STATE()
"""
function detectStatesEqSystem(syst::BDAE.EqSystem)::BDAE.EqSystem
  syst = begin
    local vars::BDAE.Variables
    local eqs::Array
    local stateCrefs = Dict{DAE.ComponentRef, Bool}()
    @match syst begin
      BDAE.EQSYSTEM(vars, eqs) => begin
        for eq in eqs
          (_, stateCrefs) = BDAEUtil.traverseEquationExpressions(eq, detectStateExpression, stateCrefs)
        end
        #= Do replacements for stateCrefs =#
        @set syst.orderedVars = updateStates(vars, stateCrefs)
        syst
      end
    end
  end
  return syst
end

"""
  Detects if equations.
  Returns new temporary variables and an array of equations
"""
function detectIfEquationsEqSystem(syst::BDAE.EqSystem)::BDAE.EqSystem
  syst = begin
    local vars::BDAE.Variables
    local eqs::Array
    local tmpVarToElement = Dict{BDAE.VAR, BDAE.IF_EQUATION}()
    @match syst begin
      BDAE.EQSYSTEM(vars, eqs) => begin
        for eq in eqs
          (_, tmpVarToElement) = BDAEUtil.traverseEquationExpressions(eq,
                                                                            replaceIfExpressionWithTmpVar,
                                                                            tmpVarToElement)
        end
        #= Append the new variables to the list of variables =#
        local newVariables = collect(keys(tmpVarToElement))
        local newEquations = collect(values(tmpVarToElement))
        @debug "After conversion. New Variables: $newVariables, New Equations: $newEquations"
        @set syst.orderedEqs = vcat(syst.orderedEqs, newEquations)
        @set syst.orderedVars.varArr = vcat(syst.orderedVars.varArr, newVariables)
        syst
      end
    end
  end
  return syst
end

let
"""
  Detects if expression.
  We replace the if expression with our temporary variable.
  These variables are assigned in newly created if equations that we add to the tmpVarToElement::Dict.
Thus we create the mapping
  tmpVar -> equation it is assigned in
"""
local tick = 0
global replaceIfExpressionWithTmpVar
  function replaceIfExpressionWithTmpVar(exp::DAE.Exp, tmpVarToElement::Dict{BDAE.VAR, BDAE.IF_EQUATION})
    (newExp, cont, tmpVarToElement) = begin
      #= All these temporary variables are REAL numbers for now =#
      local varType = DAE.T_REAL_DEFAULT
      local varName = "tmp$tick"
      local var::DAE.ComponentRef = DAE.CREF_IDENT(varName, varType, nil)
      local emptySource = DAE.emptyElementSource
      local attr = BDAE.EQ_ATTR_DEFAULT_UNKNOWN
      @match exp begin
        DAE.IFEXP(cond, expThen, expElse) => begin
          local varAsCREF::DAE.CREF = DAE.CREF(var, varType)
          local backendVar = BDAE.VAR(varName, BDAE.VARIABLE(), varType)
          tmpVarToElement[backendVar] = BDAE.IF_EQUATION(list(cond),
                                             list(BDAE.EQUATION(varAsCREF,expThen, attr, emptySource)),
                                             list(BDAE.EQUATION(varAsCREF, expElse, attr,emptySource)))
          (varAsCREF, true, tmpVarToElement)
        end
        _ => begin
          (exp, true, tmpVarToElement)
        end
      end
    end
    #= Note we replace the if expression with our temporary variable =#
    return (newExp, cont, tmpVarToElement)
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
        #= add state with boolean value that does not matter, it is later onlBDAE.BACKEND_DAE(eqs = eqs)y checked if it exists at all =#
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

"""
    kabdelhak:
    Traverses all variables and uses a hashmap to determine if a variable needs
    to be updated to be a BDAE.STATE()
"""
function updateStates(vars::BDAE.Variables, stateCrefs::Dict{DAE.ComponentRef, Bool})
  vars = begin
    local varArr::Array{BDAE.Var, 1}
    @match vars begin
      BDAE.VARIABLES(varArr = varArr) => begin
        for i in 1:arrayLength(varArr)
          varArr[i] = begin
            local cref::DAE.ComponentRef
            local var::BDAE.Var
            @match varArr[i] begin
              var && BDAE.VAR(varName = cref) where (haskey(stateCrefs, cref)) => begin
                var = @set var.varKind = BDAE.STATE(0, NONE(), true)
                (var)
              end
              _ => begin
                (varArr[i])
              end
            end
          end
        end
        @set vars.varArr = varArr
        (vars)
      end
    end
  end
end

"""
    kabdelhak:
    Residualize every equation in each system of the dae by subtracting the rhs
    from the lhs.
    (daeMode)
"""
function residualizeEveryEquation(dae::BDAE.BDAEStructure)
  dae = BDAEUtil.mapEqSystems(dae, makeResidualEquations)
end

"""
    kabdelhak:
    Traverser for daeMode() to map all equations of an equation system
"""
function makeResidualEquations(syst::BDAE.EqSystem)
  syst = BDAEUtil.mapEqSystemEquations(syst, BackendEquation.makeResidualEquation)
end

@exportAll()
end
