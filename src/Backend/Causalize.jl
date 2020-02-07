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

#= ExportAll is not good practice but it makes it so that we do not have to write export after each function :( =#
using ExportAll

import Absyn
import DAE

import BackendDAE
import BackendDAEUtil
import BackendEquation

"""
    Variable can be: Variable, Discrete, Constant and Parameters
    From this create Algebraic and State Variables (Known variable)
    Traverse all equations and locate the variables that are derived.
    These we mark as states
"""
function detectStates(dae::BackendDAE.BackendDAEStructure)
  BackendDAEUtil.mapEqSystems(dae, detectStatesEqSystem)
end

"""
    kabdelhak:
    Detects all states in the system by looking for component references in
    der() calls. Updates all variables with those component references to
    varKind BackendDAE.STATE()
"""
function detectStatesEqSystem(syst::BackendDAE.EqSystem)
  syst = begin
    local vars::BackendDAE.Variables
    local eqs::BackendDAE.EquationArray
    local stateCrefs = Dict{DAE.ComponentRef, Bool}()
    @match syst begin
      BackendDAE.EQSYSTEM(vars, eqs) => begin #= qualified access possible? =#
        for eq in eqs
          (_, stateCrefs) = BackendDAEUtil.traverseEquationExpressions(eq, detectStateExpression, stateCrefs)
        end
        #= Do replacements for stateCrefs =#
        @set syst.orderedVars = updateStates(vars, stateCrefs)
        (syst)
      end
    end
  end
  return syst
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
        #= add state with boolean value that does not matter, it is later onlBackendDAE.BACKEND_DAE(eqs = eqs)y checked if it exists at all =#
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
    to be updated to be a BackendDAE.STATE()
"""
function updateStates(vars::BackendDAE.Variables, stateCrefs::Dict{DAE.ComponentRef, Bool})
  vars = begin
    local varArr::Array{BackendDAE.Var, 1}
    @match vars begin
      BackendDAE.VARIABLES(varArr = varArr) => begin
        for i in 1:arrayLength(varArr)
          varArr[i] = begin
            local cref::DAE.ComponentRef
            local var::BackendDAE.Var
            @match varArr[i] begin
              var && BackendDAE.VAR(varName = cref) where (haskey(stateCrefs, cref)) => begin
                var = @set var.varKind = BackendDAE.STATE(0, NONE(), true)
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
"""
function daeMode(dae::BackendDAE.BackendDAEStructure)
  dae = BackendDAEUtil.mapEqSystems(dae, makeResidualEquations)
end

"""
    kabdelhak:
    Traverser for daeMode() to map all equations of an equation system
"""
function makeResidualEquations(syst::BackendDAE.EqSystem)
  syst = BackendDAEUtil.mapEqSystemEquations(syst, BackendEquation.makeResidualEquation)
end


@exportAll()
end
