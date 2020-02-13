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

module BackendDAEUtil

using MetaModelica

using ExportAll
using Setfield

import DAE
import Util

import BackendDAE
import BackendEquation

"""
This function converts an array of variables to the BackendDAE variable structure
"""
function convertVarArrayToBackendDAE_Variables(vars::Array{BackendDAE.Var})::BackendDAE.Variables
  local variables::BackendDAE.Variables = begin
    BackendDAE.VARIABLES([i for i in vars])
  end
  return variables
end

function createEqSystem(vars::BackendDAE.Variables, eqs::BackendDAE.EquationArray)
  (BackendDAE.EQSYSTEM(vars, eqs, NONE(), NONE(), NONE(),
                       BackendDAE.NO_MATCHING(), nil,
                       BackendDAE.UNKNOWN_PARTITION(),
                       BackendEquation.emptyEqns()))
end

function mapEqSystems(dae::BackendDAE.BackendDAEStructure, traversalOperation::Function)
   dae = begin
     local eqs::Array{BackendDAE.EqSystem, 1}
     @match dae begin
       BackendDAE.BACKEND_DAE(eqs = eqs) => begin
         for i in 1:arrayLength(eqs)
           eqs[i] = traversalOperation(eqs[i])
         end
         @set dae.eqs = eqs
         (dae)
       end
       _ => begin
         (dae)
       end
     end
   end
end

function mapEqSystemEquations(syst::BackendDAE.EqSystem, traversalOperation::Function)
  syst = begin
    local eqs::Array{BackendDAE.Equation,1}
    @match syst begin
      BackendDAE.EQSYSTEM(orderedEqs = eqs) => begin
        for i in 1:arrayLength(eqs)
          eqs[i] = traversalOperation(eqs[i])
        end
        @set syst.orderedEqs = eqs
        (syst)
      end
    end
  end
end


function mapEqSystemEquationsNoUpdate(syst::BackendDAE.EqSystem, traversalOperation::Function, extArg)
  extArg = begin
    local eqs::Array{BackendDAE.Equation,1}
    @match syst begin
      BackendDAE.EQSYSTEM(orderedEqs = eqs) => begin
        for i in 1:arrayLength(eqs)
          extArg = traversalOperation(eqs[i], extArg)
        end
        (extArg)
      end
    end
  end
end

function mapEqSystemVariablesNoUpdate(syst::BackendDAE.EqSystem, traversalOperation::Function, extArg)
  extArg = begin
    local varArr::Array{BackendDAE.Var,1}
    @match syst begin
      BackendDAE.EQSYSTEM(orderedVars = BackendDAE.VARIABLES(varArr = varArr)) => begin
        for i in 1:arrayLength(varArr)
          extArg = traversalOperation(varArr[i], extArg)
        end
        (extArg)
      end
    end
  end
  return extArg
end

function traverseEquationExpressions(eq::BackendDAE.Equation, traversalOperation::Function, extArg::T)::Tuple{BackendDAE.Equation,T} where{T}
   (eq, extArg) = begin
     local lhs::DAE.Exp
     local rhs::DAE.Exp
     local cref::DAE.ComponentRef
     @match eq begin
       BackendDAE.EQUATION(lhs = lhs, rhs = rhs) => begin
         (lhs, extArg) = Util.traverseExpTopDown(lhs, traversalOperation, extArg)
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @set eq.lhs = lhs
         @set eq.rhs = rhs
         (eq, extArg)
       end
       BackendDAE.SOLVED_EQUATION(componentRef = cref, exp = rhs) => begin
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @set eq.rhs = rhs;
         (eq, extArg)
       end
       BackendDAE.RESIDUAL_EQUATION(exp = rhs) => begin
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @set eq.rhs = rhs;
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
function DAE_VarKind_to_BDAE_VarKind(kind::DAE.VarKind)::BackendDAE.VarKind
  @match kind begin
    DAE.VARIABLE(__) => BackendDAE.VARIABLE()
    DAE.DISCRETE(__) => BackendDAE.DISCRETE()
    DAE.PARAM(__) => BackendDAE.PARAM()
    DAE.CONST(__) => BackendDAE.CONST()
  end
end

@exportAll()
end
