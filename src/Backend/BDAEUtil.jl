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

import Absyn
using MetaModelica
using ExportAll
using Setfield


import ..DAE
import ..Util

import ..BDAE
import ..BackendEquation

"""
This function converts an array of variables to the BDAE variable structure
"""
function convertVarArrayToBDAE_Variables(vars::Array{BDAE.Var})::BDAE.Variables
  local variables::BDAE.Variables = begin
    BDAE.VARIABLES([i for i in vars])
  end
  return variables
end

function createEqSystem(vars::BDAE.Variables, eqs::Array)
  (BDAE.EQSYSTEM(vars,
                 eqs,
                 NONE(),
                 NONE(),
                 NONE(),
                 BDAE.NO_MATCHING(),
                 nil,
                 BDAE.UNKNOWN_PARTITION(),
                 BackendEquation.emptyEqns()))
end

function mapEqSystems(dae::BDAE.BDAEStructure, traversalOperation::Function)
   dae = begin
     local eqs::Array{BDAE.EqSystem, 1}
     @match dae begin
       BDAE.BACKEND_DAE(eqs = eqs) => begin
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

function mapEqSystemEquations(syst::BDAE.EqSystem, traversalOperation::Function)
  syst = begin
    local eqs::Array{BDAE.Equation,1}
    @match syst begin
      BDAE.EQSYSTEM(orderedEqs = eqs) => begin
        for i in 1:arrayLength(eqs)
          eqs[i] = traversalOperation(eqs[i])
        end
        @set syst.orderedEqs = eqs
        (syst)
      end
    end
  end
end


function mapEqSystemEquationsNoUpdate(syst::BDAE.EqSystem, traversalOperation::Function, extArg)
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

function mapEqSystemVariablesNoUpdate(syst::BDAE.EqSystem, traversalOperation::Function, extArg)
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
       BDAE.EQUATION(lhs = lhs, rhs = rhs) => begin
         (lhs, extArg) = Util.traverseExpTopDown(lhs, traversalOperation, extArg)
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @set eq.lhs = lhs
         @set eq.rhs = rhs
         (eq, extArg)
       end
       BDAE.SOLVED_EQUATION(componentRef = cref, exp = rhs) => begin
         (rhs, extArg) = Util.traverseExpTopDown(rhs, traversalOperation, extArg)
         @set eq.rhs = rhs;
         (eq, extArg)
       end
       BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
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

include("backendDump.jl")
@exportAll()
end
