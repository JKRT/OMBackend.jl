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

module BackendEquation

using MetaModelica

#= ExportAll is not good practice but it makes it so that we do not have to write export after each function :( =#
using ExportAll

import DAE
import BackendDAE

"""
    kabdelhak:
    Create an empty equation array
"""
function emptyEqns()
  eqns::BackendDAE.EquationArray = []
  (eqns)
end

"""
    kabdelhak:
    Transform a single equation to residual form by subtracting the rhs from
    the lhs
"""
function makeResidualEquation(eqn::BackendDAE.Equation)::BackendDAE.Equation
  eqn = begin
    local lhs::DAE.Exp
    local rhs::DAE.Exp
    local source::DAE.ElementSource
    local attr::EquationAttributes
    @match eqn begin
      BackendDAE.EQUATION(lhs, rhs, source, attr) => begin
        BackendDAE.RESIDUAL_EQUATION(makeResidualExp(lhs, rhs), source, attr)
      end

      BackendDAE.IF_EQUATION(__) => begin
        (eqn) #makeResidualIfEquation(eqn)
      end

      _ => begin
        (eqn)
      end
    end
  end
end

function makeResidualIfEquation(eqn::BackendDAE.Equation)::BackendDAE.Equation
    local trueEquations::List{List{BackendDAE.Equation}} = nil
    local tmpTrue::List{DAE.Equation}
    local falseEquations::List{DAE.Equation}
    (eqn)
end

"""
   kabdelhak:
   Splits an if equation in multiple equations containing if expressions.
   INTENDED?: Currently requires correct ordering of branch equations
"""
function splitIfEquationResidual(eqn::BackendDAE.Equation)::List{BackendDAE.Equation}
  residualEqs::List{BackendDAE.Equation}=nil
  residualExps = splitIfEquationResidualTraverse(eqn.conditions, Util.transposeNestedList(eqn.eqnstrue), eqn.eqnsfalse, nil)
  for exp in residualExps
    residualEqs = BackendDAE.RESIDUAL_EQUATION(exp, eqn.source, eqn.attr) <| residualEqs
  end
  return residualEqs
end

function splitIfEquationResidualTraverse(lstCond::List{DAE.Exp}, lstlstTrue::List{List{BackendDAE.Equation}}, lstFalse::List{BackendDAE.Equation}, acc::List{DAE.Exp})::List{DAE.Exp}
  acc = begin
      local tmpTrue::List{BackendDAE.Equation}
      local restTrue::List{List{BackendDAE.Equation}}
      local tmpFalse::BackendDAE.Equation
      local restFalse::List{BackendDAE.Equation}
      @match (lstlstTrue, lstFalse) begin
      (nil, _) => begin
        (acc)
      end
      (tmpTrue <| restTrue, tmpFalse <| restFalse) => begin
         acc  = makeNestedIfExpressionResidual(lstCond, tmpTrue, tmpFalse) <| acc
         splitIfEquationResidualTraverse(lstCond, restTrue, restFalse, acc)
      end
    end
  end
  return acc
end

function makeNestedIfExpressionResidual(lstCond::List{DAE.Exp}, lstTrue::List{BackendDAE.Equation}, eqFalse::BackendDAE.Equation)
  exp = begin
    local cond::DAE.Exp
    local restCond::List{DAE.Exp}
    local eqTrue::DAE.Equation
    local restTrue::List{DAE.Equation}
    @match (lstCond, lstTrue) begin
      (nil, _) => begin
         #= fails for anything but BackendDAE.EQUATION() =#
         makeResidualExp(eqFalse.lhs, eqFalse,rhs)
      end
      (cond <| restCond, eqTrue <| restTrue) => begin
         #= fails for anything but BackendDAE.EQUATION() =#
         DAE.IFEXP(cond, makeResidualExp(eqTrue.lhs, eqTrue.rhs), makeNestedIfExpressionResidual(restCond, restTrue, eqFalse))
      end
    end
  end
end

function makeResidualExp(lhs::DAE.Exp, rhs::DAE.Exp)::DAE.Exp
  DAE.BINARY(lhs, DAE.SUB(DAE.T_REAL_DEFAULT), rhs)
end

"""
    kabdelhak:
    create a single array with all equations
"""
function concenateEquations(eqs::Array{BackendDAE.EqSystem})::Array{BackendDAE.Equation}
  eqArr::Array{BackendDAE.Equation} = []
  for eq in eqs
    eqArr = vcat(eqArr, eq.orderedEqs)
  end
  return eqArr
end


@exportAll()
end
