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
import BackendDump

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

"""
    kabdelhak:
    Creates a residual equation from an if equation.
"""
function makeResidualIfEquation(eqn::BackendDAE.Equation)::BackendDAE.Equation
    local trueEquations::List{List{BackendDAE.Equation}} = nil
    local tmpTrue::List{DAE.Equation}
    local falseEquations::List{DAE.Equation}
    # put splitting of if equations here, needs to return an array equation or
    # multiple equations
    (eqn)
end

"""
   kabdelhak:
   Splits an if equation in multiple equations containing if expressions.
   INTENDED?: Currently does not require correct ordering of branch equations
     because of residual form.
   NOTE: Not used yet
"""
function splitIfEquationResidual(eqn::BackendDAE.Equation)::List{BackendDAE.Equation}
  residualEqs::List{BackendDAE.Equation}=nil
  residualExps = splitIfEquationResidualTraverse(eqn.conditions, Util.transposeNestedList(eqn.eqnstrue), eqn.eqnsfalse, nil)
  for exp in residualExps
    residualEqs = BackendDAE.RESIDUAL_EQUATION(exp, eqn.source, eqn.attr) <| residualEqs
  end
  return residualEqs
end


"""
   kabdelhak:
   Helper function to splitIfEquationResidual. Traverses all expressions to create the array
   of if expressions.
   e.g.
     if COND_1 then
       x = TRUE_1_1;
       y = TRUE_2_1;
       ...
       n = TRUE_N_1;
     elseif COND_2 then
       x = TRUE_1_2;
       y = TRUE_2_2;
       ...
       n = TRUE_N_2;
     else
       x = ELSE_1;
       y = ELSE_2;
       ...
       n = ELSE_N;
     end if;

   becomes:

    0 = if COND_1 then TRUE_1_1 - x elseif COND_2 then TRUE_1_2 - x else ELSE_1 - x;
    0 = if COND_1 then TRUE_2_1 - y elseif COND_2 then TRUE_2_2 - y else ELSE_2 - y;
       ...
    0 = if COND_1 then TRUE_N_1 - n elseif COND_2 then TRUE_N_2 - n else ELSE_N - n;
"""
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

"""
    kabdelhak:
    Creates nested if expressions to represent elseif.
"""
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

"""
    kabdelhak:
    Creates a residual exp from two expressions lhs and rhs.
    => lhs - rhs
"""
function makeResidualExp(lhs::DAE.Exp, rhs::DAE.Exp)::DAE.Exp
  DAE.BINARY(lhs, DAE.SUB(DAE.T_REAL_DEFAULT), rhs)
end

"""
    kabdelhak:
    Creates a single array with all equations from all equation systems.
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
