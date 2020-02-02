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

#= ExportAll is not good practice but it makes it so that we do not have to write export after each function :( =#
using ExportAll

import DAE
import Expression
import BackendDAE

mapFunc = Function

function mapEqSystems(dae::BackendDAE.BackendDAE, mapFunc::mapFunc)
  dae = begin
    local eqs::List{BackendDAE.EqSystem}
    local acc::List{BackendDAE.EqSystem} = nil
    @match dae
      BackendDAE.BACKENDDAE(eqs = eqs) => begin #= qualified access possible? =#
        for syst in eqs
          acc = mapFunc(syst) <| acc
        end
        BackendDAE.BACKENDDAE(eqs = listReverse(eqs))
      end
    end
  end
end

function traveseEquationExpressions(eq::BackendDAE.Equation, mapFunc::mapFunc, extArg::T) where{T}
  (eq, extArg) = begin
    @match eq begin
      local lhs::DAE.Exp
      local rhs::DAE.Exp
      local cref::DAE.ComponentRef
      BackendDAE.EQUATION(lhs = lhs, rhs = rhs) => begin
        (lhs, extArg) = Expression.traverseExpTopDown(lhs, mapFunc, extArg)
        (rhs, extArg) = Expression.traverseExpTopDown(rhs, mapFunc, extArg)
        (BackendDAE.EQUATION(lhs = lhs, rhs = rhs), extArg)
      end
      BackendDAE.SOLVED_EQUATION(cref = cref, rhs = rhs) => begin
        (rhs, extArg) = Expression.traverseExpTopDown(rhs, mapFunc, extArg)
        (BackendDAE.SOLVED_EQUATION(cref = cref, rhs = rhs), extArg)
      end
      BackendDAE.RESIDUAL_EQUATION(exp = rhs) => begin
        (rhs, extArg) = Expression.traverseExpTopDown(rhs, mapFunc, extArg)
        (BackendDAE.SOLVED_EQUATION(exp = rhs), extArg)
      end
    end
  end
end

@exportAll()
end
