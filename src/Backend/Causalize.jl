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

#= ExportAll is not good practice but it makes it so that we do not have to write export after each function :( =#
using ExportAll

import DAE
import ComponentReference
import BackendDAE
import BackendDAEUtil

function detectStates(dae::BackendDAE.BackendDAE)
  dae = BackendDAEUtil.mapEqSystems(dae, detectStatesEqSystem)
end

function detectStatesEqSystem(syst::BackendDAE.EqSystem)
  syst = begin
    @match syst begin
      local vars::Array{BackendDAE.Var,1}
      local eqs::Array{BackendDAE.Equation,1}
      local stateCrefs::List{DAE.ComponentRef} = nil
      BackendDAE.EQSYSTEM(orderedVars = vars, orderedEqs = eqs) => begin #= qualified access possible? =#
        for eq in eqs
          stateCrefs = BackendDAEUtil.traveseEquationExpressions(eq, detectStateExpression, stateCrefs)
        end
        #= Do replacements for stateCrefs =#
        (syst)
      end
    end
  end
end

function detectStateExpression(exp::DAE.Exp, stateCrefs::List{DAE.ComponentRef})
  stateCrefs = begin
    local state::DAE.ComponentRef
    @match exp begin
      DAE.CALL(path = Absyn.IDENT("der"), list(DAE.CREF(componentRef = state))) => begin
        (state <| stateCrefs)
      end
    end
  end
end

function updateStates(vars::Array{BackendDAE.Var,1}, stateCrefs::List{DAE.ComponentRef})
  vars = begin
    local state::DAE.ComponentRef
    local rest::List{DAE.ComponentRef}
    @match stateCrefs begin
      state <| rest => begin
        for i in arrayLength(vars)
          vars[i] = begin
            local cref::DAE.ComponentRef
            @match vars[i] begin
              BackendDAE.VAR(cref = cref) where (ComponentReference.crefEqual(cref, state)) => begin
                !set vars[i].varKind = BackendDAE.STATE(0, NONE(), true)
              end
            end
          end
        end
      end
    end
    nil => begin
      (vars)
    end
  end
end

function daeMode(dae::BackendDAE.BackendDAE)
  dae = BackendDAEUtil.mapEqSystems(dae, makeResidualEquations)
end

function makeResidualEquations(syst::BackendDAE.EqSystem)
  syst = BackendDAEUtil.mapEqSystemEquations(syst, BackendEquation.makeResidualEquation)
end


@exportAll()
end
