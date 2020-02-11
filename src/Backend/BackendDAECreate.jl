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
"""
This module contain the various functions that are related to the lowering
of the DAE IR into Backend DAE IR (BDAE IR). BDAE IR is the representation we use
before code generation.
"""
module BackendDAECreate

using MetaModelica
using ExportAll

import DAE
import BackendDAE
import BackendDAEUtil

"""
  This function translates a DAE, which is the result from instantiating a
  class, into a more precise form, called BackendDAE.BackendDAE defined in this module.
  The BackendDAE.BackendDAE representation splits the DAE into equations and variables
  and further divides variables into known and unknown variables and the
  equations into simple and nonsimple equations.
  The variables are inserted into a dictonary, Backend_DAE_Map.
  The equations are put in an expandable array.
  Where adding a new equation can be done in O(1) time if space is available.
  inputs:  lst: DAE.DAE_LIST, inCache: FCore.Cache, inEnv: FCore.Graph
  outputs: BackendDAE.BackendDAE"""
function lower(lst::DAE.DAE_LIST)::BackendDAE.BackendDAEStructure
  local outBackendDAE::BackendDAE.BackendDAEStructure
  local eqSystems::Array{BackendDAE.EqSystem}
  local varArray::Array{BackendDAE.Var}
  local eqArray::Array{BackendDAE.Equation}
  local name = listHead(lst.elementLst).ident
  (varArray, eqArray) = begin
    local elementLst::List{DAE.Element}
    local variableLst::List{BackendDAE.Var}
    local equationLst::List{BackendDAE.Equation}
    @match lst begin
      DAE.DAE_LIST(elementLst) => begin
        (variableLst, equationLst) = splitEquationsAndVars(elementLst)
        (listArray(variableLst), listArray(equationLst))
      end
    end
  end
  local variables = BackendDAEUtil.convertVarArrayToBackendDAE_Variables(varArray)
  #= We start with an array of one system =#
  eqSystems = [BackendDAEUtil.createEqSystem(variables, eqArray)]
  outBackendDAE = BackendDAE.BACKEND_DAE(name, eqSystems, BackendDAE.SHARED_DUMMY())
end

"""
  Splits a given DAE.DAEList into equations and variables
"""
function splitEquationsAndVars(elementLst::List{DAE.Element})::Tuple
  local variableLst::List{BackendDAE.Var} = nil
  local equationLst::List{BackendDAE.Equation} = nil
  for elem in elementLst
    _ = begin
      local backendDAE_Var
      local backendDAE_Equation
      @match elem begin
        DAE.VAR(__) => begin
          variableLst = BackendDAE.VAR(elem.componentRef,
                                       BackendDAEUtil.DAE_VarKind_to_BDAE_VarKind(elem.kind),
                                       elem.direction,
                                       elem.ty,
                                       elem.binding,
                                       elem.dims,
                                       elem.source,
                                       elem.variableAttributesOption,
                                       NONE(), #=Tearing=#
                                       elem.connectorType,
                                       false #=We do not know if we can replace or not yet=#
                                       ) <| variableLst
        end
        DAE.EQUATION(__) => begin
          equationLst = BackendDAE.EQUATION(elem.exp,
                                            elem.scalar,
                                            elem.source,
                                            #=TODO: Below might need to be changed =#
                                            BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN) <| equationLst
        end

        DAE.WHEN_EQUATION(__) => begin
          equationLst = lowerWhenEquation(elem) <| equationLst
        end

        DAE.IF_EQUATION(__) => begin
          equationLst = lowerIfEquation(elem) <| equationLst
        end

        DAE.COMP(__) => begin
          variableLst,equationLst = splitEquationsAndVars(elem.dAElist)
        end
        _ => begin
          continue
        end
      end
    end
  end
  return (variableLst, equationLst)
end

function lowerWhenEquation(eq::DAE.Element)::BackendDAE.Equation
  local whenOperatorLst::List{BackendDAE.WhenOperator} = nil
  local whenEquation::BackendDAE.WhenEquation
  local elseOption::Option{BackendDAE.WhenEquation} = NONE()
  local elseEq::BackendDAE.Element
  whenOperatorLst = createWhenOperators(eq.equations, whenOperatorLst)
  if isSome(eq.elsewhen_)
    SOME(elseEq) = eq.elsewhen_
    elseOption = SOME(lowerWhenEquation(elseEq))
  end
  whenEquation = BackendDAE.WHEN_STMTS(eq.condition, whenOperatorLst, elseOption)
  return BackendDAE.WHEN_EQUATION(1, whenEquation, eq.source, BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN)
end

function createWhenOperators(elementLst::List{DAE.Element},lst::List{BackendDAE.WhenOperator})::List{BackendDAE.WhenOperator}
  lst = begin
    local rest::List{DAE.Element}
    local acc::List{BackendDAE.WhenOperator}
    local cref::DAE.ComponentRef
    local e1::DAE.Exp
    local e2::DAE.Exp
    local e3::DAE.Exp
    local source::DAE.ElementSource
    @match elementLst begin
      DAE.EQUATION(exp = e1, scalar = e2, source = source) <| rest => begin
        acc = BackendDAE.ASSIGN(e1, e2, source) <| lst
        createWhenOperators(rest, acc)
      end
      DAE.ASSERT(condition = e1, message = e2, level = e3, source = source) <| rest => begin
        acc = BackendDAE.ASSERT(e1, e2, e3, source) <| lst
        createWhenOperators(rest, acc)
      end
      DAE.TERMINATE(message = e1, source = source) <| rest => begin
      acc = BackendDAE.TERMINATE(e1, source) <| lst
      createWhenOperators(rest, acc)
      end
      DAE.REINIT(componentRef = cref, exp = e1, source = source) <| rest => begin
        acc = BackendDAE.REINIT(cref, e1, source) <| lst
        createWhenOperators(rest, acc)
      end
      DAE.NORETCALL(exp = e1, source = source) <| rest => begin
        acc = BackendDAE.NORETCALL(e1, source) <| lst
        createWhenOperators(rest, acc)
      end
      #= MAYBE MORE CASES NEEDED =#
      nil => begin
        (lst)
      end
      _ <| rest => begin
        createWhenOperators(rest, lst)
      end
    end
  end
end

function lowerIfEquation(eq::DAE.Element)::Backend.Equation
  local trueEquations::List{List{BackendDAE.Equation}} = nil
  local tmpTrue::List{DAE.Equation}
  local falseEquations::List{DAE.Equation}

  for lst in eq.equations2
    (_, tmpTrue) = splitEquationsAndVars(lst)
    trueEquations = tmpTrue <| trueEquations
  end

  trueEquations = listReverse(trueEquations)

  (_, falseEquations) = splitEquationsAndVars(eq.equations3)

  return BackendDAE.IF_EQUATION(eq.condition1,
                                trueEquations,
                                falseEquations,
                                eq.source,
                                BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN)
end

@exportAll()
end
