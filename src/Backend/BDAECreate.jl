#=
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
module BDAECreate

using MetaModelica
using ExportAll

import DAE
import ..BDAE
import ..BDAEUtil
import ..Causalize
import OMFrontend

"""
  This function translates a DAE, which is the result from instantiating a
  class, into a more precise form, called BDAE.BDAE defined in this module.
  The BDAE.BDAE representation splits the DAE into equations and variables
  and further divides variables into known and unknown variables and the
  equations into simple and nonsimple equations.
  inputs:  lst: DAE.DAE_LIST
  outputs: BDAE.BACKEND_DAE"""
function lower(lst::DAE.DAE_LIST)::BDAE.BACKEND_DAE
  local outBDAE::BDAE.BACKEND_DAE
  local eqSystems::Vector{BDAE.EQSYSTEM}
  local varArray::Vector{BDAE.VAR}
  local eqArray::Vector{BDAE.Equation}
  local name = listHead(lst.elementLst).ident
  (varArray, eqArray, initialEquations) = begin
    local elementLst::List{DAE.Element}
    local variableLst::List{BDAE.VAR}
    local equationLst::List{BDAE.Equation}
    @match lst begin
      DAE.DAE_LIST(elementLst) => begin
        (variableLst, equationLst, initialEquations) = splitEquationsAndVars(elementLst)
        (listArray(listReverse(variableLst)), listArray(listReverse(equationLst)), initialEquations)
      end
    end
  end
  local variables = BDAEUtil.convertVarArrayToBDAE_Variables(varArray)
  @debug "varArray:" length(variableLst)
  @debug "eqLst:" length(equationLst)
  #= We start with an array of one system =#
  eqSystems = [BDAEUtil.createEqSystem(variables, eqArray)]
  outBDAE = BDAE.BACKEND_DAE(name, eqSystems, BDAE.SHARED([], [], listArray(initialEquations)))
end

"""
  Lowers a FlatModelica defined in the new frontend into BDAE.
  1. We translate all different components of the flat model into the DAE representation.
  2. We convert this representation into the BackendDAE representation.
  3. We return backend DAE to be used in the remainder of the compilation before code generation.
"""
function lower(frontendDAE::OMFrontend.Main.FlatModel)
  local equations = [equationToBackendEquation(eq) for eq in OMFrontend.Main.convertEquations(frontendDAE.equations)]
  local variables = [variableToBackendVariable(var) for var in OMFrontend.Main.convertVariables(frontendDAE.variables, list())] 
  local algorithms = [alg for alg in frontendDAE.algorithms]
  local iAlgorithms = [iAlg for iAlg in frontendDAE.initialAlgorithms]
  local initialEquations = [equationToBackendEquation(ieq) for ieq in OMFrontend.Main.convertEquations(frontendDAE.initialEquations)]  
  eqSystems = [BDAEUtil.createEqSystem(variables, equations)]
  #= 
  The resulting backend DAE.   
  =#
  outBDAE = BDAE.BACKEND_DAE(frontendDAE.name, eqSystems, BDAE.SHARED([], [], initialEquations))
end

function convertVariableIntoBDAEVariable(var::OMFrontend.Main.Variable)
  elem = OMFrontend.Main.convertVariable(var, OMFrontend.Main.VARIABLE_CONVERSION_SETTINGS(true, false, true))
  BDAE.VAR(elem.componentRef,
           BDAEUtil.DAE_VarKind_to_BDAE_VarKind(elem.kind),
           elem.direction,
           elem.ty,
           elem.binding,
           elem.dims,
           elem.source,
           elem.variableAttributesOption,
           NONE(), #=Tearing=#
           elem.connectorType,
           false #=We do not know if we can replace or not yet=#
           )
end


                    
"""
  Splits a given DAE.DAEList and converts it into a set of BDAE equations and BDAE variables.
  In addition provides the initial equations for the system.
  TODO: Optimize by using List instead of array.
"""
function splitEquationsAndVars(elementLst::List{DAE.Element})::Tuple{List, List, List}
  local variableLst::List{BDAE.VAR} = nil
  local equationLst::List{BDAE.Equation} = nil
  local initialEquationLst::List{BDAE.Equation} = nil
  for elem in elementLst
    _ = begin
      local backendDAE_Var
      local backendDAE_Equation
      @match elem begin
        DAE.VAR(__) => begin
          variableLst = BDAE.VAR(elem.componentRef,
          BDAEUtil.DAE_VarKind_to_BDAE_VarKind(elem.kind),
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
          equationLst = BDAE.EQUATION(elem.exp,
                                      elem.scalar,
                                      elem.source,
                                      #=TODO: Below might need to be changed =#
                                      BDAE.EQ_ATTR_DEFAULT_UNKNOWN) <| equationLst
        end
        DAE.WHEN_EQUATION(__) => begin
          equationLst = lowerWhenEquation(elem) <| equationLst
        end
        DAE.IF_EQUATION(__) => begin
          equationLst = lowerIfEquation(elem) <| equationLst
        end
        DAE.INITIALEQUATION(__) => begin
          initialEquationLst = BDAE.EQUATION(elem.exp1,
          elem.exp2,
          elem.source,
          #=TODO: Below might need to be changed =#
          BDAE.EQ_ATTR_DEFAULT_UNKNOWN) <| initialEquationLst
        end
        DAE.COMP(__) => begin
          variableLst,equationLst,initialEquationLst = splitEquationsAndVars(elem.dAElist)
        end
        _ => begin
          @error "Skipped:" elem
          throw("Unsupported equation: $elem")
        end
      end
    end
  end
  return (variableLst, equationLst, initialEquationLst)
end

function equationToBackendEquation(elem::DAE.Element)
  @match elem begin
    DAE.EQUATION(__) => begin
      equationLst = BDAE.EQUATION(elem.exp,
                                  elem.scalar,
                                  elem.source,
                                  #=TODO: Below might need to be changed =#
                                  BDAE.EQ_ATTR_DEFAULT_UNKNOWN)
    end
    DAE.WHEN_EQUATION(__) => begin
      equationLst = lowerWhenEquation(elem)
    end
    DAE.IF_EQUATION(__) => begin
      equationLst = lowerIfEquation(elem)
    end
    DAE.INITIALEQUATION(__) => begin
      initialEquationLst = BDAE.EQUATION(elem.exp1,
                                         elem.exp2,
                                         elem.source,
                                         #=TODO: Below might need to be changed =#
                                         BDAE.EQ_ATTR_DEFAULT_UNKNOWN)
    end
    DAE.COMP(__) => begin
      variableLst,equationLst,initialEquationLst = splitEquationsAndVars(elem.dAElist)
    end
    _ => begin
      @error "Skipped:" elem
      throw("Unsupported equation: $elem")
    end
  end
end

function variableToBackendVariable(elem::DAE.Element)
  @match elem begin
    DAE.VAR(__) => begin
      variableLst = BDAE.VAR(elem.componentRef,
      BDAEUtil.DAE_VarKind_to_BDAE_VarKind(elem.kind),
      elem.direction,
      elem.ty,
      elem.binding,
      elem.dims,
      elem.source,
      elem.variableAttributesOption,
      NONE(), #=Tearing=#
      elem.connectorType,
      false #=We do not know if we can replace or not yet=#)
    end
  end
end


function lowerWhenEquation(eq::DAE.Element)::BDAE.WHEN_EQUATION
  local whenOperatorLst::List{BDAE.WhenOperator} = nil
  local whenEquation::BDAE.WhenEquation
  local elseOption::Option{BDAE.WhenEquation} = NONE()
  local elseEq::BDAE.Element
  whenOperatorLst = createWhenOperators(eq.equations, whenOperatorLst)
  if isSome(eq.elsewhen_)
    SOME(elseEq) = eq.elsewhen_
    elseOption = SOME(lowerWhenEquation(elseEq))
  end
  whenEquation = BDAE.WHEN_STMTS(eq.condition, whenOperatorLst, elseOption)
  return BDAE.WHEN_EQUATION(1, whenEquation, eq.source, BDAE.EQ_ATTR_DEFAULT_UNKNOWN)
end

function createWhenOperators(elementLst::List{DAE.Element},lst::List{BDAE.WhenOperator})::List{BDAE.WhenOperator}
  lst = begin
    local rest::List{DAE.Element}
    local acc::List{BDAE.WhenOperator}
    local cref::DAE.ComponentRef
    local e1::DAE.Exp
    local e2::DAE.Exp
    local e3::DAE.Exp
    local source::DAE.ElementSource
    @match elementLst begin
      DAE.EQUATION(exp = e1, scalar = e2, source = source) <| rest => begin
        acc = BDAE.ASSIGN(e1, e2, source) <| lst
        createWhenOperators(rest, acc)
      end
      DAE.ASSERT(condition = e1, message = e2, level = e3, source = source) <| rest => begin
        acc = BDAE.ASSERT(e1, e2, e3, source) <| lst
        createWhenOperators(rest, acc)
      end
      DAE.TERMINATE(message = e1, source = source) <| rest => begin
        acc = BDAE.TERMINATE(e1, source) <| lst
        createWhenOperators(rest, acc)
      end
      DAE.REINIT(componentRef = cref, exp = e1, source = source) <| rest => begin
        #= BDAE uses an exp here instead of a cref =#
        expTy = if typeof(cref.identType) == DAE.T_ARRAY
          #= If we are refeering to an array it is the content of the array that is the type of the exp. =#
          cref.identType.ty #= Note this would be wrong if we would consider other compound types. =#
        else
          cref.identType #=OK it is the type of the component reference directly=#
        end
        local crefExp = DAE.CREF(cref, expTy)
        acc = BDAE.REINIT(crefExp, e1, source) <| lst
        createWhenOperators(rest, acc)
      end
      DAE.NORETCALL(exp = e1, source = source) <| rest => begin
        acc = BDAE.NORETCALL(e1, source) <| lst
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

"""
  Transform a DAE if-equation into a BDAE if-equation
"""
function lowerIfEquation(eq::IF_EQ)::BDAE.IF_EQUATION where {IF_EQ}
  local trueEquations::List{List{BDAE.Equation}} = nil
  local tmpTrue::List{BDAE.Equation}
  local falseEquations::List
  for lst in eq.equations2
    (_, tmpTrue, _) = splitEquationsAndVars(lst)
    trueEquations = tmpTrue <| trueEquations
  end
  trueEquations = listReverse(trueEquations)
  (_, falseEquations, _) = splitEquationsAndVars(eq.equations3)
  res = BDAE.IF_EQUATION(eq.condition1,
                         trueEquations,
                         falseEquations,
                         eq.source,
                         BDAE.EQ_ATTR_DEFAULT_UNKNOWN)
  return res
end

@exportAll()
end
