
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

"""
  This function translates a DAE, which is the result from instantiating a
  class, into a more precise form, called BackendDAE.BackendDAE defined in this module.
  The BackendDAE.BackendDAE representation splits the DAE into equations and variables
  and further divides variables into known and unknown variables and the
  equations into simple and nonsimple equations.
  The variables are inserted into a dictonary, Backend_DAE_Map. 
  The equations are put in an expandable array. 
  Where adding a new equation can be done in O(1) time if space is available.
  inputs:  lst: DAE.DAElist, inCache: FCore.Cache, inEnv: FCore.Graph
  outputs: BackendDAE.BackendDAE"""
function lower(lst::DAE.DAElist)::List{BackendDAE.BackendDAEStructure}
  local outBackendDAE::BackendDAE.BackendDAEStructure
  local eqSystems::Array{BackendDAE.EqSystem}
  local varArray::Array{BackendDAE.Var, 1}
  local eqArray::Array{BackendDAE.Equation, 1}
  (varArray, eqArray) = begin
    local elementLst::List{DAE.Element}
    local variableLst::List{BackendDAE.Var} #= init empty ? =#
    local equationLst::List{BackendDAE.Equation}
    @match lst begin
      DAE.DAE_LIST(elementLst) => begin
        (variableLst, equationLst) = splitEquationsAndVars(elementLst)
        (listArray(variableLst), listArray(equationLst))
      end
    end
  end
  eqSystems = BackendDAEUtil.createEqSystem(variableArray, equationArray)
  outBackendDAE = BackendDAE.DAE(eqSystems, SHARED_DUMMY())
end

"""
  Splits a given DAE.DAEList into equations and variables
"""
function splitEquationsAndVars(elementLst::List{DAE.Element})::Tuple(List, List)
  local variableLst::List{BackendDAE.Var} = nil
  local equationLst::List{BackendDAE.Equation} = nil
  for elem in elementLst
    _ = begin
      local cref::DAE.ComponentRef
      local kind::DAE.VarKind
      local lhs::DAE.Exp
      local rhs::DAE.Exp
      @match elem begin
        DAE.VAR(componentRef = cref, kind = kind) => begin
          variableLst = BackendDAE.VAR(cref = cref, kind = kind) <| variableLst
        end
        DAE.EQUATION(exp = lhs, scalar = rhs) => begin
          equationLst = BackendDAE.EQUATION(lhs = lhs, rhs = rhs) <| equationLst
        end
      end
    end
  end
end

@exportAll()
end
