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


module BackendDump

using MetaModelica

#= ExportAll is not good practice but it makes it so that we do not have to write export after each function :( =#
using ExportAll

import BackendDAE
import BackendDAEUtil

const DOUBLE_LINE = "============================================"::String
const LINE = "---------------------------------------------"::String


function dumpBackendDAEStructure(dae::BackendDAE.BackendDAEStructure, heading::String)
  print(DOUBLE_LINE + "\n")
  print("BackendDAE: " + heading + "\n")
  print(DOUBLE_LINE + "\n")

  for eq in dae.eqs
    print("\nEqs:\n")
    print(LINE + "\n")
    BackendDAEUtil.mapEqSystemEquationsNoUpdate(eq, printEqTraverse, 0)
    print("\nVars:\n")
    print(LINE + "\n")
    BackendDAEUtil.mapEqSystemVariablesNoUpdate(eq, printVarTraverse, 0)
  end
end

function printAnyTraverse(any, extArg)
  print(any)
  print("\n")
  (extArg)
end

function printVarTraverse(var::BackendDAE.Var, extArg)
  print(var.varName.ident)
  print(" | ")
  print(var.varKind)
  print("\n")
  (extArg)
end

# kabdelhak: Very ugly, this needs improvement!
# We need a reasonable printExp() function
function printEqTraverse(eq::BackendDAE.Equation, extArg)
  _ = begin
    local lhs::DAE.Exp
    local rhs::DAE.Exp
    local cref::DAE.ComponentRef
    @match eq begin
      BackendDAE.EQUATION(lhs = lhs, rhs = rhs) => begin
        print(lhs)
        print(" = ")
        print(rhs)
        print("\n")
      end
      BackendDAE.SOLVED_EQUATION(componentRef = cref, exp = rhs) => begin
      print(cref)
      print(" = ")
      print(rhs)
      print("\n")
      end
      BackendDAE.RESIDUAL_EQUATION(exp = rhs) => begin
      print("0 = ")
      print(rhs)
      print("\n")
      end
    end
  end
end

@exportAll()
end
