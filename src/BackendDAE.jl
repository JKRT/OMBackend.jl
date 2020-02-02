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

module BackendDAE

using MetaModelica

#= ExportAll is not good practice but it makes it so that we do not have to write export after each function :( =#
using ExportAll

import DAE

@Uniontype BackendDAE begin
  @Record BACKENDAE
    eqs #= The list of equation systems =# :: List{EqSystem}
  end
end

@Uniontype EqSystem begin
  @Record EQSYSTEM begin
    orderedVars #= Ordered Variables, only states and alg. vars =# :: Array{Var,1}
    orderedEqs #= Ordered Equations =# :: Array{Equation,1}
  end
end

@Uniontype Var begin
  @Record VAR begin
    cref :: DAE.ComponentRef
    kind :: DAE.VarKind
  end
end

@Uniontype Equation begin
  @Record EQUATION begin
    lhs :: DAE.Exp
    rhs :: DAE.Exp
  end

  @Record SOLVED_EQUATION begin
    cref :: DAE.ComponentRef
    rhs :: DAE.Exp
  end

  @Record RESIDUAL_EQUATION begin
    exp :: DAE.Exp
  end
end

@exportAll()
end
