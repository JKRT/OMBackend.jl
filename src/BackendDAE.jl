#=
# This file is part of OpenModelica.
#
# Copyright (c) 1998-2020, Open Source Modelica Consortium (OSMC),
# c/o Linköpings universitet, Department of Computer and Information Science,
# SE-58183 Linköping, Sweden.
#
# All rights reserved.
#
# THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
# THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
# ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
# RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
# ACCORDING TO RECIPIENTS CHOICE.
#
# The OpenModelica software and the Open Source Modelica
# Consortium (OSMC) Public License (OSMC-PL) are obtained
# from OSMC, either from the above address,
# from the URLs: http:www.ida.liu.se/projects/OpenModelica or
# http:www.openmodelica.org, and in the OpenModelica distribution.
# GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
#
# This program is distributed WITHOUT ANY WARRANTY; without
# even the implied warranty of  MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
# IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
#
# See the full OSMC Public License conditions for more details.
#
=#

using MetaModelica



# A single expression of a equation
# Currently a String like "x+y*sin(2*pi*time)*func(42)"
Expression = String


# Equation uniontyp representing a single equation
@Uniontype Equation begin
  @Record EQUATION begin
    # General equation of form lhs = rhs
    lhs::Expression     # Left hand side of expression
    rhs::Expression     # Right hand side of expression
  end

  @Record RESIDUAL_EQUATION begin
    # Residual equations 0 = expression
    expression::Expression
  end
end

# Define a system of equations
@Uniontype EqSystem begin
  @Record EQSYSTEM begin
    equations::Array{Equation,1}
  end
end

# Main data strucutre containing all equation systems, variables and shared information
@Uniontype BackendDAE begin
  @Record DAE begin
    eqSyst::EqSystem
    #DAEVariables variables;
    # TODO: Add "shared" stuff for everything else
  end
end
