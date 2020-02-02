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

include("CodeGenerationUtil.jl")
#include("BackendDAE.jl")



function generateSingleResidualEquation(equation#=::Equation=#)


end


"""
Write function modelName_DAE_equations() to file.
"""
function generateDAEFunction(file::IOStream, backendDAE, modelName::String)

  write(file, "function $(modelName)_DAE_equations(res, dx, x, p, t)\n")

  equationSystems = backendDAE.eqs::EqSystem
  index = 1
  for eqSystem in equationSystems   # Loop over equation Systems
    for equation in eqSystem.orderedEqs # Loop over equations
      @match(equation)
      @case(eq as RESIDUAL_EQUATION(__)) then
        write(file, string("  res[$index] = ", exp, "\n"))  # TODO: Add DAE.Exp to string() function
        index += 1
    end
  end

  write(file, "end\n")
end

```

Generate a julia file containing functions to simulate the DAE
```
function generateCode(backendDAE#=::BackendDAE=#, path::Strings, modelName::String)

  # Create file
  filename = string(modelName, ".jl")
  file = open(joinpath(path,filename), "w")
  write(file, copyRightString())

  # Write DAEFunction
  generateDAEFunction(file, backendDAE, modelName)

  # Finished code generation
  close(file)
  println("Generated file succesfull")
end
