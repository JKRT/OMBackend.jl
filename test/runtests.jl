
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
*/ =#

#= Author John Tinnerholm & Andreas Heureman =#

import CSV
import Tables

# Load module for test DAEs

include("debugUtil.jl")
#= Uncomment to turn on logging =#
#ENV["JULIA_DEBUG"] = "OMBackend"

using OMBackend
using OMBackend.ExampleDAEs
using Test

import OMBackend.Runtime.OMSolution

"""
  Test translation to backend DAE.
  Also test simulation from 0.0 to 1.0 seconds. 
  Does not verify semantics, only checks that the model was simulated correctly.
"""
function testBackend(TEST_CASES::Array; mode)
  local MODEL_NAME = ""
  #= These tests does not validate the semantics. Only that we can successfully compile and simulate this set of models=#
  @testset "Sanity test. Simple translation and simulation of non-hybrid-DAE" begin
    @testset "Simple non-hybrid continuous systems" begin
      for testCase in TEST_CASES
        @testset "$testCase" begin
          @testset "Compile" begin
            frontendDAE::OMBackend.DAE.DAE_LIST = getfield(ExampleDAEs,Symbol("$(testCase)_DAE"))
            try
              #= TODO: We should check this with some reference IR =#
              (MODEL_NAME, modelCode) = OMBackend.translate(frontendDAE; BackendMode = mode)
              @debug generateFile(testCase, modelCode)
              @info "Translated: $MODEL_NAME"
              @test true
            catch e
              @info e
              throw(e)
              @test false
            end
          end
          @testset "Simulate" begin
            @info "Simulating: $MODEL_NAME"
            try
              simulationResult = getSimulationResult(OMBackend.simulateModel(MODEL_NAME; tspan=(0.0, 1.0)))
              @test simulationResult == :Success
            catch e
              @info "Simulation failure: see $(testCase)_result.csv"
              @test false
            end
          end
        end
      end
    end
  end
end

"Get the solution code for an OMSolution"
function getSimulationResult(solution::OMBackend.CodeGeneration.OMSolution)
  return solution.diffEqSol.retcode
end

"""
  Get solution from a \"regular\" solution from DifferentialEquation.jl suite
"""
function getSimulationResult(regularSol)
  return regularSol.retcode
end

"""
  Main driver to run the test cases.
"""
function runTests()
  #= DAE's to use for the sanity check. =#
  TEST_CASES_BASIC = ["helloWorld", "lotkaVolterra", "vanDerPol"]
  TEST_CASES_ADVANCED = ["simpleMech", "simpleCircuit"]
  TEST_CASES_HYBRID = ["bouncingBallReals"]
  @testset "Backend test" begin
    @testset "DifferentialEquations.jl Backend tests" begin
      testBackend(TEST_CASES_BASIC; mode= OMBackend.DAE_MODE)
      testBackend(TEST_CASES_ADVANCED; mode= OMBackend.DAE_MODE)
      #= Currently failing due to  https://github.com/SciML/Sundials.jl/issues/292 and https://github.com/SciML/Sundials.jl/issues/290=#
      #testBackend(TEST_CASES_HYBRID; mode = OMBackend.DAE_MODE)
    end
    @testset "MTK backend test" begin
      testBackend(TEST_CASES_BASIC; mode= OMBackend.MODELING_TOOLKIT_MODE)
    end
  end
end

runTests()
