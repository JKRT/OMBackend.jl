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

import OMBackend.CodeGeneration.OMSolution

"""
   DAE's to use for the sanity check. 
"""
global TEST_CASES = ["helloWorld", "lotkaVolterra", "vanDerPol", "simpleMech", "simpleCircuit"]
global MODEL_NAME = ""

#= These tests does not validate the semantics. Only that we can successfully compile and simulate this set of models=#
@testset "Sanity test. Simple translation and simulation of Hybrid-DAE" begin
  @testset "DAE-mode" begin
    for testCase in TEST_CASES
      @testset "$testCase" begin
        @testset "compile" begin
          global MODEL_NAME
          frontendDAE::OMBackend.DAE.DAE_LIST = getfield(ExampleDAEs,Symbol("$(testCase)_DAE"))
          try
            #= TODO: We should check this with some reference IR =#
            (MODEL_NAME, modelCode) = OMBackend.translate(frontendDAE)
            @debug generateFile(testCase, modelCode)
            @info "Translated: $MODEL_NAME"
            @test true
          catch e
            @info e
            throw(e)
            @test false
          end
        end
        @testset "simulate" begin
          global MODEL_NAME
          @info "Simulating: $MODEL_NAME"
          try
            simulationResults = OMBackend.simulateModel(MODEL_NAME; tspan=(0.0, 1.0))
            @test simulationResults.diffEqSol.retcode == :Success
          catch e
            @info "Simulation failure: see $(testCase)_result.csv"
            @info "Trying to simulate again"
          end
        end
      end
    end
  end
end
