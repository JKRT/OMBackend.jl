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
#=Author John Tinnerholm & Andreas Heureman=#

import CSV
import Tables

# Load module for test DAEs
using Test
using OMBackend

include("debugUtil.jl")
#= logging =#
ENV["JULIA_DEBUG"] = "all"

const CURRENT_DIRECTORY = @__DIR__
const EXAMPLE_DAE_DIRECTORY = CURRENT_DIRECTORY * "/ExampleDAE"
if ! (CURRENT_DIRECTORY in LOAD_PATH && EXAMPLE_DAE_DIRECTORY in LOAD_PATH)
  @info("Setting up loadpath..")
  push!(LOAD_PATH, CURRENT_DIRECTORY, EXAMPLE_DAE_DIRECTORY)
  @info("Done setting up loadpath: $LOAD_PATH")
end
using ExampleDAEs

"""
These currently work
"""
global TEST_CASES = ["helloWorld", "lotkaVolterra", "vanDerPol"]
#global TEST_CASES = ["bouncingBall"]

#using ExampleDAEs
global MODEL_NAME = ""
@testset "UnitTests" begin
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
            @info MODEL_NAME
            @test true
          catch e
            @info e
            throw(e)
            @test false
          end
        end
        @testset "simulate" begin
          global MODEL_NAME
          try
            simulationResults = OMBackend.simulateModel(MODEL_NAME, (0.0,10.0))
            simTable = Tables.table(simulationResults)
            CSV.write("$(testCase)_result.csv", simTable)
            @debug "Simulation results:" simulationResults
            @test true
          catch e
            @info e
            @test false
          end
        end
        #= TODO: Currently issues due to difference in csv format =#
        @testset "validate solution" begin
          local omc_table::Array = Tables.matrix(CSV.read("./OMC_Reference/$(testCase)_res.csv"))
          local julia_table::Array = Tables.matrix(CSV.read("$(testCase)_result.csv"))
          @test_broken omc_table ≈ julia_table
          #= Delete the file =#
          #rm("$(testCase)_result.csv")
        end
      end
    end
  end

  @testset "ODE-mode" begin
    for testCase in TEST_CASES
      @testset "$testCase" begin
        @testset "compile" begin
          global MODEL_NAME
          frontendDAE::OMBackend.DAE.DAE_LIST = getfield(ExampleDAEs,Symbol("$(testCase)_DAE"))
          try
            #= TODO: We should check this with some reference IR =#
            (MODEL_NAME, modelCode) = OMBackend.translate(frontendDAE)
            @debug generateFile(testCase, modelCode)
            @info MODEL_NAME
            @test true
          catch e
            @info e
            throw(e)
            @test false
          end
        end
        @testset "simulate" begin
          global MODEL_NAME
          try
            simulationResults = OMBackend.simulateModel(MODEL_NAME, (0.0,10.0))
            simTable = Tables.table(simulationResults)
            CSV.write("$(testCase)_result.csv", simTable)
            @debug "Simulation results:" simulationResults
            @test true
          catch e
            @info e
            @test false
          end
        end
        #= TODO: Currently issues due to difference in csv format =#
        @testset "validate solution" begin
          local omc_table::Array = Tables.matrix(CSV.read("./OMC_Reference/$(testCase)_res.csv"))
          local julia_table::Array = Tables.matrix(CSV.read("$(testCase)_result.csv"))
          @test_broken omc_table ≈ julia_table
          #= Delete the file =#
          #rm("$(testCase)_result.csv")
        end
      end
  end
end
