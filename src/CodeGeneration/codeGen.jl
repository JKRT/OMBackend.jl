#=
  This file contains the code generation for the DifferntialEquations.jl backend. 

TODO: 
  Add support for if equations 
  Current approach. One separate function for each branch. 

  Author: John Tinnerholm
=#

"""
  Contains the headerstring defining the OpenModelica copyright notice.
"""
const HEADER_STRING ="
  $(copyRightString())"

#= To keep track of generated callbacks. =#
let CALLBACKS = 0
  global function ADD_CALLBACK()
    CALLBACKS += 1
    return CALLBACKS
  end
  global function RESET_CALLBACKS()
    CALLBACKS = 0
    return CALLBACKS
  end
  global function COUNT_CALLBACKS()
    return CALLBACKS
  end
end

"""
    Write  modelName_DAE_equations() to file.
"""
function writeDAE_equationsToFile(fileName::String, contents::String)
  local fdesc = open(fileName, "w")
  write(fdesc, contents)
  close(fdesc)
end

"""
  Generate Julia code from SimCode.
  Returns an expression in memory representing the program and a string
"""
function generateCode(simCode::SimulationCode.SIM_CODE)::Tuple{String, Expr}
  local stateVariables::Array = []
  local algVariables::Array = []
  local stateDerivatives::Array = []
  local parameters::Array = []
  local crefToSimVarHT = simCode.crefToSimVarHT
  for varName in keys(crefToSimVarHT)
    ixAndVar = crefToSimVarHT[varName]
    local varType = ixAndVar[2].varKind
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algVariables, varName)
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end
  #= 
    Check if we have state variables
  =#
  local systemOfDifferentials = ! isempty(stateVariables)
  #= 
  Decide on what type of solver we are to generate code for
  =#
  if systemOfDifferentials
    (modelName, program) =
      createHelperFunctionsForSystemOfDifferentialEquations(simCode.name,
                                                            simCode,
                                                            stateVariables,
                                                            algVariables,
                                                            stateDerivatives,
                                                            parameters)
  else
    msg = "Systems without differential terms are not yet supported.\n Current system:\n \nAlgebraic variables: $(algVariables) \n state variables: $(stateVariables)\n Parameters: $(parameters)"
    @error msg
  end
  #= Decide what runnable to generate =#
  # Return file content
  return (modelName, program)
end

function createHelperFunctionsForSystemOfDifferentialEquations(modelName::String,
                                                               simCode::SimulationCode.SIM_CODE,
                                                               stateVariables,
                                                               algVariables,
                                                               stateDerivatives,
                                                               parameters)::Tuple{String, Expr}
  stateMarkings = createStateMarkings([], stateVariables, simCode)
  #=Start conditions of algebraic and state variables=#
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquations(algVariables, stateVariables, simCode)
  local DAE_STATE_UPDATE_VECTOR_START_COND = createRealToStateVariableMapping(stateVariables, simCode; toFrom = ("x", "reals"))
  @debug("Generating start conditions")
  local START_CONDTIONS = quote
    function $(Symbol("$(modelName)StartConditions"))(aux, t)
      local x = zeros($(length(stateVariables)))
      local dx = zeros($(length(stateVariables)))
      local p = aux[1]
      local reals = aux[2]
      $START_CONDTIONS_EQUATIONS
      $(DAE_STATE_UPDATE_VECTOR_START_COND...)
      return x, dx
    end
  end
  auxFuncSymbol = Symbol(modelName, "AuxVarsHandler")
  local DAE_EQUATION_FUNCTION = createSolverCode(Symbol(modelName, "DAE_equations"),
                                                 auxFuncSymbol,
                                                 stateVariables,
                                                 simCode.residualEquations,
                                                 simCode.ifEquations,
                                                 simCode;
                                                 eqLhsName = "res", eqRhsName = "reals")
  #= 
    Create function to handle auxilary variables. 
    That is variables that are not a part of the system that is sent to the solver. 
  =#
  local AUX_VARS_FUNCTION = createAuxEquationCode(auxFuncSymbol, algVariables,simCode.residualEquations, simCode; eqLhsName = "reals")
  local DIFFERENTIAL_VARS_FUNCTION = quote
    begin
      function $(Symbol("$(modelName)DifferentialVars"))()
        return $stateMarkings
      end
    end
  end
  local PARAMETER_FUNCTION = createParameterCode(modelName, parameters, stateVariables, algVariables, simCode)
  local CALL_BACK_EQUATIONS = createCallbackCode(modelName, simCode)
  local RUNNABLE = createDAERunnable(modelName, simCode)
  @debug("Code-generation done")
  RESET_CALLBACKS()
  #= 
    The state machine is used to select the active mode 
    0 is the default mode.
  =#
  program = quote    
    $(HEADER_STRING)
    using DiffEqBase
    using DifferentialEquations
    using Plots
    using Sundials
    import OMBackend
    Mode = [0]
    $(START_CONDTIONS)
    $(AUX_VARS_FUNCTION)
    $(DIFFERENTIAL_VARS_FUNCTION)
    $(DAE_EQUATION_FUNCTION)
    $(PARAMETER_FUNCTION)
    $(CALL_BACK_EQUATIONS)
    $(RUNNABLE)
  end
  return (modelName, program)
end

function createCallbackCode(modelName::N, simCode::S) where {N, S}
  local WHEN_EQUATIONS = createEquations(simCode.whenEquations, simCode)
  #= 
    For if equations we create zero crossing functions (Based on the conditions). 
    The body of these equations are evaluated in the main body of the solver itself.
  =#
  local IF_EQUATIONS = createIfEquationCallbacks(simCode.ifEquations, simCode)
  local SAVE_FUNCTION = createSaveFunction(modelName)
  quote
    $(Symbol("saved_values_$(modelName)")) = SavedValues(Float64, Tuple{Float64,Array})
    function $(Symbol("$(modelName)CallbackSet"))(aux)
      local p = aux[1]
      $(LineNumberNode((@__LINE__), "WHEN EQUATIONS"))
      $(WHEN_EQUATIONS...)
      $(LineNumberNode((@__LINE__), "IF EQUATIONS"))
      $(IF_EQUATIONS...)
      $(SAVE_FUNCTION)
      return $(Expr(:call, :CallbackSet, returnCallbackSet()...))
    end
  end
end

function createParameterCode(modelName, parameters, stateVariables, algVariables, simCode)::Expr
  local PARAMETER_EQUATIONS = createParameterEquations(parameters, simCode)
  quote
    function $(Symbol("$(modelName)ParameterVars"))()
      local aux = Array{Array{Float64}}(undef, $(2))
      local p = Array{Float64}(undef, $(arrayLength(parameters)))
      local reals = Array{Float64}(undef, $(arrayLength(stateVariables) + arrayLength(algVariables)))
      aux[1] = p
      aux[2] = reals
      $(PARAMETER_EQUATIONS...)
      return aux
    end
  end
end

"""
  Creates equation code from the set of residual equations and a supplied set of variables.
  This is the method used for the solver.  
  $(SIGNATURES)
"""
function createSolverCode(functionName::Symbol,
                          auxFuncSymbol::Symbol,
                          variables::Vector{V},
                          residuals::Vector{R},
                          ifEquations::Vector{IF_EQ},
                          simCode::SimulationCode.SIM_CODE;
                          eqLhsName, eqRhsName)::Expr where {V, R, IF_EQ}
  local UPDATE_VECTOR = createRealToStateVariableMapping(variables, simCode)
  #= Creates the equations =#
  local EQUATIONS = createEquations(variables, residuals, simCode; eqLhsName = eqLhsName, eqRhsName)
  local IF_EQUATIONS = createEquations(variables, ifEquations, simCode; eqLhsName = eqLhsName, eqRhsName = eqRhsName)
  local modelName = simCode.name
  quote
    function $(functionName)(res, dx, x, aux, t)
      $(auxFuncSymbol)(res, dx, x, aux, t)
      local p = aux[1]
      local reals = aux[2]
      $(EQUATIONS...)
      $(IF_EQUATIONS...)
      $(UPDATE_VECTOR...)
    end
  end
end

"""
  This function creates code for the auxilary variables. 
  These variables are not passed to the solver, rather they are updated here 
  after a solver has run it's iterations.
"""
function createAuxEquationCode(functionName::Symbol, variables::Array{V},
                               residuals::Array{R}, simCode::SimulationCode.SIM_CODE;
                               eqLhsName) where {V, R}
  local AUX_EQUATIONS = createAuxEquationCode(variables, simCode; arrayName = eqLhsName)
  quote
    function $(functionName)(res, dx, x, aux, t)
      local p = aux[1]
      local reals = aux[2]
      $(AUX_EQUATIONS...)
    end
  end  
end


function createDAERunnable(modelName::String, simCode::SimulationCode.SIM_CODE)
  quote
    function $(Symbol("$(modelName)Simulate"))(tspan = (0.0, 1.0))
      $(LineNumberNode((@__LINE__), "Auxilary variables"))   
      local aux = $(Symbol("$(modelName)ParameterVars"))()
      (x0, dx0) =$(Symbol("$(modelName)StartConditions"))(aux, tspan[1])
      local differential_vars = $(Symbol("$(modelName)DifferentialVars"))()
      #= Pass the residual equations =#
      local problem = DAEProblem($(Symbol("$(modelName)DAE_equations")), dx0, x0, 
                                 tspan, aux, differential_vars=differential_vars, 
                                 callback=$(Symbol("$(modelName)CallbackSet"))(aux))
      #= Solve with IDA =#
      local solution = solve(problem, IDA())
      #= Convert into OM compatible format =#
      local savedSol = map(collect, $(Symbol("saved_values_$(modelName)")).saveval)
      local t = [savedSol[i][1] for i in 1:length(savedSol)]
      local vars = [savedSol[i][2] for i in 1:length(savedSol)]
      local T = eltype(eltype(vars))
      local N = length(aux[2])
      local nsolution = DAESolution{Float64,N,typeof(vars),Nothing, Nothing, Nothing, typeof(t),
                                    typeof(problem),typeof(solution.alg),
                                    typeof(solution.interp),typeof(solution.destats)}(
                                      vars, nothing, nothing, nothing, t, problem, solution.alg,
                                      solution.interp, solution.dense, 0, solution.destats, solution.retcode)
      ht = $(SimulationCode.makeIndexVarNameUnorderedDict(simCode.matchOrder, simCode.crefToSimVarHT))
      omSolution = OMBackend.Runtime.OMSolution(nsolution, ht)
      return omSolution
    end
  end
end


"""
  This method creates a runnable for a linear/non-linear system of equations.
  That is a system that does not contain differential equations
"""
function createLinearRunnable(modelName::String, simCode::SimulationCode.SIM_CODE)
  quote
    import NonlinearSolve
    function $(Symbol("$(modelName)Simulate"))(tspan = (0.0, 1.0))
      $(LineNumberNode((@__LINE__), "Auxilary variables"))   
      local aux = $(Symbol("$(modelName)ParameterVars"))()
      (x0, dx0) =$(Symbol("$(modelName)StartConditions"))(aux, tspan[1])
      local differential_vars = $(Symbol("$(modelName)DifferentialVars"))()
      #= Pass the residual equations =#
      local problem = NonlinearProblem($(Symbol("$(modelName)DAE_equations")), dx0, x0, 
                                       tspan, aux, differential_vars=differential_vars, 
                                       callback=$(Symbol("$(modelName)CallbackSet"))(aux))
      #= Solve with IDA =#
      local solution = Runtime.solve(problem::NonlinearProblem, IDA())
      #= Convert into OM compatible format =#
      local savedSol = map(collect, $(Symbol("saved_values_$(modelName)")).saveval)
      local t = [savedSol[i][1] for i in 1:length(savedSol)]
      local vars = [savedSol[i][2] for i in 1:length(savedSol)]
      local T = eltype(eltype(vars))
      local N = length(aux[2])
      local nsolution = DAESolution{Float64,N,typeof(vars),Nothing, Nothing, Nothing, typeof(t),
                                    typeof(problem),typeof(solution.alg),
                                    typeof(solution.interp),typeof(solution.destats)}(
                                      vars, nothing, nothing, nothing, t, problem, solution.alg,
                                      solution.interp, solution.dense, 0, solution.destats, solution.retcode)
      ht = $(SimulationCode.makeIndexVarNameUnorderedDict(simCode.matchOrder, simCode.crefToSimVarHT))
      omSolution = OMBackend.Runtime.OMSolution(nsolution, ht)
      return omSolution
    end
  end
end



"""
  This functions creates the update equations for the auxilary variables.
  The set of auxilary variables is the set of variables of other types than state variables.
  That is booleans integers and algebraic variables.
TODO:
  Currently only being done for the algebraic variables.
"""
function createAuxEquationCode(algVariables::Array{V},
                               simCode::SimulationCode.SIM_CODE
                               ;arrayName)::Array{Expr} where {V}
  #= Sorted equations for the algebraic variables. =#
  local auxEquations::Array{Expr} = [] 
  auxEquations = vcat(createSortedEquations([algVariables...], simCode; arrayName = "reals"))
  return auxEquations
end

function createStateMarkings(algVariables::Array, stateVariables::Array, simCode::SimulationCode.SIM_CODE)::Array{Bool}
  local stateMarkings::Array = [false for i in 1:length(stateVariables) + length(algVariables)]
  for sName in stateVariables
    stateMarkings[simCode.crefToSimVarHT[sName][1]] = true
  end
  return stateMarkings
end

"""
  Creates the save-callback. 
  saved_values_\$(modelName) is provided 
  as a shared global for the specific model under compilation.
"""
function createSaveFunction(modelName)::Expr
  ADD_CALLBACK()
  local callbacks = COUNT_CALLBACKS()
  local cbSym = Symbol("cb$(callbacks)")
  return quote
    savingFunction(u, t, integrator) = let
      (t, deepcopy(integrator.p[2]))
    end
    $cbSym = SavingCallback(savingFunction, $(Symbol("saved_values_$(modelName)")))
  end
end

"""
  Returns the argument array for the callback set.
"""
function returnCallbackSet()::Array
  local cbs::Array{Symbol} = []
  for t in 1:COUNT_CALLBACKS()
    cb = Symbol("cb", t)
    push!(cbs, cb)
  end
  expr::Array = if length(cbs) < 1
    []
  else
    cbs
  end
  return cbs 
end

function createRealToStateVariableMapping(stateVariables::Array, simCode::SimulationCode.SIM_CODE; toFrom::Tuple=("reals", "x"))::Array{Expr}
  local daeStateUpdateVector::Array = []
  for svName in stateVariables
    local varIdx = simCode.crefToSimVarHT[svName][1]
    push!(daeStateUpdateVector, :($(Symbol(toFrom[1]))[$varIdx] = $(Symbol(toFrom[2]))[$varIdx]))
  end
  return daeStateUpdateVector
end

"""
      TODO: We use state for everything now...
      Algebraic and state should be initilised in a different way.
"""
function createStartConditionsEquations(algVariables::Array,
                                        stateVariables::Array,
                                        simCode::SimulationCode.SIM_CODE)::Expr
  #= Create the sorted set of equation to occur in the preamble of the method. =#
  local startEquations = createSortedEquations([algVariables..., stateVariables...], simCode; arrayName = "reals")
  return quote
    $(getStartConditions(algVariables, "reals", simCode))
    $(getStartConditions(stateVariables, "reals", simCode)) #TODO: Should be DX0?
    $(startEquations...)
  end
end

"
  Creates sorted equations for all supplied variables that are present in the strongly connected components.
  Author: johti17
  # Arguments
- `variables::Array`: Variables involved in the set of equations to be sorted
- `simCode::SimulationCode.SIM_CODE`: Simulation code
- `arrayName::String` : The name of the array used in the generated code. 
"
function createSortedEquations(variables::Array, simCode::SimulationCode.SIM_CODE; arrayName::String)::Array{Expr}
  #= Assign the variables according to sorting/matching =#
  local ht = SimulationCode.makeIndexVarNameDict(simCode.matchOrder, simCode.crefToSimVarHT)
  local vIndicesSupplied = [simCode.crefToSimVarHT[v][1] for v in variables]
  local iterationComponents = []
  local sortedEquations = []
  for i in simCode.matchOrder
    variableIdx = MetaGraphs.get_prop(simCode.equationGraph, i, :vID)
    if ! (ht[variableIdx] in variables)
      continue
    end
    push!(iterationComponents, i)
  end
  local residuals = simCode.residualEquations
  #=
    For state variables we write to dx0. Thus, we are looking for dx0[<index>].
  =#
  for i in iterationComponents
    variableIdx = MetaGraphs.get_prop(simCode.equationGraph, i, :vID)
    equationIdx = simCode.matchOrder[variableIdx]
    local equation = stripComments(expToJuliaExp(residuals[equationIdx].exp, simCode; varPrefix = arrayName))
    local rewrittenEquation::Expr = stripBeginBlocks(arrayToSymbolicVariable(equation))
    local varToSolve = simCode.crefToSimVarHT[ht[variableIdx]][2]
    local varToSolveExpr::Symbol = SimulationCode.isState(varToSolve) ? Symbol("dx_$(variableIdx)") : Symbol("$(arrayName)_$(variableIdx)")
    #= Solve equation symbolically (Or numerically) =#
    local reduceSolution::Expr = Reduce.Algebra.solve(rewrittenEquation, varToSolveExpr)[1]
    local postProcessSolution::Expr = symbolicVariableToArrayRef(reduceSolution)
    push!(sortedEquations, postProcessSolution)
  end
  return sortedEquations
end

"""
  Create equation for a set of variables V and a set of residual equations, E.
  Equations are only created if V ∈ E.
C is the context of the residual equations it should either be a branch or simcode. 
"""
function createEquations(variables::Array{V}, equations::Array{BDAE.RESIDUAL_EQUATION}, context::C;
                         eqLhsName::String, eqRhsName::String) where {V, C}
  #= Check if we have any equations =#
  if isempty(equations)
    return Expr[]
  end
  #= Generate the sorted set of equations and variables =#
  local eqs::Vector{Expr} = Expr[]
  for (equationCounter, variable) in enumerate(variables)
    local equation = SimulationCode.getEquationSolvedIn(variable, context)
    local eqJL::Expr = eqToJulia(equation, context, equationCounter; eqLhsName = eqLhsName, eqRhsName = eqRhsName)
    push!(eqs, eqJL)
  end
  return eqs
end

"""
  Creates the set of if equations
"""
function createEquations(variables::Vector{V},
                         equations::Vector{SimulationCode.IF_EQUATION},
                         simCode::SimulationCode.SIM_CODE;
                         eqLhsName::String, eqRhsName::String) where {V}
  local if_eqs::Vector{Expr} = Expr[]
  for if_eq in equations #= Of type if equation =#    
    res = createIfEquation(variables,
                           if_eq,
                           simCode;
                           eqLhsName = eqLhsName,
                           eqRhsName = eqRhsName)
    push!(if_eqs, res)
  end
  return if_eqs
end

"""
  This function creates a single if-equation.
  It does so by creating a if-elseif-else statement
  where the contents of each branch is casualised. 
  Which branch is active depends on the active mode
"""
function createIfEquation(variables::Vector{V},
                          if_eq::SimulationCode.IF_EQUATION,
                          simCode::SimulationCode.SIM_CODE;
                          eqLhsName::String, eqRhsName::String) where {V}
  local data = Expr[]
  local branches = if_eq.branches
  #= If empty opt out =#
  if isempty(branches)
    return data
  end
  local branches = if_eq.branches
  local nBranches = length(branches)
  local containsElse = branches[nBranches].identifier == SimulationCode.ELSE_BRANCH
  local haveElseIfBranches = nBranches > 2 #= Should be like this if the frontend did the job :) =#
  #= (The first branch is the if. The second branch is the else the reminder should be the number of elseifs) =#
  local nElseIfBranches = nBranches - 2
  if !(haveElseIfBranches)
    return
    if nBranches == 1
      :(if Mode[1] == $(branches.identifier)
        $(flattenExprs(createEquations(variables, branches[1].residualEquations, branches[1]; eqLhsName = eqLhsName, eqRhsName = eqRhsName)))
        end)
    else
      :(if Mode[1] == $(branches.identifier)
          $(flattenExprs(createEquations(variables, branches[1].residualEquations, branches[1]; eqLhsName = eqLhsName, eqRhsName = eqRhsName)))
          else
          $(flattenExprs(createEquations(variables, branches[2].residualEquations, branches[2]; eqLhsName = eqLhsName, eqRhsName = eqRhsName)))
        end
        )
    end
  end
  #= We have a number of elseif branches > 0 =#
  eq = flattenExprs(createEquations(variables, branches[1].residualEquations, branches[1]; eqLhsName = eqLhsName, eqRhsName = eqRhsName))
  local ifElseList::Expr = Expr(:if, :(Mode[1] == $(branches[1].identifier)), eq)
  #= Create the elseif branches for this if equation. 
  The else if branches are indexed with i + 1 since the first index is contained by the if equation =#
  iter = ifElseList
  for i in 1:nElseIfBranches
    eq = flattenExprs(createEquations(variables, branches[i + 1].residualEquations, branches[i + 1], eqLhsName = eqLhsName, eqRhsName = eqRhsName))
    push!(iter.args, Expr(:elseif, :(Mode[1] == $(branches[i + 1].identifier)), eq))
    #= Append the expression to the back =#
    iter = iter.args[end]
  end
  eq = flattenExprs(createEquations(variables, branches[end].residualEquations, branches[nBranches]; eqLhsName = eqLhsName, eqRhsName = eqRhsName))
  push!(iter.args, eq)
  #=Done!=#
  res = ifElseList
  return res
end


"""
 Create a set for all equations T. 
"""
function createEquations(equations::Array{T}, simCode::SimulationCode.SIM_CODE)::Array{Expr} where T
  local eqs::Array{Expr} = Expr[]
  for (equationCounter, eq) in enumerate(equations)
    local eqJL::Expr = eqToJulia(eq, simCode, equationCounter)
    push!(eqs, eqJL)
  end
  return eqs
end

"""
  Create equations for the parameters.
"""
function createParameterEquations(parameters::Array, simCode::SimulationCode.SimCode)
  local parameterEquations::Array = []
  local hT = simCode.crefToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    local simVarType::SimulationCode.SimVarType = simVar.varKind
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => exp
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
    push!(parameterEquations,
          quote
          $(LineNumberNode(@__LINE__, "$param"))
          p[$index] = $(expToJuliaExp(bindExp, simCode))
          end
          )
  end
  return parameterEquations
end

"""
 $(SIGNATURES)
"""

function eqToJulia(eq::BDAE.RESIDUAL_EQUATION,
                   context::C,
                   arrayIdx::Int64; eqLhsName::String, eqRhsName::String) where {C}
  quote 
    $(Symbol((eqLhsName)))[$arrayIdx] = $(expToJuliaExp(eq.exp, context; varPrefix = eqRhsName))
  end
end

"""
 Transforms a given BDAE.Equation equation into Julia code.
 $(SIGNATURES)

 If callbacks are to be generated increase the internal variable, callback counter representing
 how many callbacks we currently have.
"""
function eqToJulia(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, arrayIdx::Int64)
  throw("Unknown equation")
end

function createIfEquationCallbacks(ifEqs::Vector{SimulationCode.IF_EQUATION}, simCode::SimulationCode.SIM_CODE)
  res = []
  for eq in ifEqs
    push!(res, createIfEqCallback(eq, simCode))
  end
  if isempty(res)
    return res
  end
  return res[1]
end

"""
    Converts a SIMCODE.IF_EQUATION into a Julia representation.
    Note that this only creates the callback for the particular equation (Not the branching logic).
    The branching logic itself is generated in the main function which is passed to a solver
    along with the residuals for each branch. 
"""
function createIfEqCallback(ifEq::SimulationCode.IF_EQUATION, simCode::SimulationCode.SIM_CODE)
  exprs = Expr[]
  branches = ifEq.branches
  for branch in branches
    #= No callback to be generated for the else branch. =#
    if branch.identifier == SimulationCode.ELSE_BRANCH
      continue
    end
    ADD_CALLBACK()
    local callbacks = COUNT_CALLBACKS()
    #=TODO: 
    Continuous callback if the equation is a part of the system that is solved. 
    If not some other type should be generated.
    The type should be decided depending on the kind of variable involved in the condition.
    =#
    local cond = transformToZeroCrossingCondition(branch.condition)
    local res = quote
      function $(Symbol("condition$(callbacks)"))(x,t,integrator) 
        $(expToJuliaExp(cond, simCode))
      end
      #= Active this branc =#
      function $(Symbol("affect$(callbacks)!"))(integrator)
        @info "Branch activated"
        global Mode[1] = $(branch.identifier)
      end
      #= Deactive this branch =#
      function $(Symbol("affect_neg$(callbacks)!"))(integrator)
        @info "Branch deactivated"
        global Mode[1] = 0
      end

    #=Below zero triggers callback zero crossing from positive to negative. The reverse deactivates it =#
    $(Symbol("cb$(callbacks)")) = ContinuousCallback($(Symbol("condition$(callbacks)")), 
                                                     $(Symbol("affect_neg$(callbacks)!")),
                                                     rootfind=true, save_positions=(true, true),
                                                     affect_neg! = $(Symbol("affect$(callbacks)!")),)
    end
    push!(exprs, res)
  end
  return exprs
end

"""
  This function creates a representation of a when equation in Julia.
  This function shall be called in the process of constructing the different callbacks.
"""
function eqToJulia(eq::BDAE.WHEN_EQUATION, simCode::SimulationCode.SIM_CODE, arrayIdx::Int64)::Expr
  local wEq = eq.whenEquation
  local whenStmts = createWhenStatements(wEq.whenStmtLst, simCode, prefix="integrator.u")
  local cond = transformToZeroCrossingCondition(wEq.condition)
  ADD_CALLBACK()
  local callbacks = COUNT_CALLBACKS()
  quote
    function $(Symbol("condition$(callbacks)"))(x,t,integrator) 
      $(expToJuliaExp(cond, simCode))
    end
    function $(Symbol("affect$(callbacks)!"))(integrator)
      $(whenStmts...)
    end
    $(Symbol("cb$(callbacks)")) = ContinuousCallback($(Symbol("condition$(callbacks)")), 
                                                            $(Symbol("affect$(callbacks)!")),
                                                            rootfind=true, save_positions=(false, false),
                                                            affect_neg! = $(Symbol("affect$(callbacks)!")),)
  end
end


"""
 Creates Julia code for the set of whenStatements in the when equation.
 There are some constructs that may only occur in a when equations. 
"""
function createWhenStatements(whenStatements::List, simCode::SimulationCode.SIM_CODE; prefix="x")::Array{Expr}
  local res::Array{Expr} = []
  @debug "Calling createWhenStatements with: $whenStatements"
  for wStmt in  whenStatements
    @match wStmt begin
      BDAE.ASSIGN(__) => begin
        push!(res, quote
              $(expToJuliaExp(wStmt.left, simCode, varPrefix=prefix)) = $(expToJuliaExp(wStmt.right, simCode, varPrefix=prefix))
              end)
      end
      #= Handles reinit =#
      BDAE.REINIT(__) => begin
        (index, type) = simCode.crefToSimVarHT[BDAE.string(wStmt.stateVar)]
        push!(res, quote 
              integrator.u[$(index)] = $(expToJuliaExp(wStmt.value, simCode, varPrefix=prefix))
              end)
      end
      _ => ErrorException("$whenStatements in @__FUNCTION__ not supported")
    end
  end
  return res
end


"""
  Converts a DAE expression into a Julia expression
  $(SIGNATURES)
The context can be any type that contains a set of residual equations.

"""
function expToJuliaExp(exp::DAE.Exp, context::C, varSuffix=""; varPrefix="x")::Expr where {C}
  hashTable = context.crefToSimVarHT
  local expr::Expr = begin
    local int::Int64
    local real::Float64
    local bool::Bool
    local tmpStr::String
    local cr::DAE.ComponentRef
    local e1::DAE.Exp
    local e2::DAE.Exp
    local e3::DAE.Exp
    local expl::List{DAE.Exp}
    local lstexpl::List{List{DAE.Exp}}
    @match exp begin
      DAE.BCONST(bool) => quote $bool end
      DAE.ICONST(int) => quote $int end
      DAE.RCONST(real) => quote $real end
      DAE.SCONST(tmpStr) => quote $tmpStr end
      DAE.CREF(cr, _)  => begin
        varName = BDAE.string(cr)
        builtin = if varName == "time"
          true
        else
          false
        end
        if ! builtin
          #= If we refer to time, we  return t instead of a concrete variable =#
          indexAndVar = hashTable[varName]
          varKind::SimulationCode.SimVarType = indexAndVar[2].varKind
          @match varKind begin
            SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
            SimulationCode.STATE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName state"))
              $(Symbol(varPrefix))[$(indexAndVar[1])]
            end
            SimulationCode.PARAMETER(__) => quote
              $(LineNumberNode(@__LINE__, "$varName parameter"))
              p[$(indexAndVar[1])]
            end
            SimulationCode.ALG_VARIABLE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName, algebraic"))
              $(Symbol(varPrefix))[$(indexAndVar[1])]
            end
            SimulationCode.STATE_DERIVATIVE(__) => :(dx$(varSuffix)[$(indexAndVar[1])] #= der($varName) =#)
          end
        else #= Currently only time is a builtin variabe. Time is represented as t in the generated code =#
          quote
            t
          end
        end
      end
      DAE.UNARY(operator = op, exp = e1) => begin
        o = DAE_OP_toJuliaOperator(op)
        quote
          $(o)($(expToJuliaExp(e1, context, varPrefix=varPrefix)))
        end
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        a = expToJuliaExp(e1, context, varPrefix=varPrefix)
        b = expToJuliaExp(e2, context, varPrefix=varPrefix)
        o = DAE_OP_toJuliaOperator(op)
        quote
          $o($(a), $(b))
        end
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        lhs = expToJuliaExp(e1, context, varPrefix=varPrefix)
        o = DAE_OP_toJuliaOperator(op)
        quote
          $o($(lhs))
        end
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExp(e1, context, varPrefix=varPrefix) + " " + BDAE.string(op) + " " + expToJuliaExp(e2, simCode, varPrefix=varPrefix))
        end
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        lhs = expToJuliaExp(e1, context, varPrefix=varPrefix)
        o = DAE_OP_toJuliaOperator(op)
        rhs = expToJuliaExp(e2, context, varPrefix=varPrefix)
        quote 
          $o($(lhs), $(rhs))
        end
      end
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        throw(ErrorException("If expressions not allowed in backend code. 
                              They should have been eliminated by a previous backend  pass."))
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst)  => begin
        DAECallExpressionToJuliaCallExpression(tmpStr, explst, context, hashTable, varPrefix=varPrefix)
      end
      DAE.CAST(ty, exp)  => begin
        quote
          $(generateCastExpression(ty, exp, context, varPrefix))
        end
      end
      _ =>  throw(ErrorException("$exp not yet supported"))
    end
  end
  return expr
end

function generateCastExpression(ty, exp, simCode, varPrefix)
  return @match ty, exp begin
    (DAE.T_REAL(__), DAE.ICONST(__)) => float(eval(expToJuliaExp(exp, simCode, varPrefix=varPrefix)))
    (DAE.T_REAL(__), DAE.CREF(cref)) where typeof(cref.identType) == DAE.T_INTEGER  => float(eval(expToJuliaExp(exp, simCode, varPrefix=varPrefix)))
    _ => throw("Cast $ty: for exp: $exp not yet supported in codegen!")
  end
end

"""
  Generates the start conditions.
  All variables default to zero if they are not specified by the user.
"""
function getStartConditions(vars::Array, condName::String, simCode::SimulationCode.SimCode)::Expr
  local startExprs::Array{Expr} = []
  local residuals = simCode.residualEquations
  local ht::Dict = simCode.crefToSimVarHT
  if length(vars) == 0
    return quote
    end
  end
  for var in vars
    (index, simVar) = ht[var]
    local simVarType = simVar.varKind
    local optAttributes::Option{DAE.VariableAttributes} = simVar.attributes
    if simVar.attributes == nothing
      continue
    end
    () = @match optAttributes begin
      SOME(attributes) => begin
        () = @match (attributes.start, attributes.fixed) begin
          (SOME(start), SOME(fixed)) || (SOME(start), _)  => begin
            @debug "Start value is:" start
            push!(startExprs,
                  quote
                  $(LineNumberNode(@__LINE__, "$var"))
                  $(Symbol("$condName"))[$index] = $(expToJuliaExp(start, simCode))
                  end)
            ()
          end
          (NONE(), SOME(fixed)) => begin
            push!(startExprs, :($(condName)[$(index)] = 0.0))
            ()
          end
          (_, _) => ()
        end
        NONE() => ()
      end
    end
  end
  return quote
    $(startExprs...)
  end
end
