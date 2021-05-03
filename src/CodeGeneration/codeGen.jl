"""
"""
const HEADER_STRING ="
  $(copyRightString())"

"Counter for generating callbacks for a model"
global CALLBACK_COUNTER = 0

"""
    Write  modelName_DAE_equations() to file.
"""
function writeDAE_equationsToFile(fileName::String, contents::String)
  local fdesc = open(fileName, "w")
  write(fdesc, contents)
  close(fdesc)
end

"
    ODE-mode code generation
  "
function generateCode(simCode::SimulationCode.EXPLICIT_SIM_CODE)
  throw("Not implemented..")
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
  local stateMarkings::Array = []
  local crefToSimVarHT = simCode.crefToSimVarHT
  local modelName::String = simCode.name
  local exp::DAE.Exp
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
  stateMarkings = createStateMarkings([], stateVariables, simCode)
  #=Start conditions of algebraic and state variables=#
  #= Discrete events=#
  #=Continuous events=#
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

  local AUX_FUNC_NAME = Symbol("$(modelName)AuxVarsHandler")
  local DAE_STATE_UPDATE_VECTOR = createRealToStateVariableMapping(stateVariables, simCode)
  local DAE_EQUATIONS = createResidualEquations(stateVariables, simCode.residualEquations, simCode)
  local DAE_EQUATION_FUNCTION = quote
    function $(Symbol("$(modelName)DAE_equations"))(res, dx, x, aux, t)
      $(AUX_FUNC_NAME)(dx, x, aux, t)
      reals = aux[2]
      p = aux[1]
      $(DAE_EQUATIONS...)
      $(DAE_STATE_UPDATE_VECTOR...)
    end
  end

  local AUX_VARS_EQUATIONS = createAuxVarEquations(algVariables::Array, stateVariables::Array, simCode::SimulationCode.SIM_CODE)
  local AUX_VARS_FUNCTION = quote
    function $(AUX_FUNC_NAME)(dx, x, aux, t)
      p = aux[1]
      reals = aux[2]
      $(AUX_VARS_EQUATIONS...)
    end
  end
  
  local DIFFERENTIAL_VARS_FUNCTION = quote
    begin
      function $(Symbol("$(modelName)DifferentialVars"))()
        return $stateMarkings
      end
    end
  end

  local PARAMETER_EQUATIONS = createParameterEquations(parameters, simCode)
  local PARAMETER_VARS = quote
    function $(Symbol("$(modelName)ParameterVars"))()
      aux = Array{Array{Float64}}(undef, $(2))
      p = Array{Float64}(undef, $(arrayLength(parameters)))
      rs = Array{Float64}(undef, $(arrayLength(stateVariables) + arrayLength(algVariables)))
      aux[1] = p
      aux[2] = rs
      $(PARAMETER_EQUATIONS...)
      return aux
    end
  end
  local WHEN_EQUATIONS = createEquations(simCode.whenEquations, simCode)
  local IF_EQUATIONS = createEquations(simCode.ifEquations, simCode)
  local SAVE_FUNCTION = createSaveFunction(modelName)
  local CALL_BACK_EQUATIONS = quote
    $(Symbol("saved_values_$(modelName)")) = SavedValues(Float64, Tuple{Float64,Array})
    function $(Symbol("$(modelName)CallbackSet"))(p)
      $(LineNumberNode((@__LINE__), "WHEN EQUATIONS"))
      $(WHEN_EQUATIONS...)
      $(LineNumberNode((@__LINE__), "IF EQUATIONS"))
      $(IF_EQUATIONS...)
      $(SAVE_FUNCTION)
      return CallbackSet($(evaluateCallBackset()))
    end
  end
  local RUNNABLE = quote
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
      ht = $(SimulationCode.makeIndexVarNameDict(simCode.matchOrder, simCode.crefToSimVarHT))
      omSolution = OMBackend.CodeGeneration.OMSolution(nsolution, ht)
      return omSolution
    end
  end
  @debug("Code-generation done")
  global CALLBACK_COUNTER = 0
  program = quote    
    $(HEADER_STRING)
    using DiffEqBase
    using DifferentialEquations
    using Plots
    using Sundials
    using DataStructures
    import OMBackend
    $(START_CONDTIONS)
    $(AUX_VARS_FUNCTION)
    $(DIFFERENTIAL_VARS_FUNCTION)
    $(DAE_EQUATION_FUNCTION)
    $(PARAMETER_VARS)
    $(CALL_BACK_EQUATIONS)
    $(RUNNABLE)
  end
  # Return file content
  return (modelName, program)
end

"
  This functions creates the update equations for the auxilary variables.
  Two sets of equations are created.
  The first  handle the sorted equations of the set of algebraic variables. 
  The second set of equations handle the 
   
"
function createAuxVarEquations(algVariables::Array, stateVariables::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  #= Sorted equations for the algebraic variables. =#
  local auxEquations::Array{Expr} = [] 
  auxEquations = vcat(createSortedEquations([algVariables...], simCode; arrayName = "reals"))
  return auxEquations
end

function createStateMarkings(algVariables::Array, stateVariables::Array, simCode::SimulationCode.SIM_CODE)::Array{Bool}
  local stateMarkings::Array = [false for i in 1:length(stateVariables)+length(algVariables)]
  for sName in stateVariables
    stateMarkings[simCode.crefToSimVarHT[sName][1]] = true
  end
  return stateMarkings
end

"
  Creates the save-callback. 
  saved_values_\$(modelName) is provided 
  as a shared global for the specific model under compilation.
"
function createSaveFunction(modelName)::Expr
  global CALLBACK_COUNTER += 1
  local cbSym = Symbol("cb$(CALLBACK_COUNTER)")
  return quote
    savingFunction(u, t, integrator) = let
      (t, deepcopy(integrator.p[2]))
    end
    $cbSym = SavingCallback(savingFunction, $(Symbol("saved_values_$(modelName)")))
    return $cbSym
  end
end

function evaluateCallBackset()::Expr
  local cbs::Array{Expr} = []
  for t in 1:CALLBACK_COUNTER
    push!(cbs,
          quote
          $(Symbol("cb$(t)"))
          end)
  end
  expr::Expr = if length(cbs) < 1
    quote end
  else
    quote
      $((cbs...))
    end
  end
  return expr 
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
  local startEquations = createSortedEquations([algVariables..., stateVariables...], simCode; arrayName = "reals")
  return quote
    $(getStartConditions(algVariables, "reals", simCode))
    $(getStartConditions(stateVariables, "reals", simCode)) #Should be DX0?
    $(startEquations...)
  end
end

function createAlgebraicEquations(algVariables::Array, simCode::SimulationCode.SIM_CODE)
  local algEquations = createSortedEquations([algVariables...], simCode)
  return quote
    $(algEquations...)
  end
end

"
  Creates sorted equations for all supplied variables that are present in the strongly connected components.
  Author: johti17
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

"
 Creates the residual equations. 
TODO: only create a residual equations if an equation contains a state variable and we have atleast one state variable in the system
"
function createResidualEquations(stateVariables, equations::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  local eqs::Array{Expr} = []
  local ht = simCode.crefToSimVarHT
  for i in 1:length(stateVariables)
    local equationCounter = i
    #= Get the variable index for this statevariable =#
    local variableIdx = ht[stateVariables[i]][1]
    #= Look up the correct equation index with this index=#
    local equationIdx = simCode.matchOrder[variableIdx]
    #= Use this equation=#
    local equation = equations[equationIdx]
    #= Generate the residual =#
    push!(eqs, residualEqtoJulia(equation, simCode, equationCounter))
  end
  return eqs
end

function createEquations(equations::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  local eqs::Array{Expr} = []
  for (equationCounter, eq) in enumerate(equations)
    local eqJL::Expr = eqToJulia(eq, simCode, equationCounter)
    push!(eqs, eqJL)
  end
  return eqs
end

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
 Transforms a given equation into Julia code
"""
function residualEqtoJulia(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, resNumber)::Expr
  local rhs::DAE.Exp
  local result::Expr = @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
       quote
         res[$(resNumber)] = $(expToJuliaExp(rhs, simCode;varPrefix="reals"))
       end
    end
    _ => begin
      throw("traversalError for $eq")
    end
  end
  return result
end

"""
 Transforms a given equation into Julia code
 The optional parameter prefix specifices the prefix.
 The prefix is used for...
"""
function eqToJulia(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, resNumber)::Expr
  local rhs::DAE.Exp
  @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
      "  $(prefix)[$resNumber] = " + expToJuliaExp(rhs, simCode; varPrefix = prefix) + "\n"
    end
    BDAE.WHEN_EQUATION(_, wEq, _) => begin
      local whenStmts = createWhenStatements(wEq.whenStmtLst, simCode, prefix="integrator.u")
      local cond = prepForZeroCrossing(wEq.condition)
      global CALLBACK_COUNTER += 1
      quote
        function $(Symbol("condition$(CALLBACK_COUNTER)"))(x,t,integrator) 
          $(expToJuliaExp(cond, simCode))
        end
        function $(Symbol("affect$(CALLBACK_COUNTER)!"))(integrator)
          $(whenStmts...)
        end
        $(Symbol("cb$(CALLBACK_COUNTER)")) = ContinuousCallback($(Symbol("condition$(CALLBACK_COUNTER)")), 
                                                                $(Symbol("affect$(CALLBACK_COUNTER)!")),
                                                                rootfind=true, save_positions=(false, false),
                                                                affect_neg! = $(Symbol("affect$(CALLBACK_COUNTER)!")),)
      end
    end
    BDAE.IF_EQUATION(conds, trueEqs, falseEqs) => begin
      global CALLBACK_COUNTER += 1
      f(x) = eqToJulia(x, simCode, resNumber, prefix = "integrator.u")
      local condsJL = tuple(map((x) -> expToJuliaExp(x, simCode; varPrefix="x"), conds)...)
      local trueEqJL = map(f, trueEqs)
      local falseEqJL = map(f, falseEqs)
      quote 
        function $(Symbol("condition$(CALLBACK_COUNTER)"))(x,t,integrator)
          Bool(floor($(condsJL...)))
        end
        function $(Symbol("affect$(CALLBACK_COUNTER)!"))(integrator)
          $(trueEqJL...)
        end
        $(Symbol("cb$(CALLBACK_COUNTER)")) = ContinuousCallback(Symbol("condition$(CALLBACK_COUNTER)"),
                                                                Symbol("affect$(CALLBACK_COUNTER)!"))
        $(global CALLBACK_COUNTER += 1)
        function $(Symbol("condition$(CALLBACK_COUNTER)"))(x,t,integrator)
          ! Bool(floor($(condsJL...)))
        end
        function $(Symbol("affect$(CALLBACK_COUNTER)!"))(integrator)
          $(falseEqJL...)
        end
        createSymbol("cb$(CALLBACK_COUNTER)") = ContinuousCallback(
          Symbol("condition$(CALLBACK_COUNTER)"),
          Symbol("affect$(CALLBACK_COUNTER)!"))
      end
    end
    _ => begin
      throw("Error when generating code for $eq")
    end
  end
end


"""
 Creates Julia code for  when equation.
 E.g a set of Discrete callbacks
 TODO: Handle the other cases
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


"Converts a DAE expression into a Julia expression"
function expToJuliaExp(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE, varSuffix=""; varPrefix="x")::Expr
  hashTable = simCode.crefToSimVarHT
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
          #= If we refeer to time, we simply return t instead of a concrete variable =#
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
          $(o)($(expToJuliaExp(e1, simCode, varPrefix=varPrefix)))
        end
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        a = expToJuliaExp(e1, simCode, varPrefix=varPrefix)
        b = expToJuliaExp(e2, simCode, varPrefix=varPrefix)
        o = DAE_OP_toJuliaOperator(op)
        quote
          $o($(a), $(b))
        end
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        quote
          $("(" + BDAE.string(op) + " " + expToJuliaExp(e1, simCode, varPrefix=varPrefix) + ")")
        end
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExp(e1, simCode, varPrefix=varPrefix) + " " + BDAE.string(op) + " " + expToJuliaExp(e2, simCode, varPrefix=varPrefix))
        end
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExp(e1, simCode,
                          varPrefix=varPrefix) + " "
            + BDAE.string(op) + " " + expToJuliaExp(e2, simCode,varPrefix=varPrefix))
        end
      end
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        throw(ErrorException("If expressions not allowed in backend code"))
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst)  => begin
        DAECallExpressionToJuliaCallExpression(tmpStr, explst, simCode, hashTable, varPrefix=varPrefix)
      end
      DAE.CAST(ty, exp)  => begin
        quote
          $(generateCastExpression(ty, exp, simCode, varPrefix))
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
    _ => throw("Cast $ty: not yet supported in codegen!")
  end
end

"""
      Generates the start conditions.
      All Variables default to zero if not specified.
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
