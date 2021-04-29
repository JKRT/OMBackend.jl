using ModelingToolkit
#=
# Define some variables
@parameters t L g
@variables x(t) y(t) T(t)
D = Differential(t)
eqs2 = [0 ~ - D(D(x)) + T*x,
        D(D(y)) ~ T*y - g,
        0 ~ x^2 + y^2 - L^2]
pendulum2 = ODESystem(eqs2, t, [x, y, T], [L, g], name=:pendulum)
# Turn into a first order differential equation system
first_order_sys = ModelingToolkit.ode_order_lowering(pendulum2)
# Perform index reduction to get an Index 1 DAE
new_sys = dae_index_lowering(first_order_sys)
u0 = [
  D(x)    => 0.0,
  D(y)    => 0.0,
  x       => 1.0,
  y       => 0.0,
  T       => 0.0
]
p = [
    L => 1.0,
    g => 9.8
]
prob_auto = ODEProblem(new_sys,u0,(0.0,10.0),p)
sol = solve(prob_auto, Rodas5());
=#
"""
  Generates simulation code targetting modeling toolkit
"""
function generateMDKCode(simCode::SimulationCode.SIM_CODE)
  local crefToSimVarHT = simCode.crefToSimVarHT
  local equations::Array = []
  local exp::DAE.Exp
  local modelName::String = simCode.name
  local parameters::Array = []
  local stateDerivatives::Array = []
  local stateMarkings::Array = []
  local stateVariables::Array = []
  local algebraicVariables::Array = []
  for varName in keys(crefToSimVarHT)
    ixAndVar = crefToSimVarHT[varName]
    local varType = ixAndVar[2].varKind
    @match varType  begin
      SimulationCode.INPUT(__) => @error "INPUT not supported in CodeGen"
      SimulationCode.STATE(__) => push!(stateVariables, varName)
      SimulationCode.PARAMETER(__) => push!(parameters, varName)
      SimulationCode.ALG_VARIABLE(__) => push!(algebraicVariables #= Should be alg if more specialised scheme is used.=#, varName)
      #=TODO: Do I need to modify this?=#
      SimulationCode.STATE_DERIVATIVE(__) => push!(stateDerivatives, varName)
    end
  end
  equations = simCode.residualEquations
  equations = createResidualEquationsMDK(stateVariables, equations, simCode::SimulationCode.SIM_CODE)
  #= Symbolic names =#
  local stateVariablesSym = [:($(Symbol(v))(t)) for v in stateVariables]
  local algebraicVariablesSym = [Symbol(v) for v in algebraicVariables]
  local parVariablesSym = [Symbol(p) for p in parameters]
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquationsMDK(stateVariables, algebraicVariables, simCode)
  local PARAMETER_EQUATIONS = createParameterEquationsMDK(parameters, simCode)
  #= Formulate the problem as a DAE Problem=#
  program = quote
    using ModelingToolkit
    using DiffEqBase
    using DifferentialEquations    
    function $(Symbol("$(modelName)Simulate"))(tspan = (0.0, 1.0))
        pars = @parameters begin
          $(parVariablesSym...)
          t
        end
        vars = @variables begin
          $(stateVariablesSym...)
          $(algebraicVariablesSym...)
        end
      der = Differential(t)
      eqs = [$(equations...)]
      nonLinearSystem = ODESystem(eqs, t, vars, pars, name=:($(Symbol($modelName))))
      pars = Dict($(PARAMETER_EQUATIONS...), t => tspan[1])
      initialValues = [$(START_CONDTIONS_EQUATIONS...)]
      firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
      reducedSystem = dae_index_lowering(firstOrderSystem)
      problem = ODEProblem(reducedSystem, initialValues, tspan, pars)
      solve(problem)
    end
  end
  return (modelName, program)
end

"
 Creates the residual equations in unsorted order
"
function createResidualEquationsMDK(stateVariables, equations::Array, simCode::SimulationCode.SIM_CODE)::Array{Expr}
  local eqs::Array{Expr} = []
  local ht = simCode.crefToSimVarHT
  local lhsVariables = Set([])
  for i in 1:length(stateVariables)
    local variableIdx = ht[stateVariables[i]][1]
    local equationIdx = simCode.matchOrder[variableIdx]
    local equation = equations[equationIdx]
    push!(eqs, residualEqtoJuliaMDK(equation, simCode, equationIdx))
  end
  return eqs
end

"Converts a DAE expression into a Julia expression"
function expToJuliaExpMDK(exp::DAE.Exp, simCode::SimulationCode.SIM_CODE, varSuffix=""; varPrefix="x")::Expr
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
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.PARAMETER(__) => quote
              $(LineNumberNode(@__LINE__, "$varName parameter"))
              $(Symbol(indexAndVar[2].name))
            end
            SimulationCode.ALG_VARIABLE(__) => quote
              $(LineNumberNode(@__LINE__, "$varName, algebraic"))
              $(Symbol(indexAndVar[2].name))
            end
            _ => begin
              @error "Unsupported varKind: $(varKind)"
              throw()
            end
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
          $(o)($(expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix)))
        end
      end
      DAE.BINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        a = expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix)
        b = expToJuliaExpMDK(e2, simCode, varPrefix=varPrefix)
        o = DAE_OP_toJuliaOperator(op)
        quote
          $o($(a), $(b))
        end
      end
      DAE.LUNARY(operator = op, exp = e1)  => begin
        quote
          $("(" + BDAE.string(op) + " " + expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix) + ")")
        end
      end
      DAE.LBINARY(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExpMDK(e1, simCode, varPrefix=varPrefix) + " " + BDAE.string(op) + " " + expToJuliaExpMDK(e2, simCode, varPrefix=varPrefix))
        end
      end
      DAE.RELATION(exp1 = e1, operator = op, exp2 = e2) => begin
        quote 
          $(expToJuliaExpMDK(e1, simCode,
                          varPrefix=varPrefix) + " "
            + BDAE.string(op) + " " + expToJuliaExpMDK(e2, simCode,varPrefix=varPrefix))
        end
      end
      DAE.IFEXP(expCond = e1, expThen = e2, expElse = e3) => begin
        throw(ErrorException("If expressions not allowed in backend code"))
      end
      DAE.CALL(path = Absyn.IDENT(tmpStr), expLst = explst)  => begin
        DAECallExpressionToMDKCallExpression(tmpStr, explst, simCode, hashTable, varPrefix=varPrefix)
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


"""
  Transforms a given equation into MDK Julia code
"""
function residualEqtoJuliaMDK(eq::BDAE.Equation, simCode::SimulationCode.SIM_CODE, equationIdx::Int64)::Expr
  local rhs::DAE.Exp
  local ht = SimulationCode.makeIndexVarNameDict(simCode.matchOrder, simCode.crefToSimVarHT)
  local result::Expr = @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
      #= Using Reduce to rewrite equations...=#
      variables = Set([])
      local variableIdx = MetaGraphs.get_prop(simCode.equationGraph, equationIdx, :vID)
      local varToSolve = simCode.crefToSimVarHT[ht[variableIdx]][2]
      local varToSolveExpr::Symbol = SimulationCode.isState(varToSolve) ? Symbol("der($(varToSolve.name))") : Symbol(varToSolve.name)
      eq = quote
        0 = $(expToJuliaExpMDK(rhs, simCode))
      end
      strippedEq = stripBeginBlocks(eq)
      reducSol = Reduce.Algebra.solve(strippedEq, varToSolveExpr)[1]
      #= Replace equality with ~ =#
      quote
        $(reducSol.args[1]) ~ $(reducSol.args[2])
      end
    end
    _ => begin
      throw("traversalError for $eq")
    end
  end
  return result
end


"
    Generates the initial value for the equations
    TODO: Currently unable to generate start condition in order 
  "
function createStartConditionsEquationsMDK(states::Array,
                                        algebraics::Array,
                                        simCode::SimulationCode.SimCode)::Array{Expr}
  #  local startEquations = createSortedEquations([algVariables..., stateVariables...], simCode; arrayName = "reals")
  local algInit = getStartConditionsMDK(algebraics, simCode)
  local stateInit = getStartConditionsMDK(states, simCode)
  return vcat(algInit, stateInit)
end



function getStartConditionsMDK(vars::Array, simCode::SimulationCode.SimCode)::Array{Expr}
  local startExprs::Array{Expr} = []
  local residuals = simCode.residualEquations
  local ht::Dict = simCode.crefToSimVarHT
  if length(vars) == 0
    return []
  end
  for var in vars
    (index, simVar) = ht[var]
    varName = simVar.name
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
            try
              local exp = eval(expToJuliaExpMDK(start, simCode))
              push!(startExprs, :($(Symbol("$varName")) => $(expToJuliaExpMDK(start, simCode))))
            catch
              push!(startExprs, :($(Symbol("$varName")) => pars[$(expToJuliaExpMDK(start, simCode))]))
            end
            ()
          end
          (NONE(), SOME(fixed)) => begin
            push!(startExprs, :($(varName) => 0.0))
            ()
          end
          (_, _) => ()
        end
        NONE() => ()
      end
    end
  end
  return startExprs
end


function createParameterEquationsMDK(parameters::Array, simCode::SimulationCode.SimCode)::Array{Expr}
  local parameterEquations::Array = []
  local hT = simCode.crefToSimVarHT
  for param in parameters
    (index, simVar) = hT[param]
    local simVarType::SimulationCode.SimVarType = simVar.varKind
    bindExp = @match simVarType begin
      SimulationCode.PARAMETER(bindExp = SOME(exp)) => exp
      _ => ErrorException("Unknown SimulationCode.SimVarType for parameter.")
    end
    #= Solution for https://github.com/SciML/ModelingToolkit.jl/issues/991 =#
    push!(parameterEquations,
          quote
          $(LineNumberNode(@__LINE__, "$param"))
          $(Symbol(simVar.name)) => $(expToJuliaExpMDK(bindExp, simCode))
          end
          )
  end
  return parameterEquations
end
