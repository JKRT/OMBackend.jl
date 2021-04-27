import ModelingToolkit
#=
# Define some variables
@parameters t L g
@variables x(t) y(t) T(t)
D = Differential(t)

eqs2 = [D(D(x)) ~ T*x,
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

function generateMDKSimulationCode(simCode::SimulationCode.UNSORTED_SIM_CODE)
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
    #=I assume I can put state and alg in the same pile for modelling toolkit for now=#
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
  equations = createResidualEquations(stateVariables, equations, simCode::SimulationCode.UNSORTED_SIM_CODE)

  #= Symbolic names =#
  local stateVariablesSym = [:($(Symbol(v))(t)) for v in stateVariables]
  local algebraicVariablesSym = [Symbol(v) for v in algebraicVariables]
  local parLen = length(parameters)
  local START_CONDTIONS_EQUATIONS = createStartConditionsEquations(stateVariables, algebraicVariables, simCode)
  local PARAMETER_EQUATIONS = createParameterEquations(parameters, simCode)
  #= Formulate the problem as a DAE Problem=#
  res = quote
    using ModelingToolkit
    variables = @variables begin
      $(stateVariablesSym...)
      $(algebraicVariablesSym...)
    end
    @parameters begin
      p[1:$(parLen)]
      t
    end
    $(PARAMETER_EQUATIONS...)
    DER = Differential(t)
    initialValues = $(START_CONDTIONS_EQUATIONS)
    eqs = [$(equations...)]
    nonLinearSystem = ODESystem(eqs, t, variables, p, name=Symbol($modelName))
    firstOrderSystem = ModelingToolkit.ode_order_lowering(nonLinearSystem)
    reducedSystem = dae_index_lowering(firstOrderSystem)
    #= Get initial values and values of parameters=#

    #= End =#
    prob_auto = ODEProblem(reducedSystem, initialValues, (0.0,10.0), p)
  end
  return (modelName, res)
end

"
 Creates the residual equations in unsorted order
"
function createResidualEquations(stateVariables, equations::Array, simCode::SimulationCode.UNSORTED_SIM_CODE)::Array{Expr}
  local eqs::Array{Expr} = []
  local ht = simCode.crefToSimVarHT
  for i in 1:length(stateVariables)
    local variableIdx = ht[stateVariables[i]][1]
    local equation = equations[i]
    push!(eqs, residualEqtoJulia(equation, simCode, i))
  end
  return eqs
end

"Converts a DAE expression into a Julia expression"
function expToJuliaExp(exp::DAE.Exp, simCode::SimulationCode.UNSORTED_SIM_CODE, varSuffix=""; varPrefix="x")::Expr
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
              p[$(indexAndVar[1])]
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


"""
 Transforms a given equation into MDK Julia code
"""
function residualEqtoJulia(eq::BDAE.Equation, simCode::SimulationCode.UNSORTED_SIM_CODE, resNumber)::Expr
  local rhs::DAE.Exp
  local result::Expr = @match eq begin
    BDAE.RESIDUAL_EQUATION(exp = rhs) => begin
       quote
         0 ~ $(expToJuliaExp(rhs, simCode; varPrefix="x"))
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
function createStartConditionsEquations(states::Array,
                                        algebraics::Array,
                                        simCode::SimulationCode.SimCode)::Expr
#  local startEquations = createSortedEquations([algVariables..., stateVariables...], simCode; arrayName = "reals")
  return quote
#    $(getStartConditions(states, "x", simCode))
  end
end
