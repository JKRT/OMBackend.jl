function reconfiguration()
  #=
  Rerun OCC algorithms.
  Find the root variables of the OCC graph + the variables they assign and the root sources.
  That is the reference variables for the roots.
  =#
  @time (rootIndices, variablesToSetIdx, rootSources, variablestoReset) = returnRootIndices(cb.name,
                                                                                            cb,
                                                                                            integrator.u,
                                                                                            tspan,
                                                                                            problem)
  @info "rootsources" rootSources
  @assert length(rootIndices) == length(keys(rootSources)) "Root sources and indices must have the same length. Length was $(length(rootIndices)) == $(length(keys(rootSources)))"
  @info "Root indices" rootIndices
  @info "variablesToSetIdx" variablesToSetIdx
  local newU0::Vector{Float64} = Float64[v for v in integrator.u]
  local stateVars = states(OMBackend.LATEST_REDUCED_SYSTEM)
  #= This is bad, do not use strings this way. =#
  local stateVarsAsStr = [replace(string(s), "(t)" => "") for s in stateVars]
  local OM_NameToMTKIdx = Dict()
  local rootKeys = [k for k in keys(rootSources)]
  local rootValues = [v for v in values(rootSources)]
  local variablesToSet = collect(Iterators.flatten([v for v in values(variablestoReset)]))
  @info "variablesToSet" variablesToSet
  local rootKeysToMTKIdx = indexin(rootKeys, stateVarsAsStr)
  @info "rootKeysToMTKIdx" rootKeysToMTKIdx
  local rootValsToMTKIdx = indexin(rootValues, stateVarsAsStr)
  @info "rootValsToMTKIdx" rootValsToMTKIdx
  local variablesToResetMTKIdx = indexin(variablesToSet, stateVarsAsStr)
  local rootToEquationMap::Dict{String, Symbolics.Equation} = Dict()
  @info "variablesToResetMTKIdx" variablesToResetMTKIdx
  for (i, rk) in enumerate(rootKeys)
    OM_NameToMTKIdx[rk] = rootKeysToMTKIdx[i]
  end
  for (i, rv) in enumerate(rootValues)
    #= This case is true if there is a constant value at the end=#
    if rootValsToMTKIdx[i] != nothing
      OM_NameToMTKIdx[rv] = rootValsToMTKIdx[i]
    else
      OM_NameToMTKIdx[rv] = Meta.parse(rv)
    end
  end
  for (i, vr) in enumerate(variablesToSet)
    OM_NameToMTKIdx[vr] = variablesToResetMTKIdx[i]
  end
  @info "OM_NameToMTKIdx" OM_NameToMTKIdx
  #=
  Start by setting the root variables we got from returnRootIndices.
  Each of these variables have a reference variable.
  =#
  for k in rootKeys
    rootStart = OM_NameToMTKIdx[k]
    rootSource = OM_NameToMTKIdx[rootSources[k]]
    @info "rootStart" rootStart
    @info "rootsource" rootSource
    @info "k" k
    if ! (rootSource isa Float64)
      @info "integrator.u[$(rootStart)] = integrator.u[$(rootSource)]"
      integrator.u[rootStart] = integrator.u[rootSource]
      rootToEquationMap[k] = ~(0, stateVars[rootStart] - stateVars[rootSource])
    else
      @info "integrator.u[$(rootStart)] = $(rootSource)"
      integrator.u[rootStart] = rootSource
      rootToEquationMap[k] = ~(0, stateVars[rootStart] - rootSource)
    end
  end
  #=
  Get the variables of the system
  =#
  local equationDeps = ModelingToolkit.equation_dependencies(OMBackend.LATEST_REDUCED_SYSTEM)
  #= TODO: Not good should strive to fix the internal representations. =#
  local allEquationsAsStr = map((x)-> begin
                                  if !isempty(x)
                                    [replace(string(y), "(t)" => "") for y in x]
                                  end
                                end, equationDeps)
  local rootKeys = [string(k) for k in keys(rootSources)]
  @info "rootKeys" rootKeys
  local equationToAddMap = Dict()
  #=
  Find the equations to replace
  =#
  global ALLEQUATIONSASSTR = allEquationsAsStr
  local oldEquations = equations(OMBackend.LATEST_REDUCED_SYSTEM)
  local assignedRoots = String[]
  for (i, eqVariables) in enumerate(allEquationsAsStr)
    if eqVariables !== nothing && length(eqVariables) == 2
      firstV = first(eqVariables)
      secondV = last(eqVariables)
      check1 = firstV in rootKeys && secondV in rootValues
      check2 = firstV in rootValues && secondV in rootKeys
      check3 = length(ModelingToolkit.difference_vars(oldEquations[i])) == 2
      if (check1 || check2) && check3
        @info "Root already assigned. Root was at $(i)-->" eqVariables
        push!(assignedRoots, last(eqVariables))
      end
    end
  end
  @info "assignedRoots" assignedRoots
  for (i, eqVariables) in enumerate(allEquationsAsStr)
    if eqVariables !== nothing && length(eqVariables) == 2
      local cand = first(eqVariables)
      if cand in rootKeys && ! (cand in assignedRoots)
        @info "Eq" cand
        equationToAddMap[cand] = i
      end
    end
  end
  #= The equation indices are the locations in which we are to insert our new equations=#
  @info "equationToAddMap" equationToAddMap
  #= One assignment need to be changed. =#
  for k in rootKeys
    @info "Value of rootStart" k
    varsToSet = variablestoReset[k]
    val = OM_NameToMTKIdx[k]
    @info "varsToSet" varsToSet
    for v in varsToSet
      vIdx = OM_NameToMTKIdx[v]
      @info "integrator[$(vIdx)] = integrator[$(val)]"
      integrator.u[vIdx] = integrator.u[val]
    end
  end
  push!(solutions, integrator.sol)
  push!(oldSols, (integrator.sol, getSyms(problem), activeModeName))
  newEquations = Symbolics.Equation[]
  @info "equationToAddMap" equationToAddMap
  @info "rootToEquationMap" rootToEquationMap
  for eq in oldEquations
    println(eq)
    push!(newEquations, eq)
  end
  for eq in keys(rootToEquationMap)
    if eq in assignedRoots
      @info "Eq that was in root was" eq
      continue
    end
    local replacementIdx = equationToAddMap[eq] #Hack
    local newEquation = rootToEquationMap[eq]
    @info "Setting newEquation[$(replacementIdx)] = $(eq)"
    newEquations[replacementIdx] = newEquation
    @info  "New equations at $(replacementIdx)" newEquation
  end
  #newEquations[4] = rootToEquationMap["G2"]
  println("New equations:")
  for eq in newEquations
    println(eq)
  end
  @time newSystem = ODESystem(newEquations,
                              independent_variable(OMBackend.LATEST_REDUCED_SYSTEM),
                              states(OMBackend.LATEST_REDUCED_SYSTEM),
                              parameters(OMBackend.LATEST_REDUCED_SYSTEM);
                              name = Symbol(cb.name),
                              discrete_events = ModelingToolkit.discrete_events(OMBackend.LATEST_REDUCED_SYSTEM),
                              )
  @info "New system created with" length(newEquations)
  @info "States" length(states(OMBackend.LATEST_REDUCED_SYSTEM))
  newSystem = OMBackend.CodeGeneration.structural_simplify(newSystem)
  global NEW_SYSTEM = newSystem
  local discrete_events = newSystem.discrete_events
  @info discrete_events
  newU0 = integrator.u
  events = if length(discrete_events) > 0
    RuntimeUtil.evalDiscreteEvents(discrete_events, newU0, i.t, newSystem)
  else
    []
  end
  for e in events
    newU0[e[1]] = e[2]
  end
  @time newProblem = ModelingToolkit.ODEProblem(
    newSystem,
    newU0,
    tspan,
    problem.p,
    #=
    TODO currently only handles a single structural callback.
    =#
    callback = callbackConditions #Should be changed look at method below
  )
  @time integrator = init(newProblem,
                          alg;
                          kwargs...)
  @time reinit!(integrator,
                newU0;
                t0 = i.t,
                reset_dt = true)

end
