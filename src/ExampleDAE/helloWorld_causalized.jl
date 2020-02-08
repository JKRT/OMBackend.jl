const HelloWorld_causalized_DAE = BackendDAE.BACKEND_DAE(
  BackendDAE.EqSystem[
    BackendDAE.EQSYSTEM(
      BackendDAE.VARIABLES(
        BackendDAE.Var[
          #= Var 1 =#
          BackendDAE.VAR(DAE.CREF_IDENT("a", DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), ImmutableList.ListDef.Nil{Any}()), BackendDAE.PARAM(), DAE.BIDIR(), DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), MetaModelica.SOME{DAE.Exp}(DAE.CAST(DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), DAE.ICONST(1))), ImmutableList.ListDef.Nil{Any}(), DAE.SOURCE(MetaModelica.SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\HelloWorld.mo", true, 3, 3, 3, 23, 0.0), ImmutableList.ListDef.Cons{Absyn.Within}(Absyn.WITHIN(Absyn.IDENT("Real")), ImmutableList.ListDef.Nil{Any}()), Prefix.NOCOMPPRE(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}()), MetaModelica.SOME{DAE.VariableAttributes}(DAE.VAR_ATTR_REAL(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, MetaModelica.SOME{Bool}(false), nothing)), nothing, DAE.NON_CONNECTOR(), false),
          #= Var 2 =#
          BackendDAE.VAR(
            #= varName =# DAE.CREF_IDENT("x", DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), ImmutableList.ListDef.Nil{Any}()),
            #= varKind =# BackendDAE.VARIABLE(),
            #= varDirection =# DAE.BIDIR(),
            #= varType =# DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()),
            #= bindExp =# nothing,
            #= arryDim =# ImmutableList.ListDef.Nil{Any}(),
            #= source =# DAE.SOURCE(MetaModelica.SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\HelloWorld.mo", true, 2, 3, 2, 20, 0.0), ImmutableList.ListDef.Cons{Absyn.Within}(Absyn.WITHIN(Absyn.IDENT("Real")), ImmutableList.ListDef.Nil{Any}()), Prefix.NOCOMPPRE(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}()),
            #= values =# MetaModelica.SOME{DAE.VariableAttributes}(DAE.VAR_ATTR_REAL(nothing, nothing, nothing, nothing, nothing, MetaModelica.SOME{DAE.Exp}(DAE.CAST(DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), DAE.ICONST(1))), nothing, nothing, nothing, nothing, nothing, nothing, nothing, MetaModelica.SOME{Bool}(false), MetaModelica.SOME{DAE.Exp}(DAE.SCONST("binding")))),
            #= tearingSelectOption =# nothing,
            #= connectorType =# DAE.NON_CONNECTOR(),
            #= unreplaceable =# false),
            #= Var 2 =#
          BackendDAE.VAR(
            #= varName =# DAE.CREF_IDENT("der(x)", DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), ImmutableList.ListDef.Nil{Any}()),
            #= varKind =# BackendDAE.STATE(1, MetaModelica.SOME{DAE.ComponentRef}(DAE.CREF_IDENT("x", DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), ImmutableList.ListDef.Nil{Any}())), true),
            #= bindExp =# nothing,
            #= arryDim =# ImmutableList.ListDef.Nil{Any}(),
            #= source =# DAE.SOURCE(MetaModelica.SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\HelloWorld.mo", true, 2, 3, 2, 20, 0.0), ImmutableList.ListDef.Cons{Absyn.Within}(Absyn.WITHIN(Absyn.IDENT("Real")), ImmutableList.ListDef.Nil{Any}()), Prefix.NOCOMPPRE(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}()),
            #= values =# MetaModelica.SOME{DAE.VariableAttributes}(DAE.VAR_ATTR_REAL(nothing, nothing, nothing, nothing, nothing, MetaModelica.SOME{DAE.Exp}(DAE.CAST(DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), DAE.ICONST(1))), nothing, nothing, nothing, nothing, nothing, nothing, nothing, MetaModelica.SOME{Bool}(false), MetaModelica.SOME{DAE.Exp}(DAE.SCONST("binding")))),
            #= tearingSelectOption =# nothing,
            #= connectorType =# DAE.NON_CONNECTOR(),
            #= unreplaceable =# false)
          ]),
      BackendDAE.Equation[BackendDAE.EQUATION(DAE.CALL(Absyn.IDENT("der"), ImmutableList.ListDef.Cons{DAE.Exp}(DAE.CREF(DAE.CREF_IDENT("x", DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), ImmutableList.ListDef.Nil{Any}()), DAE.T_REAL(ImmutableList.ListDef.Nil{Any}())), ImmutableList.ListDef.Nil{Any}()), DAE.CALL_ATTR(DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), false, true, false, false, DAE.DEFAULT_INLINE(), DAE.NO_TAIL())), DAE.BINARY(DAE.UNARY(DAE.UMINUS(DAE.T_REAL(ImmutableList.ListDef.Nil{Any}())), DAE.CREF(DAE.CREF_IDENT("a", DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), ImmutableList.ListDef.Nil{Any}()), DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()))), DAE.MUL(DAE.T_REAL(ImmutableList.ListDef.Nil{Any}())), DAE.CREF(DAE.CREF_IDENT("x", DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()), ImmutableList.ListDef.Nil{Any}()), DAE.T_REAL(ImmutableList.ListDef.Nil{Any}()))), DAE.SOURCE(MetaModelica.SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\HelloWorld.mo", true, 5, 3, 5, 19, 0.0), ImmutableList.ListDef.Cons{Absyn.Within}(Absyn.WITHIN(Absyn.IDENT("HelloWorld")), ImmutableList.ListDef.Nil{Any}()), Prefix.NOCOMPPRE(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}(), ImmutableList.ListDef.Nil{Any}()), BackendDAE.EQUATION_ATTRIBUTES(false, BackendDAE.UNKNOWN_EQUATION_KIND(), BackendDAE.EVALUATION_STAGES(false, false, false, false)))], nothing, nothing, nothing, BackendDAE.NO_MATCHING(), ImmutableList.ListDef.Nil{Any}(), BackendDAE.UNKNOWN_PARTITION(), Any[])], BackendDAE.SHARED_DUMMY())
