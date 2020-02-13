const vanDerPol_DAE = DAE.DAE_LIST(Cons{DAE.Element}(DAE.COMP("VanDerPol", Cons{DAE.Element}(
  DAE.VAR(
    DAE.CREF_IDENT("x", DAE.T_REAL(Nil{Any}()), Nil{Any}()),
    DAE.VARIABLE(),
    DAE.BIDIR(),
    DAE.NON_PARALLEL(),
    DAE.PUBLIC(),
    DAE.T_REAL(Nil{Any}()),
    nothing,
    Nil{Any}(),
    DAE.NON_CONNECTOR(),
    DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\VanDerPol.mo", true, 4, 3, 4, 34, 0.0), Cons{Within}(WITHIN(IDENT("Real")), Nil{Any}()), Prefix.NOCOMPPRE(), Nil{Any}(), Nil{Any}(), Nil{Any}(), Nil{Any}()),
    SOME{DAE.VariableAttributes}(
      DAE.VAR_ATTR_REAL(
        nothing, #= quantity =#
        nothing, #= unit =#
        nothing, #= displayUnit =#
        nothing, #= min =#
        nothing, #= max =#
        SOME{DAE.Exp}(DAE.CAST(DAE.T_REAL(Nil{Any}()), DAE.ICONST(1))), #= start value =#
        SOME{DAE.Exp}(DAE.BCONST(true)),
        nothing, nothing, nothing, nothing, nothing, nothing, SOME{Bool}(false), SOME{DAE.Exp}(DAE.SCONST("binding")))), SOME{SCode.Comment}(SCode.COMMENT(nothing, nothing)), NOT_INNER_OUTER()), Cons{DAE.Element}(DAE.VAR(DAE.CREF_IDENT("y", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.VARIABLE(), DAE.BIDIR(), DAE.NON_PARALLEL(), DAE.PUBLIC(), DAE.T_REAL(Nil{Any}()), nothing, Nil{Any}(), DAE.NON_CONNECTOR(), DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\VanDerPol.mo", true, 5, 3, 5, 34, 0.0), Cons{Within}(WITHIN(IDENT("Real")), Nil{Any}()), Prefix.NOCOMPPRE(), Nil{Any}(), Nil{Any}(), Nil{Any}(), Nil{Any}()), SOME{DAE.VariableAttributes}(DAE.VAR_ATTR_REAL(nothing, nothing, nothing, nothing, nothing, SOME{DAE.Exp}(DAE.CAST(DAE.T_REAL(Nil{Any}()), DAE.ICONST(1))), SOME{DAE.Exp}(DAE.BCONST(true)), nothing, nothing, nothing, nothing, nothing, nothing, SOME{Bool}(false), SOME{DAE.Exp}(DAE.SCONST("binding")))), SOME{SCode.Comment}(SCode.COMMENT(nothing, nothing)), NOT_INNER_OUTER()), Cons{DAE.Element}(DAE.VAR(DAE.CREF_IDENT("lambda", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.PARAM(), DAE.BIDIR(), DAE.NON_PARALLEL(), DAE.PUBLIC(), DAE.T_REAL(Nil{Any}()), SOME{DAE.Exp}(DAE.RCONST(0.3)), Nil{Any}(), DAE.NON_CONNECTOR(), DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\VanDerPol.mo", true, 6, 3, 6, 30, 0.0), Cons{Within}(WITHIN(IDENT("Real")), Nil{Any}()), Prefix.NOCOMPPRE(), Nil{Any}(), Nil{Any}(), Nil{Any}(), Nil{Any}()), SOME{DAE.VariableAttributes}(DAE.VAR_ATTR_REAL(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, SOME{Bool}(false), nothing)), SOME{SCode.Comment}(SCode.COMMENT(nothing, nothing)), NOT_INNER_OUTER()), Cons{DAE.Element}(DAE.EQUATION(DAE.CALL(IDENT("der"), Cons{DAE.Exp}(DAE.CREF(DAE.CREF_IDENT("x", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.T_REAL(Nil{Any}())), Nil{Any}()), DAE.CALL_ATTR(DAE.T_REAL(Nil{Any}()), false, true, false, false, DAE.DEFAULT_INLINE(), DAE.NO_TAIL())), DAE.CREF(DAE.CREF_IDENT("y", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.T_REAL(Nil{Any}())), DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\VanDerPol.mo", true, 8, 3, 8, 13, 0.0), Cons{Within}(WITHIN(IDENT("VanDerPol")), Nil{Any}()), Prefix.NOCOMPPRE(), Nil{Any}(), Nil{Any}(), Nil{Any}(), Nil{Any}())),

  Cons{DAE.Element}(DAE.EQUATION(
    DAE.CALL(
      IDENT("der"),
      Cons{DAE.Exp}(DAE.CREF(DAE.CREF_IDENT("y", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.T_REAL(Nil{Any}())), Nil{Any}()),
      DAE.CALL_ATTR(DAE.T_REAL(Nil{Any}()), false, true, false, false, DAE.DEFAULT_INLINE(), DAE.NO_TAIL())),
    DAE.BINARY(
      DAE.BINARY(
        DAE.CREF(DAE.CREF_IDENT("lambda", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.T_REAL(Nil{Any}())),
        DAE.MUL(DAE.T_REAL(Nil{Any}())),
        DAE.BINARY(
          DAE.BINARY(
            DAE.RCONST(1.0),
            DAE.SUB(DAE.T_REAL(Nil{Any}())),
            DAE.BINARY(DAE.CREF(DAE.CREF_IDENT("x", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.T_REAL(Nil{Any}())), DAE.POW(DAE.T_REAL(Nil{Any}())), DAE.RCONST(2.0))),
          DAE.MUL(DAE.T_REAL(Nil{Any}())),
          DAE.CREF(DAE.CREF_IDENT("y", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.T_REAL(Nil{Any}()))
        )
      ),
      DAE.SUB(DAE.T_REAL(Nil{Any}())),
      DAE.CREF(DAE.CREF_IDENT("x", DAE.T_REAL(Nil{Any}()), Nil{Any}()), DAE.T_REAL(Nil{Any}()))
    ),
    DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\VanDerPol.mo", true, 9, 3, 9, 36, 0.0), Cons{Within}(WITHIN(IDENT("VanDerPol")), Nil{Any}()), Prefix.NOCOMPPRE(), Nil{Any}(), Nil{Any}(), Nil{Any}(), Nil{Any}())), Nil{Any}()))))), DAE.SOURCE(SOURCEINFO("", false, 0, 0, 0, 0, 0.0), Nil{Any}(), Prefix.NOCOMPPRE(), Nil{Any}(), Nil{Any}(), Nil{Any}(), Nil{Any}()), SOME{SCode.Comment}(SCode.COMMENT(nothing, SOME{String}("Van der Pol oscillator model")))), Nil{Any}()))
