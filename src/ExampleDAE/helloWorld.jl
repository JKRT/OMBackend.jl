dae = DAE.DAE_LIST(Cons {
	DAE.Element
}(DAE.COMP("HelloWorld", Cons {
	DAE.Element
}(DAE.VAR(DAE.CREF_IDENT("x", DAE.T_REAL(Nil {
	Any
}()), Nil {
	Any
}()), DAE.VARIABLE(), DAE.BIDIR(), DAE.NON_PARALLEL(), DAE.PUBLIC(), DAE.T_REAL(Nil {
	Any
}()), nothing, Nil {
	Any
}(), DAE.NON_CONNECTOR(), DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\HelloWorld.mo", true, 2, 3, 2, 20, 0.0), Cons {
	Within
}(WITHIN(IDENT("Real")), Nil {
	Any
}()), Prefix.NOCOMPPRE(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}()), SOME {
	DAE.VariableAttributes
}(DAE.VAR_ATTR_REAL(nothing, nothing, nothing, nothing, nothing, SOME {
	DAE.Exp
}(DAE.CAST(DAE.T_REAL(Nil {
	Any
}()), DAE.ICONST(1))), nothing, nothing, nothing, nothing, nothing, nothing, nothing, SOME {
	Bool
}(false), SOME {
	DAE.Exp
}(DAE.SCONST("binding")))), SOME {
	SCode.Comment
}(SCode.COMMENT(nothing, nothing)), NOT_INNER_OUTER()), Cons {
	DAE.Element
}(DAE.VAR(DAE.CREF_IDENT("a", DAE.T_REAL(Nil {
	Any
}()), Nil {
	Any
}()), DAE.PARAM(), DAE.BIDIR(), DAE.NON_PARALLEL(), DAE.PUBLIC(), DAE.T_REAL(Nil {
	Any
}()), SOME {
	DAE.Exp
}(DAE.CAST(DAE.T_REAL(Nil {
	Any
}()), DAE.ICONST(1))), Nil {
	Any
}(), DAE.NON_CONNECTOR(), DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\HelloWorld.mo", true, 3, 3, 3, 23, 0.0), Cons {
	Within
}(WITHIN(IDENT("Real")), Nil {
	Any
}()), Prefix.NOCOMPPRE(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}()), SOME {
	DAE.VariableAttributes
}(DAE.VAR_ATTR_REAL(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, SOME {
	Bool
}(false), nothing)), SOME {
	SCode.Comment
}(SCode.COMMENT(nothing, nothing)), NOT_INNER_OUTER()), Cons {
	DAE.Element
}(DAE.EQUATION(DAE.CALL(IDENT("der"), Cons {
	DAE.Exp
}(DAE.CREF(DAE.CREF_IDENT("x", DAE.T_REAL(Nil {
	Any
}()), Nil {
	Any
}()), DAE.T_REAL(Nil {
	Any
}())), Nil {
	Any
}()), DAE.CALL_ATTR(DAE.T_REAL(Nil {
	Any
}()), false, true, false, false, DAE.DEFAULT_INLINE(), DAE.NO_TAIL())), DAE.BINARY(DAE.UNARY(DAE.UMINUS(DAE.T_REAL(Nil {
	Any
}())), DAE.CREF(DAE.CREF_IDENT("a", DAE.T_REAL(Nil {
	Any
}()), Nil {
	Any
}()), DAE.T_REAL(Nil {
	Any
}()))), DAE.MUL(DAE.T_REAL(Nil {
	Any
}())), DAE.CREF(DAE.CREF_IDENT("x",
	DAE.T_REAL(Nil {
		Any
	}()), Nil {
		Any
	}()), DAE.T_REAL(Nil {
	Any
}()))), DAE.SOURCE(SOURCEINFO("C:\\home\\adrpo33\\dev\\julia\\OMCompiler.jl\\lib\\omc\\HelloWorld.mo", true, 5, 3, 5, 19, 0.0), Cons {
	Within
}(WITHIN(IDENT("HelloWorld")), Nil {
	Any
}()), Prefix.NOCOMPPRE(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}())), Nil {
	Any
}()))), DAE.SOURCE(SOURCEINFO("", false, 0, 0, 0, 0, 0.0), Nil {
	Any
}(), Prefix.NOCOMPPRE(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}(), Nil {
	Any
}()), SOME {
	SCode.Comment
}(SCode.COMMENT(nothing, nothing))), Nil {
	Any
}()))
