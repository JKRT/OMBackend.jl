"
Returns true if simvar is either a algebraic or a state variable
"
function isStateOrAlgebraic(simvar::SimVar)::Bool
  res = @match simvar.varKind begin
    STATE(__) || ALG_VARIABLE(__) => true
    _ => false
  end
end
