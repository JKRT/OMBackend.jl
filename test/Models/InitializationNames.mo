package InitializationNames "The initialization's parameters and pins matched by the variables' names"

  model FreeArrayParameter "free parameters that are array elements (printed var\"k[1]\"): unknowns of the initialization"
    parameter Real k[2](each fixed = false, each start = 1);
    Real x[2](start = {1, 2}, each fixed = true);
    Real y[2];
  equation
    der(x) = -k .* x;
    y = k .* x;
  initial equation
    y = {3, 8};
  end FreeArrayParameter;

end InitializationNames;
