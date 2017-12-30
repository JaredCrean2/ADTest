This repo contains a small example of using the complex step method, dual numbers,
and hand-written algorithmic differentiation to compute the jacobian of a
modified Roe flux function with respect to the left and right states.
To simplify the benchmark, the SAT term is replaced by a simple arithmetic
function.

The benchmark evaluates the jacobian at 100,000 different states.

Some timings:
```
Flux evaluation:
  0.003527 seconds
  0.003495 seconds
Dual number {8} evaluation:
  0.020388 seconds
  0.020402 seconds
Complex number evaluation:
  0.284895 seconds
  0.289854 seconds
Hand-written AD evaluation:
  0.012319 seconds
  0.011586 seconds
```

The `Flux evaluation` results time the evaluation of the Roe flux using Float64s.

The `Dual number {8}` results time the evaluation of the flux function using Dual
numbers with 8 dual parts.  8 dual parts are necessary to differentiate with
respect to the left state (the first 4 dual parts) and the right state the
remaining 4 dual parts).

The `Complex number` result time using the complex step method to perturb one
input variable at a time and evaluation the flux function.  As a result, the
flux function is evaluated 800,000 times using Complex128 datatypes.

The `Hand-written AD` result time the evaluation of hand-written AD
code that computes the jacobian directly.  This code is run using Float64 datatypes.

## Summary of files

 * `main.jl`: runs the benchmark, also verifies the hand-written AD code is
consistent with the complex-step results
 * `roe.jl`: contains the Roe flux function, as well as the Euler flux and
 pressure calculation required by the Roe flux
 * `roe_diff.jl`: contains the differentiated versions of the functions
in `roe.jl`
 * `DualMod.jl`: original dual number implementation, provided by Jianfeng Yan
 * `DualMod_static.jl`: dual number implementation modified to use `StaticArrays`.
This removes the need to heap allocate the datatype.
