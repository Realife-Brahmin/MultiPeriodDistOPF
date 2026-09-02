# Typed constraint residual vector

## Diagnosis

The post-optimization large10k profile split the remaining value update into
individual products. The value-Hessian recursion used only `31.751 MiB` and
about `0.36 s` at a warm stage. Nearly the entire remaining update cost came
from `omega' * c`: about `1975 MiB` and `4.0 s` per stage.

`omega` was a `Matrix{Float64}`, but `c` was a `Vector{Any}` because
`Vector(constraints.c(x, u))` preserved the dynamically generated callback's
abstract element type. This prevented BLAS dispatch and forced a generic
matrix-vector multiplication over the full dense equality-sensitivity map.

## Change

Construct `c` as `Vector{T}`, where `T` is the solver's numeric type. For the
network experiments this converts an already Float64-valued residual vector
from abstract storage to `Vector{Float64}`. No residual value or equation is
changed.

## Validation

- A full strict IEEE123C `T = 3` regression retained 48 iterations and
  objective `2808.924770488853`. Final residual changes were at approximately
  `1e-18`, consistent with BLAS summation order.
- On a bounded large10kC `T = 3` probe, warm-stage update allocation fell from
  `2011.138 MiB` to `33.909 MiB` (98.31%). The `omega' * c` product fell from
  about `3.98 s` and `1975.262 MiB` to `0.0101 s` and `0.008 MiB`.
- A full strict IEEE2522C `T = 12` run retained the same displayed 79-row
  iteration trace. Runtime fell from `962.181 s` to `737.043 s`, a further
  23.40% reduction. The objective changed by `2e-12`, while final residuals
  agreed to reported precision and practical `1e-6` convergence remained at
  iteration 76.

Relative to the original sparse IEEE2522C `T = 12` runtime of `1902.103 s`,
the full sequence of exact storage and dispatch improvements gives a 61.25%
single-run reduction. Timing percentages are not statistical estimates.
