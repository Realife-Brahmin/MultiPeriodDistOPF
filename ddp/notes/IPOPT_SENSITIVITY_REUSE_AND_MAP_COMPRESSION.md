# Ipopt KKT reuse and backward-map compression

This investigation follows the successful stagewise Ipopt oracle test. It asks
whether the full FilterDDP backward message can be obtained without either
forming all `nx` sensitivity columns or re-solving the nonlinear stage twice
per direction.

## Can ordinary Ipopt.jl reuse Ipopt's final factorization?

Not through its supported public interface. The installed Ipopt.jl wrapper and
Ipopt C interface expose the final primal/dual solution, iteration callbacks,
current iterates and violations, and derivative callbacks. They do not expose
the internal MUMPS factor object or a public operation equivalent to
`solve_last_factorization!(additional_rhs)`.

Reconstructing the Jacobian and Lagrangian Hessian after Ipopt returns and
factorizing their KKT matrix externally is possible, but that is essentially
the operation FilterDDP already performs. It would use Ipopt as the nonlinear
stage optimizer without automatically avoiding the expensive sensitivity
solve.

## The purpose-built Ipopt sensitivity route exists

Ipopt's documented `sIPOPT` component applies the implicit-function theorem to
the KKT conditions and computes parametric NLP sensitivities from a solved NLP.
That is conceptually the exact operation required for FilterDDP's `beta` and
`omega`. The local Ipopt binary artifact already contains:

- `bin/ipopt_sens.exe`;
- `bin/libsipopt-3.dll`;
- the `Sens*.hpp` C++ headers.

However, Ipopt.jl does not wrap this interface. The available integration paths
are therefore:

1. generate an AMPL `.nl` model with the required sIPOPT parameter suffixes and
   call `ipopt_sens.exe`;
2. write a small C++ bridge against `libsipopt` and expose it to Julia;
3. modify or extend Ipopt.jl with a sensitivity wrapper.

The first route is awkward for the current JuMP model because the sIPOPT AMPL
suffixes are not ordinary JuMP constraints or variables. The C++ bridge is the
most controlled research prototype, but it is a real integration task rather
than a one-line Ipopt option.

Official references:

- <https://coin-or.github.io/Ipopt/SPECIALS.html#sipopt-optimal-sensitivity-based-on-ipopt>
- <https://coin-or.github.io/Ipopt/INTERFACES.html>

## Is the backward map low rank?

`analyze_backward_map_compressibility.jl` analyzed the converged IEEE123C,
`T = 3`, stage-1 capture (`nx = 51`). It separately tested:

- the mathematically best truncated-SVD approximation;
- a practical approximation formed from pivot-selected state columns;
- `beta`, `omega`, and the incoming `Vxx` spectrum.

The optimal ranks required to meet relative Frobenius tolerances are:

| Relative error | `beta` rank | `omega` rank | `Vxx` rank |
|---:|---:|---:|---:|
| 10% | 38 of 51 | 7 of 51 | 14 of 51 |
| 5% | 48 of 51 | 14 of 51 | 26 of 51 |
| 1% | 51 of 51 | 33 of 51 | 44 of 51 |

The control map `beta` is therefore not usefully low rank. Even an optimal
rank-40 SVD has 8.8% error. Forty physically realizable selected state columns
produce 10.4% error, and their worst reconstructed individual column remains
42.8% wrong. At rank 48, selected columns still have 5.2% global error and a
26.4% worst-column error.

`omega` is much more compressible globally, but FilterDDP needs the control
map as well. Compressing only `omega` does not remove the dominant `nu x nx`
`beta` and bound-sensitivity work. `Vxx` is moderately compressible at loose
tolerances but nearly full rank at 1%.

## Consequence

The simple strategy “solve only a handful of representative battery
directions and reconstruct all the rest” is rejected on IEEE123. It sacrifices
too much control sensitivity before achieving a meaningful reduction in the
number of right-hand sides. This does not rule out a feeder-aware sparse,
hierarchical, or inexact representation, but there is no evidence for a generic
small global low-rank basis.

The most technically faithful next route is a small sIPOPT bridge on IEEE123.
Its first milestone should request the same three directions already validated
by complete Ipopt re-solves and compare both accuracy and marginal time. Only
if sIPOPT reuses the solved NLP information cheaply should it be extended to
all 51 directions and then tested for memory behavior. This route may improve
constant factors and avoid repeated nonlinear solves; the rank result warns
that it probably cannot remove the fundamental `O(nx)` information content of
the control sensitivity map.

