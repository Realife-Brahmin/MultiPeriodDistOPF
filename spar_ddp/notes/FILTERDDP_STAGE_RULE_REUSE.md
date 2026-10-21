# Reusing sparse FilterDDP stage-rule storage

## What was duplicated

After solving a sparse stage KKT system, the implementation copied its dense
solution into temporary `alpha_beta` and `psi_omega` blocks, copied their
sensitivity columns again into `beta` and `omega`, and then replaced the
stage's already allocated update rule. The previous rule for that stage is not
used again before the next forward pass.

The sparse path now copies the solved blocks directly into the existing stage
rule. Its vectors and sensitivity matrices are updated in place. The dense and
unconstrained paths retain their original behavior.

## Validation

- IEEE123C `T = 3` retained the exact strict 48-iteration result: objective
  `2808.924770488853`, primal infeasibility `6.227371345322e-8`, dual
  infeasibility `2.303870314635e-8`, and complementarity
  `1.624842969030e-9`.
- On a bounded large10kC `T = 3` probe, warm-stage update allocation fell from
  `3525.163 MiB` to `2011.139 MiB`, saving `1514.024 MiB` or 42.95% per stage
  construction. Persistent update-rule storage remains `757.012 MiB`; this
  change removes temporary copying rather than changing retained data.
- A full strict IEEE2522C `T = 12` cold-start run retained the exact
  79-iteration trace, objective, and residuals. Runtime was `1039.905 s`,
  versus `1092.507 s` with only in-place KKT RHS solving: a further 4.82%
  reduction. All practical `1e-6` criteria were first met at iteration 76.

The runtime percentages are single-run measurements. The exact claims are
trajectory preservation and elimination of the redundant array copies.
