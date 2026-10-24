# Reusable KKT right-hand-side workspace

## Change

The sparse backward pass previously built a fresh dense block matrix
`[-Qu -B; -c -cx]` at every stage. Besides allocating the final right-hand
side, bracket concatenation created large intermediate copies. A single dense
workspace is now allocated lazily on the first sparse stage of each backward
pass, filled block by block, overwritten by `ldiv!`, copied into the existing
stage rule, and then refilled for the next stage. Compact models never allocate
this network workspace.

The workspace also remains safe across regularization retries because every
block is overwritten before each solve. Explicit KKT capture still copies the
unsolved contents before factorization.

## Validation

- IEEE123C `T = 3` retained the exact strict 48-iteration result, objective
  `2808.924770488853`, and established final residuals.
- On a bounded large10kC `T = 3` probe, warm-stage KKT assembly allocation
  fell from `1760.525 MiB` to `68.569 MiB`, saving `1691.956 MiB` or 96.10%
  per stage. The first stage in a backward sweep allocates the reusable
  755.343-MiB workspace; following stages refill it without reallocating it.
- The full strict IEEE2522C `T = 12` trace remained byte-identical through all
  79 iterations. Runtime fell from `1039.905 s` to `962.181 s`, a further
  single-run reduction of 7.47%. The objective, residuals, and first practical
  `1e-6` milestone at iteration 76 are unchanged.

Relative to the original sparse IEEE2522C `T = 12` runtime of `1902.103 s`,
the complete sequence of exact storage and allocation improvements reaches a
single-run reduction of 49.41%. These timings are evidence of practical
effect, not statistical estimates.
