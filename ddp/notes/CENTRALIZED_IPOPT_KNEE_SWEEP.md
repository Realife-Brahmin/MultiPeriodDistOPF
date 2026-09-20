# Centralized IPOPT horizon knee sweep

This continuation extends the centralized BFM-NL IPOPT timing study beyond the
previous horizon limits. Every profiled case preserves IPOPT's full timing
statistics and records sampled peak working set, dimensions, nonzeros,
iterations, objective, validation status, and raw logs.

## med2522 T=384

The case converged in 95 iterations, passed the independent constraint
validator, and reached objective `8765.79217322525`. The uninstrumented control
required `1235.065 s` of JuMP solve time; the profiled rerun required
`1241.290 s`, only `0.50%` more. Timing statistics are therefore not a material
perturbation at this scale. The sampled control-run peak was `12288.523 MiB`
(`12.0005 GiB`).

In the profiled run, IPOPT's overall algorithm time was `1239.475 s`.
`ComputeSearchDirection` consumed `1099.842 s` (88.73%),
`StdAugSystemSolverMultiSolve` consumed `1067.848 s`, and
`LinearSystemBackSolve` consumed `355.223 s` (28.66%). Because the bundled
MUMPS interface reports zero for `LinearSystemFactorization`, the established
paper convention estimates factorization as augmented-system time minus
back-solve time: `712.625 s` (57.49%). Function evaluations consumed
`62.807 s` (5.07%). The category arithmetic is consistent up to overlapping
IPOPT timers and rounding.

Machine-readable rows are in
`ddp/results/centralized_ipopt/centralized_ipopt_knee.csv`; the uninstrumented
control is retained separately. The next med2522 point should be `T=576`, after
which the per-period time and RAM curve can be compared against `T=288/384/576`.

## med2522 T=576

The case converged and validated in 100 iterations with objective
`8774.12264117318`. JuMP solve time is `1963.617 s` and sampled peak working set
is `18590.453 MiB` (`18.155 GiB`). IPOPT reports `1960.352 s` overall,
`1742.942 s` in `ComputeSearchDirection`, `565.345 s` in triangular back-solves,
`1126.196 s` inferred factorization, and `95.912 s` in function evaluations.

There is no clear knee through this point. Solve time per period is `3.36 s` at
`T=144`, `3.23 s` at `T=384`, and `3.41 s` at `T=576`; memory per period is
also essentially linear between the last two points. A `T=768` run is therefore
justified and is projected at roughly 45 minutes and 24--25 GiB peak.
