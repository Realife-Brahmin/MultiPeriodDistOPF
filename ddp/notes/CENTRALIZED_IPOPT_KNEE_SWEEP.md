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
