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

## med2522 T=768 failure

The next point did not complete an optimization iteration. MUMPS returned
`INFO(1)=-13` during the first factorization when requesting one additional
contiguous `2147483647`-byte allocation. The model had 8,320,512 variables,
6,001,920 equalities, and 4,446,720 inequalities. This brackets the centralized
memory knee on this machine between the successful `T=576` point and failed
`T=768`; the raw failed stdout is retained under `logs/`.

## Matched periodic large10k T=96

The matched periodic-profile case converged in 90 iterations with objective
`3178495.9248074344`. JuMP solve time was `1682.259 s`, total driver wall time
was `1769.549 s`, and sampled peak working set was `10205.840 MiB`
(`9.966 GiB`). IPOPT's overall algorithm time was `1677.230 s`, including
`1465.386 s` in `ComputeSearchDirection`, `308.779 s` in triangular
back-solves, `1214.339 s` inferred factorization, and `60.720 s` in function
evaluations. The full timing log and machine-readable row are in
`ddp/results/centralized_ipopt_matched_knee/`.

## Matched periodic large10k T=144

The next horizon also converged, in 100 iterations with objective
`3178495.9315857552`. JuMP solve time was `3319.996 s`, total driver wall time
was `3460.443 s`, and sampled peak working set was `16133.801 MiB`
(`15.756 GiB`). IPOPT's overall algorithm time was `3317.672 s`, including
`2979.736 s` in `ComputeSearchDirection`, `520.931 s` in triangular
back-solves, `2540.831 s` inferred factorization, and `98.789 s` in function
evaluations. Memory remains close to linear in horizon, while solve time per
period rises from `17.52 s` at `T=96` to `23.06 s` at `T=144`; continue at
`T=192` to determine whether this is the beginning of the runtime knee.

## Matched periodic large10k T=192

The case converged in 57 iterations with objective `3178495.9383642841`.
JuMP solve time was `9164.712 s` (`2.55 h`), total driver wall time was
`9354.272 s`, and sampled peak working set was `18026.996 MiB`
(`17.604 GiB`). IPOPT's overall algorithm time was `9155.207 s`, including
`8789.629 s` in `ComputeSearchDirection`, `396.417 s` in triangular
back-solves, `8570.559 s` inferred factorization, and `68.406 s` in function
evaluations. Solve time per period more than doubled from `23.06 s` at T=144
to `47.73 s` at T=192 while peak RAM rose by only 11.7%, establishing a clear
runtime knee. One final T=288 point is authorized to determine whether the
curve collapses or merely bends; preserve complete failure logs if it cannot
finish.

## Matched periodic large10k T=288 failure

The final escalation point failed during the first MUMPS factorization.
MUMPS returned `INFO(1)=-13` while requesting an additional `6705 MB`; IPOPT
then exited from restoration after one iteration. The model had 12,771,360
variables, 9,211,104 equalities, and 2,972,160 inequalities. Sampled peak
working set before the failed allocation was `19086.961 MiB` (`18.640 GiB`).
This brackets the matched large10k centralized limit between successful
`T=192` and failed `T=288`: runtime had already bent sharply by `T=192`, and
the next standard horizon encountered the physical-memory/factorization wall.
No larger centralized run is justified on this host without changing the
linear-solver memory behavior or adding RAM.
