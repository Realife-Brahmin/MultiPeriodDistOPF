# FilterDDP per-iteration timing breakdown

## Question

The earlier matrix inventory established dimensions and sparsity, but not how
much time each operation consumes. This experiment times every stage of every
backward sweep and separately times the complete forward rollout/line search.
It covers both requested axes: IEEE2522C at `T=3` and `T=12`, and large10kC at
`T=3`.

## Accounting

One accepted FilterDDP iteration normally has one complete `T -> 1` backward
sweep and one forward rollout. A barrier reduction changes the global
complementarity target and triggers another complete backward sweep without a
forward rollout or incrementing the displayed iteration. Final convergence is
also certified by a last backward sweep. Therefore:

```
backward sweeps = accepted iterations + barrier-update sweeps + final sweep.
```

All three strict `1e-7` runs made 10 monotone barrier reductions, from `mu=1`
through the same schedule to the `1e-8` floor. IEEE2522 `T=3` used 56 accepted
iterations and 67 sweeps; IEEE2522 `T=12` used 79 and 90; large10k `T=3` used
115 and 126. Thus the usual iteration count hides 11 full backward sweeps in
each run. The barrier-only sweeps consumed 12.678 s, 49.523 s, and 272.550 s.

## Measured totals

The compact summary is
`ddp/results/network_filterddp/iteration_timing_summary.csv`; the full compact
per-sweep records retain iteration, barrier iteration, actual `mu`, outcome,
step size, backtracks, and all timing categories. Percentages below use
measured backward plus forward time rather than process setup time.

| Case | RHS solve | Factorization | KKT assembly | Derivatives | Forward rollout |
|---|---:|---:|---:|---:|---:|
| IEEE2522 `T=3` | 34.8 s (34.9%) | 26.4 s (26.5%) | 11.1 s (11.1%) | 7.5 s (7.5%) | 6.2 s (6.2%) |
| IEEE2522 `T=12` | 202.0 s (42.0%) | 155.2 s (32.3%) | 24.3 s (5.0%) | 38.7 s (8.1%) | 23.4 s (4.9%) |
| large10k `T=3` | 1140.6 s (44.4%) | 851.4 s (33.1%) | 97.1 s (3.8%) | 134.3 s (5.2%) | 170.3 s (6.6%) |

`derivative_s` contains the separately logged first- and second-order callback
times, so those subcategories must not be added again. Residual unclassified
backward overhead is 1.2--6.6%.

## Interpretation

The dominant cost is not translating equations into matrices. KKT assembly is
only 3.8--11.1% of measured algorithm time. Nor is exact derivative evaluation
dominant: it is 5.2--8.1%. Sparse factorization plus the many-right-hand-side
sensitivity solve consumes 61.4% on IEEE2522 `T=3`, 74.3% on IEEE2522 `T=12`,
and 77.5% on large10k `T=3`.

Both increasing the horizon on one feeder and increasing feeder/state size
shift a larger fraction of runtime into linear algebra. The RHS solve is the
largest category in every case. This supports the existing diagnosis:
FilterDDP is slow mainly because every stage propagates a full state-sensitivity
map, requiring a factorization and an `nx+1`-column solve, not because Hessian
callbacks are intrinsically expensive.

The large10k barrier sweeps also become more expensive as `mu` falls: the first
barrier-update sweep took 19.35 s and the last 38.12 s. This is consistent with
the measured deterioration in KKT conditioning. It is a secondary multiplier
on the repeated sensitivity workload, not the primary scaling source.

## Validation and caveats

All three convergence traces are byte-for-byte identical to their prior
optimized timing-matrix traces. Instrumented solve times were 102.767 s,
484.169 s, and 2572.413 s. IEEE2522 `T=3` is 6.3% slower than its prior
single-run 96.664-s baseline, so these totals are used for within-run shares,
not as replacement performance benchmarks. No matrix payloads were serialized.

Instrumentation is supplied by
`ddp/patches/iteration_timing_diagnostic.patch`, applied after the production
patch stack. `scripts/extract_filterddp_timing.ps1` converts stdout to one CSV
row per complete backward sweep.
