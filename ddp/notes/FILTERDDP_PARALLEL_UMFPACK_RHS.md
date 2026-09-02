# Parallel UMFPACK right-hand-side experiment

## Setup

The active-row and cached-Jacobian improvements leave the 1021-column sparse
KKT solve as the dominant large10k stage operation. One exact stage-1 system
(`96968 x 96968`, 1021 RHS) was benchmarked on the lab PC's 10-core/20-thread
Intel i9-10900X. Columns were split across 1, 2, 4, and 8 Julia workers.

The safe benchmark gives each worker an independent UMFPACK factorization.
Every configuration reproduced the same complete solution norm; representative
column residuals were between `5.3e-15` and `4.8e-14`.

## Independent factors

| workers | concurrent factorization | concurrent solve | total |
|---:|---:|---:|---:|
| 1 | 0.494 s | 3.150 s | 3.644 s |
| 2 | 0.950 s | 2.612 s | 3.562 s |
| 4 | 1.757 s | 1.611 s | 3.367 s |
| 8 | 3.385 s | 1.232 s | 4.617 s |

Four workers nearly halve the solve itself, but creating four numeric factors
is effectively serialized or memory-bandwidth limited. The resulting total is
only 7.59% below the one-worker total. Eight workers are slower overall.
Independent factors also increase native factor memory, which Julia object-size
accounting cannot measure reliably.

## Shared factor

A read-only shared-factor probe was kept isolated because UMFPACK does not
promise this usage through the Julia wrapper. Minimal two- and four-worker runs
were numerically correct, but solve times remained `3.153 s` and `3.172 s`:
no faster than serial. This indicates serialized access to the factor or its
solve workspace. A repeated four-worker stress process did not return a result,
so shared-factor concurrency is rejected on both performance and robustness
grounds.

## Decision

No production parallel path is added. A 7.6% isolated linear-algebra gain does
not justify four native factors, extra RAM, mandatory multi-threaded Julia
startup, and a second complex solve path. The result also confirms that merely
splitting RHS columns cannot cure FilterDDP scaling.

The benchmark should be revisited only with a sparse solver designed for
parallel factorization/triangular solves (for example a properly configured
multi-process MUMPS build) or after a mathematical reduction in the number of
sensitivity columns.
