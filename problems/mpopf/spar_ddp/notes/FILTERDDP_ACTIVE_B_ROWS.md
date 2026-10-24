# Active-row mixed-derivative representation

## Structure found

For large10k, FilterDDP formed `B` as a dense `54665 x 1020` matrix. A direct
stage-1 scan found exactly `1,040,400` nonzeros (`1.866%` density): all 1020
columns were active, but only 1020 rows were active. Those rows were contiguous
(`41284:42303`) and correspond to the battery-power controls through which the
battery-energy state enters the dynamics.

This is not an approximate low-rank claim. In the network formulation,
`lux`, `fux`, and `cux` are identically zero, while sparse `fu` has nonzeros
only in the battery-power columns. Consequently `B` is exactly a dense
`nx x nx` block embedded in an otherwise zero `nu x nx` matrix.

## Guarded exact implementation

The sparse backward path activates the compact representation only when `fu`
is sparse and `lux`, `fux`, and `cux` contain no nonzeros. It extracts the
active control rows from `fu`, forms only the corresponding dense `B_active`
block, writes that block into the same dense KKT RHS workspace, and computes
the value update from the matching packed rows of `beta`. Other models retain
the original dense path.

The active block matched the corresponding dense `B` rows exactly. The compact
`beta' * B` result agreed with the dense product to relative error `9.58e-16`.

## Measured effect

On a warm large10k stage before the production change, dense `B` formation
required about `0.62--0.72 s` and `450.971 MiB`; dense `beta' * B` required
about `0.33--0.36 s` and `7.938 MiB`. The active-block micro-benchmark required
about `0.259 s`/`27.773 MiB` to form the block and `0.041 s`/`16.140 MiB` for
the packed product.

In the production one-iteration profile, warm-stage algebra fell from roughly
`1.20--1.40 s` and `1001.743 MiB` to `0.134 s` and `174.839 MiB`; update time
fell from about `0.54--0.58 s` to `0.19--0.20 s`. The bounded large10k solve
phase fell from `29.987 s` to `26.589 s` (11.33%).

## Full validation

- IEEE123C `T = 3` retained 48 iterations and a byte-for-byte identical trace.
- IEEE2522C `T = 12` retained 79 iterations and a byte-for-byte identical trace;
  practical `1e-6` convergence remains at iteration 76.
- IEEE2522C solve time fell from `577.190 s` to `481.518 s`, a further 16.58%
  single-run reduction. Build time changed from `3.340 s` to `3.540 s`.
- Relative to the original sparse IEEE2522C runtime of `1902.103 s`, the full
  exact optimization sequence now gives a 74.68% single-run reduction.

The full dense `beta` and `omega` policy maps and the 1021-column KKT solve
remain. This optimization exploits exact MPOPF state-to-control incidence; it
does not alter FilterDDP's convergence mechanism.
