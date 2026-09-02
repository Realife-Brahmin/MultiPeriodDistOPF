# Optimized FilterDDP timing matrix

Date: 2026-09-02

## Purpose

The original horizon sweep established that sparse FilterDDP could solve the
network MPOPF cases, but it predated the exact memory and allocation fixes
developed during the bottleneck investigation. This sweep repeats the same
systems and horizons with the complete validated patch stack and the
factor-backed policy representation enabled. The blocked-RHS minimum-memory
mode was disabled because it is slower.

All cases were cold-started, run sequentially on the same lab PC, and used the
strict `1e-7` stopping test. The solve time excludes model construction. Peak
working set is the largest process working set sampled by the queue runner; it
is not a component-wise memory attribution. The IEEE123 `T=3` peak is missing
because the first runner process completed before its bookkeeping fallback was
fixed.

## Result

The complete machine-readable comparison is
`ddp/results/network_filterddp/optimized_timing_comparison.csv`. The largest
network results are:

- IEEE2522C improves by 58--75% across all horizons, including 1902.103 s to
  476.400 s at `T=12` and 9359.596 s to 3940.608 s at `T=96`.
- large10k improves from 7123.782 s to 2633.856 s at `T=3` (63.03%) and from
  23490.218 s to 14809.508 s at `T=6` (36.95%).
- The optimized large10k peak working sets are 3808 MiB (`T=3`) and 4216 MiB
  (`T=6`). This confirms the practical effect of replacing the original
  horizon-persistent dense policy maps with compact factor-backed actions.

The full IEEE2522 `T=12`, `T=48`, and `T=96` traces and the large10k `T=6`
trace are byte-identical to their stored pre-sweep references. large10k `T=3`
has the same iteration count and identical final row. IEEE123 remains strictly
converged but its `T=12` and larger iteration counts differ from the oldest
timing table, so those rows are timing evidence rather than claims of an
identical floating-point trajectory.

## Interpretation

These are single runs, not statistical timing estimates. Nevertheless, the
consistent large reductions across every medium/large horizon show that the
changes remove real repeated work rather than only shifting memory around.
They do not change FilterDDP's mathematical burden: the backward pass still
solves for all state-sensitivity columns needed by the value recursion. The
next algorithmic target remains an exact structured or Schur-complement
treatment of that many-right-hand-side solve.

large10k `T=12` was deliberately not launched. The original run took 19.43 h;
even extrapolating the measured `T=6` improvement leaves a potentially
half-day run, which is not justified solely to fill a timing-table cell.
