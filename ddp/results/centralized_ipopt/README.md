# Centralized IPOPT timing results

These are fresh cold-start JuMP--IPOPT timings for reconstructing the
centralized `C (s)` column of the computational comparison. The authoritative
table is `centralized_ipopt_timing.csv`; `raw/` retains the corresponding IPOPT
tail and independent validator summary for every completed case.

The paper uses `jump_solve_time_s` for `C (s)`. `ipopt_reported_s` is IPOPT's
own timer, while `solve_wall_s` includes the surrounding `optimize!` call. Do
not mix any of these with model parsing/building or the complete pipeline time.

## Completed ieee123 cases

`T = 3` converged in 34 iterations to `2808.92465122283`, with 0.280 s
inside IPOPT, 0.294 s from JuMP's solver timer, and 1.112 s solve wall time.
This fresh result supersedes the older 0.188-s timing in the TPEC table.

`T = 6` converged locally in 44 IPOPT iterations to
`2973.55334064433`. The independent validator accepted every constraint. The
fresh run recorded 0.640 s inside IPOPT, 0.654 s from JuMP's solver timer,
1.396 s solve wall time, and a 919.598 MiB sampled Julia-process working set.
The working set includes the Julia runtime and loaded packages and is not an
IPOPT-only memory measurement.

Run and publication instructions are in
`ddp/notes/CENTRALIZED_IPOPT_TIMING_SWEEP.md`.
