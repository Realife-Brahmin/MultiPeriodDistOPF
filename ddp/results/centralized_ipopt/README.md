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

`T = 12` converged locally in 41 IPOPT iterations to `2781.5319546502`. The
independent validator accepted every constraint. The fresh run recorded
0.811 s inside IPOPT, 0.828 s from JuMP's solver timer, 1.561 s solve wall
time, and a 937.926 MiB sampled Julia-process working set. Problem size
(7968 variables, 5220 equality and 4896 inequality constraints) is exactly
double the `T = 6` row and quadruple `T = 3`, consistent with the per-stage
counts of 664 variables, 435 equality and 408 inequality constraints seen at
every horizon so far. The TPEC table previously carried no centralized entry
for this row (`---`), so this is a new value rather than a supersession.

`T = 24` converged locally in 42 IPOPT iterations to `2821.36125475254`. The
independent validator accepted every constraint. The fresh run recorded
3.035 s inside IPOPT, 3.055 s from JuMP's solver timer, 4.395 s solve wall
time, and a 960.039 MiB sampled Julia-process working set. Problem size
(15936 variables, 10440 equality and 9792 inequality constraints) again
scales linearly at 664/435/408 per stage. Unlike the `T = 3` row, this
fresh 3.055 s is close to the old TPEC entry of 3.060 s (0.005 s, 0.16%
apart) rather than a large supersession — a useful cross-check that the
reconstruction agrees with the prior measurement where the prior measurement
happened to already be accurate.

`T = 48` converged locally in 43 IPOPT iterations to `2844.65680064169`. The
independent validator accepted every constraint. The fresh run recorded
6.141 s inside IPOPT, 6.173 s from JuMP's solver timer, 7.528 s solve wall
time, and a 991.754 MiB sampled Julia-process working set. Problem size
(31872 variables, 20880 equality and 19584 inequality constraints) again
scales linearly at 664/435/408 per stage. The TPEC table previously carried
no centralized entry for this row (`---`), so this is a new value rather
than a supersession.

`T = 96` converged locally in 46 IPOPT iterations to `2857.16211357218`. The
independent validator accepted every constraint. The fresh run recorded
12.242 s inside IPOPT, 12.288 s from JuMP's solver timer, 13.406 s solve wall
time, and a 1068.379 MiB sampled Julia-process working set. Problem size
(63744 variables, 41760 equality and 39168 inequality constraints) again
scales linearly at 664/435/408 per stage. The fresh 12.288 s is 2.2% below
the old TPEC entry of 12.570 s — a moderate supersession, larger than `T = 24`
(0.16%) but far smaller than `T = 3` (56%).

`T = 144` converged locally in 50 IPOPT iterations to `2861.42874574193`. The
independent validator accepted every constraint. The fresh run recorded
20.848 s inside IPOPT, 20.917 s from JuMP's solver timer, 22.179 s solve wall
time, and a 1105.289 MiB sampled Julia-process working set. Problem size
(95616 variables, 62640 equality and 58752 inequality constraints) again
scales linearly at 664/435/408 per stage. Unlike every other ieee123 row,
FilterDDP has never been run at this horizon (Table I's tested-horizons list
stops at 96), so this centralized value has no `Optimized (s, x)` counterpart
yet in the TPEC table's Table II. Recorded here per the sweep matrix; whether
and how to add a Table II row is a paper-content decision left to the user.

## Completed ieee2522 cases

`T = 3` converged locally in 52 IPOPT iterations to `8515.8800338502`. The
independent validator accepted every constraint. The fresh run recorded
5.527 s inside IPOPT, 5.576 s from JuMP's solver timer, 6.783 s solve wall
time, and a 1008.52 MiB sampled Julia-process working set. Problem size is
32511 variables, 23448 equality and 17379 inequality constraints. This
supersedes the old TPEC entry of `1.730` (Gurobi-sourced); the new IPOPT
value is about 3.2x higher, reflecting a solver change rather than a
regression, per the user's decision to switch med2522/large10k's Table II
`C` column from Gurobi to IPOPT to match ieee123.

Run and publication instructions are in
`ddp/notes/CENTRALIZED_IPOPT_TIMING_SWEEP.md`.
