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

`T = 6` converged locally in 68 IPOPT iterations to `9040.01299652051`. The
independent validator accepted every constraint. The fresh run recorded
15.037 s inside IPOPT, 15.064 s from JuMP's solver timer, 16.108 s solve wall
time, and a 1087.504 MiB sampled Julia-process working set. Problem size
(65022 variables, 46896 equality and 34758 inequality constraints) scales
linearly against `T = 3` at 10837/7816/5793 per stage. This supersedes the
old TPEC entry of `5.334` (Gurobi-sourced); the new IPOPT value is about
2.8x higher, the same direction as the `T = 3` change.

`T = 12` converged locally in 70 IPOPT iterations to `8512.51541531388`. The
independent validator accepted every constraint. The fresh run recorded
31.083 s inside IPOPT, 31.171 s from JuMP's solver timer, 32.39 s solve wall
time, and a 1225.25 MiB sampled Julia-process working set. Problem size
(130044 variables, 93792 equality and 69516 inequality constraints) scales
linearly at 10837/7816/5793 per stage. This supersedes the old TPEC entry of
`11.015` (Gurobi-sourced); the new IPOPT value is about 2.8x higher, the same
direction as `T = 3` and `T = 6`.

`T = 24` converged locally in 74 IPOPT iterations to `8632.27397397632`. The
independent validator accepted every constraint. The fresh run recorded
66.572 s inside IPOPT, 66.793 s from JuMP's solver timer, 68.152 s solve wall
time, and a 1630.133 MiB sampled Julia-process working set. Problem size
(260088 variables, 187584 equality and 139032 inequality constraints) scales
linearly at 10837/7816/5793 per stage. This supersedes the old TPEC entry of
`27.224` (Gurobi-sourced); the new IPOPT value is about 2.5x higher.

`T = 48` converged locally in 78 IPOPT iterations to `8701.14567443687`. The
independent validator accepted every constraint. The fresh run recorded
137.531 s inside IPOPT, 138.733 s from JuMP's solver timer, 140.475 s solve
wall time, and a 2311.605 MiB sampled Julia-process working set. Problem
size (520176 variables, 375168 equality and 278064 inequality constraints)
scales linearly at 10837/7816/5793 per stage. This supersedes the old TPEC
entry of `58.151` (Gurobi-sourced); the new IPOPT value is about 2.4x
higher. This centralized objective also matches FilterDDP's own reported
`T = 48` objective, `8701.148876`, to about 4e-7 relative — an independent
cross-check that the two models solve the same problem.

`T = 96` converged locally in 83 IPOPT iterations to `8737.64779735172`. The
independent validator accepted every constraint. The fresh run recorded
304.973 s inside IPOPT, 305.233 s from JuMP's solver timer, 309.564 s solve
wall time, and a 3853.211 MiB sampled Julia-process working set. Problem
size (1040352 variables, 750336 equality and 556128 inequality constraints)
scales linearly at 10837/7816/5793 per stage, completing the med2522
matrix. This supersedes the old TPEC entry of `121.305` (Gurobi-sourced);
the new IPOPT value is about 2.5x higher. The objective also matches the
`8737.6487` value already cited in the paper's own `T = 96` FilterDDP
comparison to about 1e-7 relative, confirming it is the same reference
value rather than a new one.

`T = 144` converged locally in 88 IPOPT iterations to `8750.07667926917`.
The independent validator accepted every constraint. The fresh run recorded
483.091 s inside IPOPT, 483.468 s from JuMP's solver timer, 489.479 s solve
wall time, and a 5555 MiB sampled Julia-process working set. Problem size
(1560528 variables, 1125504 equality and 834192 inequality constraints)
scales linearly at 10837/7816/5793 per stage. As with ieee123 `T = 144`,
FilterDDP has never been run at this horizon for med2522 (Table I's
tested-horizons list stops at 96, and Table II has no row for it), so per
the same user decision this value is recorded here only; the TPEC
repository is not touched for this row.

## Completed large10k cases

`T = 3` converged locally in 110 IPOPT iterations to `2985283.79326033`. The
independent validator accepted every constraint. The fresh run recorded
55.762 s inside IPOPT, 55.809 s from JuMP's solver timer, 56.922 s solve wall
time, and a 1333.93 MiB sampled Julia-process working set. Problem size is
133035 variables, 95949 equality and 71109 inequality constraints (44345/
31983/23703 per stage). This supersedes the old TPEC entry of `6.304`
(Gurobi-sourced); the new IPOPT value is about 8.85x higher — a much larger
IPOPT-vs-Gurobi gap than med2522's 2.5-3x, plausible given Gurobi's native
conic solver is far better suited to a huge SOCP than IPOPT's general
interior-point method.

`T = 6` converged locally in 99 IPOPT iterations to `2944457.00867062`. The
independent validator accepted every constraint. The fresh run recorded
104.760 s inside IPOPT, 104.843 s from JuMP's solver timer, 106.239 s solve
wall time, and a 1758.609 MiB sampled Julia-process working set. Problem
size (266070 variables, 191898 equality and 142218 inequality constraints)
scales linearly against `T = 3` at 44345/31983/23703 per stage. No prior
centralized timing was recorded for this row (`---`), so this is a new
value rather than a supersession.

`T = 12` converged locally in 88 IPOPT iterations to `2976105.08846598`. The
independent validator accepted every constraint. The fresh run recorded
209.086 s inside IPOPT, 209.447 s from JuMP's solver timer, 211.369 s solve
wall time, and a 2530.766 MiB sampled Julia-process working set. Problem
size (532140 variables, 383796 equality and 284436 inequality constraints)
scales linearly at 44345/31983/23703 per stage. This objective matches the
paper's already-cited FilterDDP `T = 12` stabilized objective,
`2,976,105.0928`, to about 1.5e-9 relative — the paper's text currently
says no independent centralized `T = 12` objective is stored for large10k;
that statement needs updating once this row is adopted into the TPEC
table.

`T = 24` converged locally in 112 IPOPT iterations to `2996960.01465613`.
The independent validator accepted every constraint. The fresh run
recorded 572.724 s inside IPOPT, 574.926 s from JuMP's solver timer,
578.422 s solve wall time, and a 4270.617 MiB sampled Julia-process working
set. Problem size (1064280 variables, 767592 equality and 568872
inequality constraints) scales linearly at 44345/31983/23703 per stage.
Unlike large10k `T = 3, 6, 12`, FilterDDP has never been run at this
horizon (Table I's large10k tested-horizons list is `3, 6; 12†` only, and
Table II has no row for it), so per the same hold-off precedent as the
ieee123/med2522 `T = 144` rows, this value is recorded here only; the TPEC
repository is not touched for this row.

Run and publication instructions are in
`ddp/notes/CENTRALIZED_IPOPT_TIMING_SWEEP.md`.
