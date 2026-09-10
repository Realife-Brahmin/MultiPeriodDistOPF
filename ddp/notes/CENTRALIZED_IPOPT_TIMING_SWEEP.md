# Centralized IPOPT timing sweep

This sweep reconstructs the centralized-time column, `C (s)`, using fresh,
cold-start JuMP--IPOPT runs of the same BFM-NL MPOPF model and profiles used by
the tADMM study.

Decided 2026-09-10: the sweep now covers all three systems and does
supersede the older native-conic Gurobi timings for `ieee2522C_1ph` and
`large10kC_1ph`, not just the already-IPOPT-attributed `ieee123C_1ph`. The
TPEC Table II caption is updated to say IPOPT for all three systems once the
first non-ieee123 row lands.

## Matrix and order

- `ieee123C_1ph`: `T = 3, 6, 12, 24, 48, 96, 144`
- `ieee2522C_1ph`: `T = 3, 6, 12, 24, 48, 96, 144`
- `large10kC_1ph`: `T = 3, 6, 12, 24, 48`

Run only one case at a time, normally in the order above. The authoritative
machine-readable table is
`ddp/results/centralized_ipopt/centralized_ipopt_timing.csv`.

## One-case command

From repository root on Windows:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_centralized_ipopt_case.ps1 `
  -System ieee123C_1ph -Horizon 6
```

The runner forces IPOPT, disables the user Julia startup file, samples the
Julia process working set, preserves the IPOPT log and validation summary, and
upserts one CSV row. Never edit `config.jl` for a sweep case.

## Required evidence per row

Do not publish a `C (s)` value unless all of the following are retained:

- system, `T`, `Delta t = 24/T`, cold-start status, and validation result;
- objective, IPOPT iterations, IPOPT-reported time, JuMP solve time, and solve
  wall time;
- variable/constraint counts, Jacobian/Hessian nonzeros, IPOPT version and
  linear solver;
- sampled peak working set and the complete final IPOPT log;
- an explicit failure reason instead of a fabricated timing if the run fails.

Use `jump_solve_time_s` as the paper's `C (s)` value, matching the solver-time
meaning of the existing computational table. Keep `ipopt_reported_s` and
`solve_wall_s` beside it so that the choice remains auditable.

## Transaction after every run

1. Inspect the raw log and require an accepted solve status.
2. Require the independent post-solve validator to report `FEASIBLE`.
3. Sanity-check objective and dimensions against prior same-system results.
4. Update this repository's CSV/README and any relevant shared context.
5. Commit and push those files on `ddp-understanding-sep02`.
6. Update the TPEC source table from the verified CSV row; compile and visually
   inspect the PDF.
7. Commit and push the TPEC source and PDF on its `main` branch.
8. Only then start the next case. Never overlap cases.

If interrupted, inspect the CSV, raw logs, both Git histories, and process list.
Resume at the first row missing either a validated CSV entry or both pushed
repository commits.

## Initial cases

The initial `ieee123C_1ph`, `T = 3` and `T = 6` runs on 2026-09-09 used IPOPT
3.14.19 with MUMPS 5.8.2. They converged in 34 and 44 iterations to objectives
`2808.92465122283` and `2973.5533406443265`; all constraints passed validation.
These fresh timing rows, rather than old paper-table values, are authoritative
for this reconstruction.

The `ieee123C_1ph`, `T = 12` run on 2026-09-10 used the same IPOPT 3.14.19 /
MUMPS 5.8.2. build. It converged in 41 iterations to `2781.5319546502`; all
constraints passed validation. No prior centralized timing existed for this
row in the TPEC table, so this is a new entry rather than a supersession.

The `ieee123C_1ph`, `T = 24` run on 2026-09-10 used the same build. It
converged in 42 iterations to `2821.36125475254`; all constraints passed
validation. The fresh `jump_solve_time_s = 3.055` is close to the old TPEC
entry of `3.060`, unlike the larger T=3 supersession.

The `ieee123C_1ph`, `T = 48` run on 2026-09-10 used the same build. It
converged in 43 iterations to `2844.65680064169`; all constraints passed
validation. No prior centralized timing existed for this row in the TPEC
table, so this is a new entry rather than a supersession.

The `ieee123C_1ph`, `T = 96` run on 2026-09-10 used the same build. It
converged in 46 iterations to `2857.16211357218`; all constraints passed
validation. The fresh `jump_solve_time_s = 12.288` is 2.2% below the old
TPEC entry of `12.570`.

The `ieee123C_1ph`, `T = 144` run on 2026-09-10 used the same build. It
converged in 50 iterations to `2861.42874574193`; all constraints passed
validation. This completes the ieee123 matrix in this file's own order, but
FilterDDP has never been run at `T = 144` for ieee123 (Table I stops at 96
and Table II has no row for it), so there is no existing TPEC table cell to
update. Do not invent a Table II row; ask the user how they want this
represented before touching the TPEC repository for this case.

Decided 2026-09-10: hold off entirely. The TPEC repository is not touched
for `ieee123C_1ph` `T = 144` until FilterDDP is actually run at that horizon.
The centralized value stays recorded only in this repository's CSV/README.

The `ieee2522C_1ph`, `T = 3` run on 2026-09-10 used the same build. It
converged in 52 iterations to `8515.8800338502`; all constraints passed
validation. This is the first IPOPT row for this system; it supersedes the
old Gurobi-sourced TPEC entry of `1.730`.

The `ieee2522C_1ph`, `T = 6` run on 2026-09-10 used the same build. It
converged in 68 iterations to `9040.01299652051`; all constraints passed
validation. Dimensions scale linearly against `T = 3` (10837/7816/5793 per
stage). Supersedes the old Gurobi-sourced TPEC entry of `5.334`.

The `ieee2522C_1ph`, `T = 12` run on 2026-09-10 used the same build. It
converged in 70 iterations to `8512.51541531388`; all constraints passed
validation. Dimensions scale linearly at 10837/7816/5793 per stage.
Supersedes the old Gurobi-sourced TPEC entry of `11.015`.

The `ieee2522C_1ph`, `T = 24` run on 2026-09-10 used the same build. It
converged in 74 iterations to `8632.27397397632`; all constraints passed
validation. Dimensions scale linearly at 10837/7816/5793 per stage.
Supersedes the old Gurobi-sourced TPEC entry of `27.224`.

The `ieee2522C_1ph`, `T = 48` run on 2026-09-10 used the same build. It
converged in 78 iterations to `8701.14567443687`; all constraints passed
validation. Dimensions scale linearly at 10837/7816/5793 per stage.
Supersedes the old Gurobi-sourced TPEC entry of `58.151`. Objective matches
FilterDDP's own stored `T = 48` objective, `8701.148876`, to ~4e-7 relative.

The `ieee2522C_1ph`, `T = 96` run on 2026-09-10 used the same build. It
converged in 83 iterations to `8737.64779735172`; all constraints passed
validation. Dimensions scale linearly at 10837/7816/5793 per stage,
completing the med2522 matrix. Supersedes the old Gurobi-sourced TPEC entry
of `121.305`. Objective matches the paper's already-cited `T = 96`
centralized value, `8737.6487`, to ~1e-7 relative.

The `ieee2522C_1ph`, `T = 144` run on 2026-09-10 used the same build. It
converged in 88 iterations to `8750.07667926917`; all constraints passed
validation. Dimensions scale linearly at 10837/7816/5793 per stage. As with
ieee123 `T = 144`, FilterDDP has never been run at this horizon for med2522
either, so per the same 2026-09-10 hold-off decision this value stays
recorded only in this repository's CSV/README; the TPEC repository is not
touched for this row.

The `large10kC_1ph`, `T = 3` run on 2026-09-10 used the same build. It
converged in 110 iterations to `2985283.79326033`; all constraints passed
validation. Dimensions are 44345/31983/23703 per stage. Supersedes the old
Gurobi-sourced TPEC entry of `6.304` with `jump_solve_time_s = 55.809`
(about 8.85x higher, a larger gap than med2522 showed).
