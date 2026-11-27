# Storage-coordinated tADMM: first checks (2026-11-15)

Verification only, on the 309 PC, serial. Iteration counts and solution quality
are the result; no wall-clock number here is a benchmark.

## What is tested

`storage_coordinator.jl` keeps the parallel step of the review-pilot tADMM (one
branch-flow OPF per period, same model, same Ipopt settings, MA57) and changes
the coordination:

- consensus variable: each battery's power at each period (was: energy copies
  over a three-period window);
- Update II: one small problem per battery over the whole horizon, the schedule
  closest to the requested powers that satisfies the energy dynamics and the
  energy and power limits (was: averaging neighbouring copies);
- the period subproblems carry no energy variables.

Information therefore crosses the whole horizon in every iteration, where
neighbour averaging moves it one period per iteration.

Instances: current periodic exports, July battery coefficient
(`TADMM_CB_PROFILE=july`), free terminal SOC, corrected losses. The reference is
the centralized solve of the same model. Near-optimal ("NO") means objective
within 0.1% of the reference and original-model violation at most 1e-6
(ieee123) or 1e-5 (med2522), evaluated on the network dispatch with energies
integrated from its battery powers. The old scheme is `run_admm` with the
harness defaults (`baseline_rounds.jl`), not a tuned configuration.

## Iterations to the NO point

ieee123 (51 batteries):

| T | old scheme | new, best fixed rho | new, decreasing rho |
|---|---|---|---|
| 24 | 407 | 50 | 50 |
| 48 | not reached in 1000 (1.1% off) | 43 | 34 |
| 96 | not reached in 1000 (3.5% off) | 63 | 50 |

med2522 (249 batteries), T = 24: 55 at rho = 3000, 65 at 10000, 78 at 30000.
ChatGPT's best run of the old scheme with the same coefficient is 126
(`results/july_weight_comparison`).

No run met the strict stop (consensus 2.5e-7, dual 1e-4, original 5e-7) within
its 120 to 200 iterations; the consensus error flattens near 1e-7 to 1e-6. Not
yet investigated.

## The penalty

`rho_unit = (price spread) * kVA_B * dt / (median battery power rating)`, the
penalty at which the day's price spread moves a median battery by its rating in
one iteration. It halves when T doubles (12000, 6000, 3000 on ieee123).

- The best fixed penalty is 0.25 `rho_unit` at T = 24, 48 and 96 (grid 0.0625,
  0.125, 0.25, 0.5).
- A geometric decrease from 2 to 0.0625 `rho_unit` over 20 iterations, then
  held, matches or beats it: 50, 34, 50.
- Increasing schedules are worse at every horizon (77 to 107, or not reached).
- A high penalty reduces the violation quickly and the objective slowly; a low
  one does the reverse (T = 24: rho = 1e5 has violation below 1e-6 at iteration
  17 and objective within 0.1% at 84; rho = 1e4 at 62 and 8).

The 0.25 factor is fitted on ieee123, not derived.

## Files

Per case: `trace_july_<tag>.csv` (one row per iteration and penalty setting),
`baseline_july_trace.csv` (old scheme), `reference_july.txt`, `schedule.log`.
`batch.log` is the first batch. Launchers: `run_storage_coordinator.sh`,
`batch_storage_coordinator_nov15.sh` (run from the repository root).
