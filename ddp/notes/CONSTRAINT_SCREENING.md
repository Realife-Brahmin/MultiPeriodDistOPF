# Constraint screening: which inequalities matter, and what the others cost

Agenda item of 2026-10-07 (constraint screening / presolve). Method, as asked:
count the active and inactive inequalities first, then act on what the data
shows. Cases: ieee123 and med2522 at `T = 24`, large10k at `T = 6` and `24`,
matched instances (periodic profile, per-system `C_B`, soft terminal SOC).

## Short answers

- **About three quarters of the inequalities are inactive at the optimum.**
  Voltage limits are the largest group and are 99.9% inactive.
- **An inactive bound costs FilterDDP no linear algebra** (a bound is one
  diagonal term). **It costs step length.** The line search halves the step
  whenever any variable would get too close to any bound.
- **One redundant bound was throttling every run.** The substation voltage is
  fixed at 1.05 pu by an equality, and its upper limit is 1.05 pu. Every full
  step lands on that bound. On ieee123 it set the step in 80 of 85
  iterations; no run ever took a full step.
- **Three exact changes** cut the iterations to near-optimality from 85 / 77 /
  101 to 19 / 39 / 11 (ieee123 `T=24`, med2522 `T=24`, large10k `T=6`).
- **Ipopt gains from two of the three as well** (53 to 12 iterations on
  large10k), so the comparison must be redone with both solvers on the same
  model and start. Timed results: Section 8.

## 1. How many inequalities, how many active

Counted on Ipopt's solution (`constraint_activity.jl`). *Active*: slack below
1e-6 of the range, or multiplier larger than slack. *Near*: within 1% of the
range.

| | ieee123 T=24 | med2522 T=24 | large10k T=6 | large10k T=24 |
|---|---|---|---|---|
| all inequalities | 19,560 | 277,896 | 284,430 | 1,137,720 |
| active | 24.7% | 24.6% | 22.8% | 23.0% |
| active, leaving out the SOC rows | 12.1% | 3.6% | 1.6% | 1.9% |

Active share by type:

| inequality | ieee123 | med2522 | large10k T=6 | large10k T=24 |
|---|---|---|---|---|
| voltage lower limit | 0 | 0.1% (80) | 0 | 0 |
| voltage upper limit | 0.8% (24) | 0.05% (30) | 1.0% | 1.0% |
| battery power, charge / discharge | 2.5% / 3.3% | 0.6% / 3.0% | 0 / 0 | 2.5% / 4.2% |
| battery energy, lower / upper | 32.5% / 17.2% | 31.0% / 11.0% | 33.3% / 16.7% | 33.1% / 19.2% |
| DER reactive, lower / upper | 0 / 93.7% | 0.03% / 61.9% | 0 / 0 | 0 / 0 |
| `ell >= 0` | 4.7% | 2.0% | 0 | 0 |
| SOC relaxation | 92.9% | 100% | 99.0% | 99.0% |
| `P_Subs >= 0` | 0 | 0 | 0 | 0 |

- The SOC rows are tight almost everywhere: they are equalities in practice.
  The loose ones are lines of near-zero resistance.
- **Voltage.** On ieee123 and large10k the only active limits are upper
  limits at buses joined to the substation by near-zero impedance: they sit
  at 1.05 pu because the substation does. Lowest voltage 1.014 pu (ieee123),
  0.984 / 0.973 pu (large10k).
- **med2522 is the one case with real voltage activity**: 80 lower limits on
  24 buses at the far end of the feeder (depth 188-275 of 275), in periods
  1, 2, 18, 19, 20. One upper limit two buses below the substation is also
  genuinely active (it caps reactive export).
- Battery energy limits and the DER reactive upper limit are the constraints
  that shape the solution. Nothing to screen there.

## 2. What an inactive bound costs: step length

FilterDDP's forward pass tries step 1 and halves it whenever a trial point
breaks the fraction-to-boundary rule (a control or a bound multiplier would
lose more than 99% of its distance to zero). Those halvings were not counted
as backtracks, so the logs never showed them. `FILTERDDP_FTB_DIAGNOSTIC=1` now
prints each one; `step_blockers.py` summarises them.

Baseline runs (flat start, all bounds, halving line search):

| | ieee123 T=24 | med2522 T=24 | large10k T=6 |
|---|---|---|---|
| iterations to near-optimality | 85 | 77 | 101 |
| full steps taken | 0 | 0 | 0 |
| steps of 0.5 / 0.25 / smaller | 48 / 35 / 2 | 34 / 31 / 12 | 27 / 1 / 73 |
| step set by the substation voltage limit | 94% | 40% | 30% |
| by the energy-window upper limit | | | 35% |
| by the multiplier of `ell >= 0` | 2% | 4% | 18% |
| by the multiplier of `P_Subs >= 0` | | | 8% |

- **Substation voltage.** Fixed at 1.05^2 by its equality row; upper limit
  1.05^2. The Newton step puts it exactly on the bound, the rule rejects the
  full step, the solver takes 0.5, and the slack halves each iteration:
  1.5e-2, 7.4e-3, 3.7e-3, ..., 2.2e-16. Once the slack reaches machine precision even 0.5 is rejected and
  the step drops to 0.25 (the last 35 iterations on ieee123).
- **large10k** spends its first 74 iterations at steps of 0.125 or less. From
  the flat start the solver has to find the power flow through the barrier:
  the substation import must grow from 1 to over 800 p.u., and the multiplier of
  `P_Subs >= 0` lets it roughly double per iteration.

## 3. Limits that cannot bind: exact rules, no OPF solution needed

`voltage_screening.jl`. In receiving-end form the voltage drop on line
`(i, j)` is

```
v_i - v_j = 2 (r P_recv + x Q_recv) + |z|^2 ell ,
```

and the flow arriving at `j` is the net load of `j`'s subtree plus the losses
inside it. With `ell >= 0`, `r, x >= 0`, every feasible point has
`P_recv >= Pmin_j`, `Q_recv >= Qmin_j`: the lossless subtree load with every
battery discharging at rating and every DER at its reactive limit. So
`v_i - v_j >= d_j := 2 (r Pmin_j + x Qmin_j)`, one tree sweep per period.

| rule | statement |
|---|---|
| substation | `v_root` is fixed by its own equality: its limits are redundant |
| upper bound | `v_j <= 1.05^2 - sum of d along the path`; if that is within the limit, the upper limit cannot bind |
| downhill | if `d_j >= 0`, voltage cannot rise from `i` to `j`: `j`'s upper limit is implied by `i`'s, and `i`'s lower limit by `j`'s |
| `ell >= 0` | implied by the SOC row `P^2 + Q^2 + s = v ell`, `s >= 0`, `v > 0` |
| `P_Subs >= 0` | `P_Subs >= sum of Pmin over the substation's lines`; if positive, cannot bind |

The "passive chain" idea (a lateral with no battery or PV: lower limit at the
bottom, upper limit at the top) is the downhill rule on a subtree with no
devices. The downhill rule also covers subtrees whose devices cannot outweigh
their load: 60% / 80% / 88% of lines are downhill, against 5% / 42% / 21% of
buses passive.

What the rules keep (share of bus-period limits):

| | ieee123 | med2522 | large10k |
|---|---|---|---|
| voltage upper limits kept | 0% | 13.3% | 0% |
| voltage lower limits kept | 50.9% | 34.6% | 16.6% |
| `ell >= 0`, `P_Subs >= 0` kept | none | none | none |
| premises broken by the solution | 0 | 0 | 0 |

- Every rule is a relaxation, never a guess at the active set: a dropped
  limit is satisfied by every point that meets the remaining constraints.
  Each run also prints `BOUND_CHECK` (all voltages, `ell`, `P_Subs` against
  the original limits at the solution): no violation in any run.
- Lower voltage limits are kept in the solver even where the rules allow
  dropping them. They are real on med2522, they keep `v > 0` during the
  iterations, and dropping them changed nothing on ieee123.
- Limits in the FilterDDP model are per control, not per period, so a limit
  is dropped only if no period needs it (med2522: 578 of 2,522 buses keep
  their upper limit).

## 4. The power-flow start

`loadflow_start.jl`, `FILTERDDP_LOADFLOW_START=1`. With batteries idle and DER
reactive power zero, the network state of each period is one
backward/forward sweep on the tree (6-10 sweeps). The solve then starts
from a point that already satisfies the network equations (largest residual
1e-6 on large10k) instead of from voltages 1.0 and flows 0.

- large10k: decisive (Section 6).
- ieee123: no effect. The run rejoins the same path after 13 iterations.
- med2522: no effect. The idle-battery power flow dips to 0.881 pu, below the
  limit, so the solver has to move it anyway.

## 5. Exact largest step (affine line search)

`affine_linesearch.jl`, `FILTERDDP_AFFINE_LINESEARCH=1`. The battery dynamics
are affine, so the whole trial trajectory is affine in the step size: one
unit-step rollout gives the direction of every control and multiplier. Then

- the largest step the fraction-to-boundary rule allows is a ratio test (as
  in Ipopt), not "halve until it passes": 0.97 is taken as 0.97, not 0.5;
- each further trial is an axpy and a constraint evaluation; the stage
  systems are solved once per iteration, not once per trial.

The acceptance tests (filter, switching, Armijo) are unchanged. This changes
the iterates. Tried and rejected: Ipopt-style separate steps for the bound
multipliers (`FILTERDDP_SPLIT_STEP=1`): 60 iterations against 31 on ieee123.

## 6. Effect on iterations

Verification runs (not timed), diagonal Hessian. *Screen* =
`substation,vupper,ell,psubs`. Iterations to near-optimality:

| configuration | ieee123 T=24 | med2522 T=24 | large10k T=6 |
|---|---|---|---|
| baseline | 85 | 77 | 101 |
| substation limit dropped | 33 | 53 | 90 |
| screen + power-flow start | 32 | 48 | 20 |
| screen + exact step | 19 | 38 | 71 |
| **screen + power-flow start + exact step** | **19** | **39** | **11** |

To the solver's own tolerance (1e-7), and Ipopt (MA57, its own 1e-8):

| | ieee123 T=24 | med2522 T=24 | large10k T=6 |
|---|---|---|---|
| FilterDDP baseline | 93 | 94 | 128 |
| FilterDDP, all three | 33 | 77 | 24 |
| Ipopt as before | 38 | 60 | 53 |
| Ipopt, screen | 29 | 58 | 40 |
| Ipopt, screen + power-flow start | 21 | 53 | 12 |

- Objectives agree with Ipopt to 1e-7 (ieee123), 2e-8 (med2522), 1e-9
  (large10k) in the strict runs.
- ieee123 at `T = 6` and `96`: 20 and 55 iterations to near-optimality
  (67 and 120 before).

## 7. Caveats

- **Near-optimality is now declared at a much looser point.** Before, primal
  feasibility arrived last, by which time the objective was within 1e-7 to
  5e-5 of the optimum. Now the iterate is feasible early and the 0.5%
  objective test decides: the gap at the stopping point is 0.19-0.28% on
  med2522 and 0.03-0.34% on large10k (large10k `T=6`: `mu = 0.04`, smallest
  `ell` 0.42 p.u., SOC slacks still loose). The like-for-like figure is the
  time to a fixed, tighter gap or to the solver's own tolerance (Section 8).
- **Timings so far included Julia's one-time compilation**: 15-30 s, all
  inside the first iteration. That is 14.7 of the 17.9 s reported for ieee123
  `T=6`, and 23 of 260 s for large10k `T=6`. Ipopt's solve time has no such
  term. `FILTERDDP_WARMUP=full` runs the whole solve once first, so the
  timed (second) solve excludes it; same iterates. Reported times are solve
  times from here on (user, 2026-10-05).
- **The same model and start help Ipopt**, most on large10k. Comparisons
  must use Ipopt with `IPOPT_SCREEN` and `IPOPT_LOADFLOW_START`.
- med2522 gains least: its remaining short steps come from genuinely active
  limits (the upper limit below the substation, the energy window).

## 8. Timed results

`run_table2_screen.sh`, 2026-10-05, 309 lab PC, Balanced power plan,
background load logged per case (0.5-1.0 core, as in earlier clean runs).
Both solvers on the screened model from the power-flow start. FilterDDP:
diagonal Hessian, ieee123 with UMFPACK on one thread, med2522 and large10k
with the tree solver on eight; time to near-optimality, compilation excluded.
Ipopt: solve time to its own tolerance; HSL is the best of MA57 and MA97.

| case | FilterDDP before | FilterDDP now | iterations | Ipopt MUMPS | Ipopt HSL | vs HSL | vs MUMPS |
|---|---|---|---|---|---|---|---|
| ieee123 T=6 | 17.9 | 1.0 | 67 -> 20 | 0.25 | 0.13 | 7.7x | 4.0x |
| ieee123 T=24 | 25.7 | 2.9 | 85 -> 19 | 0.74 | 0.45 | 6.3x | 3.9x |
| ieee123 T=96 | 70.4 | 32.4 | 120 -> 55 | 3.56 | 2.04 | 15.9x | 9.1x |
| med2522 T=6 | 57.5 | 18.0 | 68 -> 33 | 5.25 | 3.11 | 5.8x | 3.4x |
| med2522 T=24 | 177.8 | 82.5 | 77 -> 39 | 31.2 | 18.7 | 4.4x | 2.6x |
| med2522 T=96 | 811.7 | 418.8 | 96 -> 49 | 158.9 | 90.0 | 4.7x | 2.6x |
| large10k T=6 | 259.7 | 31.4 | 101 -> 11 | 12.3 | 5.0 | 6.3x | 2.6x |
| large10k T=24 | 924.2 | 105.2 | 98 -> 9 | 65.2 | 28.6 | 3.7x | 1.6x |
| large10k T=48 | 1537.0 | 279.4 | 85 -> 11 | 145.4 | 60.7 | 4.6x | 1.9x |

"Before" is the fastest earlier configuration and includes compilation
(15-30 s). The same runs with compilation inside: 17.0 / 19.2 / 43.3,
38.8 / 101.3 / 422.8, 53.6 / 124.9 / 292.2 s.

Ipopt on the same nine cells, before -> now (solve time, s):

| | MUMPS | best HSL | iterations (MA57) |
|---|---|---|---|
| ieee123 T=6 / 24 / 96 | 0.39 / 1.47 / 23.2 -> 0.25 / 0.74 / 3.56 | 0.19 / 0.74 / 3.91 -> 0.13 / 0.45 / 2.04 | 37 / 38 / 44 -> 20 / 21 / 25 |
| med2522 T=6 / 24 / 96 | 7.2 / 41.1 / 220.0 -> 5.2 / 31.2 / 158.9 | 4.4 / 24.4 / 127.3 -> 3.1 / 18.7 / 90.0 | 44 / 60 / 78 -> 39 / 53 / 64 |
| large10k T=6 / 24 / 48 | 51.2 / 282.9 / 623.7 -> 12.3 / 65.2 / 145.4 | 27.6 / 132.4 / 272.0 -> 5.0 / 28.6 / 60.7 | 53 / 66 / 75 -> 12 / 17 / 18 |

Reading:

- FilterDDP is 1.9-2.0x faster on med2522 and 5.4-8.5x on large10k, compilation
  aside (earlier runs less their first-iteration compile time). Ipopt gains 1.3-1.4x on med2522 and 4.5-5.5x on large10k from the
  model and start alone.
- Net: FilterDDP stands at 3.7-6.3x Ipopt's best HSL time and 1.6-3.4x MUMPS
  on the two larger systems (5.7-12.9x and 2.5-8.3x before, compilation included).
- Iteration counts are now comparable. What remains is the cost of one
  iteration: 8.4 s against 1.7 s at large10k `T=24`, 1.8 s against 0.35 s at
  med2522 `T=24` (FilterDDP against Ipopt-MA57).
- ieee123 `T=96` is the outlier (55 iterations): not examined yet.
- The FilterDDP figures stop at near-optimality, which is now a looser point
  (Section 7); Ipopt runs to 1e-8. Strict-tolerance FilterDDP runs of the
  nine cells are queued (`run_screen_followup.sh`, part E) to give the time
  to a fixed gap.

### Time to a fixed objective gap

From the strict-tolerance runs (`run_screen_followup.sh`, part E;
`near_opt_posthoc.py --gaps=...`), same primal thresholds. Seven of nine
cells so far. Starred times are from runs with other load on the machine
(iterations valid, times indicative); large10k `T=48` not run yet.

| case | gap 0.5% | gap 1e-4 | gap 1e-5 | own tolerance | Ipopt HSL | at 1e-5 vs HSL |
|---|---|---|---|---|---|---|
| ieee123 T=6 | 20 it, 1.1 s | 20, 1.1 | 23, 1.2 | 31, 1.4 | 0.13 | 9x |
| ieee123 T=24 | 19, 2.6 | 23, 3.0 | 26, 3.4 | 33, 4.4 | 0.45 | 7.5x |
| ieee123 T=96 | 55, 26.6 | 55, 26.6 | 55, 26.6 | 61, 29.5 | 2.04 | 13x |
| med2522 T=6 | 33, 18.1 | 48, 26.6 | 51, 28.4 | 67, 37.1 | 3.11 | 9.1x |
| med2522 T=24 | 39, 80.6 | 57, 116.2 | 62, 126.8 | 77, 157.0 | 18.7 | 6.8x |
| med2522 T=96 | 49, 409.8* | 71, 582.7* | 77, 630.5* | 94, 757.0* | 90.0 | 7.0x* |
| large10k T=6 | 11, 31.3 | 16, 46.1 | 17, 50.3 | 24, 74.0 | 5.0 | 10.1x |
| large10k T=24 | 9, 122.4* | 21, 257.5* | 24, 300.1* | 38, 445.2* | 28.6 | 10.5x* |

- The earlier near-optimal points sat at a gap of 1e-7 to 5e-5, so the
  1e-5 column is the like-for-like one. There FilterDDP is 7-10x Ipopt's
  best HSL time on the two larger systems: about where it stood before
  (5.7-12.9x). Both solvers got faster; the ratio did not move much.
- At equal accuracy the changes save med2522 about 20% of its iterations
  (62 against 77 at `T=24`) and large10k about 80% (17 against 101 at `T=6`).
  Ipopt's gain from the same model and start is of the same size on each.
- The 3.7-6.3x figure above holds only under the 0.5% rule, whose stopping
  point is now looser. Which criterion the paper uses is open (user).

## 9. Follow-up (2026-10-05 evening to 2026-10-06)

### large10k with its PV units

Every large10k number above is from the instance without PV (CLAUDE.md). On
the corrected instance (logs tagged `pv`):

- Counts are unchanged at this resolution: 23.0% of 1,137,720 inequalities
  active at `T=24`, 1.9% without the conic rows, 0.5% of the voltage limits.
  The rules now keep 189 of 247,680 upper voltage limits at `T=24` (30 at
  `T=6`), where PV can raise the voltage, and 18.6% of the lower limits.
- Iterations to near-optimality at `T=6`: earlier configuration 100,
  substation limits dropped 82, S+P 21, S+E 59, S+P+E 10. To the solver's own
  tolerance: 120 against 23. (`ddp/results/screen_ablation/`)
- Step rejections of the earlier configuration: no full step in 100
  iterations, 66 at 1/4 or less; the step was set by the energy window in 42%
  of iterations, the substation voltage limit in 33%, the multiplier of
  `P_Subs >= 0` in 13% and of `ell >= 0` in 9%.

| large10k (PV) | `T=6` | `T=24` | `T=48` |
|---|---|---|---|
| Ipopt unscreened, iterations | 256 | 65 | 69 |
| ... MA57 / MA97 / MUMPS (s) | 115.2 / 321.0 / 206.5 | 130.6 / 113.2 / 238.0 | 362.6 / 224.2 / 506.6 |
| Ipopt screened + start, iterations | 12 | 18 | 19 |
| ... MA57 / MA97 / MUMPS (s) | 5.0 / 5.2 / 12.3 | 30.2 / 31.8 / 67.2 | 65.3 / 75.2 / 147.8 |
| FilterDDP to near-optimality | 10 it, 28.8 s | 11 it, 128.9 s | 11 it, 272.1 s |
| FilterDDP to its own tolerance | 23 it, 82.5 s | 37 it, 465.6 s | 32 it, 838.0 s |
| FilterDDP to a gap of 1e-5 | 17 it, 58.5 s | 23 it, 300.9 s | 21 it, 584.1 s |

MA97 needs 670 iterations on the unscreened `T=6` instance (321 s is its time
there). The FilterDDP times in this table are those of 2026-10-05; the faster
sequential times of 2026-10-06 are in `PARALLEL_IN_TIME.md`, Section 9.

### Ipopt was not throttled the way FilterDDP was

Its step is already a ratio test and it relaxes every bound by a relative
`1e-8`. One bound dropped at a time, large10k `T=6`, MA57:

| dropped | nothing | substation limits | upper voltage limits | `ell >= 0` | `P_Subs >= 0` | all four | all four + start |
|---|---|---|---|---|---|---|---|
| iterations | 256 | 95 | 79 | 30 | 1126 | 39 | 12 |

Its gain from screening is real but erratic and not tied to one bound.

### Longer horizons, current configuration

FilterDDP to near-optimality against the faster of MA57 and MA97 on the same
screened model and start:

| case | FilterDDP | Ipopt HSL | ratio | before |
|---|---|---|---|---|
| med2522 `T=192` | 51 it, 1102 s | 225 s | 4.9 | 108 it, 1905 s, 6.9 |
| med2522 `T=384` | 64 it, 2897 s | 474 s | 6.1 | 114 it, 4745 s, 8.0 |
| large10k `T=96` | 12 it, 646 s | 139 s | 4.6 | |
| large10k `T=192` | 13 it, 1249 s | 328 s | 3.8 | |

Peak memory: 7.3 against 8.3 GiB at med2522 `T=384`, 13.4 against 17.8 GiB at
large10k `T=192`. med2522 `T=1152` and `T=1536` were cancelled (user,
2026-10-06: horizons that long will not be reported).

### `C_B = 0` with the diagonal Hessian

Both cases that failed now reach near-optimality: large10k `T=6` in 17
iterations (48.2 s; Ipopt-MA57 on the screened model 22 iterations, 8.5 s;
806 iterations and 376 s unscreened) and med2522 `T=96` in 55 iterations
(443.7 s).

### Exact against diagonal Hessian

At the 0.5% stop the diagonal Hessian is faster on the two larger feeders, but
it stops at a gap of 0.1-0.3% while the exact-Hessian runs, whose feasibility
arrives last, stop at `1e-5` to `1e-7`. At a gap of `1e-5`:

| case | diagonal | exact |
|---|---|---|
| med2522 `T=6` | 51 it, 28.4 s | 42 it, 24.4 s |
| med2522 `T=24` | 62 it, 126.8 s | 52 it, 122.4 s |
| med2522 `T=96` | 77 it, 787 s | 66 it, 488 s (gap 2.7e-5) |
| large10k `T=6` | 17 it, 58.5 s | 12 it, 58.9 s |
| large10k `T=24` | 23 it, 300.9 s | 12 it, 225.8 s |

So which Hessian is faster depends on the accuracy asked for. (Single runs of
the overnight queue.)

### Gurobi 13, for reference

Unscreened model, default threads / one thread (s): ieee123 `T=24` 0.80 / 0.62,
`T=96` 9.1 / 3.5; med2522 `T=6` 4.5 / 3.4, `T=24` 24.6 / 28.5, `T=96`
219 / 102.5; large10k `T=6` 17.5 / 27.0, `T=24` 49.5 / 57.4. Its presolve
removes almost nothing here (12 rows on large10k `T=6`). Not yet run on the
screened model, which needs the cone stated explicitly once `ell >= 0` is
dropped (Gurobi recognises a rotated cone only with nonnegative variables).

### Check against OpenDSS

33 dispatches of the current configuration, horizons up to `T=384`
(`run_opendss_check_all.sh`, `ddp/results/opendss_check/`). OpenDSS converges
on all and no voltage limit is violated in any.

| stopping point | dispatches | voltages agree to | import overstated by |
|---|---|---|---|
| solver's own tolerance | 9 | 9.2e-8 pu | 0 |
| 0.5% stop, diagonal Hessian | 13 | 1.1e-3 pu | 0.0004% to 0.37% |
| 0.5% stop, exact Hessian | 8 | 8.3e-6 pu | up to 0.003% |

At the 0.5% stop the conic relaxation is not yet tight, so part of the
optimizer's loss is not physical and the dispatch costs less than the
objective reports.

### Where the four findings stand in the literature

From the review of 2026-10-06 (DOIs checked in Crossref; the paper's Section
on iteration reduction cites them):

- The no-interior bound: Waechter and Biegler (2006), Sec. 3.5, name the
  structure and relax bounds against it; presolvers remove the fixed variable
  (Gondzio 1997). The step throttling it causes under a halving line search
  was not found reported.
- The screening rules: the voltage bound is Lemma 1 of Gan, Li, Topcu and Low
  (2015), used there and in Low (2014) as a condition for exactness of the
  relaxation, not as a screening rule. `ell >= 0` from the cone is elementary;
  the substation import bound is standard presolve.
- The power-flow start: reported for interior-point AC OPF on transmission
  systems (Kardos et al. 2022). Not found for branch-flow multi-period OPF or
  interior-point DDP.
- The exact step: the ratio test is standard for all-at-once interior-point
  methods and for Riccati-based MPC with linear dynamics (Rao, Wright and
  Rawlings 1998); the closest DDP work applies it to a linearized slack
  direction (Prabhu, Rangarajan and Kothare 2025). That linear dynamics make
  the DDP rollout affine in the step, hence the test exact, was not found
  stated.

## Files

| | |
|---|---|
| `constraint_activity.jl`, `run_constraint_activity.sh` | the count; logs in `ddp/results/constraint_activity/` |
| `voltage_screening.jl` | the rules |
| `loadflow_start.jl` | power-flow start |
| `DDP4OPF.jl/src/affine_linesearch.jl` | exact-step line search |
| `DDP4OPF.jl/src/forward_pass.jl`, `step_blockers.py` | step-rejection diagnostic and its summary |
| `ieee123c_filterddp.jl` | `FILTERDDP_SCREEN`, `FILTERDDP_LOADFLOW_START`, `FILTERDDP_WARMUP`, `BOUND_CHECK` |
| `centralized_ipopt_matched.jl` | `IPOPT_SCREEN`, `IPOPT_LOADFLOW_START` |
