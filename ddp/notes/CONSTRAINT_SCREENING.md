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
