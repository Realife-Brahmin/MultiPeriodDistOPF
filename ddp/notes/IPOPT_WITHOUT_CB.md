# Centralized Ipopt with and without the battery term C_B

Question from the meeting of 2026-09-25: does the second battery term in the
objective, `C_B * pbase^2 * dt * sum P_B^2`, make the problem easier for the
optimizers, and how much does it restrain battery use? Centralized Ipopt only;
FilterDDP not yet.

## Setup

Every instance of the paper's Table II, solved twice with
`centralized_ipopt_matched.jl`, identical except for `C_B`:

- `C_B = 0` (`run_ipopt_no_cb.sh`, logs in `ddp/results/ipopt_no_cb/logs/`)
- `C_B = 1e-3`, the matched value (`run_ipopt_cb_metrics.sh`, logs in
  `ddp/results/ipopt_cb_1e-3/logs/`). The iteration counts reproduce the
  matched race exactly; this set was re-run only to record battery use, since
  the race logs (`matched_ipopt_race/logs/`, untouched) predate the metric.

Periodic profile, soft terminal SOC with the per-system gamma (kept in both
sets), same exported instances and Ipopt options. Median background load
0.50-0.82 cores in every run. Every solve ends `LOCALLY_SOLVED`.

The driver now prints `CENTRAL_IPOPT_BATTERY` after the timed solve:

- **mean power use**: mean of `|P_B| / P_B^rated` over batteries and periods
- **energy window used**: per battery, `(max - min)` of its SOC trajectory
  (including `B0`) over the usable window `(soc_max - soc_min) * B_R`, averaged
- **at power limit**: share of battery-periods with `|P_B| >= 0.99 P_B^rated`
- **throughput**: `sum |P_B| dt` (kWh)
- the objective split into energy, `C_B` and terminal terms

The two ieee123 no-`C_B` cells first run before the window metric existed were
re-run (identical iterations and objective); their first logs are in
`ipopt_no_cb/logs/superseded/`.

## Result

| System | T | iterations (C_B=1e-3 -> 0) | solve s | mean power use | energy window | at power limit | throughput |
|---|---:|---|---|---|---|---|---|
| ieee123 | 6 | 37 -> 38 | 0.4 -> 0.3 | 20% -> 20% | 100% -> 100% | 0% -> 0% | 2,498 -> 2,501 kWh |
| ieee123 | 24 | 38 -> 42 | 1.2 -> 1.5 | 34% -> 34% | 100% -> 100% | 6% -> 20% | 4,066 -> 4,116 kWh |
| ieee123 | 96 | 97 -> 77 | 22.5 -> 19.5 | 34% -> 34% | 100% -> 100% | 6% -> 30% | 4,107 -> 4,120 kWh |
| med2522 | 6 | 46 -> 47 | 7.5 -> 7.9 | 20% -> 20% | 100% -> 100% | 0% -> 0% | 6,505 -> 6,510 kWh |
| med2522 | 24 | 61 -> 64 | 42.9 -> 44.1 | 34% -> 34% | 100% -> 100% | 6% -> 17% | 10,772 -> 10,892 kWh |
| med2522 | 96 | 78 -> 82 | 218.1 -> 233.2 | 34% -> 34% | 100% -> 100% | 6% -> 23% | 10,792 -> 10,939 kWh |
| large10k | 6 | 51 -> 79 | 48.5 -> 76.7 | 5% -> 21% | 30% -> 100% | 0% -> 0% | 0.48 -> 2.10 GWh |
| large10k | 24 | 55 -> 77 | 255.9 -> 306.2 | 5% -> 34% | 30% -> 100% | 0% -> 17% | 0.49 -> 3.45 GWh |
| large10k | 48 | 66 -> 87 | 563.8 -> 687.3 | 5% -> 34% | 30% -> 100% | 0% -> 31% | 0.49 -> 3.45 GWh |

## Reading

- **ieee123 and med2522: C_B = 1e-3 barely restrains the batteries.** They
  already sweep their whole energy window with it, so they are energy-limited,
  not penalty-limited. Without `C_B` the same energy is moved in a more
  bang-bang way (more periods at the power limit). Ipopt's effort is about the
  same (T=96 on ieee123 is even faster without the term).
- **large10k: C_B = 1e-3 largely suppresses the batteries.** With it they use
  5% of rated power and 30% of their energy window, moving 0.49 GWh; without it
  they move 4-7x more and use the full window. The centralized problem is also
  harder without the term: 32-55% more iterations and 20-59% more solve time.
  The large10k batteries are large (median rating 1,667 kWh against 40 and 21
  kWh on the other feeders), so a single `C_B` value weighs very differently
  across the three systems.
- So the meeting's concern holds, and more strongly than "the term helps the
  solver": the matched large10k instance used in Table II is a mostly-idle-
  battery problem because of `C_B`.

Not yet done: FilterDDP without `C_B`. There, the diagonal `2 C_B S^2 dt`
contribution to the stage Hessian (see CLAUDE.md, reduced-space section) is
exactly what vanishes, so its behaviour may differ more than Ipopt's.

## Sizing C_B (2026-09-26)

**Rule.** Away from its bounds, a battery's power in the Ipopt objective
(`c_t pbase dt P_Subs + C_B pbase^2 dt P_B^2`) follows
`P_B(t) = (c_t - λ)/(2 C_B)` kW, with `λ` near the mean price for an
energy-neutral cycle. Its peak swing is therefore
`P* = (c_max - c̄)/(2 C_B)`, independent of battery size. Asking
`P* = α P_rated` gives

    C_B = (c_max - c̄) / (2 α P_rated)

With `c_max - c̄ = 0.06 $/kWh` (the periodic profile at full resolution,
identical on all three feeders; T=6 sampling gives 0.052, so the rule uses the
unsampled value and `C_B` does not vary with `T`) and each system's median
rating:

| System | median P_rated | α at C_B = 1e-3 | C_B for α = 1 | C_B for α = 3 |
|---|---:|---:|---:|---:|
| ieee123 | 10.0 kW | 3.0 | 3.0e-3 | 1.0e-3 |
| med2522 | 5.3 kW | 5.7 | 5.7e-3 | 1.9e-3 |
| large10k | 416.7 kW | 0.072 | 7.2e-5 | 2.4e-5 |

Ratings vary about ±33% within each system, so one value per system puts every
battery between roughly 0.7α and 1.5α.

**The ceiling is energy, not power.** All batteries are 4 h (E/P = 4) with a
30-95% SOC window, i.e. 2.6 h of usable full-power energy, and the price has
about two cycles a day. Without `C_B`, every feeder moves the same share of
rated power (34%) and uses 100% of the energy window. So 34% mean power use is
the most these instances allow; more battery action would need a different
instance (more energy per kW), not a smaller `C_B`.

**Check on large10k T=24** (`CB=<value> CELLS="large10kC_1ph:24"
run_ipopt_cb_metrics.sh`, logs in `ipopt_cb_<value>/logs/`):

| C_B | α | iterations | solve s | mean power use | energy window | at power limit | throughput (MWh) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1e-3 | 0.07 | 55 | 255.9 | 4.9% | 30% | 0% | 487 |
| 7.2e-5 | 1 | 73 | 351.0 | 22.0% | 86% | 0% | 2,205 |
| 2.4e-5 | 3 | 66 | 283.4 | 33.7% | 100% | 7% | 3,436 |
| 0 | inf | 77 | 306.2 | 33.8% | 100% | 17% | 3,450 |

`α = 3` moves 99.6% of the no-penalty energy, fills the window, and keeps the
profile smooth (7% of battery-periods at the power limit against 17% without
the term), the same behaviour ieee123 shows at its current `α = 3` (99% of the
energy, 6% at the limit). `α = 1` is not enough: 64% of the energy, 86% of
the window. Ipopt takes 66 iterations at `α = 3`, between the 55 with the old
value and 77 without the term.

**Recommendation.** `α = 3`, one `C_B` per system: ieee123 1.0e-3 (unchanged),
large10k 2.4e-5, med2522 1.9e-3 for a strictly uniform rule, or 1e-3 (`α` =
5.7) kept as is, which already gives 99% of its no-penalty energy. Keeping
med2522 at 1e-3 means only large10k has to be re-run.

## Table II re-run at the per-system C_B (2026-09-26)

Adopted by the user: one `C_B` per system from the median rating, alpha = 3
(`CB_BY_SYSTEM` in `terminal_soc_penalty.jl`, `REDUCED_CB=system`). All nine
Table II cells re-run with Ipopt (`run_ipopt_cb_metrics.sh` with `CB=system`,
logs `ipopt_cb_system/logs/`) and FilterDDP in Table II's configuration
(diagonal Hessian, exact rewrites, factor-backed policy, blocked solve; eight
Julia threads for med2522/large10k, one for ieee123; logs
`kkt_ordering/fullrun_blocked/*_cbsys_*`). Every FilterDDP run reaches
near-optimality; background load below 1.1 other cores throughout.

| System | T | C_B | Ipopt it / s | FilterDDP it / s | ratio | mean power use | energy window | at power limit |
|---|---:|---|---|---|---:|---:|---:|---:|
| ieee123 | 6 | 1.0e-3 | 37 / 0.32 | 67 / 20.7 | 64.2x | 20% | 100% | 0% |
| ieee123 | 24 | 1.0e-3 | 38 / 1.25 | 85 / 32.1 | 25.6x | 34% | 100% | 6% |
| ieee123 | 96 | 1.0e-3 | 97 / 22.8 | 120 / 95.3 | 4.2x | 34% | 100% | 6% |
| med2522 | 6 | 1.9e-3 | 44 / 6.9 | 68 / 102.0 | 14.8x | 20% | 100% | 0% |
| med2522 | 24 | 1.9e-3 | 60 / 41.3 | 77 / 371.4 | 9.0x | 33% | 100% | 4% |
| med2522 | 96 | 1.9e-3 | 78 / 217.5 | 96 / 1751.2 | 8.1x | 34% | 100% | 4% |
| large10k | 6 | 2.4e-5 | 53 / 51.7 | 101 / 1222.9 | 23.6x | 21% | 100% | 0% |
| large10k | 24 | 2.4e-5 | 66 / 275.1 | 98 / 3983.5 | 14.5x | 34% | 100% | 7% |
| large10k | 48 | 2.4e-5 | 75 / 618.8 | 85 / 6360.8 | 10.3x | 34% | 100% | 6% |

(Battery columns are from the centralized solutions.) ieee123 is the same
problem as before (`C_B` unchanged) and reproduces its earlier Ipopt objective
exactly. The batteries now behave the same way on all three feeders. At
large10k Ipopt needs more iterations than at the old `C_B` (53/66/75 against
51/55/66) and FilterDDP about the same (101/98/85 against 103/96/81), so the
large10k ratios improve (23.6/14.5/10.3x against 28.0/15.7/11.0x at the old
`C_B` with the same solver settings). The paper's Table II, settings table
and text now use these numbers.

## FilterDDP without C_B (2026-09-26)

FilterDDP with `C_B = 0` on all nine Table II cells, in Table II's
configuration, stopped against the `C_B = 0` Ipopt runs above (`CB=0` in
`run_blocked_solve_fullrun.sh`, logs `fullrun_blocked/*_cb0_*`). Background
load below 1.5 other cores. Seconds (iterations); "reg" counts iterations where
FilterDDP's KKT regularization had to fire (it never fires at the per-system
`C_B`).

| System | T | Ipopt, C_B = 0 | FilterDDP, per-system C_B | FilterDDP, C_B = 0 | ratio, C_B = 0 | reg |
|---|---:|---|---|---|---:|---:|
| ieee123 | 6 | 0.34 (38) | 20.7 (67) | 21.6 (70) | 64.0x | 0 |
| ieee123 | 24 | 1.50 (42) | 32.1 (85) | 33.4 (88) | 22.2x | 0 |
| ieee123 | 96 | 19.5 (77) | 95.3 (120) | 139.0 (132) | 7.1x | 34 |
| med2522 | 6 | 7.9 (47) | 102.0 (68) | 117.1 (78) | 14.9x | 0 |
| med2522 | 24 | 44.1 (64) | 371.4 (77) | 397.3 (83) | 9.0x | 0 |
| med2522 | 96 | 233.2 (82) | 1751.2 (96) | **fails** at 93: primal 1.2e-4 > 1e-5, objective within 0.03% | -- | 4 |
| large10k | 6 | 76.7 (79) | 1222.9 (101) | **fails** at 2: line search | -- | 0 |
| large10k | 24 | 306.2 (77) | 3983.5 (98) | 4511.1 (107) | 14.7x | 5 |
| large10k | 48 | 687.3 (87) | 6360.8 (85) | 7197.5 (95) | 10.5x | 4 |

The large10k `T=6` failure was re-run with every speedup off (no blocked
solve, no factor-backed policy, one Julia thread; `*_plain_cb0_*`) and fails
identically at iteration 2, so it belongs to the formulation, not the
implementation.

**Reading.** The penalty matters more to FilterDDP than to Ipopt. Without it
Ipopt solves all nine cells, in 15% less to 48% more time than at the
per-system `C_B`. FilterDDP fails two of nine:
large10k `T=6` at once, and med2522 `T=96` stalling just short of the primal
threshold. Where it succeeds it is 4-46% slower, with 3-12 more iterations,
and its KKT regularization fires on five cells. With `C_B = 0` the only
curvature on the battery powers is the barrier term plus the value-function
term, so the stage KKT is close to singular in those directions; `C_B` supplies
a well-conditioned diagonal. The per-system `C_B` keeps that benefit while
letting the batteries use their full energy window.

## The failure belongs to the diagonal Hessian (2026-09-26)

Question (user): are diagonalization and `C_B` correlated? large10k `T=6`
re-run with the EXACT stage Hessian, otherwise in Table II's configuration
(`fullrun_blocked/*_T6_exact_tableII_jt8_*`):

| large10k T=6 | Hessian | C_B | iterations to near-opt. | time (s) | s/iteration | KKT regularization |
|---|---|---|---:|---:|---:|---:|
| Table II | diagonal | 2.4e-5 | 101 | 1222.9 | 12.1 | never |
| control | exact | 2.4e-5 | 97 | 5633.1 | 58.1 | never |
| no C_B | diagonal | 0 | fails at 2 (line search) | -- | -- | -- |
| no C_B | exact | 0 | 108 | 7690.2 | 71.2 | never |

With the exact Hessian `C_B = 0` converges without any regularization, so the
`C_B = 0` failure is a property of the diagonal approximation. The diagonal
approximation keeps only `dt^2 diag(V_xx)` of the value-function curvature on
the battery powers. With `C_B > 0` the constant `2 C_B S^2 dt` on those
diagonals keeps the stage well conditioned; with `C_B = 0` the dropped
off-diagonal `V_xx` coupling is what held it together. So diagonalizing is
safe only while `C_B` (or the barrier) dominates the battery-power block --
the diagonal-dominance measurement (`run_hessian_dominance.sh`, built but not
yet run, user's decision 2026-09-26) would quantify that. As before, the exact
Hessian needs fewer iterations (97 vs 101) but costs about 4.8x more per
iteration here.
