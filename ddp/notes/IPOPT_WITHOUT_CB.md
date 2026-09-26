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
