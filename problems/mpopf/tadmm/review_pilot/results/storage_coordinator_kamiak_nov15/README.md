# Storage-coordinated tADMM on Kamiak: first measured runs (2026-11-15)

One run per case. Each job runs on one node, times the centralized solve
(Ipopt, MA57, one core) first, then tADMM with one pinned process per period,
each after a full warm-up pass. Logs, traces and OpenDSS replay summaries are
copied here; dispatch files and replay tables stay on Kamiak under
`~/kamiak_jobs/<folder>`.

Clocks: tADMM `algorithm_s` covers Update I, II and III of every iteration with
the messages to and from the workers, plus the residual test; a reported time
adds the wall time of the one repair solve whose dispatch is reported.
Per-iteration quality checks, model building and the start computation are
outside both solvers' clocks. tADMM uses T+1 cores, centralized one.

Every reported dispatch is replayed chronologically in OpenDSS
(`sequential_dss.jl`: live Storage elements, energy carried through the horizon).

## Folders

- `tadmm_sc_nov15`: first protocol. Stops after the near-optimal points, which
  are identified with the centralized cost in hindsight. Job 31817171 failed
  before solving: node cn185 lacks the system OpenBLAS that `libma57.so` links
  to; later jobs carry a private copy.
- `tadmm_sc_nov15_v2`: both solvers to convergence. tADMM stops on its own
  residuals (primal rms 1e-5, dual rms 1e-3 of the multipliers); the centralized
  0.1% point is read from the Ipopt iteration trace.
- `tadmm_sc_nov15_v3`: as v2, but a repair solve that stops at Ipopt's
  acceptable level was accepted on its measured violation.

## Results

Whole 64-core Xeon 6438Y+ node except where noted. Times in seconds.

| case | centralized, converged | tADMM, converged (iterations) | centralized, 0.1% | tADMM repaired, 0.1% (iterations) |
|---|---|---|---|---|
| ieee123 T24 | 0.434 | 0.488 (21) | 0.327 | 0.243 (10) |
| med2522 T24 | 14.2 | 9.03 (24), replay fails, see below | 11.1 | 3.34 (9) |
| large10k T24 | 33.1 | 33.9 (39) | 24.9 | 7.89 (9) |
| large10k T48 | 84.3 | 46.7 (43) | 66.3 | 10.8 (9) |

All centralized, converged and repaired 0.1% dispatches pass the replay, with
one exception.

**med2522 converged dispatch (v3) fails the replay's state-agreement test.**
With the batteries fixed at the coordinator's schedule, period 1 (voltages at
both limits, an inverter at its reactive limit) is infeasible by 3e-7 in a
squared-voltage limit; Ipopt stops at its acceptable level with a cone that is
not tight, and the OPF's line flows then differ from OpenDSS by 0.09 kW
(tolerance 0.01 kW). Batteries execute exactly and no voltage limit is violated
at the replay's 1e-4 p.u. tolerance. v3's acceptance rule is therefore
withdrawn: from `tadmm_sc_nov15_v4` on, a repair counts only if every period
solves to the strict tolerance. On this case the repair is strict at iterations
1-10, 13, 20-21, 27-36 and 52-59 (v1 trace), so the stop moves from 24 to 27.

Without repair the first-protocol runs reach the 0.1% point at iteration 63
(med2522 T24, 24.8 s), 92 (large10k T24, 85.2 s) and 78 (large10k T48, 83.7 s),
no faster than centralized, and the large10k T48 no-repair dispatch fails the
replay.

Queued when this was written (v4): med2522 T24, and ieee123, med2522 and
large10k at T = 96 on 128-core nodes.

## Later the same day

- `tadmm_sc_nov15_v4`: strict repair rule. med2522 T24 converges at iteration 27
  in 9.85 s against 14.1 s centralized; 3.38 s against 11.0 s to within 0.1% of
  the centralized optimum. All three dispatches pass the replay. Its T = 96 jobs
  were cancelled before starting and resubmitted under v6.
- `tadmm_sc_nov15_v5_explore`: penalty schedules at T = 24 on shared nodes
  (26 CPUs), so use the iteration counts, not the times. Iterations to
  convergence (med2522 / large10k): 2 -> 0.0625 over 20 then held, 27 / 39; the
  same with residual balancing, 32 / 24; 2 -> 0.25 over 10 then held, 23 / 20;
  the same with balancing, 20 / 20; 2 -> 0.125 over 14, 20 / 21. The point
  within 0.1% of the centralized optimum stays at iteration 8 or 9 throughout.
  Every converged dispatch whose replay had finished when the logs were copied
  passes it.

Chosen schedule from here on: 2 -> 0.25 `rho_unit` over 10 iterations, then
held (0.25 is the best fixed value found on ieee123), without balancing. On
ieee123 it converges at iteration 18 (T = 24) and 22 (T = 96), 309 PC.
`tadmm_sc_nov15_v6` holds the timed whole-node runs with it.
