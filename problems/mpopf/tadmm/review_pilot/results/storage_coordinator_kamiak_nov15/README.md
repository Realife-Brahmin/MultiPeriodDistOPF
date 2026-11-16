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

## Timed runs with the 2 -> 0.25 schedule (`tadmm_sc_nov15_v6`)

Whole 64-core nodes, two runs per case (a / b). Every centralized, converged
and 0.1% dispatch passes the replay. "0.1%" is the time to a dispatch within
0.1% of the centralized optimum, for both solvers.

| case | centralized, converged | tADMM, converged (iterations) | centralized, 0.1% | tADMM repaired, 0.1% (iterations) |
|---|---|---|---|---|
| ieee123 T24 | 0.433 / 0.435 | 0.425 / 0.426 (18) | 0.326 / 0.327 | 0.225 / 0.229 (9) |
| med2522 T24 | 14.5 / 14.1 | 8.12 / 7.96 (23) | 11.4 / 11.0 | 3.01 / 2.95 (8) |
| large10k T24 | 31.6 / 31.6 | 16.0 / 16.1 (20) | 23.7 / 23.7 | 7.65 / 7.40 (9) |
| large10k T48 | 85.4 / 80.9 | 61.1 / 59.6 (54) | 67.2 / 63.6 | 9.89 / 9.66 (8) |

Against the 2 -> 0.0625 schedule, convergence is faster at T = 24 (large10k 39
-> 20 iterations) and slower at T = 48 (43 -> 54). The cost is within 0.003% by
iteration 20 in both; the remaining iterations shrink the primal residual from
8e-4 to 1e-5. `tadmm_sc_nov15_v7_explore` tests over-relaxation and residual
balancing for that tail. The T = 96 runs were still queued.

## Tail of the convergence: what was tried (`tadmm_sc_nov15_v7_explore`, `tadmm_sc_nov15_v8`)

Iterations to convergence (primal rms 1e-5, relative dual 1e-3, strict repair).
v7 jobs ran on shared nodes: counts only. The point within 0.1% of the
centralized optimum is at iteration 7 to 10 under every setting.

| setting | ieee123 T24 | med2522 T24 | large10k T24 | large10k T48 |
|---|---|---|---|---|
| 2 -> 0.0625 over 20, held ("base") | 21 | 27 | 39 | 43 |
| base + residual balancing | 21 | 32 | 24 | 28 |
| 2 -> 0.25 over 10, held | 18 | 23 | 20 | 54 |
| 2 -> 0.25 over 10 + balancing | - | 20 | 20 | 118 |
| 2 -> 0.25 over 10 + over-relaxation 1.6 | 18 | 22 | 105 | 72 |
| 2 -> 0.25 over 10 + balancing + over-relaxation | - | 21 | 20 | running |
| base + over-relaxation 1.6 | - | - | - | 54 |

No setting is best everywhere. Over-relaxation is dropped. The two candidates
are "2 -> 0.25 over 10" (best at T = 24, poor at T = 48) and "base + balancing"
(24 to 32 everywhere tested). These were tuned on the test cases themselves.

Timed whole-node runs of "base + balancing" (`tadmm_sc_nov15_v8`), seconds,
centralized / tADMM: ieee123 T24 converged 0.429 / 0.491 (21), 0.1% 0.323 /
0.244 (10); med2522 T24 converged 14.1 / 12.3 (32), 0.1% 11.0 / 3.29 (9);
large10k T24 converged 33.1 / 20.9 (24), 0.1% 24.8 / 7.45 (9). All pass the replay.

T = 96: whole 128-core nodes were not expected before 2026-11-17, so two
variants were added: 98 cores on a shared 128-core node (`*_shared98`, times
indicative) and 48 workers, two periods each, on a whole 64-core node (`*_w48`).
The ieee123 T = 96 instance was not on Kamiak; the copy exported on the 309 PC
(sha256 69cc566b...) was uploaded to `~/kamiak_jobs/tadmm_sc_inputs_extra`.

## T = 48 and T = 96 with "base + balancing" (`tadmm_sc_nov15_v8`)

Seconds, centralized / tADMM; iterations in brackets. "0.1%" is the time to a
dispatch within 0.1% of the centralized optimum, for both solvers. Every
centralized, converged and 0.1% dispatch below passes the chronological replay.

| case | node, workers | converged | ratio | 0.1% | ratio |
|---|---|---|---|---|---|
| large10k T48 | whole 64-core, 48 | 80.9 / 31.4 (28) | 2.6 | 63.6 / 10.7 (9) | 6.0 |
| ieee123 T96 | whole 64-core, 48 | 1.88 / 1.25 (23) | 1.5 | 1.44 / 0.51 (10) | 2.8 |
| med2522 T96 | whole 64-core, 48 | 71.0 / 21.8 (26) | 3.3 | 53.6 / 6.85 (8) | 7.8 |
| large10k T96 | whole 64-core, 48 | 176.7 / 58.7 (29) | 3.0 | 139.6 / 21.7 (10) | 6.4 |
| ieee123 T96 | shared 128-core, 96 | 1.44 / 1.03 (23) | 1.4 | 1.11 / 0.44 (10) | 2.5 |
| med2522 T96 | shared 128-core, 96 | 60.2 / 14.0 (26) | 4.3 | 45.4 / 3.70 (8) | 12.3 |
| large10k T96 | shared 128-core, 96 | 165.5 / 36.9 (29) | 4.5 | 131.1 / 13.7 (10) | 9.6 |

One run each. The 48-worker rows put two periods on each worker (Xeon 6438Y+);
the 96-worker rows had 98 of a 128-core Xeon 6760P node with other users' jobs
on the rest, so their times are indicative. Iterations do not depend on the
number of workers. The whole-128-core T = 96 jobs (31818612-14) were still
queued, expected 2026-11-17.

Iterations to convergence with this setting over everything run: 21 to 32 at
T = 24, 28 at T = 48, 23 to 29 at T = 96; to the 0.1% point, 8 to 10.

Peak job memory (Slurm MaxRSS, whole job): large10k T96 149 GB with 96 workers
and 105 GB with 48; med2522 T96 78 and 47 GB; ieee123 T96 62 and 35 GB. Most of
it is one Julia process per worker; the centralized solve's own peak was not
separated in these jobs.

The last exploratory job (2 -> 0.25 with balancing and over-relaxation,
large10k T48) did not converge in 150 iterations. The three T = 96 jobs of the
2 -> 0.25 schedule in `tadmm_sc_nov15_v6` were cancelled unstarted.

## Lean driver and measured memory (`tadmm_sc_nov15_v9_lean`)

`storage_coordinator_lean.jl`: each worker holds only its own periods' data and
one model per period; the repair is attempted only when the residual test is
met, and every attempt is inside the tADMM clock. Memory is the sum of Pss over
the process tree, sampled every 2 s (`sample_memory.py`); shared pages count
once. Whole 64-core nodes, "base + balancing", one run each. Seconds and GB.

| case | workers | centralized: time, memory | tADMM: time (iterations, repair attempts), memory |
|---|---|---|---|
| ieee123 T96 | 48 | 1.90, 0.74 | 1.14 (23, 1), 21.4 |
| med2522 T24 | 24 | 14.6, 1.17 | 30.7 (32, 11), 13.3 |
| med2522 T96 | 48 | 70.9, 2.77 | 25.9 (26, 4), 28.4 |
| large10k T24 | 24 | 33.0, 2.79 | 19.0 (24, 1), 19.1 |
| large10k T96 | 48 | 176.9, 9.48 | 56.3 (29, 1), 48.3 |

Every centralized and converged dispatch passes the replay.

- **Memory.** tADMM needs 5 to 30 times the centralized memory. About 0.45 GB
  per worker is the Julia process itself (ieee123: 21.4 GB for 48 workers of a
  tiny problem). The earlier Slurm MaxRSS figures were roughly twice the Pss:
  med2522 T96 with 48 workers read 47 GB there and 33.2 GB here with the same
  (diagnostic) driver, job 31819485. The lean driver brings that to 28.5 GB.
- **Failed repair attempts are expensive on med2522.** At T24 the residual test
  is met from iteration 22 but the repair is strict only at 32; the ten failed
  attempts run long (Ipopt stops at its acceptable level after 100+ iterations)
  and tADMM ends at 30.7 s against 14.6 s centralized. The diagnostic driver
  reported 12.3 s for the same 32 iterations because it charged only the
  accepted repair. Its converged times for med2522 are therefore optimistic.
  large10k and ieee123 accept the first attempt and are unaffected.

The diagnostic "before" run for large10k T96 was still running.

## Memory before and after, same sampler

large10k T96, 48 workers: 89.2 GB with the diagnostic driver (job 31819486),
48.3 GB with the lean one. med2522 T96, 48 workers: 33.2 GB and 28.5 GB.
Centralized: 9.5 GB and 2.8 GB.

## Correction: why the repair is refused on med2522

Earlier entries say that with the batteries fixed, period 1 of med2522 is
"infeasible by 3e-7". That figure is Ipopt's largest single constraint
violation at the point where it gives up, and it is spread evenly over 5,314
equalities (the root-voltage equation and every voltage-drop equation). The
shortfall itself is larger: the period solves strictly only when the voltage
limits are relaxed by between 1e-4 and 3e-4 p.u. squared (`diagnose_repair.jl`,
iteration 24, where the two battery schedules differ by 2e-4 p.u. in that
period). Relaxing them by 1e-6, 1e-5 or 3e-5, or using Ipopt's default bound
relaxation, does not help. Period 1 has voltages at both limits and needs its
batteries where its own solve put them. The strict rule was the right response;
the size quoted for the infeasibility was wrong.

## Pinned repair (`tadmm_sc_nov15_v10_lean`)

When the repair is refused in a few periods, those periods keep the battery
powers their own solve asked for, the battery step is re-solved with them
pinned, and the repair is run once more with the resulting schedule. On the
309 PC med2522 T24 then converges at iteration 22 with one pinned period, in
half the time of the refused-attempt version; the accepted dispatch is 0.012%
above the centralized optimum (3e-4% before) with violation 7e-10. Timed runs
of all cases were queued when this was written.

## Timed runs with the pinned repair (`tadmm_sc_nov15_v10_lean`)

Lean driver, every repair attempt inside the clock, "base + balancing", whole
64-core nodes, one run each. Seconds and GB (Pss over the process tree). Every
centralized and converged dispatch passes the replay.

| case | workers | centralized: time, memory | tADMM: time (iterations), memory | ratio |
|---|---|---|---|---|
| ieee123 T24 | 24 | 0.435, - | 0.482 (21), 11.0 | 0.90 |
| med2522 T24 | 24 | 14.5, 1.09 | 9.66 (22, one pinned period), 13.4 | 1.50 |
| large10k T24 | 24 | 31.6, 2.78 | 19.2 (24), 19.3 | 1.64 |
| large10k T48 | 48 | 85.6, 4.95 | 37.1 (28), 36.7 | 2.30 |
| ieee123 T96 | 48 | 1.88, 0.83 | 1.13 (23), 21.5 | 1.66 |
| med2522 T96 | 48 | 67.3, 2.92 | 18.6 (23, pinned), 28.1 | 3.61 |
| large10k T96 | 48 | 219.7, 9.85 | 68.4 (29), 49.6 | 3.21 |

med2522 T24 went from 30.7 s (ten refused repairs) to 9.66 s. The centralized
large10k T96 solve took 219.7 s here against 176.9 s and 176.7 s in the two
earlier runs on other nodes of the same type; the tADMM time moved the same way
(68.4 s against 56.3 s), so the ratio is unchanged.

T = 144 and T = 288 (48 workers) and the threaded driver
(`tadmm_sc_nov15_v11_threaded`) were queued when this was written. Instances for
T = 144 and 288 were exported on the 309 PC and uploaded; their sha256 prefixes
are 4455a857 / 3d920ab1 / 32d26723 (T144) and ff53ad25 / cafc8dca / 9a49bdaf
(T288) for ieee123 / med2522 / large10k.

Thread probe on the 309 PC (`probe_threads.jl`, `storage_coordinator_threaded.jl`,
8 threads): same iterations, gap and violation as the process version on ieee123
and med2522 T24; the whole process peaks at 1.2 GB and 2.4 GB, centralized phase
included; about 28% slower than 8 worker processes on med2522.

## T = 144 and T = 288 (`tadmm_sc_nov15_v10_lean`, 48 workers, whole 64-core nodes)

Lean driver with the pinned repair; seconds and GB (Pss). With 48 workers each
worker solves three (T = 144) or six (T = 288) periods one after another.

| case | centralized: time, memory | tADMM: time (iterations), memory | ratio |
|---|---|---|---|
| ieee123 T144 | 3.41, 0.91 | 1.62 (23), 23.2 (whole job) | 2.11 |
| med2522 T144 | 111.8, 3.83 | 30.1 (23), 32.2 | 3.71 |
| large10k T144 | 267.9, 14.8 | 85.6 (29), 63.7 | 3.13 |
| ieee123 T288 | 6.44, 1.07 | 3.29 (23), 25.2 | 1.96 |
| med2522 T288 | 257.2, 7.10 | 59.7 (23), 40.9 | 4.31 |

Centralized and tADMM dispatches pass the replay for the first four rows; for
med2522 T288 the centralized one had passed and the tADMM replay was still
running when the logs were copied. large10k T288 was still running. Iterations
stay at 23 to 29 at every horizon from 24 to 288.

## Conic baseline on Kamiak (`tadmm_sc_nov15_v12_conic`)

Ipopt/MA57 and Clarabel in the same job, one thread each, on a Xeon 6438Y+
node shared with the other conic jobs (4 cores per job). Seconds.

| case | Ipopt/MA57 | Clarabel |
|---|---|---|
| ieee123 T24 | 0.45 | 0.64 |
| ieee123 T96 | 2.02 | 3.74, almost optimal |
| med2522 T24 | 15.1 | 21.6, almost optimal, violation 1.1e-6 |
| med2522 T96 | 74.7 | 107.1, almost optimal |
| large10k T24 | 35.2 | 110.2, almost optimal, violation 1.1e-3, cone gap 8e-2 |

large10k T48 and T96 were still running. See also
`../conic_baseline_nov15/README.md` for the 309 PC runs with Gurobi.

## Binding constraints and relaxation gap (`tadmm_sc_nov15_v12_conic/results/binding_stats.csv`)

`binding_stats.jl` on every saved centralized and tADMM dispatch of
`tadmm_sc_nov15_v10_lean`, threshold 1e-6. Voltage limits bind on med2522 only
(10 lower and 8 upper bus-periods of 60,504 at T24; 142 and 97 of 726,048 at
T288); ieee123 has one bus per period at the upper limit; large10k has none.
Batteries sit at a power limit in 21 to 35% of battery-periods. The largest cone
gap is 3e-5 p.u. squared and the fictitious losses are below 1.2e-7 of the real
losses in every case.

## Start computation (309 PC, `time_start_computation.jl`)

The power-flow start for all periods in turn: 0.32 s (ieee123 T96), 7.2 s
(med2522 T96), 39.6 s (large10k T144); per period 0.003, 0.075 and 0.27 s. It
is inside the model build and outside both solvers' clocks.

## Evening of 2026-11-15

- **large10k T288, 48 workers** (`tadmm_sc_nov15_v10_lean`, job 31819922):
  centralized 618.9 s and 28.4 GB; tADMM 248.6 s (29 iterations, first repair
  accepted) and 108.2 GB; ratio 2.49. Both dispatches pass the replay. med2522
  T288's tADMM dispatch also passed (59.7 s against 257.2 s).
- **med2522 T96, one worker per period** on 98 cores of a shared 128-core node
  (job 31820119, indicative): 10.5 s against 58.5 s, 46.8 GB against 2.7 GB,
  replay passed.
- **Conic baseline, large10k** (`tadmm_sc_nov15_v12_conic`): Ipopt/MA57 against
  Clarabel 88.2 / 264.1 s at T48 and 189.4 / 556.1 s at T96; Clarabel stops at
  "almost optimal" with violations of 8e-3 and 6e-4 and cone gaps near 0.1.
- **Threads** (`tadmm_sc_nov15_v11_threaded`, `tadmm_sc_nov15_v13_threaded_gc`):
  same iterations and answers as the process version and about a tenth of the
  memory (ieee123 T24: 1.2 GB against 11.0 GB; med2522 T24: 2.7 GB against
  13.4 GB), but far slower: 3.5 s against 0.48 s on ieee123 T24 and 101 s
  against 9.7 s on med2522 T24, with 24 threads. A serial OpenBLAS in place of
  the pthread one, and running the Julia collector only between iterations,
  changed nothing (3.0 to 3.1 s). Slurm accounting shows where the time goes:
  med2522 T24 used 61 min of system CPU time against 12 min of user time (the
  process version: under 1 min against 14 min). The threads contend in the
  kernel; memory mapping is the suspect. Tests with glibc allocator settings,
  tcmalloc and MUMPS were queued.

- **large10k T96, one worker per period** on 98 cores of a shared 128-core
  node (job 31820120, indicative): 38.4 s against 172.7 s (ratio 4.5), 70.6 GB
  against 9.5 GB, both dispatches pass the replay.
- **Repeat of med2522 T96, 48 workers** (job 31820123): 18.3 s against 67.2 s;
  the first run was 18.6 s against 67.3 s.
- **MOSEK** against Ipopt/MA57 on Kamiak, one thread: ieee123 T24 0.83 / 0.44 s,
  ieee123 T96 4.57 / 1.98 s, med2522 T24 24.9 / 14.8 s; MOSEK reports optimal
  and the objectives agree to 2e-7 relative. The larger cases were running.
