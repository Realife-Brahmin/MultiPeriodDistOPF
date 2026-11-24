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

## Late evening of 2026-11-15

**One worker per period, shared 128-core Xeon 6760P nodes (indicative).**
large10k T144, 96 workers (job 31820121): 78.8 s against 372.2 s centralized
(ratio 4.7), 83.4 GB against 14.3 GB, both dispatches pass the replay.
large10k T288 with 96 workers was running.

**Repeats** (lean driver, whole 64-core nodes), tADMM / centralized seconds,
first run in brackets: med2522 T24 9.08 / 14.06 (9.66 / 14.51); large10k T24
19.22 / 31.55 (19.24 / 31.63); large10k T96 55.2 / 176.7 (68.4 / 219.7, and
56.3 / 176.9 in `tadmm_sc_nov15_v9_lean`); med2522 T96 18.3 / 67.2
(18.6 / 67.3). Iterations identical; ratios within 4%. The one slow node run
(220 s centralized) slowed both solvers alike.

**MUMPS for reference** (`mumps_*` folders in `tadmm_sc_nov15_v10_lean`): both
the centralized solve and tADMM's solves use Ipopt's MUMPS instead of MA57.
tADMM / centralized seconds: ieee123 T24 0.74 / 0.77, med2522 T24 14.9 / 23.3,
large10k T24 39.3 / 66.3. Both solvers are 1.6 to 2.1 times slower than with
MA57 and the ratios barely move (1.04, 1.56, 1.69 against 0.90, 1.50, 1.64).
Same iterations; dispatches pass the replay. The T96 runs were not finished.

**MOSEK on Kamiak** (`tadmm_sc_nov15_v12_conic`, logs `mosek-*`), one thread,
Ipopt/MA57 / MOSEK seconds: med2522 T96 73.3 / 97.7; large10k T24 34.8 / 58.6,
T48 88.6 / 135.8, T96 188.1 / 285.2. MOSEK reports optimal throughout, agrees
with Ipopt's objective to 1e-8 relative on large10k, and leaves a violation of
at most 1.5e-5. It is 1.3 to 1.7 times slower.

**Thread slowdown, more tests** (`tadmm_sc_nov15_v13_threaded_gc`, ieee123
T24, 24 threads, 0.48 s as processes): glibc allocator settings (large mmap
and trim thresholds) 2.9 s; passive OpenMP waiting 4.6 s; MUMPS instead of
MA57 7.4 s. None helps, and it is not specific to MA57. The tcmalloc test
failed at start-up (its libstdc++ is older than the one Julia needs). Cause
still unknown; the process version remains the one to use.

- **large10k T288, 96 workers** (three periods each; 98 cores of a shared
  128-core node, job 31820122, indicative): 160.0 s against 743.5 s centralized
  (ratio 4.65), 29 iterations, 125.9 GB against 30.3 GB. The centralized
  dispatch passed the replay; the tADMM replay was still running when the logs
  were copied.
- **MUMPS, med2522 T96, 48 workers** (job 31820399): 33.1 s against 110.4 s
  (ratio 3.34; with MA57 18.6 s against 67.3 s, ratio 3.61). Both dispatches
  pass the replay. large10k T96 with MUMPS was running.

- **large10k T288 with 96 workers:** the tADMM dispatch has now passed the
  replay as well.
- **MUMPS, large10k T96, 48 workers** (job 31820400): 110.9 s against 392.0 s
  centralized (ratio 3.53; with MA57 55 to 68 s against 177 to 220 s, ratio
  3.2). Both dispatches pass the replay.
- **Memory per worker** (from `memory_samples.csv`, Pss; all workers divided by
  their number): ieee123 0.43 to 0.51 GB; med2522 0.51 to 0.71 GB; large10k
  0.66 to 1.69 GB (0.64 GB with one period per worker, 1.69 GB with six),
  against 0.7 to 1.1, 1.1 to 7.1 and 2.8 to 29 GB for the centralized solve.
  The coordinator's own memory is NOT measured by these runs: it is the process
  that did the centralized solve first, and its Pss stays at that level
  (27.3 GB at large10k T288 against 28.8 GB centralized). So the tADMM totals
  quoted above include that leftover. `tadmm_sc_nov15_v14_coordinator` reruns
  five cases without the centralized solve to measure it.

## Coordinator memory (`tadmm_sc_nov15_v14_coordinator`)

Lean driver without the centralized solve (`SC_REFERENCE`), whole 64-core
nodes, so the coordinator process holds nothing left over from it. Pss, GB.

| case | workers | coordinator | all processes | all processes when the centralized solve ran first |
|---|---|---|---|---|
| ieee123 T96 | 48 | 0.44 | 21.2 | 21.5 |
| med2522 T96 | 48 | 0.50 | 26.1 | 28.1 |
| large10k T24 | 24 | 0.92 | 17.0 | 19.3 |
| large10k T96 | 48 | 1.45 | 40.6 | 49.6 |

The coordinator needs about as much as one worker (a Julia process and the
instance data). Iterations, times and gaps are the same as in the runs with
the centralized solve, and the tADMM dispatches pass the replay. large10k T288
was still running.

large10k T288 without the centralized solve (job 31820757, 48 workers):
coordinator 2.46 GB, all processes 83.1 GB (108.2 GB when the centralized solve
ran first), tADMM 207.9 s in 29 iterations, replay passed. This job has no
centralized time of its own; the matched pair is 248.6 s against 618.9 s in
`tadmm_sc_nov15_v10_lean`. No Kamiak jobs were left after this one.

## Workers over several nodes (`tadmm_sc_nov15_v15_workers`)

`start_pool!` has a Slurm backend (`SC_CLUSTER=slurm`): one `srun` step with one
task per worker, each pinned to its own CPU by `worker_slurm.sh`;
`kamiak_storage_coordinator_multinode.sh` launches the lean driver that way.
Gate (job 31822345, ieee123 T96, 96 workers on two 128-core nodes): 23
iterations, gap 2.24e-06, 0.78 s, both dispatches pass the replay; the
single-node run gives the same iterations and gap.

### T = 288 by number of workers (128-core Xeon 6760P nodes)

All through the multi-node launcher, the centralized solve skipped
(`SC_REFERENCE`) except where a centralized time is given. Kamiak does not allow
exclusive nodes, so these nodes are shared with other users' jobs. The memory
sampler sees the head node only (the workers are started by `srun`), so worker
memory is not measured in these runs.

| system | workers | nodes | periods per worker | tADMM (s) | iterations | gap | replay | job |
|---|---|---|---|---|---|---|---|---|
| large10k | 48 | 1 | 6 | 222.8 | 29 | 2.14e-06 | pass | 31822491 |
| large10k | 72 | 1 | 4 | 157.9 | 29 | 2.14e-06 | pass | 31822490 |
| large10k | 96 | 1 | 3 | 151.6 | 29 | 2.14e-06 | pass | 31822489 |
| large10k | 144 | 2 | 2 | 74.2 | 29 | 2.14e-06 | pass | 31822488 |
| med2522 | 48 | 1 | 6 | 57.2 | 23 | 9.28e-06 | pass | 31822496 |
| med2522 | 72 | 1 | 4 | 30.5 | 23 | 9.28e-06 | pass | 31822495 |
| med2522 | 96 | 1 | 3 | 27.5 | 23 | 9.33e-06 | pass | 31822494 |
| med2522 | 144 | 2 | 2 | 20.3 | 23 | 9.32e-06 | pass | 31822493 |

Centralized large10k T288 on the same node type, in the 96-worker job: 552.2 s
(383.2 s to within 0.1%), replay passed. Against it tADMM is 2.5, 3.5, 3.6 and
7.4 times faster with 48, 72, 96 and 144 workers. On one node the time stops
falling between 72 and 96 workers (157.9 to 151.6 s), while 144 workers on two
nodes, 72 on each, halve it: a node saturates at about 72 workers of this size,
and more nodes help where more workers on the same node do not. Shared nodes,
so single runs carry some noise; the saturation point is one reading, not a law.

Centralized med2522 T288 on the same node type, in the 96-worker job: 211.1 s
(144.4 s to within 0.1%), replay passed. Against it tADMM is 3.7, 6.9, 7.7 and
10.4 times faster with 48, 72, 96 and 144 workers. Here the gain flattens
earlier (30.5 to 27.5 to 20.3 s). The likely reason, not measured here: the
periods are small, so the battery step and the exchange with the coordinator
take a growing share of each iteration.

The two 288-worker jobs (three nodes, 98 CPUs each) had not started after five
hours in the queue (31822487, 31822492).

## A conic solver inside tADMM (`tadmm_sc_nov15_v16_conic`)

`SC_STAGE_SOLVER=mosek|clarabel|gurobi` (`storage_coordinator_conic_stage.jl`)
runs tADMM's own solves with a conic solver: each period's OPF with the cone
stated explicitly and all original bounds, the battery step as the same
quadratic program, the repair and the pinned repair unchanged. A conic
interior-point solve cannot start from the previous iterate, so every iteration
solves every period from scratch. Tolerances are those of the centralized conic
baseline (1e-8). The job times, on the same node, the centralized problem with
that conic solver (`conic_baseline.jl`), the centralized problem with
Ipopt/MA57, and tADMM. Gap is against the Ipopt/MA57 objective.

| case | workers | solver | centralized, that solver (s) | centralized, Ipopt/MA57 (s) | tADMM, that solver (s) | iterations | gap | violation | replay | job |
|---|---|---|---|---|---|---|---|---|---|---|
| ieee123 T24 | 24 | MOSEK | 0.82 | 0.43 | 0.52 | 21 | 5.5e-07 | 9.2e-07 | pass | 31822503 |
| med2522 T24 | 24 | MOSEK | 24.6 | 14.5 | 10.8 | 22 | 9.5e-05 | 2.3e-06 | pass | 31822504 |
| ieee123 T96 | 48 | MOSEK | 4.40 | 1.87 | 1.24 | 23 | 2.3e-06 | 7.0e-07 | pass | 31822625 |
| ieee123 T24 | 24 | Clarabel | 0.66 | 0.44 | 0.57 | 21 | 5.7e-07 | 5.1e-07 | pass | 31822629 |
| med2522 T96 | 48 | MOSEK | 93.5 | 67.2 | 19.1 | 23 | 7.3e-05 | 1.2e-06 | pass | 31822626 |

The T24 cases ran on 26 CPUs of a 64-core node, ieee123 T96 on a whole one. med2522 needed one pinned
period, as with Ipopt. With Ipopt/MA57 inside, tADMM took 9.7 s on med2522 T24,
so the solver inside changes tADMM's time little, while the centralized solve
is 1.7 times slower with MOSEK. Memory with MOSEK inside: 13.5 GB (ieee123) and
15.9 GB (med2522) for 25 processes, about 0.1 GB per process more than with Ipopt.

med2522 T96 with MOSEK inside used four pinned periods. With Ipopt/MA57 inside
it took 18.6 s, so again the solver inside matters little to tADMM (19.1 s),
and with MOSEK on both sides tADMM is 4.9 times faster than centralized.

Three cases failed because a conic solve stalled short of its tolerance in the
first iteration, and the code then stopped: large10k T24 (MOSEK, a period OPF,
`SLOW_PROGRESS`, job 31822627), large10k T96 (MOSEK, the battery step, job
31822628) and med2522 T24 (Clarabel, a period OPF, `NUMERICAL_ERROR`, job
31822630). The centralized conic times of those jobs stand: MOSEK 57.0 s
(large10k T24) and 275.0 s (large10k T96), Clarabel 20.8 s (med2522 T24,
`ALMOST_OPTIMAL`). The code now uses the stalled solve's last point inside the
iterations and prints `SC_CONIC_INEXACT`; the repair still accepts only a solve
that met its tolerance. The three cases were resubmitted at tolerance 1e-8 and
1e-6 (`CONIC_TOL`, applied to the centralized conic solve as well), jobs
31824404 to 31824409.

### Resubmitted cases (read 2026-11-16 02:52; the Kamiak login expired at 03:15)

| case | workers | solver, tolerance | centralized, that solver (s) | centralized, Ipopt/MA57 (s) | tADMM, that solver (s) | iterations | gap | violation | replay | job |
|---|---|---|---|---|---|---|---|---|---|---|
| large10k T24 | 24 | MOSEK, 1e-8 | 55.3 | 32.0 | 45.5 | 24 | 3.3e-07 | 7.6e-07 | pass | 31824404 |

large10k T24 is the case where the solver inside does matter: with Ipopt/MA57
inside tADMM took 19.2 s, with MOSEK 45.5 s, because every iteration solves
every period from scratch while Ipopt restarts from the previous iterate. So
with MOSEK on both sides tADMM is 1.2 times faster than centralized, and it is
slower than centralized Ipopt/MA57 (32.0 s). Three of the 576 period solves of
the timed pass stalled (`SLOW_PROGRESS`, periods 1 and 6) and their last point
was used; the repair solves all met the tolerance and the dispatch passes the
replay. The earlier sentence that the solver inside changes tADMM's time little
holds for ieee123 and med2522 only.

Clarabel inside tADMM on med2522 T24 at 1e-8 failed again (job 31824408): the
period-18 solve ends in `NUMERICAL_ERROR` with no point to use. Clarabel's
centralized time at 1e-6 is 20.4 s (`OPTIMAL`).

Not read before the login expired: large10k T24 with MOSEK at 1e-6 (31824405),
med2522 T24 with Clarabel at 1e-6 (31824409), large10k T96 with MOSEK at 1e-8
and 1e-6 (31824406, 31824407), and the two 288-worker jobs (31822487,
31822492), still waiting in the queue after five hours.

### Read 2026-11-16 09:20 (overnight jobs)

| case | workers | solver, tolerance | centralized, that solver (s) | centralized, Ipopt/MA57 (s) | tADMM, that solver (s) | iterations | gap | violation | replay | job |
|---|---|---|---|---|---|---|---|---|---|---|
| large10k T96 | 48 | MOSEK, 1e-8 | 277.0 | 177.2 | 96.9 | 29 | 4.4e-07 | 6.3e-06 | pass | 31824406 |
| large10k T24 | 24 | MOSEK, 1e-6 | 51.0 | 34.2 | not converged in 150 | - | - | - | - | 31824405 |
| large10k T96 | 48 | MOSEK, 1e-6 | 255.4 | 218.9 | not converged in 150 | - | - | - | - | 31824407 |
| med2522 T24 | 24 | Clarabel, 1e-6 | 20.4 | 14.5 | not converged in 150 | - | - | - | - | 31824409 |

- large10k T96, MOSEK at 1e-8 on both sides: tADMM is 2.9 times faster than
  centralized MOSEK and 1.8 times faster than centralized Ipopt/MA57 on the same
  node. With Ipopt/MA57 inside tADMM took 55 to 68 s, so MOSEK inside costs
  about 1.5 times here (2.4 times at T24). Two period solves of the timed pass
  stalled and their last point was used; the dispatch passes the replay.
- At tolerance 1e-6 the conic solves are too loose for tADMM's own stopping
  rule (primal residual 1e-5): the MOSEK runs never met the residuals in 150
  iterations (no repair was attempted), and the Clarabel run met them but every
  one of its 129 repairs was refused (violation 1.2e-04). So the conic solver
  inside needs its tight tolerance.
- The two 288-worker jobs were still waiting in the queue (over 11 hours).

### Gurobi on Kamiak

Kamiak has Gurobi 13.0.1 (`/opt/apps/gurobi/13.0.1`, module under `Licensed`)
with a token-server licence that answers from compute nodes only (the login
node gets "Server must be on the same subnet"). Use it with
`GRB_LICENSE_FILE=/opt/apps/gurobi/13.0.1/licenses/gurobi.lic`; Gurobi.jl 1.9.3
was added to `conic_env` on Kamiak.

| case | workers | solver | centralized, that solver (s) | centralized, Ipopt/MA57 (s) | tADMM, that solver (s) | iterations | gap | violation | replay | job |
|---|---|---|---|---|---|---|---|---|---|---|
| ieee123 T24 | 24 | Gurobi | 0.59 | 0.43 | 1.82 | 21 | 5.8e-07 | 3.0e-08 | pass | 31825855 |

Gurobi inside tADMM works (24 workers each took a licence token), but on this
small case it is three times slower than centralized Gurobi: 21 iterations of
from-scratch solves at about 0.08 s each, where Ipopt or MOSEK inside take about
0.5 s in total. Larger cases were submitted (jobs 31825876 to 31825880).

### Why the 288-worker jobs never started

The `kamiak` partition allows one user 180 CPUs at a time (QOS `kamiak`,
`MaxTRESPU cpu=180`), and each 288-worker job asks for 294. The `vcea`
partition has no such cap but only three 64-core nodes (192 CPUs). So 288
single-CPU workers cannot run under this account as it stands; both jobs were
cancelled. The largest worker count that fits is 144 (148 CPUs, two nodes).
One worker per period is therefore measured at T = 144 instead: jobs 31825882
to 31825887 (large10k and med2522, 48, 72 and 144 workers).

### Gurobi inside tADMM, larger cases (read 2026-11-16 09:41)

| case | workers | centralized, Gurobi (s) | centralized, Ipopt/MA57 (s) | tADMM, Gurobi inside (s) | tADMM, Ipopt inside, earlier (s) | iterations | gap | violation | replay | job |
|---|---|---|---|---|---|---|---|---|---|---|
| ieee123 T96 | 48 | 3.08 | 1.87 | 3.75 | about 1 | 23 | 2.3e-06 | 1.4e-08 | pass | 31825876 |
| med2522 T24 | 24 | 15.9 | 14.5 | 18.3 | 9.7 | 22 | 7.4e-05 | 7.1e-07 | pass | 31825877 |
| med2522 T96 | 48 | 74.8 | 71.2 | 34.8 | 18.6 | 23 | 5.9e-05 | 4.1e-06 | pass | 31825878 |
| large10k T24 | 24 | 82.7 | 31.7 | 87.1 | 19.2 | 26 | 3.1e-07 | 6.1e-06 | pass | 31825879 |
| large10k T96 | 48 | 283.5 | 186.3 | 194.3 | 55 to 68 | 32 | 4.0e-07 | 2.9e-10 | pass | 31825880 |

No stalled solves. Gurobi inside is the slowest of the three inner solvers:
about twice Ipopt or MOSEK inside on med2522 and more on ieee123. With Gurobi on
both sides tADMM wins only at med2522 T96 (2.2 times); at ieee123 T96 and
med2522 T24 it is 0.8 to 0.9 of centralized Gurobi. Centralized Gurobi itself is
close to centralized Ipopt/MA57 on med2522. On large10k T24 Gurobi inside takes
87.1 s in 26 iterations (24 with Ipopt or MOSEK inside), no better than
centralized Gurobi (82.7 s) and far behind centralized Ipopt/MA57 (31.7 s).
On large10k T96 it takes 194.3 s in 32 iterations (29 with Ipopt or MOSEK
inside): 1.5 times faster than centralized Gurobi, level with centralized
Ipopt/MA57 (186.3 s), three times slower than tADMM with Ipopt inside.

### T = 144 by number of workers (128-core nodes), first rows

| system | workers | nodes | periods per worker | tADMM (s) | iterations | gap | replay | job |
|---|---|---|---|---|---|---|---|---|
| large10k | 48 | 1 | 3 | 95.2 | 29 | 1.33e-06 | pass | 31825884 |
| large10k | 72 | 1 | 2 | 61.5 | 29 | 1.33e-06 | pass | 31825883 |
| large10k | 144 | 2 | 1 | 41.0 | 29 | 1.33e-06 | pass | 31825882 |
| med2522 | 48 | 1 | 3 | 25.3 | 23 | 9.11e-06 | pass | 31825887 |
| med2522 | 72 | 1 | 2 | 15.9 | 23 | 9.06e-06 | pass | 31825886 |
| med2522 | 144 | 2 | 1 | 9.0 | 23 | 9.14e-06 | pass | 31825885 |

Centralized large10k T144 on the same node type, in the 72-worker job: 245.9 s
(207.5 s to within 0.1%), replay passed; tADMM is 2.6, 4.0 and 6.0 times
faster with 48, 72 and 144 workers. Centralized med2522 T144, in its 72-worker
job: 91.0 s (67.3 s to within 0.1%), replay passed; tADMM is 3.6, 5.7 and 10.1
times faster. The 144-worker rows are one worker per period, on two nodes.
Shared nodes, single runs. No Kamiak jobs were left after these (10:21).

## Penalty rules compared (2026-11-17, T = 96, 48 workers, shared 128-core nodes)

Centralized skipped (`SC_REFERENCE`). `current` = 2 rho_unit decreasing to rho_unit/16 over
20 iterations, then balancing (2 / 1.5); the others start at rho_unit and use no decrease.

| rule | med2522 | large10k |
|---|---|---|
| current | 23 iterations, 15.4 s | 29 iterations, 61.9 s |
| balancing 2 / 1.5 from iteration 1 | 23 iterations, 16.5 s | not converged in 150 |
| balancing 10 / 2 (Boyd) from iteration 1 | 46 iterations, 27.1 s | not converged in 150 |
| fixed rho_unit | 52 iterations, 36.3 s | 124 iterations, 248.9 s |

Every converged run passes the OpenDSS replay. Jobs 31849402-31849409. The decreasing start is
what makes large10k converge; none of the simpler rules matches it.

## Final large10k table (128-core Xeon 6760P, centralized in the same job)

| T | workers | nodes | centralized (s) | tADMM (s) | iterations | ratio | job |
|---|---|---|---|---|---|---|---|
| 96 | 96 | 1 | 183.9 | 55.2 | 29 | 3.33 | 31849399 |
| 144 | 144 | 2 | 248.3 | 38.1 | 29 | 6.51 | 31849400 |
| 288 | 144 | 2 | 557.4 | 77.5 | 29 | 7.19 | 31851704 |

T = 96: both dispatches pass the replay; all processes 71.2 GB at peak, the coordinator holding
the centralized solve's 9.2 GB. T = 144 and 288 with 144 workers on two nodes were still queued.

### Whole-number schedules (2026-11-17, same setting; `SC_RHO_STEP`)

Start 2 rho_unit, end rho_unit/16 (five halvings), then the same balancing (2 / 1.5).

| schedule | med2522 | large10k |
|---|---|---|
| halve every 4 iterations (20) | 23 iterations, 14.1 s | 30 iterations, 49.0 s |
| halve every 2 iterations (10) | 26 iterations, 21.8 s | 44 iterations, 89.2 s |
| halve every iteration (5) | 19 iterations, 13.8 s | 77 iterations, 124.3 s |
| current (smooth, 20) | 23 iterations, 15.4 s | 29 iterations, 61.9 s |

All pass the replay (jobs 31850334-31850339; the first submission failed at start-up because the
worker files on Kamiak were older than the driver). "Halve every four iterations" matches the
smooth schedule to within one iteration, so it can replace it. The length of the decrease matters:
on large10k a shorter decrease (5 or 10 iterations) costs 1.5 to 2.7 times the iterations.

T = 144, 144 workers: both dispatches pass the replay (run separately, job 31851708, because
Kamiak's srun has no `--overlap` and the per-node sampler then stopped the job script before its
replays; fixed, the sampler is now started by the first worker task on each node). Memory, from a
rerun without the centralized solve (job 31851705, 40.1 s, 29 iterations): 47.7 GB on the head node
(72 workers and the 2.4 GB coordinator) and 45.3 GB on the other (72 workers), so 90.7 GB for 144
workers, 0.63 GB each.

T = 288, 144 workers (job 31851704, per-node sampling working): both dispatches pass the replay.
Memory in the timed pass: 85.7 GB on the head node (72 workers and the coordinator, which still
held the centralized solve's memory) and 56.4 GB on the other (72 workers, 0.78 GB each). These
three rows replace the 48-worker large10k rows at T = 96, 144 and 288 in the paper's Table V.

Penalty comparison, conclusion: the decreasing start (2 rho_unit to rho_unit/16 over 20 iterations,
then balancing) is needed. Fixed rho, or balancing from the first iteration (standard 10/2 or mild
2/1.5), costs 2-4 times the iterations or does not converge on large10k. The decrease can be applied
as halving every 4 iterations with the same result; shortening it to 10 or 5 iterations costs
1.5-2.7 times the iterations on large10k. Not yet tested: other start and end values, and other
networks (med3980, sf6937).

## New grids and lossy batteries (2026-11-18)

### med3980 and sf6937 (single-substation decks from the DDP branch; capacitors modelled)

First runs (`tadmm_sc_nov18_newgrids`, T = 96, 48 workers, current penalty rule): the residuals met
their tolerances (med3980 from iteration 28) but every repair was refused (41 and 32 attempts), and
the OpenDSS replay stopped on a bus outside the modelled network. Centralized: 184.9 s (med3980,
123 Ipopt iterations) and 318.4 s (sf6937, 124).

Diagnosis (jobs 31863298-31863300, refusals logged): the repair's Ipopt cap of 60 iterations.
med3980 with the cap at 3000: 28 iterations, one repair, accepted, 44.3 s, gap 3.3e-5. With the cap
at 60: one period hit the cap at iteration 28 ("Maximum_Iterations_Exceeded"), accepted at 31. The
first runs also lost the squared-current screening (a comma in `--export` split the setting), which
made the repairs harder still. The launcher now sets the screening and a cap of 3000 for these two
grids. The replay now skips deck lines the model leaves out and buses outside the network.

### Lossy batteries (SC_ETA, eta each way)

ieee123, T = 24, 24 workers (jobs 31863308-09): eta = 1.0 reproduces the lossless result (21
iterations, 0.44 s); eta = 0.95: centralized 0.38 s, tADMM 22 iterations, 0.42 s, gap 6.7e-7; both
dispatches pass the replay with the storage efficiencies set to 95% in OpenDSS.

### Results (2026-11-18 evening, 128-core nodes, 48 workers, T = 96)

| case | centralized (s) | tADMM (s) | iterations | repair attempts | ratio | job |
|---|---|---|---|---|---|---|
| med3980 | 140.3 | 44.1 | 28 | 1 | 3.2 | 31863440 |
| sf6937 (cap 3000, no centralized in job) | 318.4 (earlier job) | 350.0 | 72 | 6 (16 periods refused each time, 8 pinned) | 0.9 | 31863300 |
| med2522, eta = 0.95 | 60.6 | 64.2 | 43 | 17 | 0.94 | 31863441 |
| large10k, eta = 0.95 | 141.2 | 69.3 | 33 | 1 | 2.0 | 31863442 |

- sf6937: residuals were met by iteration ~54 but 16 periods refuse the repair; after six attempts
  the pinned repair succeeds at iteration 72. tADMM is not faster here.
- med2522 with losses: the residuals settle late (iteration 26 or so) and 17 repair attempts are
  needed (2 periods pinned), against 23 iterations and 1 attempt lossless. large10k with losses is
  fine (33 iterations, one attempt). Lossy replays pass on both, efficiencies set in OpenDSS.
- Replays of med3980 and sf6937 failed (voltage difference about 1.05 p.u. at some buses, for the
  centralized dispatch too): the replay placed PV injections at bus numbers instead of the deck's bus
  names. Fixed; replays re-run (job 31863542).
- Replays after the bus-name fix (job 31863542): med3980 centralized and tADMM dispatches pass,
  sf6937 tADMM dispatch passes (the diagnostic job had no centralized dispatch).

### Lossless against lossy on the same node type (2026-11-19, 128-core, 48 workers, T = 96)

| case | centralized (s) | tADMM (s) | iterations | repair attempts | ratio | job |
|---|---|---|---|---|---|---|
| med2522, lossless | 79.1 | 16.8 | 23 | 1 | 4.7 | 31874007 |
| med2522, eta 0.95, repair cap 3000 | 82.4 | 91.8 | 43 | 21 | 0.90 | 31874006 |

Losses barely change the centralized solve (79.1 to 82.4 s) but cost tADMM a factor of 5, all in
the repair: periods 1-4 refuse it at every attempt from iteration 23 ("Solved_To_Acceptable_Level",
some "Infeasible_Problem_Detected"), so the cap was not the cause. Both lossy dispatches pass the
replay with 95% efficiencies. sf6937 (lossless) fails the same way in periods 2-8 and 69-76
("Infeasible_Problem_Detected"): fixing the batteries at the agreed schedule leaves no feasible
network in periods at their voltage limits.

## Lossy table (2026-11-23/24, `tadmm_sc_nov15_v15_workers`, eta = 0.95 each way)

Timing counts the numerical optimization only (no build, worker start-up or
compilation). Both methods start from the same per-period load flow with idle
batteries. Repair: back off 3 iterations after a failed check; `SC_REFUSE_RHO=4`
with `SC_REFUSE_RHO_MAX=16`. Whole 64-core nodes; 144 workers = 2 nodes x 74 CPUs.
Every converged dispatch passes the OpenDSS replay with battery efficiency.

Ratio = centralized (Ipopt/MA57) / tADMM, both to convergence.

| System | T | Workers | Ratio | Iterations | Repair attempts |
|---|---|---|---|---|---|
| ieee123 | 24 | 24 | 0.87 | 22 | 1 |
| ieee123 | 96 | 48 | 1.87 | 24 | 1 |
| ieee123 | 96 | 96 | 1.62 | 24 | 1 |
| ieee123 | 144 | 48 | 1.89 | 24 | 1 |
| ieee123 | 144 | 144 | 1.97 | 24 | 1 |
| ieee123 | 288 | 48 | 2.03 | 23 | 1 |
| med2522 | 24 | 24 | 1.18-1.33 | 26 | 2 |
| med2522 | 96 | 48 | 3.13 | 26 | 2 |
| med2522 | 96 | 96 | 5.42 | 26 | 2 |
| med2522 | 144 | 48 | 2.65 | 29 | 3 |
| med2522 | 144 | 144 | 7.36 | 26 | 2 |
| med2522 | 288 | 48 | 1.70-1.83 | 44 | 8 |
| large10k | 24 | 24 | 1.37 | 24 | 1 |
| large10k | 48 | 48 | 1.77 | 33 | 1 |
| large10k | 96 | 96 | 3.50 | 33 | 1 |
| large10k | 144 | 144 | 4.04 | 34 | 1 |
| large10k | 288 | 144 | 6.06 | 29 | 1 |

### Max-norm repair gate (`SC_EPS_MAX`), med2522, 48 workers: rejected

| T | Gate | Iterations | Attempts | Ratio |
|---|---|---|---|---|
| 96 | RMS <= 1e-5 | 26 | 2 | 3.13 |
| 96 | RMS and max <= 1e-5 | 114 | 4 | 0.64 |
| 288 | RMS <= 1e-5 | 44 | 8 | 1.70-1.83 |
| 288 | RMS and max <= 1e-5 | 92 | 1 | 1.14 |

At T = 288, iteration 24 has RMS 4e-6 but a worst battery-period mismatch of
1.8e-4, and the worst period moves (6, 239, 195-199). The refused repairs come
from these few lagging periods, but waiting for the max norm costs far more
iterations than the refusals. Next candidate: a per-period penalty raised only
where a period lags.

### Table IX statistics, large10k

No binding voltage limits at any horizon. 101 loose cone branches per period,
identical in the centralized and tADMM solutions (a model property, not a tADMM
artifact). The maximum cone gap grows with T (1.6e-5 at T = 24, 9.4e-5 at
T = 288); fictitious loss share 1e-9 to 1.3e-8.

## Net power in the pinned battery step (`SC_PIN_NET=1`, 2026-11-24)

The pinned repair on med2522 kept failing the original-model check on the SOC
bound by 1.1-1.7e-5. Cause: with the early (undervoltage, charging) periods
pinned to their own powers, some batteries would overfill by ~1e-5 with no
earlier period to compensate, and the battery QP stayed feasible by charging
and discharging at once (up to 2.7e-5 p.u.; the throughput term costs only
1e-6). The check and OpenDSS use net power, so they saw the overfill.
`SC_PIN_NET=1` replaces each such pair by the net power with the same energy
change. Only med2522 runs ever used the pinned step; ieee123 and large10k rows
are unaffected. The 1e-9 threshold also rewrites interior-point noise in almost
every entry (harmless, ~1e-10 p.u.); 1e-7 would be cleaner.

Per-period penalties (`SC_PERIOD_RHO=f`): f = 10 never fires (lagging periods
are 4-7x the RMS mismatch), f = 5 and 3 do not help (1.78, 1.82 at T = 288).
`SC_REFUSE_PERIOD=1` alone: 41 iterations, 7 attempts, 2.09. Not needed with
`SC_PIN_NET`.

med2522 with `SC_PIN_NET=1` (all pass OpenDSS):

| T | Workers | Ratio | Iterations | Attempts | Gap |
|---|---|---|---|---|---|
| 24 | 24 | 1.50 | 23 | 1 | 1.75e-4 |
| 96 | 48 | 3.88 | 23 | 1 | 9.1e-5 |
| 96 | 96 | 6.48 | 23 | 1 | 9.9e-5 |
| 144 | 48 | 3.92 | 23 | 1 | 8.2e-5 |
| 144 | 144 | 10.18 | 23 | 1 | 9.0e-5 |
| 288 | 48 | 5.25 | 23 | 1 | 6.9e-5 |
| 24 (AC) | 24 | 2.19 | 23 | 1 | 1.75e-4 |

The gap is larger than before (~1e-5) because the repair is accepted about 20
iterations earlier with the pinned periods at their own powers; all are well
inside 0.1%.
