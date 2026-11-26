# MA97 OpenMP sweep: retain one thread for the tested case

Completed2026-10-16 Pacific. Four sequential independent jobs31251874--31251877
on the same Gold6438Y+64-core nodecn157; all64 CPUs requested,300GiB,
occupancy guard passed. Same matched periodic large10kT96 input, costs,
physical initial values, MA97OpenMP binary, AMD ordering, tol1e-6 and full
optimizer warm-up. BLAS1,Julia1; only MA97 OpenMP setting varied1/4/16/64.
Loaded-library omp_get_max_threads returned the requested values. This verifies
the runtime limit, not continuous utilization of every requested thread.

## Result

Optimizer seconds were93.711,129.649,167.734,191.846 at1,4,16,64 threads.
The multithreaded settings are substantially slower, so retain OMP1/BLAS1 as
the centralized default for this case. One observation per setting; large
gaps do not warrant another campaign merely to decide the current default.
Do not extrapolate a universal MA97 threading rule to other matrices/BLAS
implementations or hardware.

The internal timing explains the loss. Numerical factorization wall time is
30.625s at1 thread and improves to20.515s at4, but triangular backsolve time
grows from29.247 to75.496s. At16/64 threads, factorization takes25.810/49.605s
and backsolve103.713/106.182s. Thus the four-thread factorization benefit is
overwhelmed by solve overhead. This is measured category behavior; the exact
lower-level cause (synchronization, memory behavior, runtime policies) was not
isolated. Timing categories are nested/overlapping and must not all be summed.
These are directly reported factorization timings, not subtraction estimates.

## Solution and quality checks

All four jobs LOCALLY_SOLVED,18 iterations, objective2704806.123823306058.
All recorded iteration/objective/callback-primal/callback-dual strings match
the one-thread trajectory exactly, apart from their timing. Final original
primal is9.84e-13 at printed precision; callback primal5.928169e-8, dual
1.861526e-11, complementarity3.2054874829955406e-7. Do not confuse callback
and original residual definitions. Every dispatch passed OpenDSS with zero
reported lower/upper voltage violations. Logging validators and completion
markers passed; full stdout/stderr inspected for runtime errors.

Posthoc shared-quality time at objective gap0.5% and original primal<=1e-5:
54.298492,76.286688,97.993307,109.374302s. At0.1%:
65.820465,92.246736,118.901139,134.022926s. Final one-thread objective is
the common reference. Original primal values retain Ipopt table precision.

Process historical peak RSS:10.2543,10.3041,10.4124,11.4265GiB. These are
Sys.maxrss, not node-wide simultaneous memory. Model-build times remain
93.760,94.153,94.301,93.899s. Central runner wall times (including setup and
warm-up, excluding some wrapper/module startup and subsequent OpenDSS) are
296.195,366.111,444.751,499.139s. Slurm allocation durations are6m16,7m26,
8m45,9m39 and additionally include validation/launch overhead.

## Evidence and instrumentation caveats

comparison.csv/json contain generated scalar comparisons; summarize.py
reproduces them from evidence_v2/logs. Raw text_evidence.tgz and unpacked logs
preserve internal Ipopt timing, original-constraint iteration tables,
callbacks, library/source/input hashes, metadata, memory, OpenDSS and statuses.
Numerical sources are copied from the finished96-case campaign and unchanged.
Snapshot /home/aryan.jha/kamiak_jobs/ma97_thread_sweep_oct16_v2; binary dispatch
files remain on Kamiak and are excluded from the text archive.

The first attempted runner failed before Julia because /usr/bin/time is absent
(31251633,exit127); blocked dependents were cancelled. Its evidence is retained
under failed_v1. v2 uses built-in RSS instead. A SECOND instrumentation defect
was found: the wrapper's global wall_start name is overwritten by the included
central runner, so MA97_COMPLETE_WALL seconds is INVALID. Do not use that field.
It does not affect separately measured CENTRAL_IPOPT solve_time_s/build_s/wall_s
or internal Ipopt timers, and its maxrss_bytes remains valid. No optimizer rerun
is needed to recover authoritative timings. Use an isolated local wrapper
timer before future reuse. Raw flawed values remain intact for auditability.

No paper was edited: cloud workspace owns the paper overhaul. Completed
96-case campaign evidence is separately available at
ddp/results/proper_logging_campaign_oct16 on the same branch (cc420d20).
