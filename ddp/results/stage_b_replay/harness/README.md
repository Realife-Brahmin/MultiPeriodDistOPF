# Stage B replay campaign, 2026-10-15

Purpose: isolate pairwise value-function composition from network stage work
and optimization progress. This is a diagnostic, not a new solve-time result.

Cluster snapshot: `/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15`.
Base: frozen `matched_benchmark_oct15_large96_threeway_v3/repo`.
Capture job 31187159 records all 96 real large10k Stage A outputs at the first
backward pass, before any pairwise composition; it then exits. The capture
launcher disables optimizer warm-up, because no optimizer time is reported.

Replay jobs 31187170 and 31187171 run sequentially on Gold 6438Y+ 64-core nodes.
They request all 64 CPUs per node and check for other running jobs. The two-node
allocation requests 128 CPUs but launches 96 workers plus a coordinator.
Smaller worker subsets alternate between the two nodes. Every replay uses
the same stage matrices, round-robin stage ownership, one BLAS thread per
worker, and the production pairwise round schedule. Garbage collection occurs
between replays and is held off within each replay, matching once-per-pass GC.

One-node shared-memory counts: 7,16,32,48,63. One-node explicit-fetch controls:
32,48,63. Two-node explicit-fetch counts: 32,48,63,96. Each configuration has
one untimed-for-summary warm-up replay and three measured replays. Every replay
must match the ordinary sequential value recursion at every stage within
relative infinity-norm error 1e-7. Validation and resetting inputs are outside
the Stage B timer. Warm-up, setup and reference construction remain visible
in job elapsed time but are excluded from replay time.

`REPLAY_ROUND` records round wall time, summed compute/fetch/GC time,
the busiest worker's summed task time, explicit fetch bytes, and cross-node
task counts. `scheduling_wait_est_s` is wall time minus the largest worker task
sum, clipped at zero: it is an estimate including dispatch and scheduling
overhead, not a direct measurement of every worker's barrier wait. Fetch time
includes serialization and waiting for the source worker. Concurrent summed
times must not be added to wall time. SharedArray access consumes memory
bandwidth even though explicit fetch bytes are zero.

No mathematical approximation is introduced. The capture source is diagnostic
only and must not replace the production solver: it intentionally exits after
Stage A. Source and input checksums are kept in the cluster snapshot.
