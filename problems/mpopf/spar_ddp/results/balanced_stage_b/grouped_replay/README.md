# Grouped Phase B replay, November 4, 2026

The (T=96), 1,020-battery pilot compares the existing balanced doubling
schedule with the grouped schedule at (b=2). Both use the original frozen
Stage A capture and the same dense in-place merge, SharedArray layout,
63 worker processes, one coordinator, one BLAS thread per process, and
task-index round-robin dispatch. Job 31554338 ran on Gold 6438Y+ node cn155
with 64 allocated CPUs and 300 GiB. The startup guard found no other
running job on that node. Julia 1.12.7 came from the historical runner's
configured executable.

## Timing

- Balanced: **1.319297 s** median (range 1.311554–1.335228 s), 368 merges
  in seven rounds.
- Grouped: **1.222001 s** median (range 1.212331–1.225087 s), 256 merges
  in eight rounds.
- The grouped median is **7.37% lower** (1.080× speedup), with 30.43% fewer
  merges. The grouped schedule was faster in each of the five paired repeats.

Round medians in execution order:

- balanced: 0.1819, 0.1834, 0.1845, 0.1824, 0.1816, 0.2467, 0.1596 s. Merge counts: 48, 48, 48, 48, 48, 64, 64.
- grouped: 0.1816, 0.1745, 0.1747, 0.1726, 0.1727, 0.1720, 0.0841, 0.0874 s. Merge counts: 48, 24, 24, 24, 24, 32, 32, 48.

The first five balanced rounds assign one merge each to 48 workers.
Rounds six and seven assign 64 merges across 63 workers, with one worker
assigned two. Grouped rounds use 24–48 workers, each assigned one merge.
`passes.csv` records run order and `rounds.csv` retains the full 63-worker
assignment vector for every round, indexed in ascending worker PID order.

Each schedule has one warm-up, followed by five timed pairs in alternating
order. Phase B includes dispatch, waits, and round bookkeeping. Slot reset,
worker GC before each pass, logging, and reference checks are outside the
timer. Round clocks cover assignment and dispatch through completion.
These are frozen-replay times and are separate from optimizer benchmarks.

## Accuracy and memory

All 96 capture files passed SHA-256 verification against the archived
manifest (manifest SHA-256
`0e89abe1962bf2eb6721527c58f036463e865d7410a7ed1c48af3b064325d26a`).
Every round was checked for distinct left writes, disjoint right operands,
and correct interval adjacency and suffix coverage. All 96 gradients and
Hessians were checked after every warm and timed pass against the ordinary
dense recursion. Maximum scaled infinity-norm errors were 2.48e-15 for
gradients and 2.98e-14 for Hessians (gate <1e-11).

Both schedules used the same 2,398,648,320 bytes (2.23 GiB) of shared slots.
Lifetime peak RSS was 1.28 GiB for the largest worker and 5.92 GiB for the
coordinator. Slurm batch MaxRSS was 70.13 GiB for the combined allocation.
Those measurements cover both schedules and setup, and cannot isolate a
per-schedule peak. The worker RSS sum in the raw log includes shared mappings.

## Late frozen correctness check

Job 31554339 applied the same schedules and dense merge serially to the
existing 24-period, 1,020-battery late capture (iteration 14, barrier iteration
8). Julia 1.12.4 and the saved capture's numerical source were used. The
leaves and all four saved stage packets passed SHA-256 checks. This job
validates reassociation and policy recovery on a late state, with no timing
comparison.

- All 24 value gradients and Hessians passed. Maximum scaled error versus
  the production balanced replay was 1.29e-15, and versus ordinary recursion
  was 5.35e-11 (late gate <1e-7).
- The balanced reference reproduced the saved control/equality probes exactly.
  Across the four saved stages and two probes per stage, grouped control,
  equality and bound-multiplier direction errors were at most 4.79e-13 scaled
  and 1.40e-11 with the relative floor (gates <1e-7 and <1e-5).
- All 24 battery state and power directions passed the same gates, with a
  maximum scaled error of 7.28e-15.

## Provenance and scope

Original T96 capture and numerical repo:
`/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/{capture,repo}`.
The timed and late snapshots are
`/home/aryan.jha/kamiak_jobs/grouped_phase_b_nov04.q7hdoX` and
`/home/aryan.jha/kamiak_jobs/grouped_phase_b_late_nov04.vkDMWI`, respectively. Source hashes, stdout, stderr,
and accounting are in `evidence/`. Both jobs completed with exit code 0.
Late stderr contains package precompilation messages. The initial submission
with Slurm's exclusive flag was rejected before a job was created. The
successful submission used the full 64-CPU request and startup occupancy
guard.

The (T=288), (b=4) schedule passed symbolic checks at 708 merges in
11 rounds versus 1,392 in nine. A frozen (T=288) capture was not found in
the inspected replay locations, so no (T=288) timing was run. No optimizer
run or production integration was performed.
