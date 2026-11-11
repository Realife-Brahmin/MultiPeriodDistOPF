# Large10k interim horizon scaling

Profile job 31680905 completed in 13 s. The periodic exporter preserved the
24-hour day, topology, device inventory and initial battery states. T60 uses
24-minute periods and T72 uses 20-minute periods. The load and cost profiles
agree with T48 at shared sample times. The PV profile follows the existing
exporter's daylight sampling. The input manifest pins both new files.

The T60 cap20 and exact cases run on fully allocated 64-core Xeon 6438Y+
nodes. T72 cap20 and exact plus a T48 cap20 control use fully allocated
76-core Xeon 8368 CPU nodes. Each case has one worker per period and a full
MA57 centralized reference. Both local-solve policies use the exploratory
T48 penalty transfer (rho0=15000, factor5, floor24). The T48 control enables
a same-CPU comparison with T72. T60 compares to the existing 6438Y+ T48.

The run saves `central_trace.csv`, the full Ipopt warmup and timed logs,
central barrier iteration counts, solver and wall clocks, tADMM round and
worker traces, and explicit physical and memory audits. Current T24 references
record 25 centralized iterations and 32.5–33.8 s wall. The 13 archived
T48 runs in the tuning directory all record 30 and 86.6–88.7 s wall.
The extra five iterations alone do not explain the elapsed-time increase.

The original source (054b39c1) is retained at
`/home/aryan.jha/kamiak_jobs/tadmm_review_nov11_large10k_interim`.
The 76-core jobs 31680908 and 31680910 exited before solving because the
older frozen inner wrapper hard-coded a 6438Y+ check. Their audits had no
completed solve to inspect. Obsolete pending job 31680914 was cancelled
after the replacement was submitted. No T60 job was interrupted.

The corrected 76-core source (9ea64e7e) is a separate frozen copy at
`/home/aryan.jha/kamiak_jobs/tadmm_review_nov11_large10k_interim_ice`.
It reuses the hashed profiles and checks the actual 8368 hardware.
`ice-submissions.tsv` lists the corrected jobs. At verification both T60
jobs were running and all three corrected 76-core jobs were queued.

## Same-node T48 to T72 comparison

Jobs 31680924 and 31680922 used cap20 on cn130 with the same penalty
schedule. Centralized went from 93.2275 s (30 Ipopt iterations) to
136.6908 s (29 iterations). tADMM core time went from 133.2674 to
280.1445 s, with 123 versus 224 rounds and 1.0835 versus 1.2506 s per
round on average. More rounds account for most of the tADMM time increase.
At the observed T72 round cost, about 109 rounds would tie centralized
and about 55 would reach 2x on core time.

The archived `paper_cost_1e-3` milestone uses owned-copy RMS at most
1e-3 and relative objective gap at most 0.5%. It occurs at round 45 and
55.5853 s core for T48, and round 47 and 69.6200 s for T72. The maximum
original-model violations at these points are 0.04973 and 0.06083 in the
model's per-unit convention. These milestone times describe the older
stopping rule. A passive chronological replay of the hashed saved
dispatches, job 31681180, measures the battery and cone discrepancies
without rerunning either optimizer.
