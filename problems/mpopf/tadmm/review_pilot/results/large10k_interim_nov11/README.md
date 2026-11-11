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
