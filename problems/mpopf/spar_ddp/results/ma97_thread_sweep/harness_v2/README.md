# Bounded MA97 OpenMP sweep

Authorized after proper_logging_campaign_oct16 completed. large10k T96 only,
OMP1/4/16/64, BLAS1 and Julia1; existing matched AMD ordering, physical start,
full warm-up, proper original-constraint/per-iteration/internal timing logging
and OpenDSS validation. Four sequential jobs on Gold6438Y+64-core nodes,300GiB,
full CPU allocation and occupancy guard. Existing campaign is not modified.
No paper updates are authorized while the cloud workspace overhauls it.

Uses an independent copied frozen repo. The only numerical execution wrapper
sets BLAS1 and verifies the loaded MA97 library's omp_get_max_threads equals
the requested count before including the unchanged centralized runner.
Record model/preprocess/optimizer/full wall times separately. GNU time reports
process maximum RSS; cgroup memory logs are retained by the established runner.
Correction: original job31251633 exited127 before Julia because /usr/bin/time
is absent on this cluster. v2 removes that optional executable and uses Julia
Sys.maxrss plus the unchanged cgroup log. No numerical work occurred in v1;
its logs and blocked dependent job IDs are preserved.
One observation per setting initially; repeat only if consequential noise makes
the choice unclear. Allocation may move among matching nodes, so report CPU
identity and masks. Do not infer universal best threading from one case.
