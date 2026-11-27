# Completed Kamiak checks and first benchmark repeats

All three T=6 systems completed full convergence, matching warm-up/measured
traces and independent OpenDSS checks with both MA57 and MA97. System gates
31643454 and 31643455 completed successfully. These shared-node verification
times are not benchmark results.

Backend sensitivity remains substantial: med2522 required 77 ADMM iterations
with MA57 and 1703 with MA97; large10k required 132 with MA57 and 568 with MA97.
The original frozen configuration eventually converged without adopting any
diagnostic solver-setting changes.

Benchmark array 31643458 tasks 0 and 1 completed both backends and all validation
checks on fully allocated 64-CPU Gold nodes. These are the two IEEE123 T=6
repeats, with six period workers plus one coordinator and the matched centralized
baseline allowed the same seven-CPU budget. The faster two-repeat median backend
was MA97 for both algorithms. Actual full-convergence loop medians, excluding
builds and warm-up, are 3.4420957805 s for parallel tADMM and 0.2641107235 s for
the centralized solve. tADMM loop time includes diagnostics, checkpointing,
output and real communication. The hypothetical ideal parallel estimate is not
used. This small-case result favors centralized solving; it is not a conclusion
about T=24, larger systems, near-optimal stopping or matched SPAR-DDP performance.

The separate adaptive candidate job 31643493 completed MA57 but failed MA97
during warm-up with ALMOST_LOCALLY_SOLVED on a local OPF. The production driver
rejected it. It is preserved as calibration evidence, not a validated result,
and has not replaced the frozen benchmark protocol.

The archive includes completed verification and benchmark traces, independent
OpenDSS logs, captured sources, resource samples, orchestration and all three
saved-subproblem diagnostic outputs. Julia serialized instances, dispatches and
checkpoints are omitted from this backup and remain on Kamiak. No in-progress
T=24 benchmark output is included.

Archive SHA256:
`e7fb268a452316f3889932fb721817fdb9788949aa56cf62b1418be9762be8b3`
