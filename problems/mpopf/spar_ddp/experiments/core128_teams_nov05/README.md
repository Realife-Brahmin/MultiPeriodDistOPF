# T96 worker/BLAS teams:48x2 versus24x4

Authorized November5 after the completed128-core T96 study. Original large10k,
1020 individual batteries,T96. Compare groupedb2/right-sided dense kernel at
48workers x2BLAS threads versus24workers x4threads, both using97physical cores
including the coordinator, inside one128-core Xeon6760P allocation.

Two full-warmup complete solves each, order48x2/24x4 then24x4/48x2.
Same physical start, frozen input, Cb2.4e-5/gamma420, screening, passGC,
tolerances, solver source and hardware for every arm. No fresh centralized
runs. This tests complete optimizer performance, not only Phase B.

Reuse exact immutable variants/source from dense_merge_core128_nov05.
Only the launcher chooses the different worker/thread settings. The same
opt-in runtime assertions, right-kernel algebra, upper-triangle/fallback
semantics, within-socket disjoint physical-core pinning and before/after
OS-thread audits remain active. Extend offline affinity regression to24x4.
Julia/OMP1, worker BLAS target2or4 verified AFTER driver setup and everymerge;
temporary _blas1 regions retained and must restore target.

Require original numerical/physical/bound and independent OpenDSS checks,
finite objective/feasibility traces, same physical starting array hashes,
source/input/harness hashes, Slurm0 and actual thread/mask/backend audits.
Use repaired compare_solutions.jl for24x4 versus48x2 SAMErepeat and reference
repeat agreement. Same-iteration and scaled state/control<=1e-10 gates are
unchanged; preserve any failure or trajectory/convergence difference.
Independent numerical/physical1e-5 and finalobjective0.1% gates also apply.

Report two-repeat complete warmed optimizer medians and ranges, A/B/C,
iterations/passes/muupdates/backtracks, trajectories, GC, worker/coordinator
peak RSS and setup/fullwarmup/processwall separately. No claim that summed
worker RSS is simultaneous private memory. Existing historical48x2 median
26.686s is context; judge against the FRESH paired baseline.
Add one third pair only if timings substantially overlap before naming a
winner, using a NEW immutable launcher and r3 paths; never rerun job.sh over
existing r1/r2. No broader tuning, central sweep or horizon expansion.

No expensive overlap, no queued/running source edits, safe repairs only new
immutable attempts with failure logs preserved. Production defaults, paper,
published benchmark rows, MSOPF and every physical battery/constraint unchanged.
Collect raw logs/hashes/audits/hardware/accounting and strict comparisons;
binary solutions stay cluster. Commit/push evidence and concise findings,
verifyremote and delete completion heartbeat only after final assessment.

