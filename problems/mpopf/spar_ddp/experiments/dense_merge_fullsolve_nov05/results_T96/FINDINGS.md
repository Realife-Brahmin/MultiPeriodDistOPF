# T96 complete-solve result, November 5, 2026

The grouped right-sided kernel with31 workers and2 BLAS threads wins the matched full-solve comparison. Median optimizer time is41.3375s versus48.2925s for the original balanced/current63x1 baseline, a14.40% reduction. Both pairs improve, and the timing ranges are disjoint from the reference and the63x1 right-kernel candidate. No third repeat is required before the authorized T288 extension.

| Arm | Configuration | Full, s | A, s | B, s | C, s | Full range, s |
|---|---|---:|---:|---:|---:|---:|
| c1 | Balanced / current / 63x1 | 48.2925 | 4.3797 | 35.1426 | 5.7205 | 48.236-48.349 |
| c2 | Grouped b=2 / current / 63x1 | 45.9345 | 4.3994 | 32.8122 | 5.6776 | 45.867-46.002 |
| c3 | Grouped b=2 / right / 63x1 | 43.7375 | 4.3795 | 30.6127 | 5.6765 | 43.682-43.793 |
| c4 | Grouped b=2 / right / 31x2 | 41.3375 | 7.2066 | 23.1692 | 6.7826 | 41.303-41.372 |

Phase B saves11.9734s (34.07%). A/C together cost an extra3.8890s and the remaining optimizer work costs an extra1.1294s, leaving6.9550s of total saving. The narrower worker team therefore slows other phases, but the measured B gain still outweighs them.

All eight main runs have16 accepted steps/solver iterations,25 backward passes,8 barrier-parameter reductions and0 backtracks. All independent physical objectives equal2704807.066934211645 at printed precision; the independent original-model maximum residual is1.136868377216e-12. Numerical residual/dual/complementarity/bound checks and OpenDSS checks pass for each run. The four preceding T3 smoke arms also passed all validators.

All six candidate-versus-same-repeat-reference comparisons pass the unchanged same-iteration/state+dispatch1e-10 scaled-infinity gate. The reference repeats are also explicitly compared. Objective/feasibility trajectories and discrepancies at common iteration/barrier keys are retained in trajectory_differences.csv. No convergence or iteration differences were hidden. Physical x0/ubar raw-array hashes agree across all eight arms.

The after-setup and after-solve audits confirm63x1 for c1/c2/c3 and31x2 for c4, with a one-thread coordinator. Worker checks at the A/C boundaries and before every dense merge remain active. The original temporary one-thread tree-factor regions are retained;31x2 denotes the worker BLAS target and its large-block setting rather than two active threads in every small operation.

One otherwise unoccupied Gold6438Y+ node cn154,64 allocated CPUs,300GiB, Julia1.12.4, full warm-up, pass-level GC and the original input/settings are used throughout. Run order is c1/c2/c3/c4 then c4/c3/c2/c1. Job31574994 COMPLETED0:0 in1h37m36s. Source/harness hashes were checked before and after the batch. Output directories and solver sources were not overwritten. Peak RSS and physical-core/thread masks are retained per run; RSS includes shared mappings and is not simultaneous private-memory consumption.

The standalone post-processing helper initially failed Julia soft-scope handling. Its failure log is preserved, and an isolated posthoc_v2 helper wraps the unchanged checks in a function. This repair did not modify any queued/running numerical source or relax a gate.

## Decision

Advance c4 against c1 at original T288, two repeats in forward/reversed order. Keep grouped b=2, right-kernel algebra, worker/thread settings and scientific controls fixed. Report T96 before submitting the extension. No further tuning is part of the extension; add a third paired repeat only if its timing ranges overlap substantially.

These are DDP implementation ablations with full-convergence timing. They do not establish new centralized speedup ratios or shared-reference near-optimality times. Production defaults and main-paper benchmark rows are untouched.

## Evidence

Raw text/source/manifests/audits/physical and OpenDSS logs, posthoc failure and corrected comparisons are in dense_fullsolve_T96_evidence.tgz. Raw binary solutions remain on Kamiak and solution.sha256 records their hashes. Per-run and median CSVs, solution-comparison CSVs, trajectories and strict validation.json accompany this report. The original immutable snapshot is /home/aryan.jha/kamiak_jobs/dense_merge_fullsolve_nov05.

To reproduce collection, extract the archive into a temporary directory and run the committed collect.py EXTRACTED/dense_merge_fullsolve_nov05 OUTPUT. The archive retains the exact original comparison helper and the separate repaired posthoc_v2 version, with their source hashes.
