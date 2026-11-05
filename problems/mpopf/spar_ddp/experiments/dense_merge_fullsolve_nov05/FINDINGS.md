# Complete dense-merge ablation result, November 5, 2026

The exact dense-merge improvements survive the tested complete original large10k solves at both horizons. The selected grouped b=2/right-sided kernel/31-workers x2-BLAS-thread configuration lowers full optimizer time by14.40% atT96 and24.03% atT288 against the original balanced/current-kernel63x1 configuration. All physical batteries and model constraints are retained.

## T288 extension

| Configuration | Full, s | A, s | B, s | C, s | Full range, s |
|---|---:|---:|---:|---:|---:|
| Original balanced / current / 63x1 | 187.5305 | 12.2257 | 150.0921 | 14.8131 | 187.217-187.844 |
| Grouped b=2 / right / 31x2 | 142.4715 | 19.9423 | 91.1305 | 18.1959 | 142.338-142.605 |

Two-repeat medians use forward/reversed order on the same otherwise unoccupied node. Phase B saves58.9616s (39.28%). A/C cost an extra11.0994s, and other optimizer work an extra2.8032s, yielding45.0590s net saving (24.03%, 1.316x). The candidate wins both repeat pairs; its slowest run is below the reference's fastest run. No third pair is needed.

All four T288 runs have18 accepted steps/solver iterations,26 backward passes,7 reductions of the barrier parameter and1 backtrack. Each passes the original numerical residual/dual/complementarity/bound1e-5 gates, independent physical residual validation and OpenDSS checks over2,972,160 bus-periods with zero voltage violations under the existing1e-4 tolerance. Independent original primal maxima are1.136868377216e-12. Independent objectives agree within about4.66e-10 absolute, or1.72e-16 relative.

Both candidate solutions pass the unchanged same-iteration and1e-10 scaled state/control comparison against the original reference from the same repeat. Maximum scaled state difference is7.47e-15 and control difference1.33e-15. The reference repeats also match. Objective/original-feasibility trajectories are retained and compared at common iteration/barrier keys; all differences remain visible in the CSV. The physical initial arrays are byte-hash identical within this horizon.

Actual after-setup and after-solve audits confirm63x1 for the reference and31x2 for the candidate; the coordinator uses one thread/core. All OS-thread masks match disjoint physical-core groups. Assertions at phase boundaries and every merge verify the chosen settings. Existing temporary _blas1 regions for small tree blocks are retained. The grouped candidate uses912 merges/10 rounds per pass versus1392/9 for the baseline.

## T96 batch

| Arm | Configuration | Full optimizer median, s |
|---|---|---:|
| c1 | balanced/current63x1 |48.2925|
| c2 | grouped/current63x1 |45.9345|
| c3 | grouped/right63x1 |43.7375|
| c4 | grouped/right31x2 |41.3375|

All eight main T96 solves and four preliminaryT3 smoke solves pass independent validation. The T96 arms have identical16 iterations/25 passes/8 barrier reductions/zero backtracks and strict final-solution agreement. Its selected candidate reduces Phase B35.142625->23.169177s, with the A/C and other-cost tradeoff detailed in results_T96/FINDINGS.md. The T96 result was reported and pushed before the T288 extension.

## Timing scope and overhead

The reported full time is the complete optimizer loop after a full warm-up, including A/B/C, directions, line search, GC and ordinary coordinator work. Model construction, warm-up and independent physical/OpenDSS validation are outside that clock. The separately retained process-wall, coordinator model-build and warm-up columns show those substantial costs; no whole-process wall-time speedup is claimed from these optimizer medians. Build time is not a clean compilation-free preprocessing benchmark. Existing benchmark rows and centralized speedup claims are not changed.

Worker and coordinator historical peak RSS, backend, affinity and all thread settings are recorded per run. Slurm batch maxRSS is254,019,824KiB over this T288 job; it is an overall high-water mark and cannot be assigned to one configuration. Worker sums include shared mappings and are not simultaneous private RAM.

## Decision

Advance the grouped/right31x2 configuration as an experimental option for the original1020-battery large10k workloads tested here. The Phase B gains outweigh the measured A/C penalties at bothT96 andT288. The grouped/right63x1 arm remains a useful control for isolating the kernel change, especially where fewer workers could slow other phases. Broader network/battery/horizon generalization is not established by these two horizons. Production-default integration requires a separate change; this campaign did not modify those defaults.

## Preservation and reproduction

T96 job31574994 completed0:0 in1h37m36s; T288 job31575079 completed0:0 in1h51m49s. Both ran on otherwise unoccupied Gold6438Y+cn154 with64 allocated CPUs/300GiB, Julia1.12.4, same initialization/model/profile/screening/Cb2.4e-5/gamma420/fullwarm-up/pass-GC/stopping policy. Input, harness and executed source manifests passed before and after each batch. Raw binary solutions remain on Kamiak; hashes accompany archived logs.

The only failed attempt was the standalone post-processing helper's Julia soft-scope error. Its original source/log are preserved. An isolated posthoc_v2 helper fixed function scope with unchanged checks, leaving active numerical files untouched. No numerical/convergence failure, gate relaxation or solver-source repair was needed.

results_T288/dense_fullsolve_T288_evidence.tgz contains raw logs, audits, masks/topology, scientific checks, source/input manifests, solution hashes and Slurm accounting. Extract it and run collect.py EXTRACTED/dense_merge_fullsolve_nov05 OUTPUT --T288 to reproduce strict validation and tables. The separate T96 evidence archive and its published comparison are retained intact.

Final assessment is complete. No extra tuning, centralized solves, expanded campaigns or paper edits were performed.
