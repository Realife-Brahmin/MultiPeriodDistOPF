# T24: exact barrier reuse benchmark

Full warmed optimizer median falls **13.909000 -> 12.414500s**,
**10.74% lower (including population and retirement/GC)**. Phase B is
31.46% lower. Two alternating matched pairs, original1020-battery large10k,
grouped b2/right, 24workers x2BLAS, otherwise-unoccupied cn179
Xeon6760P;128physicalCPU600GiB allocation. Only the reuse toggle differs.

| Metric | Fresh baseline | Cache |
|---|---:|---:|
| Full optimizer, s | 13.909000 | 12.414500
| Phase A, s | 1.689910 | 1.679737 |
| Phase B, s | 9.576299 | 6.563754 |
| Phase C, s | 1.533463 | 1.500521 |
| Legacy front-loaded target estimate (not direct elapsed), s | 7.355229 | 7.805897 |


**Quality-clock correction:** the saved legacy extractor computes target time
as logged iteration prefix PLUS the ENTIRE full-solve unallocated remainder.
That remainder includes post-loop cache retirement and cleanup, which happens
AFTER the target. Its old target-time comparison is therefore not a direct
elapsed-time measurement and must not be interpreted as a measured gain or
regression. Raw logs/CSVs and original estimates remain preserved. See
analysis/quality_clock_audit.csv/json for logged-prefix values and conservative
wall-clock bounds; exact allocation of the unlogged remainder is unavailable.
The full-convergence optimizer times and all scientific checks above are unchanged.

Full ranges: baseline13.893-13.925s,
cache12.328-12.501s. Both paired comparisons
improve with separated ranges; no third pair is needed.

All four runs pass original numerical primal/dual/complementarity and original
bounds1e-5, independent original physical objective/primal and OpenDSS with zero
voltage-limit violations. All have14iterations,
23passes,8barrier reductions and0backtracks.
Both final state/control comparisons and baseline-repeat comparison are bit
exact (unchanged1e-10 gate). All recorded objective/original-primal/dual/
complementarity trajectories and accepted steps/mu/backtracking match exactly.
Original physical starting-array SHA matches across all arms/repeats.

Cache:8eligible requests/8hits/0eligible misses/0factorization fallbacks per
solve;15fresh populations. Live numeric payload
1,298,419,200bytes (1.209GiB), never gradients. Warmup validity and
payloads are cleared; first measured pass populates. Both warm and measured
exits exercise the same retirement method. All current vectors/network records,
policies and full forward work remain charged. Posthoc matched_pass_costs.csv
separates every actual hit/population and all phase/GC fields.

Worker pass-GC summed times (not additive wall time) are
baseline4.153-4.275s, cache8.000-8.968s.
These pass counters exclude separately uninstrumented final worker retirement
GC; its wall cost IS inside full_optimizer_s. No allocation-free/detection-free
claim is made. Worker RSS sums are historical per-process high-water marks
including shared mappings/JIT/warmup, not simultaneous private RAM or a measured
private-cache increment. Whole allocation batch MaxRSS116745324KiB includes the
separate expensive diagnostic plus BOTH horizons; do not assign it to one arm.

Timing scope: coordinator model build baseline55.133-55.302s,
cache54.820-55.353s; full warmup baseline77.670-81.861s,
cache83.708-85.143s; process wall including setup/warmup/independent
validation baseline239.000-242.000s,
cache243.000-244.000s. These are
outside the optimizer clock, not clean preprocessing benchmarks or whole-process
speedup claims. The339.584s T24 shadow-validation solve is excluded entirely.

Actual backend/thread/default/OS-task masks pass before/after each run: disjoint
same-socket physical2-core teams plus coordinator1core/BLAS1. Source/input/
harness/runtime binary+libjulia+OpenBLAS checks before/after pass. Job31612890
COMPLETED0:0 in1h08m07. Full source, raw CSV/logs/audits/masks/accounting and
solution hashes are archived; binary solutions stay at the immutable cluster
snapshot. Strict reproduction: collect.py raw analysis --horizon 24; no partial
mode was used. Supplemental full trajectory comparison also passes unchanged.

Archive SHA256: common.tgz 8adb41223777a6fbf04d78cf0995a7280284ce660a3fc45809323682526dc20e; results.tgz 9607a6f4ec56a7a50fd9d731646e0c782f5ed6190fdab3c803fcbd3017891849.
Production defaults, paper and main benchmark rows are unchanged. This case is
committed separately as requested; no broader team or parameter tuning follows.
