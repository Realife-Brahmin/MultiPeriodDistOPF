# T96: exact barrier reuse benchmark

Full warmed optimizer median falls **27.008500 -> 24.358500s**,
**9.81% lower (including population and retirement/GC)**. Phase B is
25.11% lower. Two alternating matched pairs, original1020-battery large10k,
grouped b2/right, 48workers x2BLAS, otherwise-unoccupied cn179
Xeon6760P;128physicalCPU600GiB allocation. Only the reuse toggle differs.

| Metric | Fresh baseline | Cache |
|---|---:|---:|
| Full optimizer, s | 27.008500 | 24.358500
| Phase A, s | 4.155452 | 3.970219 |
| Phase B, s | 16.519537 | 12.371232 |
| Phase C, s | 3.363053 | 3.431842 |
| Legacy front-loaded target estimate (not direct elapsed), s | 14.564940 | 14.105083 |


**Quality-clock correction:** the saved legacy extractor computes target time
as logged iteration prefix PLUS the ENTIRE full-solve unallocated remainder.
That remainder includes post-loop cache retirement and cleanup, which happens
AFTER the target. Its old target-time comparison is therefore not a direct
elapsed-time measurement and must not be interpreted as a measured gain or
regression. Raw logs/CSVs and original estimates remain preserved. See
analysis/quality_clock_audit.csv/json for logged-prefix values and conservative
wall-clock bounds; exact allocation of the unlogged remainder is unavailable.
The full-convergence optimizer times and all scientific checks above are unchanged.

Full ranges: baseline26.936-27.081s,
cache24.335-24.382s. Both paired comparisons
improve with separated ranges; no third pair is needed.

All four runs pass original numerical primal/dual/complementarity and original
bounds1e-5, independent original physical objective/primal and OpenDSS with zero
voltage-limit violations. All have16iterations,
25passes,8barrier reductions and0backtracks.
Both final state/control comparisons and baseline-repeat comparison are bit
exact (unchanged1e-10 gate). All recorded objective/original-primal/dual/
complementarity trajectories and accepted steps/mu/backtracking match exactly.
Original physical starting-array SHA matches across all arms/repeats.

Cache:8eligible requests/8hits/0eligible misses/0factorization fallbacks per
solve;17fresh populations. Live numeric payload
6,392,217,600bytes (5.953GiB), never gradients. Warmup validity and
payloads are cleared; first measured pass populates. Both warm and measured
exits exercise the same retirement method. All current vectors/network records,
policies and full forward work remain charged. Posthoc matched_pass_costs.csv
separates every actual hit/population and all phase/GC fields.

Worker pass-GC summed times (not additive wall time) are
baseline14.542-14.574s, cache28.177-29.998s.
These pass counters exclude separately uninstrumented final worker retirement
GC; its wall cost IS inside full_optimizer_s. No allocation-free/detection-free
claim is made. Worker RSS sums are historical per-process high-water marks
including shared mappings/JIT/warmup, not simultaneous private RAM or a measured
private-cache increment. Whole allocation batch MaxRSS116745324KiB includes the
separate expensive diagnostic plus BOTH horizons; do not assign it to one arm.

Timing scope: coordinator model build baseline210.160-210.268s,
cache209.631-210.274s; full warmup baseline125.342-126.191s,
cache121.381-122.567s; process wall including setup/warmup/independent
validation baseline524.000-524.000s,
cache516.000-518.000s. These are
outside the optimizer clock, not clean preprocessing benchmarks or whole-process
speedup claims. The339.584s T24 shadow-validation solve is excluded entirely.

Actual backend/thread/default/OS-task masks pass before/after each run: disjoint
same-socket physical2-core teams plus coordinator1core/BLAS1. Source/input/
harness/runtime binary+libjulia+OpenBLAS checks before/after pass. Job31612890
COMPLETED0:0 in1h08m07. Full source, raw CSV/logs/audits/masks/accounting and
solution hashes are archived; binary solutions stay at the immutable cluster
snapshot. Strict reproduction: collect.py raw analysis --horizon 96; no partial
mode was used. Supplemental full trajectory comparison also passes unchanged.

Archive SHA256: common.tgz 8adb41223777a6fbf04d78cf0995a7280284ce660a3fc45809323682526dc20e; results.tgz f3eda26358f3339836f5ed793b0e91667d4ce9b1c3faeaa64e2e6f566b5ddaa9.
Production defaults, paper and main benchmark rows are unchanged. This case is
committed separately as requested; no broader team or parameter tuning follows.
