# T144 exact barrier reuse benchmark

Full warmed optimizer median **46.465 -> 41.363s**, **10.98% lower**, including cache population and retirement/GC. Phase B **29.950566 -> 22.912266s**, **23.50% lower**. Phase A6.276030->6.452266s; C5.431060->5.364063s. Baseline range46.249-46.681s and cache41.334-41.392s are separated; no third pair needed.

Two full-warmup matched pairs in ABBA order, original1020-battery large10k, grouped b2/right,48workers x2BLAS plus one-core coordinator on otherwise-unoccupied cn179 Xeon6760P,128physicalCPU600GiB. Only reuse toggle differs. This carries the established T96 team policy; no optimality claim for longer horizons.

Near-optimal target uses the HISTORICAL October16 original1020 central OBJECTIVE2707471.919112729, within0.1% and original primal1e-5, not a fresh same-node centralized timing comparison. All arms first meet target at iteration9. Direct elapsed-target timestamps were not recorded. Logged iteration prefix20.545541->18.643545s (9.26% less logged work); conservative median elapsed bounds baseline[20.545541,21.011609] and cache[18.643545,20.015601] establish earlier attainment, with reduction bounded2.58%-11.27%. Do not interpret legacy front-loaded estimates21.010609->20.014601s as actual target timestamps: they allocate post-loop cleanup before every quality point. Raw estimates retained; see quality_clock_audit.csv/json.

All four runs pass original numerical primal/dual/complementarity/bounds1e-5, independent original physical residual1.0232e-12 and OpenDSS largest difference4.87e-8pu with zero voltage violations. Both same-repeat and baseline-repeat final state/control differences0 under unchanged1e-10 gate. Full objective/original/dual/complementarity trajectories, accepted step sizes/mu/backtracks are exact. Identical all starting x0/ubar hashes. All18accepted steps,27passes,8barrier reductions,0backtracks.

Cache8eligible/8hits/0miss/0factorization fallback,19fresh populations per solve;7980populated and3360reused merges. Persistent numeric payload10,487,232,000bytes (9.767GiB). Warm payloads cleared; measured population, current vector/policy/forward work, synchronization and retirement charged. Worker pass-GC sums baseline29.435-29.642s vs cache51.630-53.085s overlap across workers and cannot be added to wall. Pass counters exclude separately uninstrumented retirement GC, whose wall cost is included. Worker RSS sums are historical high-water marks including shared mappings/JIT/warmup, not simultaneous private memory/cache increment. Batch MaxRSS154790684KiB spans all four processes.

Separate timing scopes: model build311.791-313.112s; full warmup baseline145.856-149.672s, cache153.996-155.998s; total process wall baseline710-745s, cache711-713s including setup/warmup/validation. No whole-process speedup claim.

Job31613805 COMPLETED0:0 in48m08. Actual worker count/backend/default BLAS2/OS masks and disjoint same-socket physical core teams verified before/after. Source/input/harness/runtime binary/libjulia/OpenBLAS hashes pass before/after. Empty job stderr. All raw logs/CSV/masks/hardware/code/accounting/hashes archived, actual source in runtime-source.tgz, binary solutions retained immutable on cluster. Strict collect.py --horizon144 and supplementary exact trajectory comparisons pass; gates unchanged.

Archive SHA256 common.tgz eb7aca28f490e6d0cee5bbda311da1db643770d44b896c2337d9f2a198a795da; results.tgz 7930b556f17c56ea422ef6c1e46c501f5f5a7023b8b728736e04e963bff836db.

No defaults, paper or main result rows changed. Next authorized experiment is the prepared T288 matched benchmark with passive direct elapsed markers, without repeated cache diagnostics or tuning.
