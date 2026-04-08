---
name: Results Table
description: All BF vs tADMM timing comparisons across systems, battery counts, and T values
type: project
---

## Complete Results (as of 2026-03-31)

| System | Buses | Batteries | PVs | T | BF time | tADMM time | Iters | Speedup |
|--------|-------|-----------|-----|---|---------|------------|-------|---------|
| ieee123A_1ph | 123 | 26 | ? | 24 | ? | ? | ? | **2.15x** |
| ieee2552_1ph | 2,552 | 1-10 | ? | 24 | ? | ? | ? | 0.89x |
| large10k_1ph | 10,321 | 102 | 1,022 | 4 | 289s | 144s | 9 | **2.0x** |
| large10k_1ph | 10,321 | 102 | 1,022 | 6 | 366s | 399s | 15 | 0.92x |
| large10kB_1ph | 10,321 | 1,020 | 1,022 | 4 | 417s | 189s | 9 | **2.21x** |
| large10kB_1ph | 10,321 | 1,020 | 1,022 | 6 | 868s | 497s | ? | **1.75x** |
| large10kB_1ph | 10,321 | 1,020 | 1,022 | 12 | 1,071s | 1,330s | 37 | 0.81x |
| large10kC_1ph | 10,321 | 1,020 (5x) | 1,022 (5x) | 4 | 230s | 297s | 18 | 0.77x |
| large10kC_1ph | 10,321 | 1,020 (5x) | 1,022 (5x) | 6 | 629s | 585s | 25 | **1.07x** |

## Key Findings
1. **Battery count matters more than system size** for tADMM benefit
2. **Small T (4-6) favors tADMM**, large T (12+) erases the advantage
3. **More batteries shift crossover point** — 1020 batteries keeps tADMM beneficial up to T=6 (was T=4 with 102)
4. **Iteration count is the bottleneck at high T** — T=12 needs 37 iters vs 9 at T=4
5. **Subproblem times are consistent after warmup** (~10-14 Ipopt iters per ADMM iter), so reducing ADMM iterations is a clean win
6. **k=1 cold start is the biggest bottleneck** — on large10kC T=6, k=1 took 145.6s (25% of total tADMM time), k=2 only 11s
7. **Early stopping at 0.5% gap** would improve large10kC T=6 from 1.07x to 1.38x (456s vs 585s at k=17)
