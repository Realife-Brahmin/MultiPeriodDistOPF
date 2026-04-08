---
name: large10kC System
description: 10321-bus system with 5x-scaled PVs (350kW) and batteries (400-600 kVA), results at T=4 and T=6
type: project
---

Created 2026-03-31.

- **Buses**: 10,321 (same topology as large10k_1ph)
- **Batteries**: 1,020 (same bus locations as large10kB, seed=42)
  - kVA: 400, 450, 500, 550, 600 (cycling) — 5x of large10kB
  - Avg: 417 kW rated, total ~425 MW (41% of load)
- **PVs**: 1,022 (same bus locations as large10k_1ph)
  - 350 kW / 420 kVA each — 5x of original 70 kW
  - Total: 358 MW (35% of load)
- **Loads**: 10,320 × 100 kW = 1,032 MW (unchanged)
- **Files**: `rawData/large10kC_1ph/`

## T=4 Results
- BF: $2,905,284 in 230s (311 Ipopt iters)
- tADMM: $2,905,284 in 297s (18 ADMM iters), adaptive ρ 1224.7→153.1
- Speedup: **0.77x** (tADMM slower)

## T=6 Results
- BF: $2,919,076 in 629s
- tADMM: $2,919,076 in 585s (25 ADMM iters), adaptive ρ 1500→46.9
- Speedup: **1.07x** (marginal tADMM win)
- Objective match: $0.006 difference (~0.0000%)
- k=1 cold start: 145.6s (max subproblem, 359 Ipopt iters) — 25% of total
- k=2 warm-started: ~11s per subproblem — warm-start gives 13x speedup per iteration
- Convergence data CSV: `envs/tadmm/processedData/large10kC_1ph_T6/convergence_data.csv`
- At 0.5% gap (k=17): cumulative 456s → would give 1.38x speedup

## Key Observations
- T=4→T=6 improved speedup from 0.77x to 1.07x (BF scales worse with T)
- LinDistFlow pre-solve attempted but infeasible on this network (too heavily loaded)
- Loose-tolerance SOCP pre-solve being tested as alternative k=1 warm-start strategy
