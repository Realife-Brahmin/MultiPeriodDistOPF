---
name: Active Investigation Mar 31
description: k=1 warm-start experiment in progress, loose-tol SOCP pre-solve approach, convergence CSV auto-export added
type: project
---

## Completed Today (2026-03-31)

### large10kC_1ph T=6 run — DONE
- BF: $2,919,076 in 629s
- tADMM: $2,919,076 in 585s (25 ADMM iters), ρ 1500→46.9
- Speedup: **1.07x** — improved from T=4's 0.77x
- Committed as `2776bf2b`

### Convergence CSV auto-export — DONE
- Added automatic `convergence_data.csv` export to `tadmm_socp.jl` (objfun, r, s, rho per iteration)
- No longer need to run `plot_convergence.jl` separately

### k=1 cold-start analysis
- k=1: 145.6s (359 Ipopt iters on worst subproblem) — 25% of total tADMM time
- k=2: 11s (7 Ipopt iters) — warm-start is 13x faster
- If k=1 could be as fast as k=2, total tADMM would be ~450s → 1.40x speedup

## In Progress: k=1 Warm-Start via Pre-Solve

### Approach 1: LinDistFlow pre-solve — FAILED
- Drops ℓ variable and SOC constraint, simplifies voltage drop to v_j = v_i - 2(rP+xQ)
- **LOCALLY_INFEASIBLE** even with 2% voltage bound relaxation
- LinDistFlow approximation too inaccurate for this heavily-loaded 10k-bus network with 1020 large batteries

### Approach 2: Loose-tol SOCP pre-solve — RUNNING
- Same SOCP formulation, but Ipopt limited to 30 iterations with tol=1e-2
- Accepts ITERATION_LIMIT and ALMOST_LOCALLY_SOLVED as valid statuses
- Should get a rough-but-feasible solution quickly
- Run `btdkgvent` in progress

## Current Config State (tadmm_socp.jl)
- `systemName = "large10kC_1ph"` 
- `T = 6`
- `lindistflow_warmstart = true` (currently using loose-tol SOCP approach, not actual LinDistFlow)
- `rho_base = 3000.0`, adaptive, sqrt scaling
- `enable_warm_start = true`
- LinDistFlow code path still exists in `primal_update_tadmm_socp!` (lindistflow param) but not used

## Future Tasks (user mentioned, not yet started)
- Skip BF when results already exist (reuse .jls)
- Store results as CSVs instead of JLS for offline access
- User's PhD focus: solving nonlinear problems without approximation — LinDistFlow only useful as warm-start tool
