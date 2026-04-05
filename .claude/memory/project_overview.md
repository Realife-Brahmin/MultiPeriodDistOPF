---
name: Project Overview
description: Julia MPOPF project comparing brute-force vs tADMM temporal decomposition
type: project
---

Julia project for Multi-Period Distribution Optimal Power Flow (MPOPF).
Compares **brute-force (BF)** centralized solve vs **tADMM** (temporal ADMM) decomposed solve.

- Repo: `c:\repos_addendum\MultiPeriodDistOPF`
- Solver: **Ipopt only** (no Gurobi license)
- Threads: 16 available, use `JULIA_NUM_THREADS=16`
- Entry point: `tadmm_socp.jl`
- Branch: `continued-testing-mar31` (created 2026-03-31 from master after merging `pushing-bf-mar18` via PR #147)
- GitHub: `Realife-Brahmin/MultiPeriodDistOPF`
- `gh` CLI installed and authenticated as `Realife-Brahmin`

**Why:** Investigating when temporal ADMM decomposition outperforms centralized solve for distribution OPF.
**How to apply:** Always use `JULIA_NUM_THREADS=16`, Ipopt solver, run from repo root.
