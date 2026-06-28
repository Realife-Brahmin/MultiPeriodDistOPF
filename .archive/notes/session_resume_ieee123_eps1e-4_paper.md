# Session resume: ieee123 multi-T sweep at eps_pri=1e-4 for paper

Last touched: 2026-05-22. This file exists so the work can be resumed from
any machine; the live plan file at `~/.claude/plans/glimmering-snacking-sprout.md`
is local-only.

## Why we're doing this

Paper `IAS-Trans-2025-Scaling-MPOPF-Computation-via-Temporal-Decomposition`
currently reports ieee123 as a single solution-quality point (T=96, rho=4500).
We're upgrading ieee123 to a multi-T scaling reference (T = 6, 12, 24, 48, 96,
144) using `eps_pri=1e-4` uniformly. The 1e-3 tolerance was discarded after
T=144 produced false convergence (rho=6000 ended at $1325 vs BF $2811, -53%).

The honest framing: ieee123 is small enough that tADMM never beats BF on
wall-clock — speedup < 1× at every T. We report it as a solution-quality
reference, not a speedup case, plus a methodology-notes subsection covering
the eps_pri tightening and the "sufficiently high rho_0" heuristic.

## What's done (committed on branch `theory-simulations-results-may18`)

| T   | rhos run                  | winner (NO eff time, gap)         | commits |
|-----|---------------------------|-----------------------------------|---------|
| 6   | 250, 1000, 4000           | rho=4000 (2.4s, 0.48%)            | 3       |
| 12  | 500, 2000, 8000           | rho=8000 (5.6s, 0.46%)            | 3       |
| 24  | 1000, 16000 (+4000 CRASH) | rho=16000 (22.5s, 0.31%)          | 3       |
| 144 | 6000, 24000, 96000        | rho=96000 (gap 0.36%)             | 3       |

T=24 rho=4000 hits Ipopt `NUMERICAL_ERROR` reproducibly with this hyperparam
combo; documented in `envs/tadmm/processedData/ieee123C_1ph_T24/rho_sweep_eps1e-4/rho_4000/CRASH.txt`.
Not worth chasing — rho=16000 won T=24 anyway.

Each rho lives in `envs/tadmm/processedData/ieee123C_1ph_T{T}/rho_sweep_eps1e-4/rho_{R}/`
with `near_optimal_summary.csv`, `convergence_data.csv`, `results_socp_tadmm.txt`,
`subproblem_timing_details.csv`, `tadmm_run.log`, `run_stdout.log`, `params.txt`,
`sol_socp_tadmm.jls`, `convergence/tadmm_convergence_socp.png`.

## What's pending

### Phase 1 (sims) — 2 of 5 sweeps remaining

T=48 rhos = 2000, 8000, 32000 (est ~75 min)
T=96 rhos = 4000, 16000, 64000 (est ~2 hr)

Command template:
```bash
SYSTEM_OVERRIDE=ieee123C_1ph \
SWEEP_T=<T> SWEEP_RHOS="<r1>,<r2>,<r3>" \
SWEEP_DIR_NAME="rho_sweep_eps1e-4" \
ADAPTIVE_RHO_OVERRIDE=true \
EPS_PRI_OVERRIDE=1e-4 EPS_DUAL_OVERRIDE=1e-3 \
TAU_INCR_OVERRIDE=2 TAU_DECR_OVERRIDE=2 \
julia --project=envs/tadmm --threads=16 run_rho_sweep.jl
```

After each T finishes:
1. `julia --project=envs/tadmm plot_convergence_from_csv.jl --sweep-dir envs/tadmm/processedData/ieee123C_1ph_T{T}/rho_sweep_eps1e-4`
2. `cp .../convergence/tadmm_convergence_socp.png pngs_review/ieee123C_1ph/T{TTT}_rho{RRRRRR}_eps1e-4_convergence.png` (one per rho)
3. Commit each rho atomically — message format matches the existing 12 commits

Watchpoints:
- Disable AC sleep before long runs: `powercfg -change -standby-timeout-ac 0`
  (set this session, but re-confirm if on a different machine).
- The sweep-script silently swallows Ipopt `NUMERICAL_ERROR` crashes and
  snapshots stale CSVs from the previous rho. After each sweep, audit
  per-rho `convergence_data.csv`: first iter's `rho` column must equal
  the dir name. Mismatched dirs need a direct `run_tadmm.jl` rerun
  with `RHO_OVERRIDE` or a `CRASH.txt`.

### Phase 2 (paper edits) — not started

Nine specific edits in `c:/repos_addendum/IAS-Trans-2025-Scaling-MPOPF-Computation-via-Temporal-Decomposition/sections/{results,simulation}.tex`,
plus a new figure generator `plot_ieee123_scaling.jl`. Full edit list
with file:line targets is in the local plan file. **Do not commit paper
edits — standing instruction: user reviews and commits paper repo
themselves.**

## Key context the work depends on

- `eps_pri=1e-4` is mandatory at ieee123 scale; 1e-3 produces false
  convergence below BF (memory: `feedback-eps-pri-small-systems`).
- Fastest tADMM convergence happens when `rho_0` is high enough that
  adaptive ρ only DECREASES, never needing Phase 1 ramp-up (memory:
  `feedback-rho-init-heuristic`).
- Standing instructions on the paper repo, PNG-per-run, atomic
  story-based commits, never-kill-inflight-sims, and BF-baseline
  validation all apply.

## Files of interest

- [run_rho_sweep.jl](run_rho_sweep.jl) — sweep orchestrator (supports
  `SWEEP_DIR_NAME`)
- [plot_convergence_from_csv.jl](plot_convergence_from_csv.jl) —
  per-rho convergence PNG generator (`--sweep-dir <path>` mode)
- [envs/tadmm/processedData/ieee123C_1ph_T144/rho_sweep_eps1e-4/](envs/tadmm/processedData/ieee123C_1ph_T144/rho_sweep_eps1e-4/) —
  template directory layout for Phase 1 outputs
