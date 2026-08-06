# Repo-specific context for Claude

This file is committed so any machine's Claude Code session starts from the
same ground truth. Session-local memory (`~/.claude/.../memory/`) does not
travel between machines — this file does. Keep it updated when a session
establishes something a future session, on any machine, would need.

## "DDP" means two different things in this repo — do not conflate them

- `ddp/` (repo root) is the **FilterDDP** evaluation: a reproducibility study
  of the authors' *Differential* Dynamic Programming solver
  (github.com/mingu6/FilterDDP.jl, arXiv 2504.08278 / 2606.01487), done
  2026-08-03 to present. Full status, findings, and the "poor fit for
  BFM-scale MPOPF" conclusion are in
  [ddp/README_FILTERDDP_EXPERIMENT.md](ddp/README_FILTERDDP_EXPERIMENT.md).
  The paper is at `ddp/resources/Xu_2026_FilterDDP.pdf` (and a duplicate at
  `envs/ddp/resources/Xu_2026_FilterDDP.pdf`, see manifests below).

- `envs/ddp/` is the user's own, independently-derived decomposition scheme
  for copper-plate MPOPF, titled **"Distributed Dynamic Programming"**
  (`envs/ddp/tex/ddp_copperplate_formulation.tex`, committed 2025-11-12,
  predates any Claude involvement in this repo — verified by git blame, not
  AI-authored). "Distributed" here is the user's own word choice from that
  session, not a Claude artifact.

**These are not the same algorithm** despite sharing the "DDP" initialism:

| | `envs/ddp/` (user's own) | `ddp/` (FilterDDP, the paper) |
|---|---|---|
| Backward information | Passes `μ[t]`, the dynamics-constraint dual, backward one stage per **outer forward sweep** (`envs/ddp/root_level/ddp_copperplate.jl:539-550`, `mu_prev`/`mu_coupling`) | Backward pass sweeps `t=N→1` **within one iteration**, building both `V_x` and `V_xx` (Riccati recursion, `backward_pass.jl` in the FilterDDP clone) |
| Order of approximation | First-order only: the coupling term `μ[t+1]·(B[t+1]−B[t]+Δt·P_B[t+1])` is linear in `B[t]` — no curvature crosses a stage boundary | Second-order: full local quadratic model of the cost-to-go |
| Per-stage solve | Calls an external solver (Gurobi/Ipopt) per stage per sweep | Dense in-house `nu×nu` KKT solve inside the backward pass, with filter line-search globalization |
| Propagation speed | Information from stage `T` reaches stage `1` after roughly `T` outer iterations (one stage per sweep) | Full-horizon propagation in one backward sweep per Newton-type iteration |

Important nuance established 2026-08-05: `μ[t]` in the user's method **is**
legitimately the same object as `V_x[t]` (the costate / value-gradient) —
this is not a naive approach, it's a real first-order relative of DP. The gap
is the missing curvature (`V_xx`) and the cross-iteration staleness, not the
absence of backward-looking information altogether.

## Pending task (do not start until asked)

Write a side-by-side workflow comparison — the user's exact DDP algorithm
vs. FilterDDP's algorithm — **both grounded in the exact problem statement of
the "dummy paper"**, `ddp/paper/sections/copper_plate_model.tex` (user has
approved this problem statement as the shared reference point). Requirements:

- Use **the user's own notation throughout**: `P_Subs^t`, `P_B^t`, `B^t`,
  `c^t`, `C_B`, `p_L^t`, `w`, `B_0`, `P_B^{min/max}`, `B^{min/max}` (from
  `copper_plate_model.tex`), plus `μ[t]`, `λ_Bmin[t]`, `λ_Bmax[t]` (from
  `envs/ddp/tex/ddp_copperplate_formulation.tex`). **No unqualified generic
  control-theory notation** (`L`, `l`, `Q`, `V_x`, `V_xx` unless explicitly
  mapped to one of the user's own symbols first).
- Every dual/auxiliary variable that plays the same role in both workflows
  should be given the **same symbol** in both descriptions, with the
  correspondence stated explicitly (e.g. `μ[t] ↔ V_x[t]`, the box-constraint
  duals `λ_Bmin[t]`/`λ_Bmax[t]` vs. FilterDDP's interior-point bound
  multipliers).
- Where a concept exists in FilterDDP but has no analog in the user's method
  (e.g. `V_xx`/curvature), say so explicitly rather than inventing a
  correspondence.

## Reference-paper bookkeeping

Every `resources/` folder (`ddp/resources/`, `envs/ddp/resources/`,
`envs/multi_poi/resources/`) has a `MANIFEST.txt` (`filename | url |
description`). Run `bash scripts/fetch_resources.sh` from repo root to
fetch anything missing — safe to re-run, present files are left alone. Blank
`url` fields mean the source isn't tracked down yet; fill in when found.
