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

## Julia environment: use `envs/ddp2026` for everything DDP

Since 2026-08-07 there is **one** environment for all DDP work — FilterDDP, the
centralized JuMP/Ipopt reference, and any further formulation:
`envs/ddp2026` ([README](envs/ddp2026/README.md)). FilterDDP and JuMP coexist
without conflict. Stages 5, 6 and 8 all reproduce under it.

- `envs/ddp/Project.toml` (2025-11) is **superseded**. Its `Plots`/`Crayons`/
  `Revise` deps existed because verification then meant a human reading formatted
  output. Treat `envs/ddp/` as a **read-only reference** for what the user's own
  algorithm did — not expected to be run again.
- `ddp/env` (the FilterDDP-only env) was **removed** 2026-08-07. Its `Manifest`
  was never tracked, so it preserved nothing `envs/ddp2026/Project.toml` doesn't,
  and Stages 5/6 reproduce identically under the new env. It's in git history if
  ever needed.
- Pkg **strips all comments** from `Project.toml` on every `add`/`resolve`, so put
  rationale in a README beside it, never in the file.
- Julia via the **Bash** tool, not PowerShell: PowerShell mangles quotes in
  `julia -e '...'` and has silently corrupted `Project.toml` this way.

## Current copper-plate formulation

`P_Subs^t` has **no upper bound** — only `P_Subs^t ≥ 0` (no export upstream).
This is the latest formulation, confirmed by the user 2026-08-07. Older notes
referring to an active `Psub[2] ≤ 1.35` or `≤ 1.45` are stale; the scripts and
logs were already correct and the experiment README has been fixed to match.

The paper section `ddp/paper/sections/copper_plate_model.tex` was cut down to
equations only on 2026-08-07 at the user's request — the modeling rationale it
used to carry (why no `η`, why the terminal target is a penalty, why `C_B = 0.5`
here vs. `≈10⁻⁶·min c^t` in the tADMM paper) now survives only in
`ddp/README_FILTERDDP_EXPERIMENT.md`.

**Open question the user raised, not yet decided:** demand and price in
`tab:cp_data` are nearly collinear (`r = 0.9968`). See the "Known weakness"
note in the experiment README. Changing the price series invalidates every
recorded number, so it is the user's call.

## What is and isn't verified (as of 2026-08-07)

FilterDDP is cross-checked against a **centralized JuMP/Ipopt** solve of
`eq:cp_all` on all 8 copper-plate cases (Stage 8 of the experiment README):
`2e-16` agreement on the base instances, worst gap `1.1e-10` with bounds, and
Ipopt independently certifies the 6g bound set infeasible. On the base instances
there are three mutually independent references (closed form, active-set QP,
Ipopt). **Do not redo this.**

Not verified, so don't overclaim: anything with a network (no LinDistFlow, no
BFM); horizons beyond `T = 6`; and **the user's own DDP algorithm, which has not
been compared at all** — that is the pending task below.

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
