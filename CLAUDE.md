# Repo-specific context for Claude

This file is committed so any machine's Claude Code session starts from the
same ground truth. Session-local memory (`~/.claude/.../memory/`) does not
travel between machines — this file does. Keep it updated when a session
establishes something a future session, on any machine, would need.

## Working branch

Current branch as of 2026-09-15 is `ddp-understanding-sep15` (`sep14` was merged
to `master` as PR #158 and deleted). Create a new dated branch for new work
rather than reusing it.

## The solver is `DDP4OPF`, committed at `ddp/DDP4OPF.jl`

Since 2026-09-15 the DDP solver is **our own MIT fork of FilterDDP.jl**, committed
as source, named `DDP4OPF` with a fresh UUID (`005a218a-...`). `using DDP4OPF`
everywhere; `FilterDDP` is no longer resolvable in `envs/ddp2026`. Mingda Xu's
MIT notice is retained in its `LICENSE`, and `NOTICE.md` records the upstream
repo, the fork point `513a104`, and what changed. Credit the original papers.

Why it was forked rather than kept as a clone plus patches: the working solver
existed only as uncommitted edits in a gitignored clone, and the documented patch
recipe no longer rebuilt it (3 of 9 patches failed on a clean `513a104`; the
result differed in 4 source files). The fork was verified **bit-identical** to
that clone -- 0 of 2526 full-space and 0 of 459 reduced-space solution entries
differed on ieee123 T=3.

- A new name and UUID were **required**, not cosmetic: FilterDDP is registered in
  the Julia General registry, and Julia identifies packages by UUID.
- **`ddp/patches/` is history, not a build recipe.** Do not try to rebuild from it.
- **Environment variables still use the `FILTERDDP_` prefix** on purpose: sweep and
  queue scripts set them. Renaming them is a separate decision.
- `ddp/external/FilterDDP.jl` (gitignored) is now only needed to reproduce the
  authors' own shipped example (Stage 3).
- Historical notes in `ddp/notes/` still cite `ddp/external/FilterDDP.jl` paths and
  line numbers. They record the clone as it was at the time; leave them.

## Centralized IPOPT timing sweep (complete)

That sweep is finished across all three systems and all matrix horizons, with
the TPEC table/caption/PDF updated. The procedure below is kept for reference. Reconstruct the IAS-style
centralized `C (s)` column with fresh JuMP--IPOPT runs by following
`ddp/notes/CENTRALIZED_IPOPT_TIMING_SWEEP.md`. Use
`scripts/run_centralized_ipopt_case.ps1`, one case at a time. A case is not
finished until its validated row and raw log are pushed here and the matching
TPEC table/PDF update is pushed to the TPEC repository. Resume from the first
incomplete row; never replace missing data with an older Gurobi timing.

## Two DDP codebases here — both *Differential* Dynamic Programming

**Naming: the user's method is DIFFERENTIAL Dynamic Programming. It is never
"Distributed."** Corrected throughout on 2026-08-14, including the two source
files that originally carried the wrong expansion
(`envs/ddp/tex/ddp_copperplate_formulation.tex`,
`envs/ddp/root_level/ddp_copperplate.jl`). An earlier version of this file
recorded "Distributed" as a deliberate word choice; that was wrong. Do not
reintroduce it anywhere, and do not describe the two codebases as different
algorithms that merely share an initialism — they are both DDP, differing in
**order**. (Unrelated: "Distributed Energy Resources" in `envs/multi_poi/` and
`envs/tadmm/` is correct and must be left alone; and the repo name's "Dist"
means *distribution* networks.)

- `ddp/` (repo root) is the **FilterDDP** evaluation: a reproducibility study
  of the authors' solver (github.com/mingu6/FilterDDP.jl, arXiv 2504.08278 /
  2606.01487), done 2026-08-03 to present. Full status, findings, and the
  MPOPF-fit conclusion are in
  [ddp/README_FILTERDDP_EXPERIMENT.md](ddp/README_FILTERDDP_EXPERIMENT.md).
  The paper is at `ddp/resources/Xu_2026_FilterDDP.pdf` (and a duplicate at
  `envs/ddp/resources/Xu_2026_FilterDDP.pdf`, see manifests below).

- `envs/ddp/` is the user's own, independently-derived **first-order** DDP
  scheme for copper-plate MPOPF (committed 2025-11-12, predates any Claude
  involvement — verified by git blame, not AI-authored). A clean minimal
  restatement of it lives at
  [ddp/examples/power_system/user_ddp.jl](ddp/examples/power_system/user_ddp.jl).

**The difference is the order of the cost-to-go model, not the family:**

| | `envs/ddp/` (user's own) | `ddp/` (FilterDDP, the paper) |
|---|---|---|
| Backward information | Passes `μ[t]`, the dynamics-constraint dual, backward one stage per **outer forward sweep** (`envs/ddp/root_level/ddp_copperplate.jl:539-550`, `mu_prev`/`mu_coupling`) | Backward pass sweeps `t=N→1` **within one iteration**, building both `V_x` and `V_xx` (Riccati recursion, `backward_pass.jl` in the FilterDDP clone) |
| Order of approximation | First-order only: the coupling term `μ[t+1]·(B[t+1]−B[t]+Δt·P_B[t+1])` is linear in `B[t]` — no curvature crosses a stage boundary | Second-order: full local quadratic model of the cost-to-go |
| Per-stage solve | Calls an external solver (Gurobi/Ipopt) per stage per sweep | In-house **sparse** `nu×nu` KKT solve (UMFPACK) inside the backward pass, with filter line-search globalization. An earlier version of this row said "dense"; measured at large10k a dense factorisation would need ~540 s and 22 GB per stage against 2.25 s actual |
| Propagation speed | Information from stage `T` reaches stage `1` after roughly `T` outer iterations (one stage per sweep) | Full-horizon propagation in one backward sweep per Newton-type iteration |

Important nuance established 2026-08-05: `μ[t]` in the user's method **is**
legitimately the same object as `V_x[t]` (the costate / value-gradient) —
this is not a naive approach, it's a genuine first-order DDP. The gap
is the missing curvature (`V_xx`) and the cross-iteration staleness, not the
absence of backward-looking information altogether.

**Measured 2026-08-14** on the shared `T = 6` instance, via `user_ddp.jl`:
the first-order scheme's *fixed point is correct* — fed the true `μ*` and `B*`,
one sweep reproduces the centralized optimum to `1.7e-08` — but it does not
converge to it from a cold start, settling into a period-two limit cycle
`2.1e-03` short in objective and `121` kW out in dispatch. Damping shrinks the
cycle without closing it (`a = 0.1` still `3.2e-05` short after 2000 sweeps),
because started *at* the optimum the sweep returns `μ` off by `9.2e-04`: the
optimum is not exactly a fixed point of the linearised map. FilterDDP reaches
`2.2e-10` on the same instance in 16 iterations. That is the `V_xx` gap,
quantified.

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

**The `ddp/paper/` write-up (including `copper_plate_model.tex`) was deleted
2026-09-08** — it encoded a superseded formulation. The authoritative
formulation now lives in the **TPEC repo**, not here. The modeling rationale
that section used to carry (why no `η`, why the terminal target is a penalty,
why `C_B = 0.5` here vs. `≈10⁻⁶·min c^t` in the tADMM paper) survives only in
`ddp/README_FILTERDDP_EXPERIMENT.md`.

**Instance data is now the tADMM profiles** (changed 2026-08-07 at the user's
request). Demand and price come from `envs/tadmm/root_level/config.jl`, shared via
`ddp/examples/power_system/tadmm_profiles.jl`. Three things to know:

- They **resample with T** — `tadmm_cost(3)` is NOT `tadmm_cost(6)[1:3]`. Never
  slice a fixed vector.
- **Committed results are T = 6 only** (T = 3 removed 2026-08-07 at the user's
  request; regenerating another `T` when analysis needs it is expected). But
  **inspect the price before adopting a new `T`**: at T = 3 the samples land
  where `sin` vanishes (`0, π, 2π`), so the price comes out constant and the
  instance carries no arbitrage signal at all. That, plus binding bounds only a
  few kW wide, made T = 3 both the weakest benchmark and the worst behaved
  numerically (26 iterations, trajectory agreement only `1.9e-05`).
- This **closed** the earlier collinearity concern: `r` went from `0.9968` to
  `0.644` at T = 6, thanks to the `−0.8` rad phase offset on the load.

`C_B = 0.05` since 2026-08-07 (was `0.5`). **This is the COPPER-PLATE value and
does not apply to the network cases:** the exported `network_data_*.jls` files
carry `C_B = 1.4e-07` with `dt = 8` (the tADMM regime). Established 2026-09-15 --
the difference is not cosmetic, see the reduced-space section below. `C_B` sets
how far the battery moves:
`P_B^k = (c^k - 2wΔt·s)/(2C_B)`, so the swing scales as `1/C_B` — the bounds are not
the knob. At `0.05` the T = 6 battery swings ±570 kW against a 1623-1998 kW load and
cycles `B` from 2000 down to 1070 kWh: 46% depth of discharge, matching the
`ieee123C_1ph` battery/load ratio of 44%. `cond(Q) = 601`, so the closed-form and
active-set references stay exact to ~1e-14.

Still **not** tADMM's own `C_B = 1e-6·min(c) ≈ 8e-8`: there `cond(Q) ≈ 3.6e8` and the
battery goes purely bound-limited (±340 GW absent bounds) — the tADMM regime, but
useless as a precision reference. An earlier note that this makes `Q` rank-1 was an
overstatement: `Q = 2C_B Δt I + 2w Δt² 11ᵀ` is PD for any `C_B > 0`, rank-1 only in
the limit.

**Plot battery quantities in kW and kWh, never p.u.** (user preference, 2026-08-07).
Base is tADMM's `kVA_B = 1000`: 1 p.u. = 1000 kW, 1 p.u.h = 1000 kWh. The model and
Table I stay in p.u.; the figures convert. Reference asset sizes, for sanity checks:
`ads10A_1ph` 87 kW load / 4.7 kW batt; `ieee123C_1ph` 1163 kW / 507 kW / 2027 kWh;
`ieee123_5poi_1ph` 1163 kW / 1318 kW / 5273 kWh. SOC runs 30%-95% with B_0 at 62.5%.

**Figures are generated, never hand-written** was the rule for the (now-deleted)
`ddp/paper/figures/` — `make_figure_data.jl` read the verified centralized
reference and emitted `balance.csv`, `schedule_interval.csv`, `schedule_soc.csv`
for the `.tex` files to read. Keep the same discipline (generated data, no
inlined coordinates) if/when figures are rebuilt against the TPEC repo's
formulation.

## Reduced-space MPOPF (2026-09-14/15): works, and `C_B` governs whether it is usable

Full write-up and raw data:
[ddp/notes/REDUCED_SPACE_INNER_OPF_FEASIBILITY.md](ddp/notes/REDUCED_SPACE_INNER_OPF_FEASIBILITY.md),
`ddp/results/reduced_space/`. The reduced-space work did not modify the solver
(now `ddp/DDP4OPF.jl`) -- the OCP is assembled from hand-written closures, which
is forced anyway since an Ipopt solve is not automatically differentiable.

The network can be eliminated exactly: FilterDDP optimises battery quantities
only while an inner single-period OPF recovers every network variable, with the
real-power balance duals supplying `dPhi/dP_B` for free. Verified against the
stored full-space FilterDDP solutions to `6.0e-09` (ieee123) and `1.3e-08`
(ieee2522) relative objective.

> **RETRACTED 2026-09-15 -- read before trusting the table below.** Every `T = 3`
> instance has **exactly zero price spread**: the profile generator sampled `sin`
> at `0, pi, 2pi`. With no arbitrage signal the batteries barely move (5% of
> rating at large10k), so the table measures an idle-battery problem. It is a
> valid timing on that instance and **not** evidence the decomposition works on a
> real scheduling problem. On the first instance with a real price signal
> (ieee2522 `T = 12`, 147% spread, batteries up to 97% of rating) the reduced
> method **fails to converge** while full space converges. Why is still open. Use
> `PROFILE_PERIODIC=1` exports (118% spread at `T = 3`); the exporter now warns on
> flat profiles. Full account: `ddp/notes/FINDINGS_2026-09-14_15.md`.

**Measured crossover at `T = 3`, `C_B = 1e-3`, both formulations re-run at the
same `C_B` (degenerate instances -- see retraction above):**

| | full-space | reduced (low-rank) | outcome |
|---|---|---|---|
| ieee123 (`nu` 791 -> 102) | 14.4 s / 46 it | 21.7 s / 26 it | 1.5x slower |
| ieee2522 (`nu` 13358 -> 500) | 107.3 s / 56 it | 107.1 s / 33 it | parity |
| large10k (`nu` 54665 -> 2040) | 1658.2 s / 100 it | 802.5 s / 13 it | **2.07x faster** |

That was read as "the decomposition is a LARGE-system technique"; given the
retraction, that conclusion is **unsupported** until it is re-measured on
non-degenerate instances.
Outer iteration counts run the other way with size (full-space 46/56/100,
reduced 26/33/13), and at large10k inner solves are only 25% of wall -- the
bottleneck has moved to FilterDDP's own outer cost.

**`C_B` decides whether any of this works.** It enters the stage Hessian only as
a perfectly-conditioned diagonal `2*C_B*S^2*dt`, and on ieee2522 `d2Phi` spans
`-0.148 .. 21194`. At the exported `C_B = 1.4e-07` that diagonal is **2.24**, so
there is effectively no damping (cond ~10131) and EVERY cheap curvature model
fails -- battery-term-only, low-rank Nystrom, and floored Nystrom all stall.
Only an exact per-stage Hessian converges, and it is too expensive to win. At
`C_B = 1e-3` the diagonal is 16000 (cond 2.32) and the low-rank Hessian
converges, which is what produces the table above. Size `C_B` from measurement:
`2*C_B*S^2*dt > -lambda_min(d2Phi)` is mandatory for positive definiteness, and
`~ lambda_max/kappa` sets conditioning; `lambda_max` costs ~10 inner solves by
power iteration since `H*d` is one solve.

Two further results worth not re-deriving: the battery-power box is **not**
recourse-feasible in general (true on ieee123, false on ieee2522 -- 8 of 111
dispatches rejected, all undervoltage), and `F_t` is horizon-independent, so `T`
only selects which load/PV snapshots get tested.

## What is and isn't verified (as of 2026-08-07)

FilterDDP is cross-checked against a **centralized JuMP/Ipopt** solve of
`eq:cp_all` on all six T = 6 cases (Stage 8 of the experiment README): exact on
the base instance, worst objective gap `2.2e-10` and worst trajectory gap
`6.9e-09` with bounds, and Ipopt independently certifies the 6g bound set
infeasible. The base instance carries three mutually independent references
(closed form, dense KKT, Ipopt). Stage 9 records the per-iteration trace — 17
iterations, regularisation never firing. **Do not redo any of this.**

Not verified, so don't overclaim: anything with a network (no LinDistFlow, no
BFM); any horizon other than `T = 6`; and **the user's own DDP algorithm, which
has not been compared at all** — that is the pending task below.

## Pending task (do not start until asked)

**Blocked as of 2026-09-08**: this task's shared reference point,
`ddp/paper/sections/copper_plate_model.tex`, was deleted along with the rest
of `ddp/paper/` — that formulation is superseded and the authoritative one now
lives in the **TPEC repo**. Before restarting this task, get the equivalent
problem statement from the TPEC repo and re-approve it as the shared reference
point; do not resurrect the deleted file from git history as a substitute.

Write a side-by-side workflow comparison — the user's exact DDP algorithm
vs. FilterDDP's algorithm — **both grounded in the exact problem statement of
the "dummy paper"** (previously `ddp/paper/sections/copper_plate_model.tex`,
see blocker above). Requirements:

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
