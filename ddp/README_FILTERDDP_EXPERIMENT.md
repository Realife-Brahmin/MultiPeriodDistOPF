# FilterDDP experiment — reproducibility record

**Question this set out to answer:** *can the authors' actual FilterDDP implementation
reliably solve a fully understood equality-constrained copper-plate battery problem
and reproduce a known reference solution?*

**Answer: yes.** It reproduces two independent references to machine precision on the
equality-only problem, from three different starting trajectories, with exact dynamics
feasibility. With bounds added it stays accurate to `~1e-10`, identifies the correct
active set every time, and certifies convergence at tolerances down to `1e-12`. Two
genuine limitations surfaced and are documented below rather than smoothed over.

> **Objective note.** Stages 5-6 use the tADMM paper's objective — linear energy cost
> `Σ c^t·P_Subs^t·Δt` plus a quadratic battery-throughput charge `C_B·(P_B^t)²·Δt` — not a
> quadratic cost on grid import. An earlier revision used `a^t(P_Subs^t)² + b^t P_Subs^t`;
> under that objective Stage 6d exited with a line-search failure and could not certify
> below `~1e-9`. **That was an artefact of the objective choice, not of FilterDDP**: with
> `C_B` carrying the reduced-Hessian curvature instead, the same case converges cleanly to
> `1e-12`. Both are recorded because the contrast is itself the useful result.

Stages 1-7 were run on 2026-08-03. Stage 8 was added, and every result regenerated
against the tADMM instance data, on 2026-08-07. The authors' clone was **unmodified through Stage 10**. On 2026-08-12 it received
exactly one patch, [patches/per_stage_data.patch](patches/per_stage_data.patch),
which lets the package accept a distinct objective and constraint function per
stage. See [notes/CHANGES.md](notes/CHANGES.md) for what it changes and the
regression evidence that the default path is untouched.

---

## Final status table

| Stage | Status | Evidence | Blocking issue |
| --- | --- | --- | --- |
| 1 — inspect repository, map paper↔code | **PASS** | [notes/PAPER_CODE_MAP.md](notes/PAPER_CODE_MAP.md); clone at `513a104`; both papers fetched from arXiv and text-extracted | none |
| 2 — instantiate & precompile | **PASS** | [logs/environment_setup.log](logs/environment_setup.log): both envs exit 0, 85 + 142 deps precompiled, Julia 1.12.6 | none |
| 2 — official test suite | **NOT ATTEMPTED** | [logs/tests.log](logs/tests.log): no `test/` directory, no test target; `Pkg.test()` → *"Package FilterDDP did not provide a `test/runtests.jl` file"* | **the repository ships no tests** |
| 3 — smallest official example | **PASS** | [logs/stage3_official_example.log](logs/stage3_official_example.log): 51 iters, obj `1.26574863e+00`, primal `8.09e-08` — **identical to the authors' shipped file to all 9 printed digits** | none |
| 4 — document the API | **PASS** | [notes/FILTERDDP_API.md](notes/FILTERDDP_API.md), every claim sourced to a file:line | none |
| 5 — copper-plate battery, equality only | **PASS** | [results/copper_plate/stage5_output.log](results/copper_plate/stage5_output.log): objective gap **exactly `0.0`**, the two independent references agree to `2.0e-15`, converged in 1 iteration from all 3 starts | none |
| 6a — grid-power bound | **PASS** | [results/copper_plate/stage6_output.log](results/copper_plate/stage6_output.log): 10 iters, status 0, `\|ΔJ\| 0.0`, reference active set **empty** — substation import has no upper limit, and `P_Subs ≥ 0` never binds | none |
| 6b — battery-power bounds | **PASS** | same log: 12 iters, status 0, `\|ΔJ\| 7.3e-11`, active set `P_B[2]≤0.45`, `P_B[5]≥-0.45` matches | none |
| 6c — energy bounds (slack reformulation) | **PASS** | same log: 14 iters, status 0, `\|ΔJ\| 1.5e-10`, active bounds on both ends of `1.20 ≤ B ≤ 1.95` matched | none |
| 6d — all bounds together | **PASS** | same log: 16 iters, status 0, `\|ΔJ\| 2.2e-10`, six active bounds (`B[2]≤1.95`, `P_B[2]≤0.45`, `B[4]≥1.20`, `P_B[5]≥-0.45`, `B[6]≤1.95`, `B[7]≤1.95`) matched | none |
| 6f — tolerance sweep on 6d | **PASS** | same log: status 0 at every tolerance from `1e-6` down to `1e-12`, `\|ΔJ\|` falling monotonically (10 → 20 iterations) | none |
| 6g — infeasible bound set (probe) | **FAIL, as designed** | same log: returns status 7 with primal residual `0.638`; independent enumeration confirms the feasible set is empty | **no infeasibility diagnosis — status 7 means both "infeasible" and "hard"** |
| 7 — MPOPF applicability assessment | **PASS** | [notes/MPOPF_APPLICABILITY.md](notes/MPOPF_APPLICABILITY.md) | none |
| 8 — centralized JuMP/Ipopt cross-check | **PASS** | all 7 feasible cases agree with a third, fully independent solver: worst objective gap `2.2e-10`, worst trajectory disagreement `6.9e-09`, exact on the base instance. Added 2026-08-07 | none |
| 9 — per-iteration convergence trace | **PASS** | [results/copper_plate/stage9_convergence.log](results/copper_plate/stage9_convergence.log): case 6d converges in **17** iterations; regularisation never fires. Added 2026-08-07 | none |

Status labels are used strictly: no stage is PASS without a number behind it.

---

## What is here

```
ddp/
├── README_FILTERDDP_EXPERIMENT.md   this file
├── notes/
│   ├── PAPER_CODE_MAP.md            paper equations -> source lines; 9 paper/code divergences
│   ├── FILTERDDP_API.md             the minimum API, every claim sourced
│   ├── CHANGES.md                   what was changed (nothing upstream) and my own corrections
│   └── MPOPF_APPLICABILITY.md       the skeptical assessment
├── logs/                            unedited terminal output for every stage
├── results/
│   ├── official_example/            our run vs the authors' shipped numbers
│   └── copper_plate/                Stage 5, 6, 8 and 9 output
├── examples/power_system/
│   ├── tadmm_profiles.jl            the shared tADMM load/price series
│   ├── copper_plate_battery.jl      Stage 5 model + 2 independent references
│   ├── copper_plate_battery_bounds.jl  Stage 6 bounds + active-set QP reference
│   ├── copper_plate_centralized.jl  Stage 8: eq:cp_all in JuMP/Ipopt, no reduction
│   ├── verify_against_centralized.jl   Stage 8: FilterDDP vs that reference
│   └── convergence_trace.jl         Stage 9: per-iteration trace + figure data
├── papers/                          the two arXiv papers + text extractions
└── external/FilterDDP.jl/           unmodified upstream clone @ 513a104

(the Julia environment lives at envs/ddp2026, outside this tree)
```

## Reproducing

```powershell
# upstream clone (gitignored -- carries its own .git)
git clone https://github.com/mingu6/FilterDDP.jl ddp/external/FilterDDP.jl
git -C ddp/external/FilterDDP.jl checkout 513a104
git -C ddp/external/FilterDDP.jl apply ../../patches/dynamic_network_scaling.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/factor_bound_sensitivities.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/no_copy_update_rule.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/in_place_kkt_rhs.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/reuse_stage_rule_buffers.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/reuse_kkt_rhs_workspace.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/type_constraint_residual_vector.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/in_place_B_assembly.patch
git -C ddp/external/FilterDDP.jl apply ../../patches/active_B_rows.patch

# the two papers (gitignored -- not ours to redistribute)
mkdir ddp/papers
curl -L -o ddp/papers/Xu_2026_FilterDDP.pdf https://arxiv.org/pdf/2504.08278
curl -L -o ddp/papers/Xu_2026_FilterDDP_Global_Convergence.pdf https://arxiv.org/pdf/2606.01487
# text extractions used for the quotes in notes/:
#   pdftotext -layout ddp/papers/<file>.pdf ddp/papers/<file>.txt

# environment -- envs/ddp2026 is now the single environment for all DDP work,
# FilterDDP and centralized alike. See envs/ddp2026/README.md.
julia --startup-file=no --project=envs/ddp2026 -e 'using Pkg; Pkg.instantiate()'

# Stage 3 — the authors' example (writes into their results/ dir; back it up first)
cd ddp/external/FilterDDP.jl/experiments/filterddp
julia --startup-file=no --project=. double_integrator.jl

# Stages 5 and 6
julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/copper_plate_battery.jl
julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/copper_plate_battery_bounds.jl

# Stage 8 — centralized cross-check (produce the reference, then diff against it)
julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/copper_plate_centralized.jl
julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/verify_against_centralized.jl

# Stage 9 — per-iteration convergence trace (also writes the paper's figure data)
julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/convergence_trace.jl
```

All logs were regenerated on 2026-08-07 after the instance data was switched to
the tADMM profiles (see below), so the numbers in the table above are current.
The original `ddp/env` was removed the same day: its `Manifest.toml` was never
tracked, so it recorded nothing `envs/ddp2026/Project.toml` does not, and keeping
two environments for one stack invites drift. It remains in git history.

Environment: Windows 11 (build 22621), Intel Core Ultra i7-1255U, **Julia 1.12.6**,
Symbolics 7.35.0, StaticArrays 1.9.18. The paper benchmarked on Julia 1.13.0-beta3, but
`Project.toml` only requires `julia = "1.12.4"` — 1.12.6 satisfies it and nothing failed,
so the beta was never installed.

## Stage 3 reproduction, in detail

| | iterations | status | objective | primal infeasibility | wall (ms) |
| --- | --- | --- | --- | --- | --- |
| authors' shipped file | 51 | true | `1.26574863e+00` | `8.09173759e-08` | 1.4 |
| our run | 51 | true | `1.26574863e+00` | `8.09173759e-08` | 1.6 |

This is a real convergence, not just a clean exit: the final iterate has
`du_inf = 6.5e-14`, `cs_inf = 1.1e-08`, `pr_inf = 8.1e-08` against a tolerance of `1e-7`,
and `solve!` only returns status 0 from the explicit KKT test at `src/solve.jl:30`.
The wall time is post-compilation (the script's `@benchmark` median over repeated solves);
total process time including Julia startup and symbolic model construction was 50.8 s.

Two incidental observations from the verbose trace: the regularisation fired on only 3 of
51 iterations (`lg(reg) = -4`), and the printed `alpha`/`ls` columns are **one iteration
stale**, because `iteration_status` is called before `forward_pass!` (`src/solve.jl:38,40`).

## The Stage 5 model in one paragraph

`N = 3` stages, state `x = (e, τ)` — battery energy plus an explicit time index — and
control `u = (pg, pb)` with `pb > 0` meaning discharge. Dynamics `e_{t+1} = e_t - η·pb_t`,
`τ_{t+1} = τ_t + 1`. One equality per stage, `pg_t + pb_t - d(τ_t) = 0`. Quadratic
generation cost. The time index is in the state because **FilterDDP holds one constraint
function and one stage cost for the whole horizon**, so time-varying demand and prices have
nowhere else to live; `d(τ)`, `a(τ)`, `b(τ)` are Lagrange interpolants through the nodes,
exact at `τ = 1,2,3` and smooth. The terminal energy target is a *cost*, not a constraint —
see the next section. Verified against a closed-form solution **and** an independently
constructed dense-KKT solve of the same QP; the two references agree with each other to
`2.0e-15`.

## Two limitations worth knowing before going further

**1. Terminal equality constraints are not representable.** `c` is the same function at
every stage (`src/ocp/ocp.jl:9`), so a terminal-only equality cannot be written; and a
constraint engineered to vanish at non-terminal stages would have `∇_u c = 0` there, which
makes the backward pass's `lu(AY)` singular — silently, since there is no rank check.
There is a clean partial workaround: the **terminal objective receives `u_N`**, so the
state one step past the horizon, `e_{T+1} = e_N - η·pb_N`, is exactly expressible as a
terminal *penalty*. That is what the Stage 5 model does, and it is exact algebra rather
than an approximation. But it is a penalty, not a constraint, and the reported terminal
energy residual (`2.8e-01`) reflects that honestly.

**2. There is no feasibility restoration phase**, even though **both** convergence theorems
assume one. P2's `Assumptions G` and `Assumptions B` each open by assuming *"the feasibility
restoration phase in step 9 always terminates successfully"*; the code simply sets status 7
and stops (`src/forward_pass.jl:40`). Stage 6g shows the practical consequence: a provably
infeasible problem and a merely difficult one return the same status code, and the caller
has to audit residuals to tell them apart.

Two smaller ones: `Options.s_L` defaults to `2.3` while P2 Theorem 2 requires `s_L = 1`, so
the interior-point path under default options is not the algorithm that theorem covers; and
seven `Options` fields (`quasi_newton`, `reset_cache`, `ineq_dual_init`, `κ_c`, `δ_c`,
`κ_Σ`, `linsolve_tol`) are declared but never read anywhere in `src/`. Full list of nine
paper/code divergences in [notes/PAPER_CODE_MAP.md](notes/PAPER_CODE_MAP.md) §4.

## Stage 8 — the centralized cross-check

Added 2026-08-07. `copper_plate_centralized.jl` solves `eq:cp_all` of
[paper/sections/copper_plate_model.tex](paper/sections/copper_plate_model.tex)
directly in JuMP + Ipopt, in the paper's own variables, with no reduction. It is
independent of the DDP model in three ways that could each hide a bug: it indexes
the true per-period `c^t` and `p_L^t` instead of recovering them from Lagrange
interpolants, it imposes the energy bound directly on `B` instead of routing it
through the slack control of `eq:ocp_slack`, and it carries `B^{T+1}` as an
explicit variable. Agreement therefore tests the *reformulation*, not just the
solver.

| case | iters | J (FilterDDP) | J (Ipopt) | `\|ΔJ\|` | max `\|ΔPsub\|` | max `\|ΔB\|` |
| --- | --- | --- | --- | --- | --- | --- |
| base | 1 | 1.46683686 | 1.46683686 | `0.0` | `9.5e-15` | `1.0e-14` |
| 6a | 10 | 1.46683686 | 1.46683686 | `2.2e-16` | `6.6e-11` | `8.4e-11` |
| 6b | 12 | 1.46829293 | 1.46829293 | `7.3e-11` | `3.4e-09` | `6.3e-09` |
| 6c | 14 | 1.47494196 | 1.47494196 | `1.5e-10` | `4.1e-09` | `4.2e-09` |
| 6d | 16 | 1.47517291 | 1.47517291 | `2.2e-10` | `6.8e-09` | `6.9e-09` |
| 6g | 12 | — | **LOCALLY_INFEASIBLE** | — | — | — |

Objectives agree to `~1e-10` and trajectories to `~7e-09` across every case, and
the agreement tightens monotonically as bounds are layered on rather than
degrading. An earlier revision also carried a T = 3 variant of each case; it was
removed on 2026-08-07 in favour of T = 6 alone. That variant was both the
weakest benchmark (a constant price, so no arbitrage signal) and the worst
behaved numerically (bounds a few kW wide, 26 iterations, trajectory agreement
only `1.9e-05`) -- dropping it is why the worst disagreement above is now three
orders of magnitude smaller.

> **Instance data: the tADMM profiles.** As of 2026-08-07 demand and price are the
> same series the tADMM experiments use, taken from
> `envs/tadmm/root_level/config.jl` and shared through
> [examples/power_system/tadmm_profiles.jl](examples/power_system/tadmm_profiles.jl):
> `p_L^t = 2·λ^t` with `λ^t = 0.8 + 0.1(sin(θ−0.8)+1)` and
> `c^t = 0.08 + 0.06(sin θ + 1)`, sampled at `θ ∈ range(0, 2π, T)`.
>
> This **resolves** the collinearity weakness of the previous hand-picked series
> (Pearson `r = 0.9968`). The `−0.8` rad phase offset on the load gives `r = 0.644`
> at `T = 6` and `0.682` at `T = 24` — the range real day-ahead price/load
> correlation actually occupies, and enough separation that transposing `c` and
> `p_L` would now show up.
>
> Two properties of the sampling matter. Both endpoints of `[0, 2π]` are
> included, so `c^1 = c^T` for every `T`. And **the price should be inspected
> before adopting any new `T`**: at `T = 3` the samples land on `0, π, 2π` where
> `sin` vanishes, making the price constant and the instance free of any
> arbitrage signal. That is why T = 3 was dropped.
>
> The profiles **resample** with `T`: `tadmm_cost(3)` is not `tadmm_cost(6)[1:3]`.
> Any code that slices a fixed vector is wrong.
>
> One thing was *not* adopted: tADMM sets `C_B = 1e-6·min(c) ≈ 8e-8`, which here
> would be fatal. `C_B` is exactly the reduced-Hessian curvature, and at that
> magnitude `Q = 2·C_B·Δt·I + 2·w·Δt²·11ᵀ` degenerates to a rank-1 matrix, so the
> closed-form and active-set references stop being unique and the benchmark loses
> the very thing that makes it a reference. `C_B = 0.5` is kept. The consequence
> is that the battery here is cost-limited rather than bound-limited as it is in
> tADMM, so its economic swing is small (`P_B ∈ [-0.055, 0.059]` at `T = 6`) and
> the bounds chosen to bind are correspondingly tight.

## Stage 9 — how the iterations actually go

`convergence_trace.jl` records what happens inside each iteration. FilterDDP keeps
**no iteration history** — `SolverData` (`src/solver_data.jl:8-32`) holds only the
current iterate — so the trace is captured from the verbose stream and parsed
back. The clone stays unmodified.

| | 6d, all bounds |
| --- | --- |
| iterations | **17** (`k = 0..16`) |
| `pr_inf` | `1.0e+00` → `4.4e-16` |
| `du_inf` | `9.2e-01` → `5.6e-13` |
| `cs_inf` | `1.0e+00` → `3.7e-11` |
| objective | `0.840000` → `1.4751729` |
| centralized `J*` | `1.4751729098` |
| final `\|J - J*\|` | `2.3e-10` |
| regularisation fired | **0** of 17 |
| iterations that backtracked | 8 |

The shape is worth stating: **feasibility is essentially free, and the iteration
count is set by the barrier schedule.** Both `pr_inf` and `du_inf` fall to
machine precision by `k = 4` and simply stay there; every remaining iteration is
`μ` stepping down and dragging `cs_inf` with it.

The objective is plotted against the centralized JuMP/Ipopt baseline `J*`, which
the tracer reads out of `centralized_reference.csv` rather than carrying as a
literal, so the two cannot drift. `J_k` is **not monotone**: at T = 6 it
overshoots `J*` near `k = 5` and returns from above, which is expected — the
filter line search decreases the barrier objective, not `J`. On a linear axis the
trace is visually converged by `k ≈ 5`; the log-scaled gap is what shows the
remaining eight orders of magnitude being earned one barrier step at a time.

Two things the raw output gets wrong, both handled in the tracer rather than
passed along:

1. **The printed halves are one iteration out of step.** `iteration_status` is
   called before `forward_pass!` (`src/solve.jl:38,40`), so the `alpha` and `ls`
   on row `k` belong to iteration `k-1`. A table that pairs a backward pass with
   the wrong forward pass is simply wrong, so the parser de-staggers them.

2. **`ls` is not a line-search trial count**, despite the header. `forward_pass!`
   backtracks from three places; the one at `src/forward_pass.jl:17` (rollout
   failure, i.e. the fraction-to-boundary rule) halves the step **without**
   incrementing `data.l`, while lines 25 and 37 increment. Every backtrack in
   both traces came through line 17, which is why the raw output shows
   `alpha = 0.5` alongside `ls = 0` — a contradiction on its face. The reliable
   quantity is `alpha`; the tracer reports the true backtrack count as
   `log2(1/alpha)` and relabels `l` as what it is, filter/decrease rejections only.

## Where this leaves the MPOPF question

Detail in [notes/MPOPF_APPLICABILITY.md](notes/MPOPF_APPLICABILITY.md). The compressed
version: the state dimension for an MPOPF is just the number of storage devices, which is
genuinely attractive — but every algebraic network variable must become a *control*, and
the backward pass factorises **dense `nu × nu` `StaticArrays`** at every stage. The authors
themselves recommend switching implementations above `nu = 15`. A BFM period has hundreds
to tens of thousands of algebraic variables. DDP decomposes in time, but a distribution
MPOPF coupled only through SOC is *already* nearly time-separable — that is why temporal
ADMM works on it — so this would trade one temporal decomposition for another while making
each per-period subproblem strictly harder.

Recommended order if pursued: (1) per-stage costs and constraints, (2) sparse/LAPACK stage
algebra starting from `FilterDDPLAPACK.jl`, (3) a feasibility restoration phase, (4) only
then a LinDistFlow feeder.

**Stopped here, before any distribution-network architecture, as instructed.**
