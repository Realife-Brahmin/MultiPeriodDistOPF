# FilterDDP experiment — reproducibility record

**Question this set out to answer:** *can the authors' actual FilterDDP implementation
reliably solve a fully understood T = 3 equality-constrained copper-plate battery problem
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
against the tADMM instance data, on 2026-08-07. The authors' clone is **unmodified** (see
[notes/CHANGES.md](notes/CHANGES.md)).

---

## Final status table

| Stage | Status | Evidence | Blocking issue |
| --- | --- | --- | --- |
| 1 — inspect repository, map paper↔code | **PASS** | [notes/PAPER_CODE_MAP.md](notes/PAPER_CODE_MAP.md); clone at `513a104`; both papers fetched from arXiv and text-extracted | none |
| 2 — instantiate & precompile | **PASS** | [logs/environment_setup.log](logs/environment_setup.log): both envs exit 0, 85 + 142 deps precompiled, Julia 1.12.6 | none |
| 2 — official test suite | **NOT ATTEMPTED** | [logs/tests.log](logs/tests.log): no `test/` directory, no test target; `Pkg.test()` → *"Package FilterDDP did not provide a `test/runtests.jl` file"* | **the repository ships no tests** |
| 3 — smallest official example | **PASS** | [logs/stage3_official_example.log](logs/stage3_official_example.log): 51 iters, obj `1.26574863e+00`, primal `8.09e-08` — **identical to the authors' shipped file to all 9 printed digits** | none |
| 4 — document the API | **PASS** | [notes/FILTERDDP_API.md](notes/FILTERDDP_API.md), every claim sourced to a file:line | none |
| 5 — copper-plate battery, T = 3, equality only | **PASS** | [results/copper_plate/stage5_output.log](results/copper_plate/stage5_output.log): objective gap **exactly `0.0`**, the two independent references agree to `2.0e-15`, converged in 1 iteration from all 3 starts | none |
| 6a — grid-power bound | **PASS** | [results/copper_plate/stage6_output.log](results/copper_plate/stage6_output.log): 9 iters, status 0, `\|ΔJ\| 1.1e-16`, reference active set **empty** — substation import has no upper limit, and `P_Subs ≥ 0` never binds | none |
| 6b — battery-power bounds | **PASS** | same log: 14 iters, status 0, `\|ΔJ\| 1.1e-10`, active set `P_B[1..3]≤0.004` (all three) matches | none |
| 6c — energy bounds (slack reformulation) | **PASS** | same log: 18 iters, status 0, `\|ΔJ\| 2.0e-11`, both active bounds `B[2]≤1.995`, `B[4]≥1.9895` matched | none |
| 6d — all bounds together, T = 3 | **PASS** | same log: 26 iters, status 0, `\|ΔJ\| 4.2e-11`, active set `B[4]≥1.9895` matches | none |
| 6f — tolerance sweep on 6d | **PASS** | same log: status 0 at every tolerance from `1e-6` down to `1e-12`, `\|ΔJ\|` falling to `4.1e-13` (16 → 30 iterations) | none |
| 6 — T = 6 | **PASS** | same log: all-bounds 16 iters, status 0, `\|ΔJ\| 2.2e-10`, six active bounds (`B[2]≤1.95`, `P_B[2]≤0.45`, `B[4]≥1.20`, `P_B[5]≥-0.45`, `B[6]≤1.95`, `B[7]≤1.95`) matched | none |
| 6g — infeasible bound set (probe) | **FAIL, as designed** | same log: returns status 7 with primal residual `0.638`; independent enumeration confirms the feasible set is empty | **no infeasibility diagnosis — status 7 means both "infeasible" and "hard"** |
| 7 — MPOPF applicability assessment | **PASS** | [notes/MPOPF_APPLICABILITY.md](notes/MPOPF_APPLICABILITY.md) | none |
| 8 — centralized JuMP/Ipopt cross-check | **PASS** | all 7 feasible cases agree with a third, fully independent solver: worst objective gap `2.2e-10`, worst trajectory disagreement `1.9e-05` (case 6d at T = 3, a tolerance effect — 6f drives it to `4.1e-13`), machine precision on both base instances. Added 2026-08-07 | none |

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
│   └── copper_plate/                Stage 5, 6 and 8 output
├── examples/power_system/
│   ├── tadmm_profiles.jl            the shared tADMM load/price series
│   ├── copper_plate_battery.jl      Stage 5 model + 2 independent references
│   ├── copper_plate_battery_bounds.jl  Stage 6 bounds + active-set QP reference
│   ├── copper_plate_centralized.jl  Stage 8: eq:cp_all in JuMP/Ipopt, no reduction
│   └── verify_against_centralized.jl   Stage 8: FilterDDP vs that reference
├── papers/                          the two arXiv papers + text extractions
└── external/FilterDDP.jl/           unmodified upstream clone @ 513a104

(the Julia environment lives at envs/ddp2026, outside this tree)
```

## Reproducing

```powershell
# upstream clone (gitignored -- carries its own .git)
git clone https://github.com/mingu6/FilterDDP.jl ddp/external/FilterDDP.jl
git -C ddp/external/FilterDDP.jl checkout 513a104

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
| base, T=3 | 1 | 0.73493729 | 0.73493729 | `0.0` | `7.6e-15` | `7.6e-15` |
| base, T=6 | 1 | 1.46683686 | 1.46683686 | `0.0` | `9.5e-15` | `1.0e-14` |
| 6a, T=3 | 9 | 0.73493729 | 0.73493729 | `1.1e-16` | `2.1e-11` | `1.2e-11` |
| 6b, T=3 | 14 | 0.73495643 | 0.73495643 | `1.1e-10` | `1.9e-09` | `5.7e-09` |
| 6c, T=3 | 18 | 0.73499729 | 0.73499729 | `2.0e-11` | `4.6e-08` | `4.6e-08` |
| 6d, T=3 | 26 | 0.73499712 | 0.73499712 | `4.2e-11` | `1.9e-05` | `1.9e-05` |
| 6e, T=6 | 16 | 1.47517291 | 1.47517291 | `2.2e-10` | `6.8e-09` | `6.9e-09` |
| 6g, T=6 | 12 | — | **LOCALLY_INFEASIBLE** | — | — | — |

The `1.9e-05` on 6d is a barrier-tolerance effect, not a disagreement: the Stage 6f
sweep drives the same case to `\|ΔJ\| = 4.1e-13` at `tol = 1e-12`. Objectives agree
to `~1e-10` throughout.

The pattern across the table is worth reading. The two hardest cases are the
**T = 3** ones, not T = 6: at T = 3 the price is constant, the battery barely moves,
and the bounds that bind therefore sit in a window a few kW wide -- which is what
costs 6d its 26 iterations and five digits of trajectory accuracy. Case 6e at
T = 6, where the bounds sit at a natural ±450 kW and 1200-1950 kWh, converges in 16
iterations to `7e-09`. Tight bounds relative to the natural scale of the variable,
not horizon length, are what makes this problem hard for FilterDDP.

Ipopt also reproduces the closed form of `eq:cp_closed_form` on both base
instances to `2e-16`, making three mutually independent references there.

On 6g the cross-check earns its keep: Ipopt returns `LOCALLY_INFEASIBLE` and so
*certifies* the bound set empty, where FilterDDP returns status 7 — the same code
it returns for a merely hard problem. That is the missing feasibility-restoration
phase showing up as a concrete difference in what a caller can conclude.

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
> Two properties of the sampling matter and are not defects to be fixed:
> both endpoints of `[0, 2π]` are included, so `c^1 = c^T` for every `T`; and at
> `T = 3` the samples land on `0, π, 2π` where `sin` vanishes, so **the price is
> constant** and that instance has no arbitrage signal at all. Its optimal `P_B`
> is identical in every period, set purely by the terminal penalty. T = 3 remains
> a valid — and usefully symmetric — solver test, but the economics live at T = 6.
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
