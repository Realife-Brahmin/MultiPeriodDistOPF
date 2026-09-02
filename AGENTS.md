# Repo-specific context for Codex

This file is committed so any machine's Codex session starts from the
same ground truth. Session-local memory (`~/.Codex/.../memory/`) does not
travel between machines — this file does. Keep it updated when a session
establishes something a future session, on any machine, would need.

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
  scheme for copper-plate MPOPF (committed 2025-11-12, predates any Codex
  involvement — verified by git blame, not AI-authored). A clean minimal
  restatement of it lives at
  [ddp/examples/power_system/user_ddp.jl](ddp/examples/power_system/user_ddp.jl).

**The difference is the order of the cost-to-go model, not the family:**

| | `envs/ddp/` (user's own) | `ddp/` (FilterDDP, the paper) |
|---|---|---|
| Backward information | Passes `μ[t]`, the dynamics-constraint dual, backward one stage per **outer forward sweep** (`envs/ddp/root_level/ddp_copperplate.jl:539-550`, `mu_prev`/`mu_coupling`) | Backward pass sweeps `t=N→1` **within one iteration**, building both `V_x` and `V_xx` (Riccati recursion, `backward_pass.jl` in the FilterDDP clone) |
| Order of approximation | First-order only: the coupling term `μ[t+1]·(B[t+1]−B[t]+Δt·P_B[t+1])` is linear in `B[t]` — no curvature crosses a stage boundary | Second-order: full local quadratic model of the cost-to-go |
| Per-stage solve | Calls an external solver (Gurobi/Ipopt) per stage per sweep | In-house KKT solve inside the backward pass, with filter line-search globalization; compact models retain the dense null-space route, while network models use the sparse saddle-point route added 2026-08-22 |
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

**Three-dual rolling mean tested 2026-08-23** at Rahul's suggestion. The
optional `dual_average_window=3` update in `user_ddp.jl` averages the raw dual
vectors produced by the current sweep and the preceding two sweeps; it does not
average the state trajectory. This materially reduces the cold-start error but
does not converge. The original update has a period-2 cycle, objective gap
`6.35e-4`--`2.09e-3`, and maximum battery-power error `52.7`--`120.6` kW. The
three-dual mean settles into a period-4 cycle, objective gap
`1.25e-4`--`2.94e-4`, and maximum battery-power error `15.0`--`34.0` kW. Its
result is identical at 300 and 2000 sweeps. Reproduction script and data are
`ddp/examples/power_system/first_order_dual_average_experiment.jl` and
`ddp/results/copper_plate/first_order_dual_average.csv`.

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

`C_B = 0.05` since 2026-08-07 (was `0.5`). `C_B` sets how far the battery moves:
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

**Figures are generated, never hand-written.** `ddp/paper/figures/make_figure_data.jl`
reads the verified centralized reference and emits `balance.csv`,
`schedule_interval.csv`, `schedule_soc.csv`; the `.tex` files read those tables.
Do not inline coordinates — `schedule_fig.tex` did, and went silently stale when the
instance data changed.

## What is and isn't verified (as of 2026-08-22)

FilterDDP is cross-checked against a **centralized JuMP/Ipopt** solve of
`eq:cp_all` on all six T = 6 cases (Stage 8 of the experiment README): exact on
the base instance, worst objective gap `2.2e-10` and worst trajectory gap
`6.9e-09` with bounds, and Ipopt independently certifies the 6g bound set
infeasible. The base instance carries three mutually independent references
(closed form, dense KKT, Ipopt). Stage 9 records the per-iteration trace — 17
iterations, regularisation never firing. **Do not redo any of this.**

The loss-aware BFM-SOCP network transcription is also verified against stored
centralized references. The sparse-KKT FilterDDP path solves ADS10, IEEE123C,
and corrected IEEE2522C at `T = 3`. IEEE2522C now converges in 56 iterations
and 233.098 s with objective gap `7.342e-4` and equality residual `5.631e-10`;
this supersedes the earlier zero-iteration dense-QR failure without erasing it
from the chronology. Full results are in
`ddp/results/network_filterddp/README.md` and
`sparse_kkt_filterddp_T3.csv`. IEEE123C has separately exported horizon data
through `T = 96`; do not conflate exported/centralized horizon results with a
sparse FilterDDP solve at every horizon.

The user's own first-order DDP has now been compared, as quantified above.
Sparse FilterDDP on IEEE2522C is now verified through `T = 96`: `T = 6`
converges in 82 iterations/836.933 s, `T = 12` in 79 iterations/1902.103 s,
`T = 24` in 84 iterations/3068.228 s, `T = 48` in 84 iterations/4539.491 s,
and `T = 96` in 93 iterations/9359.596 s. The `T = 96` objective gap against
the stored centralized reference is `6.084e-3` (`6.96e-7` relative), and its
maximum equality residual is `9.892e-8`; all horizons terminate strictly and
have complete iteration traces. A symbolic-factorization cache
was tested and rejected because it slowed the identical `T = 3` run from
233.098 s to 295.296 s. Large10kC sparse FilterDDP is now tested through
`T = 12`: `T = 3` converged strictly in 115 iterations/7123.782 s, and `T = 6`
in 128 iterations/23490.218 s. `T = 12` reached the 200-iteration limit after
69951.536 s with objective 2976105.092790690251, equality/primal residual
`6.821e-13`, complementarity `1.053e-9`, and dual infeasibility `5.200e-7`;
therefore it misses FilterDDP's strict `1e-7` stationarity tolerance but meets
the separately reported practical `1e-6` large-network tolerance. Complete
iteration traces are committed under `ddp/results/network_filterddp/`.

The large10k matrix path is only partly sparse. The `96968 x 96968` KKT
coefficient is sparse and uses generic sparse LU, but its `96968 x 1021`
right-hand side is explicitly dense. Per-stage `beta` (`54665 x 1020`),
`omega` (`42303 x 1020`), and both bound-dual sensitivity matrices are also
dense and persist across the horizon (about 1.57 GiB/stage), while the value
update includes dense `beta' * B`. This explains the near-linear memory growth
and is at least as important a scaling target as sparse-LU fill-in. Still
unverified: large10kC at `T >= 24`, and whether a
symmetric-indefinite, feeder-structured, selected-RHS, or matrix-free solve
materially reduces factorization, dense-sensitivity, and storage costs.

**Hessian bottleneck diagnosis (2026-08-31):** warm `T = 3` timing shows that
the expensive "Hessian work" is predominantly propagation of dense
second-order sensitivities, not evaluation of second-derivative callbacks. On
IEEE2522C the second-order callbacks totalled about 0.046 s, versus 0.791 s for
the `nx+1`-column solve and 1.585 s for policy/value updates. On large10k the
dense RHS solve and value update together consumed about 61% of the measured
backward pass, while sparse LU factorization was only 3.4%. Single-process
MUMPS was slower than UMFPACK on captured IEEE123 and IEEE2522 stage systems.
See `ddp/notes/FILTERDDP_HESSIAN_BOTTLENECK_PROFILE.md`; do not describe a
factorizer swap alone as the likely scalability solution.

**Allocation profile (2026-08-31):** the large10k sparse KKT coefficient is
only 7--23 MiB at `T = 3`, while its dense RHS and solution are about 755 MiB
each. A retained stage update is about 1607 MiB: `beta` 425 MiB, `omega`
329 MiB, and the two bound-sensitivity matrices 851 MiB. Constructing the
update creates about 6408 MiB of cumulative temporary allocation traffic per
stage (not simultaneous peak RAM). `Vxx` itself is only 7.9 MiB; its importance
is that it generates dense sensitivities, not that its own storage dominates.
Raw rows are in `backward_pass_memory_profile_T3.csv`. The scoped next
experiment is `ddp/notes/STAGEWISE_IPOPT_ORACLE_INTERFACE.md`: first verify a
single IEEE123 stage and selected sensitivity columns, rather than attempting
a full hybrid solver immediately.

**Stagewise Ipopt oracle test (2026-08-31):** this IEEE123 `T = 3`, stage-1
test now passes. Ipopt was given the same entering battery state and captured
quadratic future-cost message as FilterDDP. Central differences of complete
Ipopt stage re-solves reproduce FilterDDP's `beta*d` within 0.99--2.74% and
`omega*d` within 0.04--0.17% for two individual-battery directions and one
aggregate direction; results are stable for steps `1e-4` through `1e-6`.
Raw equality multiplier levels are not comparable because the slack model is
dual-degenerate/ill-scaled, but their directional sensitivities agree after
mapping JuMP's opposite sign convention. This validates Ipopt as a local
oracle; it does **not** solve the scaling problem, because finite-differencing
all columns would require `2*nx` nonlinear solves per stage. The next target is
reuse of Ipopt's linearized KKT system or a structured selected-direction map.

**Sensitivity reuse/compression test (2026-08-31):** ordinary Ipopt.jl does
not expose Ipopt's internal factorization for extra right-hand sides. The local
Ipopt artifact does include `ipopt_sens.exe`, `libsipopt-3.dll`, and sIPOPT C++
headers; a Julia-facing sIPOPT bridge is therefore possible but requires an
AMPL-suffix workflow, a C++ bridge, or a new Ipopt.jl wrapper. IEEE123 stage-1
compression rejects a simple low-rank/selected-battery shortcut for `beta`:
the optimal ranks for 10%, 5%, and 1% Frobenius error are 38, 48, and 51 out
of 51. Forty pivot-selected columns still have 10.4% global error and 42.8%
worst-column error. `omega` is more compressible, but that alone does not
remove the dominant control and bound sensitivities. See
`ddp/notes/IPOPT_SENSITIVITY_REUSE_AND_MAP_COMPRESSION.md`; the next faithful
prototype is a small sIPOPT bridge for the three already-validated directions.

**sIPOPT bridge spike (2026-08-31):** the official C++ example confirms that
sIPOPT retains Ipopt's algorithm/KKT objects and runs sensitivity steps after
one TNLP solve. An MPOPF bridge should introduce the entering battery state as
`nx` fixed parameter variables/constraints and attach `sens_init_constr`,
`sens_state_1`, and `sens_state_value_1` metadata. The bundled headers/import
libraries compile, but the lab PC's only compiler is Strawberry GCC 8.3 and
the resulting executable crashes at the boundary with Julia 1.12's much newer
MinGW C++ runtime after loader dependencies are resolved. Treat this as a
Windows ABI/toolchain blocker, not a failed sensitivity algorithm. Retry with
a BinaryBuilder-compatible MinGW toolchain or on Linux; first reproduce the
three validated IEEE123 directions before scaling.

**Exact bound-sensitivity storage reduction (2026-08-31):** FilterDDP's two
retained `nu x nx` bound-dual maps were merely row-scaled copies of `beta`:
`ζl = -Σ_L .* beta` and `ζu = Σ_U .* beta`. Store only the two length-`nu`
scale vectors and evaluate their action from the already-computed `beta*δx`
during rollout. On large10k this reduces bound data from 850.804 MiB to
0.834 MiB and the complete update rule from 1606.982 MiB to 757.012 MiB per
stage, exactly preserving the IEEE123 `T = 3` result. See
`ddp/notes/FILTERDDP_FACTORED_BOUND_SENSITIVITIES.md` and apply
`ddp/patches/factor_bound_sensitivities.patch` after the dynamic-network patch.
Full strict regressions preserve iteration counts and solutions. large10k
`T = 3` improved from 7123.782 s to 6660.908 s (6.50%), and IEEE2522 `T = 12`
from 1902.103 s to 1183.039 s (37.80%); the `T = 3` small/medium cases were
slightly slower. IEEE2522 `T = 24` then improved from 3068.228 s to 2330.122 s
(24.06%) with the same 84 iterations, while removing about 1.19 GiB across the
horizon. Treat the exact storage reduction as established and these single-run
runtime percentages as promising rather than statistical estimates.

**No-copy update-rule construction (2026-08-31):** after the exact bound-map
factorization, the backward pass's owned `beta`, `omega`, and vector arrays
were copied a second time by the converting `UpdateRule` constructor. A
type-specific ownership constructor removes that copy while retaining the
generic compatibility fallback. On large10k this reduces warm-stage update
allocation from 4282.176 MiB to 3525.163 MiB, exactly one 757.013-MiB update
rule (17.68%), with unchanged persistent storage and an identical IEEE123
full regression. Apply `ddp/patches/no_copy_update_rule.patch` last; see
`ddp/notes/FILTERDDP_NO_COPY_UPDATE_RULE.md`.
Full IEEE2522 strict regressions isolate a consistent modest runtime gain over
the factored-only implementation: 3.29% at `T = 3`, 4.40% at `T = 12`, and
3.47% at `T = 24`, with identical stored traces. This establishes the likely
effect size; do not spend another two hours on large10k unless a large-system
runtime number is specifically required.

**In-place sparse-KKT RHS solve (2026-09-01):** the sparse backward pass used
to retain its dense right-hand side while allocating a second dense solution
matrix. It now uses `ldiv!` to overwrite the RHS, while explicit KKT-capture
mode alone preserves a pre-solve copy. On large10k `T = 3`, warm-stage solve
allocation falls exactly from 1510.686 MiB to 755.343 MiB (50%), removing one
755.343-MiB buffer; update allocation is unchanged. The IEEE123 strict full
trajectory is identical, and a captured system reproduces its solution with
residual `3.19e-10`. Apply `ddp/patches/in_place_kkt_rhs.patch` after the
no-copy patch. See `ddp/notes/FILTERDDP_IN_PLACE_KKT_RHS.md`. This was a
bounded memory probe, not a full-run runtime benchmark.

A subsequent strict cold-start IEEE2522 `T = 12` regression preserves the
complete 79-iteration trace exactly and solves in 1092.507 s, versus
1145.392 s for no-copy alone (4.62% faster) and 1183.039 s for factored-only
(7.65% faster). Objective and final residuals are identical; all practical
`1e-6` criteria are first met at iteration 76. See
`ddp/results/network_filterddp/in_place_kkt_rhs_runtime.csv`. Treat the timing
percentages as single-run measurements, not statistical estimates.

**Sparse stage-rule buffer reuse (2026-09-01):** after the in-place KKT solve,
the sparse path still copied its solved blocks through large intermediate
matrices and replaced an already allocated stage rule. It now copies directly
into that rule's existing arrays. On large10k `T = 3`, warm-stage update
allocation falls from 3525.163 MiB to 2011.139 MiB, saving 1514.024 MiB
(42.95%) per construction while retained storage remains 757.012 MiB/stage.
IEEE123 and IEEE2522 `T = 12` preserve their exact traces. IEEE2522 `T = 12`
runtime is 1039.905 s, a further 4.82% below the in-place-RHS run. Apply
`ddp/patches/reuse_stage_rule_buffers.patch` last; see
`ddp/notes/FILTERDDP_STAGE_RULE_REUSE.md`. Timing is a single-run measurement.

**Reusable KKT RHS workspace (2026-09-01):** the sparse path now lazily
allocates one dense RHS per backward sweep, fills its four blocks in place,
and reuses it across stages and regularization retries. This replaces per-stage
matrix concatenation. On large10k `T = 3`, warm-stage KKT assembly allocation
falls from 1760.525 MiB to 68.569 MiB, saving 1691.956 MiB (96.10%). IEEE123
and IEEE2522 `T = 12` preserve exact traces. IEEE2522 `T = 12` runtime is
962.181 s, 7.47% below stage-rule reuse and 49.41% below the original sparse
1902.103-s run. Apply `ddp/patches/reuse_kkt_rhs_workspace.patch` last; see
`ddp/notes/FILTERDDP_REUSABLE_KKT_RHS_WORKSPACE.md`. Timings are single runs.

**Typed constraint residual dispatch fix (2026-09-01):** post-optimization
profiling found that network constraint residual `c` was `Vector{Any}` while
`omega` was `Matrix{Float64}`. Consequently, `omega' * c` used a generic path,
allocating about 1975 MiB and taking about 4.0 s per warm large10k stage.
Constructing `c` as `Vector{T}` restores BLAS without changing its values. The
warm-stage update falls from 2011.138 MiB to 33.909 MiB allocation (98.31%),
and that product to about 0.010 s/0.008 MiB. IEEE123 retains 48 iterations;
IEEE2522 `T = 12` retains its displayed 79-row trace and solves in 737.043 s,
23.40% below the preceding version and 61.25% below the original sparse run.
Apply `ddp/patches/type_constraint_residual_vector.patch` last; see
`ddp/notes/FILTERDDP_TYPED_CONSTRAINT_VECTOR.md`. Timings are single runs.

**In-place mixed-derivative assembly (2026-09-01):** after the dispatch fix,
the expression forming `B = lux + ux_tmp*fx + fux + cux` still created a new
dense `nu x nx` matrix for every `+`, although the three added derivative
matrices are sparse. The code now forms the dense product once and accumulates
the sparse terms into it. Warm large10k stage algebra allocation falls from
2277.949 MiB to 1001.743 MiB (56.02%); the `B` phase falls from 1276.207 MiB
to 425.402 MiB and its constraint addition from 450.971 MiB to 25.569 MiB.
The strict IEEE123 `T = 3` trace is byte-for-byte identical. Apply
`ddp/patches/in_place_B_assembly.patch` last; see
`ddp/notes/FILTERDDP_IN_PLACE_B_ASSEMBLY.md`. The strict IEEE2522 `T = 12`
trace is also byte-for-byte identical; runtime falls from 737.043 s to
701.245 s (4.86%), or 63.13% below the original 1902.103-s sparse run.

**Cached network constraint Jacobian (2026-09-01):** post-seven-fix profiling
showed that first-order construction, not Hessian callbacks, still consumed
about 2.24 s per warm large10k stage. About 1.87--1.91 s was rebuilding the
sparse `cu` matrix whose pattern and nearly all values are constant. The BFM
driver now constructs one `cu` per stage and updates only four SOC derivatives
per branch. Warm callback time falls to 0.002--0.003 s. IEEE123 `T = 3` and
IEEE2522 `T = 12` retain byte-for-byte identical traces; IEEE2522 solve time
falls from 701.245 s to 577.190 s (17.69%), or 69.66% below the original
1902.103-s run. Build time rises from 1.415 s to 3.340 s. Source tracing also
confirms that full `beta`/`omega` maps feed both the backward value recursion
and every line-search rollout, so removing them requires a faithful new action
representation rather than deleting stored columns. See
`ddp/notes/FILTERDDP_CACHED_CONSTRAINT_JACOBIAN.md`.

**KKT RHS blocking rejected for speed (2026-09-01):** one captured large10k
stage (`96968 x 96968`, 1021 RHS) was solved with UMFPACK block widths 1--1021.
The existing all-at-once solve is fastest at 3.127 s; the best smaller block,
width 128, takes 3.192 s (2.08% slower). Every layout allocates about 755.34
MiB cumulatively because all columns remain necessary. Width 128 could reduce
peak temporary RHS storage from about 755 MiB to about 95 MiB if implemented
as a direct-to-policy low-memory mode, but it is not a runtime improvement and
is not enabled. See `ddp/notes/FILTERDDP_KKT_RHS_BLOCK_BENCHMARK.md`.

**Active-row `B` representation (2026-09-01):** large10k `B` is stored as
`54665 x 1020` dense but has exactly 1020 active rows (the battery-power
controls), 1.866% total density. When sparse `fu` and zero `lux`/`fux`/`cux`
prove this structure, the backward pass now forms only the dense `1020 x 1020`
active block, writes it into the existing KKT RHS, and evaluates `beta'B` from
the matching packed beta rows. The active block is exact and its value update
agrees to `9.58e-16`. Warm large10k algebra falls from about 1.2--1.4 s and
1001.743 MiB to 0.134 s and 174.839 MiB. IEEE123 and IEEE2522 traces are
byte-for-byte identical; IEEE2522 `T = 12` falls from 577.190 s to 481.518 s
(16.58%), or 74.68% below the original 1902.103-s run. Apply
`ddp/patches/active_B_rows.patch` last; see
`ddp/notes/FILTERDDP_ACTIVE_B_ROWS.md`.

**Parallel UMFPACK RHS solve rejected (2026-09-01):** on the captured large10k
stage, four independent factors reduce the 1021-column solve from 3.150 s to
1.611 s, but concurrent factor time rises from 0.494 s to 1.757 s, so total
linear-algebra time improves only 7.59%. Eight workers are slower overall and
independent factors add native memory. A shared factor leaves two- and
four-worker solves effectively serial (~3.15--3.17 s), and a repeated stress
run did not complete cleanly. Do not add a threaded UMFPACK production path;
revisit parallelism only with a solver designed for it or after reducing the
number of sensitivity columns. See
`ddp/notes/FILTERDDP_PARALLEL_UMFPACK_RHS.md`.

**Factor-backed policy actions (2026-09-01):** the rollout needs dense `beta`
and `omega` only through their products with the current state displacement.
Those products can instead be recovered exactly from one solve with the
retained stage KKT factor and compact right-hand side.  IEEE123 `T = 3` and
IEEE2522 `T = 12` retain byte-for-byte identical traces.  Large10k retained
stage policy storage falls from 757.012 MiB to 34.888 MiB (95.39%).  A matched
idle-machine IEEE2522 `T = 12` comparison measured 474.091 s factor-backed
versus 523.357 s dense (9.41% faster), with byte-identical traces; this
supersedes the misleading busy-PC comparison of 507.278 s versus an earlier
481.518-s dense run.  Prefer `FILTERDDP_FACTOR_BACKED_POLICY=1` for network
problems while retaining dense maps as a reference.  It does not remove the
backward pass's `nx+1`-column solve.  Apply
`ddp/patches/factor_backed_policy_actions.patch` after
`active_B_rows.patch`; see
`ddp/notes/FILTERDDP_FACTOR_BACKED_POLICY_ACTIONS.md`.

**Blocked value-update RHS (2026-09-01):** combined with factor-backed policy,
the backward pass can solve sensitivity columns in blocks and immediately
accumulate their exact `Vxx`/`Vx` contribution.  At width 128, large10k's RHS
workspace falls from 755.343 MiB to 94.695 MiB (87.46%), retained policy stays
34.888 MiB/stage, and first-stage RSS growth falls from 1589.223 MiB to
242.652 MiB.  The isolated bounded solve is 32.27% slower (37.179 s versus
28.109 s), but a matched full IEEE2522 `T = 12` run took 516.285 s: 8.90%
slower than factor-backed alone and 1.35% faster than dense.  Keep
`FILTERDDP_BLOCKED_VALUE_RHS=1` as the minimum-memory option.  IEEE123 `T = 3`,
IEEE2522 `T = 3`, and the matched T = 12 modes retain byte-identical traces.
Apply `ddp/patches/blocked_value_rhs.patch` after the factor-backed patch; see
`ddp/notes/FILTERDDP_BLOCKED_VALUE_RHS.md`.

**Sensitivity locality rejected as a naive shortcut (2026-09-01):** an exact
large10k stage-1 analysis grouped every `beta`/`omega` column by feeder distance
from its perturbed battery.  Battery-power feedback is local (90% radius zero),
but median 90%-energy radii are 18 edges for real branch power, 31 for voltage,
34 for current, and 53.5 for real-power-balance multipliers; worst columns
reach 82--125 edges.  Radius 8 retains only 0.18% of the median real-balance
multiplier energy.  Do not implement a fixed-radius truncation or claim the
network sensitivities are local.  Spatial decomposition remains possible only
with explicit interface/boundary information.  See
`ddp/notes/FILTERDDP_SENSITIVITY_LOCALITY.md`.

**Optimized horizon sweep (2026-09-02):** the complete validated patch stack,
with factor-backed policies enabled and blocked RHS disabled, was rerun cold
and sequentially over IEEE123C/IEEE2522C `T=3--96` and large10kC `T=3,6`.
IEEE2522C improves by 58--75% across the sweep: `T=12` falls from 1902.103 s
to 476.400 s and `T=96` from 9359.596 s to 3940.608 s. large10kC falls from
7123.782 s to 2633.856 s at `T=3` (63.03%) and from 23490.218 s to 14809.508 s
at `T=6` (36.95%), with sampled peaks of 3808 and 4216 MiB. IEEE2522
`T=12/48/96` and large10k `T=6` preserve byte-identical traces. Some IEEE123
iteration counts differ from the oldest baseline despite strict convergence;
do not claim identical trajectories for those rows. See
`ddp/results/network_filterddp/optimized_timing_comparison.csv` and
`ddp/notes/FILTERDDP_OPTIMIZED_TIMING_MATRIX.md`. large10k `T=12` remains
pending because the extrapolated run is still potentially half a day.

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
