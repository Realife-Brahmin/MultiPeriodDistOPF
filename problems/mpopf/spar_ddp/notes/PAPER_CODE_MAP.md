# Paper ↔ Code Map — FilterDDP.jl

Repository: <https://github.com/mingu6/FilterDDP.jl> @ `513a104ccd67a881c95460736a4a3486d719804f`
(local clone: `ddp/external/FilterDDP.jl`, unmodified)

Papers (downloaded from arXiv on 2026-08-03, they were not present locally):

| Short name | File | arXiv | Version string in PDF |
| --- | --- | --- | --- |
| **P1** (main) | `ddp/papers/Xu_2026_FilterDDP.pdf` | 2504.08278 | v7 [math.OC] 31 May 2026 — "accepted to ICRA 2026" |
| **P2** (global conv.) | `ddp/papers/Xu_2026_FilterDDP_Global_Convergence.pdf` | 2606.01487 | A PREPRINT — June 2, 2026 |

Plain-text extractions (`pdftotext -layout`) sit alongside as `.txt`; line references below are to
those `.txt` files.

---

## 1. Problem class

P1 §II-A eq. (1):

```
min_{x,u}  J = Σ_{t=1..N} ℓ(x_t, u_t)
s.t.       x_1 = x̂_1
           x_{t+1} = f(x_t, u_t),  t ∈ [N-1]
           c(x_t, u_t) = 0,        t ∈ [N]
```

with the explicit standing requirement **`c : R^nx × R^nu → R^nc` where `nc ≤ nu`**
(P1 txt line 88-89), and ℓ, f, c twice continuously differentiable.

P1 §VI eq. (37) and P2 §5 eq. (49) add `u_t ≥ 0`. The **code generalises this** to two-sided
`b^L ≤ u_t ≤ b^U` (README lines 11-21, `src/ocp/control_limits.jl`). Both papers only ever write
the one-sided `u ≥ 0` form.

Code counterpart: `OCP` struct and `build_ocp` in `ddp/external/FilterDDP.jl/src/ocp/ocp.jl:4-27`.

## 2. Component-by-component mapping

| Paper object | Paper location | Source file / lines |
| --- | --- | --- |
| ℓ (stage cost) | P1 (1) | `src/ocp/objective.jl:10-37` → `Objective(l, nx, nu)` |
| ℓ^F (terminal cost) | *not in P1 (1)* — code-only split | `src/backward_pass.jl:36-50`, `src/forward_pass.jl:113-117` |
| f (discrete dynamics) | P1 (1) | `src/ocp/dynamics.jl:12-40` → `Dynamics(f, nx, nu)` |
| c (stagewise equality) | P1 (1) | `src/ocp/constraints.jl:12-52` → `EqualityConstraints(c, nx, nu)` |
| b^L, b^U (control limits) | P1 (37) as `u ≥ 0` | `src/ocp/control_limits.jl:10-18` |
| Lagrangian `L(x,u,ϕ) = ℓ + ϕᵀc` | P1 (2) | `src/backward_pass.jl:104` (`dot(c, ϕ)`) |
| KKT conditions (3a)-(3d) | P1 §II-B | checked via `data.dual_inf`, `data.primal_inf` |
| Optimality error `E(w,λ)` (7) | P1 §IV-A | `src/solve.jl:27-28` (`opt_err_0`, `opt_err_μ`) |
| Termination `E < ϵ_tol` | P1 Alg. 1 step 3 | `src/solve.jl:30` |
| Boundary conditions (8) | P1 §IV-B | `src/backward_pass.jl:16-18` (`V̂x`, `V̂xx`, `λ` ← 0) |
| Newton step (10), Q̄ blocks (11) | P1 §IV-B a) | superseded by (12)/(13) in code |
| **Perturbed Newton step (12)**, matrix `K_t` | P1 §IV-B a) | `src/backward_pass.jl:112-160` |
| `H_t`, `B_t`, `C_t` (13a-c) | P1 §IV-B | `Ĥ` `src/backward_pass.jl:98,107,114`; `B` `:101,108`; `C` `:92,106` |
| `λ̄_{t+1} · f''` tensor contraction | P1 (13) | `src/backward_pass.jl:67-69` (passes `λ` when `nc>0`) |
| `V̄_x^t`, `λ̄_t` update (14a) | P1 §IV-B b) | `src/backward_pass.jl:186-187, 190-191` |
| `P_t` update (14b) | P1 §IV-B b) | `src/backward_pass.jl:185, 189` (`V̂xx`) |
| Forward-pass trial point (15) | P1 §IV-C | `src/forward_pass.jl:72-73`, rollout `:44-123` |
| θ(w), L(w) filter criteria (16) | P1 §IV-C b) | `data.primal_1_*`, `data.barrier_lagrangian_*` |
| Sufficient decrease (17a)/(17b) | P1 §IV-C b) | `src/forward_pass.jl:34` |
| Switching condition (18) | P1 §IV-C b) | `src/forward_pass.jl:28-29` |
| `m` (19), expected decrease | P1 §IV-C b) | `data.expected_change_L`, `src/backward_pass.jl:195-196` |
| Armijo condition (20) | P1 §IV-C b) | `src/forward_pass.jl:30` |
| Filter augmentation (21) | P1 §IV-C b) | `src/solve.jl:57-61` (`update_filter!`) |
| Filter init `θ ≥ θ_max` | P1 Alg. 1 step 1 | `src/solve.jl:63-67` (`reset_filter!`), `max_primal_1 = 1e6` |
| Sub-optimal value function `V^t` (22), Thm 1 | P1 §IV-E | theory only — no direct code object |
| Barrier cost φ_μ (38) | P1 §VI | `src/backward_pass.jl:75`, `src/forward_pass.jl:111-112` |
| Perturbed error `E_μ` (39) | P1 §VI | `data.cs_inf_μ`, `src/backward_pass.jl:81,83` |
| Barrier update rule (40) | P1 §VI | `src/solve.jl:31-35` |
| Primal-dual backward pass (41) | P1 §VI-A a) | `Σ_L`,`Σ_U` `src/backward_pass.jl:96-98`; `Q̂_u` `:89` |
| Dual update `z⁺` (42) | P1 §VI-A b) | `src/forward_pass.jl:74-75`; `χ`,`ζ` from `src/backward_pass.jl:164-168` |
| Fraction-to-boundary (43) | P1 §VI-A b) | `src/forward_pass.jl:93-108`; `τ` at `:9` |
| Barrier convergence loop (Alg. 1 step 3.1) | P1 §VI-A c) | `src/solve.jl:31-36` |
| Assumption (1B) `σ_min(A_t)` bounded away from 0 | P1 §V | implicitly required by `lu(AY)`, `src/backward_pass.jl:127` |
| Assumption (1C) `ξᵀH_tξ > 0` on null(A_t) | P1 §V | `cholesky(Symmetric(Z'ĤZ); check=false)`, `src/backward_pass.jl:129-134` |
| Assumptions G1-G5, Theorem 1 (global) | P2 §4.1, §4 | see §4 below |
| Assumptions B1-B8, Theorem 2 (barrier) | P2 §5 | see §4 below |
| Nullspace decomposition `ζ = q + p`, (23)-(26) | P2 Remark 6 | **this is what the code actually implements**, `src/backward_pass.jl:117-153` |

## 3. Which examples correspond to which paper experiments

P1 §VII evaluates three task classes, 100 OCPs each. The repo ships one script per task under
`experiments/filterddp/`, plus IPOPT (`experiments/ipopt/`) and ProxDDP (`experiments/proxddp/`,
Python) baselines, and stores the authors' own outputs in `experiments/*/results/*.txt`.

| Paper experiment | FilterDDP script | states / controls (paper) | In paper? |
| --- | --- | --- | --- |
| §VII-C Cartpole swing-up w/ Coulomb friction | `experiments/filterddp/cartpole_friction.jl` | nx=4, nu=15 | yes |
| §VII-D Acrobot swing-up w/ joint limits | `experiments/filterddp/acrobot_contact.jl` | nx=4, nu=7 | yes |
| §VII-E Non-prehensile pushing | `experiments/filterddp/manip_pushing.jl` | nx=4, nu=6 | yes |
| — | `experiments/filterddp/double_integrator.jl` | nx=2, nu=3 | **no** — this is the README Quick Start ("block push", SIAM tutorial §8.1) |
| — | `experiments/filterddp/concar.jl`, `concar_quad.jl` | car-like | no |
| — | `experiments/filterddp/acrobot.jl` | unconstrained acrobot | no |

**Smallest / fastest self-contained example: `double_integrator.jl`** (N=101, nx=2, nu=3, nc=1,
one OCP, no external model files). It is the same problem as the README Quick Start and it does
exercise every part of the algorithm we care about: a nonlinear equality constraint
(`u₂ - u₃ - u₁x₂ = 0`), control limits (so the interior-point/barrier path is live), the
nullspace backward pass, the filter line search and the convergence test.

`acrobot.jl` and `cartpole_friction.jl` depend on `experiments/models/*.jl`;
`manip_pushing.jl` and `acrobot_contact.jl` are the large contact-implicit problems.
Note the caveat in README line 27: for cartpole the *paper* used the sister repo
`FilterDDPLAPACK.jl`, not this one — so `cartpole_friction.jl` here is not literally the code
that produced the paper's cartpole numbers.

## 4. Divergences between paper and code

These are things the code does **differently from**, or **not at all compared to**, what the
papers describe. Each is a direct source reading, not an inference about intent.

1. **No feasibility restoration phase.** Both P2 `Assumptions G` (txt 425-426) and
   `Assumptions B` (txt 787-789) open with *"we assume that the feasibility restoration phase in
   step 9 always terminates successfully"*, and (G5)/(B5) refer to iterations where restoration is
   invoked. **The code has no restoration phase at all.** When the line search runs out of step
   size, `src/forward_pass.jl:40` sets `status = 7` and `src/solve.jl:41` breaks out. So the
   hypothesis of P2 Theorem 1 is not met by this implementation on any run where the line search
   fails. This is the single most important gap for judging whether the global-convergence result
   covers a given problem.

2. **`δ_c` is dead.** P1 (12) and (41) put `-δ_c I` in the (2,2) block of `K_t` and P1 §IV-B says
   `δ_w, δ_c ≥ 0` are chosen "following Alg. IC in [27] so that `K_t` has an inertia of
   `(nu, nc, 0)`". In the code `δ_c = 0.` is assigned once at `src/backward_pass.jl:6` and then
   **never read**; `Options.δ_c = 1e-8` (`src/options.jl:21`) is **never referenced anywhere**.
   Only the primal regularisation `δ_w` (code: `reg`, `src/backward_pass.jl:112-114`) is live.

3. **Different linear algebra than the paper describes.** P1 §IV-B says the inertia of `K_t` is
   recovered with an `LDLᵀ` factorization (refs [32] Bunch-Kaufman, [33] rook pivoting). The code
   instead uses the **nullspace method**: QR of `[A z]` (`src/backward_pass.jl:120`), an `lu` of
   `AᵀY` for the range-space part (`:127`), and a `cholesky(...; check=false)` of the reduced
   Hessian `ZᵀĤZ` (`:129-130`) whose failure (`ck.info != 0`) triggers the regularisation bump.
   This is the scheme of P2 Remark 6/(23)-(26), and it is *equivalent in effect* — Cholesky success
   on `ZᵀĤZ` plus nonsingular `AᵀY` gives the right inertia — but it means failure is detected as
   "reduced Hessian not positive definite", and **a rank-deficient `c_u` shows up as a silently
   singular `lu(AY)`, not as an inertia failure**. There is no rank check on `AY`.

4. **`θ_min` is effectively 0.** P1 (17)/(18) and IPOPT use a threshold `θ_min` below which the
   Armijo/switching branch is taken. `solver_data(T)` sets `min_primal_1 = 1e-5`
   (`src/solver_data.jl:36`), but `reset!` — called at the top of every `solve!`
   (`src/solve.jl:15`) — overwrites it with `0.0` (`src/solver_data.jl:67`), and nothing ever
   raises it again. Since `θ ≥ 0` always, the test `θ <= data.min_primal_1`
   (`src/forward_pass.jl:31`) can only pass when `θ` is *exactly* `0.0`. For an equality-
   constrained problem the L-type/Armijo branch is therefore essentially unreachable and
   acceptance is decided by the filter branch at `src/forward_pass.jl:34`. (For `nc == 0`,
   `θ ≡ 0.0` and the Armijo branch is always taken, which is correct.)

5. **`s_L` default is outside what P2 Theorem 2 covers.** P2 §5 (txt 766-768) states *"it is
   necessary to assume that `s_L` is chosen to be 1 to establish this result"* (Theorem 2, iterates
   bounded away from their bounds — the result that in turn lets Theorem 1 apply to the barrier
   algorithm). `Options.s_L = 2.3` (`src/options.jl:30`), which is IPOPT's default `s_φ`. So with
   default options, the *interior-point* branch of the code is not the one P2 Theorem 2 proves
   things about. P1 §VII-B does say parameters match IPOPT's except `θ_μ` and `μ_init`.

6. **Split stage/terminal objective.** P1 (1) has a single `ℓ` summed over `t ∈ [N]`; the code
   requires two `Objective`s and applies `term_objective` at `t = N` only
   (`src/backward_pass.jl:36-50`). Cosmetic, but it changes how you write a problem.

7. **`ℓ^F` in the README is mis-indexed.** README line 14 writes
   `Σ_{t=1}^{N-1} ℓ(x_t,u_t) + ℓ^F(x_t,u_t)` — the terminal term should be `ℓ^F(x_N, u_N)`.
   The code does the right thing.

8. **Options declared but never used anywhere in `src/`:** `quasi_newton`, `reset_cache`,
   `ineq_dual_init`, `κ_c`, `δ_c`, `κ_Σ`, `linsolve_tol`. Verified by grep over `src/`.
   In particular `quasi_newton = true` does **nothing** — consistent with README line 7 warning
   that FilterDDP needs exact Hessians, but the option is misleading.

9. **Dual variables `z` are initialised to exactly `1.0`**, hard-coded
   (`src/solver.jl:96-98`), ignoring `Options.ineq_dual_init`.

## 5. What each paper actually proves (for later use in Stage 7)

- **P1 Theorem 2**: `w̄` feasible is a fixed point of FilterDDP ⟺ it satisfies the KKT conditions (3).
- **P1 Theorem 3 / §V**: *local quadratic* convergence in a neighbourhood of `w⋆`, under
  Assumptions 1 (1A smoothness+Lipschitz derivatives, 1B `σ_min(A_t)` uniformly bounded away from
  zero, 1C `Hₜ` positive definite on `null(A_t)`), and **explicitly assuming every full step
  `γ = 1` is accepted by the filter**.
- **P2 Theorem 1**: under Assumptions G, `lim θ(w_k) = 0` and `lim inf χ(w_k) = 0` — all limit
  points feasible, and if the iterates are bounded there is a first-order optimal limit point.
  Requires the restoration phase (not implemented — see §4.1).
- **P2 Theorem 2**: under Assumptions B (which include B7, an LICQ-type condition on
  `∇_u c` *together with* the active bound gradients `e_i`), barrier iterates stay bounded away
  from their bounds, given `s_L = 1`.

Neither paper proves convergence to a *global* minimum; "global convergence" here means
convergence to a first-order point from an arbitrary starting point, in the standard NLP sense.
