# The minimum FilterDDP.jl API

All line references are to the unmodified clone at
`ddp/external/FilterDDP.jl` @ `513a104`. Nothing here is inferred from documentation
alone; every claim names the source that establishes it.

Throughout, three labels are used:

- **[code]** — established by reading the source.
- **[paper]** — stated/proved in P1 (arXiv 2504.08278) or P2 (arXiv 2606.01487).
- **[inferred]** — my conclusion from reading, not stated anywhere. Treat as a hypothesis.

The exported surface is exactly nine names (`src/FilterDDP.jl:25-34`):
`Objective`, `EqualityConstraints`, `Dynamics`, `ControlLimits`, `build_ocp`,
`Solver`, `Options`, `solve!`, `get_trajectory`, `get_feedback`.

---

## 1. States and controls

**[code]** States and controls are `StaticArrays.SVector`s of compile-time-fixed length,
`x::SVector{nx,T}` and `u::SVector{nu,T}` (`src/solver.jl:1-7`, `TrajectoryElement`).
`nx`, `nu`, `nc` are *type parameters* of `OCP`, `Solver` and `TrajectoryElement`, so a
solver instance is specialised to one problem size and every stage has the same
dimensions.

**[code]** The horizon is `N` stages, `t = 1..N`. There are `N` states **and `N` controls** —
`u_N` exists (`src/solver.jl:56-58` allocates `N` trajectory elements; `solve!` is passed
`u::Vector{SVector{nu,T}}` of length `N`). The dynamics only link `t` to `t+1` for
`t < N` (`src/forward_pass.jl:119-121`), so `u_N` influences the solution *only* through
the terminal objective and the stage-`N` constraint.

**[code]** In user-supplied functions `x` and `u` are plain symbolic vectors (see §6), and
`x[i]` / `u[i]` indexing is the intended style — every shipped example does this.

## 2. Discrete dynamics

```julia
f = (x, u) -> x + Δ * [x[2], u[1]]     # must return a length-nx vector
dyn = Dynamics(f, nx, nu)
```

**[code]** `src/ocp/dynamics.jl:12-40`. `f` must be a pure function of `(x, u)` returning
`nx` components. There is no time argument and no per-stage dynamics — one `f` for the
whole horizon.

## 3. Stage and terminal costs

```julia
stage_obj = Objective(l,  nx, nu)      # l(x,u)  -> scalar, used for t = 1..N-1
term_obj  = Objective(lN, nx, nu)      # lN(x,u) -> scalar, used for t = N only
```

**[code]** `src/ocp/objective.jl:10-37`; the stage/terminal split is applied in
`src/backward_pass.jl:36-50` and `src/forward_pass.jl:113-117`. Both must be scalar-valued.
**Note the terminal cost takes `u` as well as `x`** — this is genuinely useful (see §7).

**[code]** There is no per-stage cost. Time-varying cost data must be carried in the state.

## 4. Nonlinear equality constraints

```julia
constraints = EqualityConstraints((x, u) -> [u[2] - u[3] - u[1]*x[2]], nx, nu)
constraints = EqualityConstraints(nx, nu)      # nc = 0, no constraints
```

**[code]** `src/ocp/constraints.jl:12-56`. `c` returns a vector; `nc = length(c(x,u))` is
determined symbolically at construction. **The same `c` is imposed at every stage
`t = 1..N`**, including `t = N` (`src/backward_pass.jl:52-63`, loop `for t = ocp.N:-1:1`).

Two hard structural requirements:

- **`nc ≤ nu`.** **[paper]** stated in P1 §II-A. **[code]** enforced implicitly and
  brutally: `src/backward_pass.jl:119` builds `@SMatrix zeros(T, nu, nu-nc)`, so `nc > nu`
  is a construction error on a negative dimension.
- **`∇_u c` must have full row rank `nc` at every stage and every iterate.**
  **[code]** `src/backward_pass.jl:118-127` forms `A = cu'`, takes a QR of `[A z]`, then
  `fk = lu(AY)` with `AY = A'*Y` (`nc × nc`). **There is no rank or conditioning check on
  `AY`** — `lu` is called without `check=false` and without inspecting the pivots, and
  `Options.linsolve_tol` (`src/options.jl:39`) is declared but never referenced anywhere in
  `src/`. **[inferred]** A rank-deficient `∇_u c` therefore does not produce a clean error;
  it produces `Inf`/`NaN` in the update rule, which then surfaces as a line-search failure
  or a nonsense answer. This is the sharpest practical trap in the library.

**[code]** A constraint that depends on `x` only (`∇_u c = 0`) violates the second
requirement identically. **Pure state constraints, including pure terminal state
constraints, are not representable.**

## 5. Inequality constraints and bounds

**[code]** The only native inequality mechanism is **bounds on controls**:

```julia
cl = ControlLimits(SVector{nu,T}([-10.0, 0.0, 0.0]),
                   SVector{nu,T}([ 10.0, Inf, Inf]))
```

`src/ocp/control_limits.jl:10-18`. `-Inf` / `Inf` (or `±floatmax`) mark a component as
unbounded on that side and it is masked out of the barrier entirely (`maskl`, `masku`).
Passing `±Inf` on every component gives `nl = nu = 0`, which disables the interior-point
machinery completely and reduces the algorithm to pure equality-constrained FilterDDP
(`src/solve.jl:19,31` gate on `ni > 0`).

**[code]** General nonlinear inequalities `b^L ≤ g(x,u) ≤ b^U` are handled by the user, not
the library, via the slack reformulation documented in README line 31: append a slack `s`
to the control vector, add the equality `g(x,u) - s = 0`, and bound `s`. Every extra
inequality therefore costs one control **and** one equality — and since `nc ≤ nu` must
hold, **the slack it adds is what pays for the constraint it adds.** `experiments/filterddp/concar.jl:89-105`
shows the shipped idiom for obstacle avoidance.

**[code]** State bounds are not native either, but they *are* representable the same way:
`e - e_min - s = 0` with `s ∈ [0, e_max - e_min]` turns a two-sided state bound into one
extra control and one extra equality, and the resulting `∇_u c` row is `[0 … 0 -1]`, which
keeps full row rank. Confirmed working in Stage 6.

**[code]** Bounds are enforced by a primal-dual interior point method with a
fraction-to-the-boundary rule (`src/forward_pass.jl:87-108`) and IPOPT's barrier update
(`src/solve.jl:31-35`). **[paper]** P1 §VI, P2 §5, following Pavlov et al. (2021). The
papers only ever write the one-sided `u ≥ 0`; the two-sided form is a code-only
generalisation.

## 6. How derivatives are computed

**[code]** Symbolically, at construction time, by `Symbolics.jl`, then compiled to native
code with `RuntimeGeneratedFunctions`. `Dynamics`, `Objective` and `EqualityConstraints`
each build first *and second* derivatives eagerly
(`src/ocp/dynamics.jl:19-35`, `objective.jl:15-32`, `constraints.jl:17-39`).

Consequences that matter in practice:

- **[code]** User functions are traced with `Symbolics.variables`, so they must be
  symbolically differentiable: **no branching on numeric values of `x`/`u`, no `if`,
  no `abs`, no `min`/`max`, no iteration counts that depend on values.** `Symbolics.simplify`
  is called on every expression and on every derivative.
- **[code]** Second derivatives are contracted against an adjoint (`λ` for dynamics,
  `ϕ` for constraints), so what is generated is a Hessian-vector product, not a 3-tensor
  (`dynamics.jl:27-35`, `constraints.jl:29-39`).
- **[inferred]** Cost scales with the *symbolic* size of the expressions. Model
  construction, not solving, is what will dominate for a large MPOPF: `Symbolics.simplify`
  on a Hessian of a several-hundred-variable expression is the bottleneck. This is a
  build-time property; it says nothing about solve speed.
- **[paper]** P1 Remark 3 and README line 7: exact Hessians are required; Gauss-Newton /
  quasi-Newton "does not perform well". **[code]** `Options.quasi_newton`
  (`src/options.jl:2`) is **never read anywhere** — setting it does nothing.

## 7. Terminal equality constraints — the honest answer

**[code]** **Not supported as constraints.** `c` is shared by all stages (§4), so a
terminal-only equality cannot be written; and a constraint made to vanish at
non-terminal stages would have `∇_u c = 0` there, breaking the rank requirement.

**[code]** But there is a clean workaround for the common case, and it is worth stating
precisely because it is easy to miss: **the terminal objective receives `u_N`**, and the
dynamics are known in closed form, so any quantity of the form `f(x_N, u_N)` — i.e. the
state *one step past the horizon* — is directly expressible in the terminal cost. In the
Stage 5 battery model, `e_{T+1} = e_N - η·pb_N` is written exactly in `lN(x,u)`, so a
terminal-energy *penalty* is exact rather than approximate. **[inferred]** Turning that
penalty into a hard equality would require either a per-stage constraint function (which
the API does not have) or an outer loop the library does not provide.

## 8. Initialisation

```julia
x1 = SVector{nx,T}([...])
ū  = [SVector{nu,T}([...]) for k = 1:N]     # length N, not N-1
solve!(solver, x1, ū)
```

**[code]** `src/solver.jl:71-106`. Only the **controls** are initialised by the user; the
state trajectory is generated by rolling the dynamics forward from `x1`
(`src/solver.jl:102-104`). **There is no way to supply an initial state trajectory** —
the method is single-shooting, not a collocation method with an infeasible start.

**[code]** Supplied controls are silently **projected into the interior** of the bounds by
the `κ_1`/`κ_2` rule (`src/solver.jl:80-93`), exactly as IPOPT does. A `ū` on or outside a
bound will not be the `ū` the solver starts from.

**[code]** Constraint multipliers start at `ϕ = 0`; bound multipliers start at exactly
`1.0` for each bounded side (`src/solver.jl:95-98`), ignoring `Options.ineq_dual_init`.

## 9. Solver options and tolerances

```julia
options = Options{Float64}(verbose=true, optimality_tolerance=1e-8, max_iterations=1000)
solver  = Solver(ocp; options=options)
```

**[code]** `src/options.jl`. `Options` is a mutable kwdef struct; fields can also be set
after construction (`solver.options.verbose = false`, as the examples do). The ones that
actually matter: `optimality_tolerance` (default `1e-8`), `max_iterations` (`1000`),
`verbose`, `print_frequency`, `μ_init` (`1.0`), `reg_1`/`reg_min`/`reg_max` and the
`κ_w_*` regularisation ladder, and the filter constants `η_L, s_L, δ, s_θ, γ_α, γ_θ, γ_L`.

**[code]** Declared but **never referenced in `src/`**: `quasi_newton`, `reset_cache`,
`ineq_dual_init`, `κ_c`, `δ_c`, `κ_Σ`, `linsolve_tol`. Verified by grep. Do not expect
them to do anything.

## 10. Calling the solve and reading the result

```julia
status = solve!(solver, x1, ū)          # or solve!(solver) to continue
x_sol, u_sol = get_trajectory(solver)   # Vector{SVector} of length N each
fb = get_feedback(solver, t)            # x -> ū_t + α_t + β_t (x - x̄_t)
```

**[code]** `src/solve.jl:1-6`, `src/solver.jl:64-69`, `:116-124`.
`get_trajectory` returns the **nominal** (accepted) trajectory.

### Convergence and failure reporting

**[code]** `solve!` returns `solver.data.status::Int`, also readable as `solver.data.status`.
Meaningful exit values (`src/solver_data.jl:5-7`, `src/print.jl:30-45`):

| status | meaning |
| --- | --- |
| `0` | optimal solution found |
| `1` | backward pass failed — no iteration matrix with correct inertia |
| `7` | line search failed to find an acceptable iterate |
| `8` | maximum iterations reached |

Values `2`–`6` are transient line-search states and are reset each iteration; if one ever
reached `on_exit` it would print `"DEBUG: This message should not display."`
(`src/print.jl:42`).

**[code]** `status == 0` is *not* merely "did not throw". The loop only breaks with status
`0` at `src/solve.jl:30`, i.e. when
`max(dual_inf, cs_inf_0, primal_inf) < optimality_tolerance`, evaluated on the current
nominal iterate in the backward pass. So `status == 0` does mean a genuine scaled-KKT test
was passed. **[inferred]** The one caveat is that these are the *scaled* residuals —
`dual_inf` and `cs_inf` are divided by IPOPT-style scaling factors
(`src/backward_pass.jl:198-201`) that are `≥ 1`, so the unscaled residuals can be larger
when the multipliers are large. Always re-audit residuals from `get_trajectory` output
rather than trusting `data.primal_inf` alone.

Other useful post-solve fields on `solver.data`: `k` (iterations), `objective`,
`primal_inf`, `dual_inf`, `cs_inf_0`, `step_size` (last accepted `γ`), `l` (line-search
backtracks on the last iteration), `reg_last`, `μ`, `j` (barrier subproblem index).

**[code]** Per-iteration history is **not stored**. `verbose=true` prints
`iter, objective, pr_inf, du_inf, cs_inf, lg(μ), lg(reg), alpha, ls`
(`src/print.jl:12-28`) and that printout is the only record of step sizes and
regularisation over time.

## 11. Dynamics feasibility during the forward pass

**[code]** **Yes, maintained exactly.** `rollout!` (`src/forward_pass.jl:44-123`) starts at
`x = solver.nominal[1].x` and sets `x = f(x, u)` at each step (`:120`), so every trial
trajectory — accepted or rejected — satisfies the dynamics to machine precision by
construction. Only the *equality constraints* `c` are violated during iteration; `θ(w)` in
the filter is `Σ_t ‖c_t‖_1` (`:83`) and does **not** include a dynamics residual, because
there is never one. This is the classic DDP property and it is genuinely present here.

## 12. Assumptions the implementation makes

| Topic | What is assumed | Basis |
| --- | --- | --- |
| Dimensions | `nc ≤ nu`, same `nx, nu, nc` at every stage | **[paper]** P1 §II-A; **[code]** `backward_pass.jl:119` |
| Smoothness | `ℓ, f, c` twice continuously differentiable **and symbolically traceable** | **[paper]** P1 §II-A; **[code]** `Symbolics` tracing |
| Hessians | exact second derivatives; no quasi-Newton fallback exists | **[paper]** P1 Remark 3; **[code]** `quasi_newton` unused |
| LICQ | `σ_min(A_t)` uniformly bounded away from 0, i.e. `∇_u c` full row rank at every stage — **a per-stage condition on the control Jacobian alone, strictly stronger than LICQ for the full NLP** | **[paper]** P1 Assumption (1B), P2 (G4)/(B4); **[code]** unchecked `lu(AY)` |
| Second-order sufficiency | `Zᵀ H Z ≻ 0` on the null space of `A_tᵀ` at every stage | **[paper]** P1 (1C), P2 (G3)/(B3); **[code]** `cholesky(...; check=false)` at `backward_pass.jl:129-134` |
| Inertia correction | failure of that Cholesky triggers a primal regularisation ladder `reg_1 → κ_w_p·reg …` up to `reg_max = 1e40`, restarting the whole backward pass | **[code]** `backward_pass.jl:136-143` |
| Initialisation | controls only; states always dynamically feasible; controls projected inside bounds | **[code]** `solver.jl:71-106` |
| Strict complementarity | not checked or required anywhere in the code | **[code]** |
| Feasibility restoration | **assumed to exist by both convergence theorems; not implemented** | **[paper]** P2 Assumptions G/B; **[code]** absent — `forward_pass.jl:40` just fails with status 7 |
| Global optimality | never claimed; "global convergence" = convergence to a first-order point from an arbitrary start | **[paper]** P2 Theorem 1 |

## 13. Minimal working template

```julia
using FilterDDP, StaticArrays

const nx, nu, N = 2, 2, 3
T = Float64

dyn        = Dynamics((x, u) -> [x[1] - u[2], x[2] + 1.0], nx, nu)
stage_obj  = Objective((x, u) -> u[1]^2, nx, nu)
term_obj   = Objective((x, u) -> u[1]^2 + 5.0 * (x[1] - u[2] - 2.0)^2, nx, nu)
constraints = EqualityConstraints((x, u) -> [u[1] + u[2] - 1.0], nx, nu)
cl         = ControlLimits(SVector{nu,T}(fill(-Inf, nu)), SVector{nu,T}(fill(Inf, nu)))

ocp    = build_ocp(N, stage_obj, term_obj, dyn, constraints, cl)
solver = Solver(ocp; options = Options{T}(verbose = true, optimality_tolerance = 1e-10))

x1 = SVector{nx,T}([2.0, 1.0])
ū  = [SVector{nu,T}([0.0, 0.0]) for _ = 1:N]
solve!(solver, x1, ū)

solver.data.status == 0 || error("did not converge: status $(solver.data.status)")
x_sol, u_sol = get_trajectory(solver)
```
