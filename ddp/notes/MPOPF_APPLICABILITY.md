# Is FilterDDP applicable to distribution-system MPOPF?

Short answer: **for copper-plate / LinDistFlow-scale problems, yes with caveats; for a
BFM-SOCP or nonconvex-BFM MPOPF on a real feeder, no — not because the mathematics
forbids it, but because the implementation's per-stage linear algebra is dense in the
number of controls, which is exactly where a distribution MPOPF is large.**

Everything below is grounded in the Stage 1-6 evidence in this directory. Where I am
extrapolating rather than reporting, I say so.

---

## 1. What maps onto what

In a branch-flow-model MPOPF over periods `t = 1..T`, the variables split cleanly:

| MPOPF quantity | DDP role | Why |
| --- | --- | --- |
| Battery state of charge `E_{i,t}` | **state** `x_t` | the only genuinely inter-temporal variable |
| Thermal/hydro commitment or ramp memory, tap-changer position | **state** | if modelled with inter-temporal limits |
| Real/reactive injections `p_{i,t}, q_{i,t}` | **control** `u_t` | free within a period |
| Branch flows `P_{ij,t}, Q_{ij,t}`, squared current `ℓ_{ij,t}`, squared voltage `v_{i,t}` | **control** `u_t` | algebraic, per-period |
| Slack variables for every inequality | **control** `u_t` | forced by the API (§5) |
| Nodal power balance, Ohm/voltage-drop, `ℓ` definition | **`c(x_t,u_t) = 0`** | per-period equalities |
| SOC recursion `E_{t+1} = E_t + η·(...)` | **`f(x_t,u_t)`** | the dynamics |
| Voltage/current/injection limits | control bounds + slacks | §5 |

**The state dimension is genuinely tiny** — `nx` = number of storage devices, not number of
buses. That is the one structurally attractive fact about DDP here: the value function
`V_t` lives on SOC space only, and the backward pass propagates an `nx × nx` Hessian.
For a feeder with 50 batteries, `nx = 50` regardless of whether the feeder has 123 or
10,000 buses.

**But `nu` is the whole per-period network.** Every algebraic variable has to be a control,
because there is nowhere else to put it.

## 2. Should the algebraic network variables be controls, eliminated, or constraints?

Three options, and the API decides for you:

1. **As controls, with the network equations as `c(x_t,u_t)=0`.** This is the only
   formulation that fits FilterDDP as written. It is also what makes `nu` explode.
2. **Eliminated through an embedded per-period OPF solve.** This would keep `nu` small
   (just the storage set-points), but it requires `f` and `ℓ` to be *implicit* functions
   defined by an inner solve — and FilterDDP builds all derivatives symbolically with
   `Symbolics.jl` at construction time (`src/ocp/*.jl`, see FILTERDDP_API.md §6). **You
   cannot hand it a function that calls a solver.** This option is closed without
   replacing the entire derivative layer.
3. **As pure algebraic constraints with no associated control.** Closed: a constraint row
   with `∇_u c = 0` makes `lu(AY)` singular in the backward pass
   (`src/backward_pass.jl:124-127`), and there is no rank check to catch it.

So option 1 is forced, and with it the `nc ≤ nu` requirement (P1 §II-A;
`src/backward_pass.jl:119`). That requirement is *satisfiable* — declaring every network
variable a control gives you roughly as many controls as equations — but it is a real
modelling constraint: **every equality you add must be paid for with a control.**

## 3. The decisive obstacle: dense `nu × nu` algebra per stage

The backward pass, for each stage, forms and factorises

- `Ĥ` — an `nu × nu` dense matrix (`src/backward_pass.jl:98`),
- a QR of an `nu × nu` matrix (`:120`),
- an `lu` of an `nc × nc` matrix (`:127`),
- a Cholesky of the `(nu-nc) × (nu-nc)` reduced Hessian (`:129-130`),

all in `StaticArrays`, i.e. **stack-allocated, compile-time-sized, dense**. There is no
sparsity anywhere.

The authors' own guidance calibrates the limit: README line 27 recommends switching to the
sister repo `FilterDDPLAPACK.jl` for "larger size OCPs … such as the contact dynamics
cartpole experiment" — and that experiment has **`nu = 15`** (P1 §VII-C). Their statically
sized implementation is already at its useful limit at fifteen controls.

A single period of a BFM MPOPF on the IEEE 123-bus feeder has on the order of several
hundred algebraic variables; on the 10k-bus system it is tens of thousands. Dense
`O(nu³)` per stage per iteration at `nu ~ 10⁴` is ~10¹² flops *per stage per iteration*,
and `SMatrix{10000,10000}` is not a thing that compiles. **[inferred, but not marginally]**
This is not a tuning problem; it is a mismatch between the algorithm's implementation and
the problem's structure.

The deeper point: **DDP buys you decomposition in time. A distribution MPOPF whose only
inter-temporal coupling is SOC is already nearly separable in time — that is precisely why
temporal ADMM works on it. DDP would trade one temporal decomposition for another while
making each per-period subproblem strictly harder**, replacing a sparse per-period OPF
with a dense per-period KKT factorisation.

## 4. Which assumptions hold, per problem class

| | copper-plate battery | LinDistFlow MPOPF | BFM-SOCP MPOPF | nonconvex BFM MPOPF |
| --- | --- | --- | --- | --- |
| Smoothness (C², symbolically traceable) | ✅ verified Stages 5-6 | ✅ linear | ✅ smooth | ⚠️ smooth only away from `v = 0` |
| `nc ≤ nu` | ✅ verified | ✅ by construction | ✅ by construction | ✅ by construction |
| `∇_u c` full row rank each stage (LICQ-like, **G4/B4**) | ✅ verified | ✅ plausible | ⚠️ see §5 | ⚠️ fails near voltage collapse |
| `ZᵀHZ ≻ 0` each stage (**G3/B3**) | ✅ verified, `reg` never fired | ✅ convex | ⚠️ convex problem, but the *stagewise* reduced Hessian is not automatically PD | ❌ not guaranteed |
| Strict complementarity | ✅ (Stage 6 active sets clean) | ⚠️ | ❌ SOC constraints are generically tight on **every** branch | ❌ |
| Feasible initialisation | ✅ trivial | ✅ | ⚠️ | ⚠️ |
| Bounded iterates (**B6**) | ✅ | ✅ | ⚠️ | ❌ |
| Tractable `nu` | ✅ `nu = 2-3` | ❌ | ❌ | ❌ |

Two entries deserve emphasis because they are counter-intuitive:

- **Convexity of the problem does not give you `ZᵀHZ ≻ 0` stagewise.** FilterDDP's `H_t`
  is not the Hessian of the full NLP; it is `L̄ᵘᵘ + f̄ᵘᵀP_{t+1}f̄ᵘ + λ̄_{t+1}·f̄ᵘᵘ` (P1 (13a)),
  which involves the *propagated* value-function Hessian `P_{t+1}` and the dual-weighted
  dynamics Hessian. A globally convex MPOPF can still produce an indefinite `H_t` at an
  intermediate iterate, which triggers the regularisation ladder
  (`src/backward_pass.jl:136-143`). This is not hypothetical — it is why the ladder exists.
- **Strict complementarity is the one I would bet against for BFM-SOCP.** In an exact SOC
  relaxation the constraint `P² + Q² ≤ v·ℓ` is *tight on every branch at the optimum*.
  That is a large set of simultaneously active inequalities, all handled by the interior
  point path, with the multipliers of a degenerate active set to resolve. Note this is a
  *conjecture*, not something Stage 6 established: with a well-conditioned objective the
  solver handled every active set I gave it cleanly down to `1e-12`. What Stage 6 does show
  is that conditioning of the reduced Hessian is what governs whether it certifies at all
  (see §7), and a fully-active conic set is exactly where that conditioning is worst.

## 5. Conic constraints `P² + Q² ≤ v·ℓ`

**Representable, but only through the slack reformulation, and it costs you.** There is no
conic or inequality primitive; the only mechanism is README line 31 / `concar.jl:89-105`:

```
v·ℓ - P² - Q² - s = 0,     s ≥ 0
```

I verified this pattern works — Stage 6c used exactly it to turn a *state* bound
(`e_min ≤ e ≤ e_max`) into one slack control plus one equality, and it reproduced the
reference active set exactly. The `∇_u` row of the slack equation is `[0 … 0 -1]`, which
preserves the rank condition. So the mechanism is sound.

The cost is arithmetic: **one slack control and one equality per conic constraint, i.e.
per branch per period.** For `n_branch` branches this adds `n_branch` to both `nu` and
`nc` in every period — compounding the §3 problem rather than relieving it. And because
these constraints are generically active, `s → 0` and the barrier must resolve a large
active set right at the boundary.

## 6. Does the global convergence theorem cover this? No.

P2 Theorem 1 concludes `lim θ(w_k) = 0` and `lim inf χ(w_k) = 0` **under Assumptions G**.
Paraphrasing what those require, and checking each against the shipped code:

- **Preamble to Assumptions G** (P2 §4.1): the analysis assumes *"the feasibility
  restoration phase in step 9 always terminates successfully"*. **The implementation has no
  feasibility restoration phase at all.** When the line search exhausts its step size it
  sets status 7 and stops (`src/forward_pass.jl:40`, `src/solve.jl:41`). Assumptions (G5)
  and (B5) are about *when restoration is invoked* and are vacuous here. **This alone puts
  the shipped solver outside the theorem's hypotheses on any run where the line search
  fails** — which Stage 6d and 6g both did.
- **(G1)** ℓ, f, c and their derivatives bounded and Lipschitz on an open set containing all
  iterates. For BFM with `ℓ = (P²+Q²)/v`, derivatives blow up as `v → 0`; boundedness is an
  assumption about the iterate path, not something you can check in advance.
- **(G3)** reduced Hessian uniformly positive definite — see §4.
- **(G4)** `σ_min(A_t) ≥ M_A > 0` uniformly — the per-stage control-Jacobian rank condition.
- **(B6)** iterates uniformly bounded; **(B7)** an LICQ-type condition on `∇_u c` *together
  with* the gradients of the active bounds.
- **P2 Theorem 2** (the barrier result that lets Theorem 1 apply to the interior-point
  extension) explicitly needs **`s_L = 1`** (P2 §5: *"it is necessary to assume that `s_L`
  is chosen to be 1 to establish this result"*). **`Options.s_L` defaults to 2.3**
  (`src/options.jl:30`). So with default options the inequality-constrained path is not the
  algorithm Theorem 2 covers.

And the standard caveat: "global convergence" here means *convergence to a first-order
point from an arbitrary start*, not convergence to a global minimum. For nonconvex BFM
MPOPF there are multiple local solutions and neither theorem says anything about which one
you land on.

## 7. Does the recursive structure survive per-period network constraints?

**The recursion itself: yes. The conditioning: unverified and I would not assume it.**

- The value-function recursion is over SOC only, so `V̂xx` stays `nx × nx` (§1). The Riccati-
  like propagation `V̂xx = C + βᵀB + ωᵀc_x` (`src/backward_pass.jl:185-189`) does not grow
  with network size. Structurally the recursion is fine.
- The per-stage KKT system, however, becomes the entire per-period OPF KKT system, solved
  densely by nullspace + Cholesky. Its conditioning is the conditioning of an OPF KKT
  matrix, which is well known to be poor near binding voltage limits — and the code has
  **no conditioning safeguard on the range-space solve** (`lu(AY)` unchecked,
  `linsolve_tol` unused; FILTERDDP_API.md §4). The only safeguard is the Cholesky on the
  reduced Hessian, which catches indefiniteness but not ill-conditioning.
- Evidence from Stage 6, and it cuts both ways. Under a badly-scaled objective (a quadratic
  cost on grid import, with the battery direction carrying almost no curvature) a 3-period,
  3-control QP stalled at `dual_inf ≈ 1.6e-8`, took 174 iterations at `tol = 1e-9` and
  failed outright at `1e-10`. Under the tADMM paper's objective, where the battery
  throughput coefficient `C_B` *is* the reduced-Hessian curvature, the identical case
  converges in 16 iterations and certifies down to `1e-12`. **[inferred]** So the barrier
  path is not hitting a fixed accuracy floor — it is highly sensitive to the conditioning
  of `Zᵀ H Z`. That is the pertinent warning for MPOPF: the tADMM paper deliberately sets
  `C_B ≈ 10⁻⁶·min_t c^t` so throughput cost cannot compete with energy cost, which is
  precisely the regime that produced the stall.

## 8. Two practical blockers I hit that are not about theory at all

1. **Time-varying data has nowhere to live.** `OCP` holds exactly one `Objective` for stages
   `1..N-1`, one for stage `N`, and one `EqualityConstraints` for all stages
   (`src/ocp/ocp.jl:4-11`). An MPOPF has a different demand, a different price and possibly
   a different topology every period. In Stages 5-6 I smuggled this in by adding a time
   index `τ` to the state and interpolating `d_t, a_t, b_t` with a Lagrange polynomial
   through the nodes — exact at the nodes, smooth, and it worked at `T = 3` and `T = 6`.
   **It does not scale**: `T = 24` or `T = 48` would need a degree-23 or degree-47
   interpolant, which is numerically hopeless and symbolically enormous. A Fourier
   representation would be better conditioned, **[inferred]** but the real fix is a
   per-stage cost/constraint API — a genuinely small change to `OCP`, `build_ocp`,
   `backward_pass!` and `rollout!`, and by far the highest-value modification if this line
   is pursued.
2. **No infeasibility diagnosis.** Stage 6g fed the solver a provably infeasible bound set;
   it returned status 7 ("line search failed") with a `0.225` power-balance residual — the
   same status it returns for a merely difficult problem (Stage 6d, which was *at* the
   optimum). **The caller cannot distinguish "infeasible" from "hard" from "converged but
   couldn't certify" without auditing residuals itself.** For MPOPF, where binding network
   limits routinely make subproblems infeasible, that is a serious operational gap.

## 9. Verdict

- **Copper-plate battery scheduling: FilterDDP works and is accurate.** Verified against
  independent references to `1e-15` (unconstrained, Stage 5) and `1e-10` (bounded,
  Stage 6), with exact active-set identification and machine-precision dynamics feasibility.
- **LinDistFlow MPOPF: plausible for a small feeder, untested here.** The assumptions are
  benign (linear, convex); the binding issue is `nu` and the missing per-stage API.
- **BFM-SOCP MPOPF: I would not expect this implementation to be competitive**, on grounds
  of dense `nu × nu` stage algebra (§3), a generically fully-active conic constraint set
  (§4), and the observed `~1e-9` accuracy floor (§7).
- **Nonconvex BFM MPOPF: no basis for confidence.** Assumptions G1/G3/B6 are unverifiable
  in advance, the covering theorem does not apply to the shipped code (§6), and there is no
  restoration phase to recover from the failures that will occur.

**An MPOPF can be written as a discrete-time OCP — that is necessary, not sufficient.** The
thing that makes FilterDDP fast on a 15-control cartpole (dense static arrays, exact
symbolic Hessians, tiny stage systems) is exactly the thing that makes it unsuitable for a
network with thousands of algebraic variables per period.

If this line is pursued, the order of work that would actually de-risk it is:
(1) per-stage costs/constraints, (2) sparse or LAPACK-backed stage algebra
(start from `FilterDDPLAPACK.jl`), (3) a feasibility restoration phase, (4) only then a
LinDistFlow feeder as the first real test.
