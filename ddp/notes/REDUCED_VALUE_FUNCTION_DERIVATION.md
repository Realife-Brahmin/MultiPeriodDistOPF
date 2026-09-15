# Why the reduced-space method stalls: a derivation

Written 2026-09-15 after three empirical fixes failed in a row (fixed penalty
coefficient, Hessian refresh, ruling out the L1 kink). The failure is
structural, not a tuning problem, and the derivation below says why.

Companion to [REDUCED_SPACE_INNER_OPF_FEASIBILITY.md](REDUCED_SPACE_INNER_OPF_FEASIBILITY.md),
which carries the measurements.

## 1. The reduced problem

Full problem, in the repository's notation:

```
min   sum_t [ c^t S dt P_Subs^t  +  C_B S^2 dt ||P_B^t||^2 ]
s.t.  B^t = B^{t-1} - dt P_B^t                        (dynamics)
      g_t(y^t, P_B^t) = 0                              (balance, voltage drop, SOC)
      y^t in Y_t                                       (v box, ell >= 0, |q_norm| <= 1,
                                                        P_Subs >= 0, soc_slack >= 0)
      P_B^t in [-P_B_R, P_B_R],   B^t in [B_min, B_max]
```

`y^t` collects every algebraic network variable. Define the **reduced stage
value** by partial minimisation over `y` at fixed `P_B`:

```
Phi_t(P_B) = min_y { c^t S dt P_Subs : g_t(y, P_B) = 0, y in Y_t }
```

giving the reduced problem

```
min   sum_t [ Phi_t(P_B^t) + C_B S^2 dt ||P_B^t||^2 ]
s.t.  B^t = B^{t-1} - dt P_B^t,  box constraints on P_B^t and B^t.
```

Partial minimisation is exact, so the reduced and full problems share optimal
values and `(B, P_B)` solutions wherever `Phi_t` is defined, i.e. for
`P_B^t` in the projected feasible set `F_t`. **This part is not in question** --
it is verified to `1e-9` against the full-space solutions.

## 2. The gradient, and the hypotheses it rests on

Write the inner Lagrangian `L(y, lambda, nu; P_B) = c S dt P_Subs +
lambda' g_t(y,P_B) + nu' (bound terms)`. At the inner optimum, the envelope
theorem gives

```
dPhi_t/dP_B = ∂L/∂P_B = lambda' ∂g_t/∂P_B
```

and since `P_B,b` appears only in the real-power balance row of battery bus `b`
with unit coefficient, `dPhi_t/dP_B,b = lambda_b`, the balance-row dual --
available free from every inner solve, and confirmed against central differences
to `1e-9`.

**The hypotheses matter, and were never checked.** The envelope theorem in this
differentiable form needs LICQ, strict complementarity, and second-order
sufficiency at the inner solution. Under those, `y*(P_B)` and `lambda(P_B)` are
`C^1` and `Phi_t` is `C^2`.

## 3. Where the regularity actually fails

Let `P_B` move so that some bound in `Y_t` changes status -- a bus voltage
reaching its limit, a DER inverter reaching its reactive limit. Strict
complementarity fails at that point, and:

* `Phi_t` stays **convex** (partial minimisation of a convex function over a
  convex set -- the BFM SOCP is convex), and
* `grad Phi_t` stays **continuous** (`lambda` is unique under LICQ), but
* `hess Phi_t` is **discontinuous across the active-set boundary**.

So

> **`Phi_t` is `C^1` but only piecewise `C^2`.**

That is exactly the wrong regularity for a second-order method. The gradient
FilterDDP consumes is sound; the Hessian it builds its local quadratic model
from is valid only on one side of a surface the step is free to cross.

**The active set here is large and `P_B`-dependent.** Measured on ieee123,
`|q_norm| = 1.0` at 33-51 of the 51 DERs in every single solve: reactive
capability is saturated because it buys loss reduction, and *which* inverters
saturate depends on the dispatch. Every outer step reshuffles the active set.

## 4. The predicted signature, and the observed one

Predicted: smooth, fast convergence while the iterates stay inside one
active-set region; a sudden loss of dual progress when a step crosses a
boundary; then oscillation between regions with the objective already essentially
correct.

Observed (ieee2522, `T = 12`, `C_B = 1e-3`, verbose trace):

```
iter   objective        du_inf       cs_inf     lg(mu)   alpha
   0   9.16353349e+03   4.8089e+02   8.3002e-03   0.00
  ...                    (clean convergence)
  13   8.57830272e+03   9.1360e-02   8.0418e-03  -3.49    1.0
  14   8.57829614e+03   9.1888e-03   8.0001e-03  -3.49    1.0   <- best point
  15   8.56980773e+03   3.2424e+00   4.6208e-03  -3.49    0.5   <- breaks
  ...
  21   8.56076846e+03   2.8367e+01   3.2000e-04  -3.49    1.0
  23   8.56076308e+03   2.6451e+01   3.2000e-04  -3.49    1.0
```

`du_inf` falls four orders of magnitude in 14 iterations, then jumps by three
and oscillates between 1 and 28. The objective is flat in its fifth digit from
iteration 17 on, and the final objective gap to full space is `2.7e-05` -- the
iterate is essentially AT the solution and cannot certify it.

**Why it is fatal rather than merely slow.** `lg(mu)` is pinned at `-3.49` from
iteration 14 onward and `cs_inf` sits exactly at `mu = 3.2e-04`. The outer
interior-point method only tightens its barrier once dual infeasibility is small
enough; once `du_inf` bounces back above that gate it can never tighten again,
so every remaining iteration burns at one barrier level. That is the whole
failure in one line.

## 5. What this explains that was previously unexplained

| observation | explanation |
|---|---|
| fixed `rho` changed almost nothing | the nonsmoothness is in `Phi` itself, not in the penalty |
| Hessian refresh changed almost nothing | re-sketching yields the *correct* Hessian of the *wrong side* |
| exact vs low-rank curvature barely differed | both are exact within a region; neither is valid across a boundary |
| refresh drove infeasible evaluations to 0 yet still failed | feasibility was never the blocker; the L1 kink is not the kink that matters |
| `lambda_min(d2Phi) = -0.148`, "within FD noise" | a finite difference taken ACROSS an active-set boundary is not a second derivative at all |
| `T = 3` converged everywhere, `T >= 12` never does | shallow cycling keeps the path inside one region; deep cycling sweeps through many |

The last row is the one that caused the wasted effort: `T = 3` success was read
as evidence the method worked, when it was evidence the test was too easy.

## 6. The fix the derivation implies

Keep the inner barrier **open** instead of resolving the active set. Define

```
Phi_t^mu(P_B) = min_y { c^t S dt P_Subs - mu sum log(slack) : g_t(y,P_B) = 0 }
```

For `mu > 0` no inequality is ever active, so `Phi^mu` is `C^infinity` in `P_B`,
its Hessian exists and varies smoothly, and `Phi^mu -> Phi` as `mu -> 0`. Tying
`mu_inner` to the outer barrier lets both tighten together, which is the
standard construction for nested interior-point problems.

Mechanically this is Ipopt's `mu_target`.

**One trap worth recording.** If `Phi^mu` is the barrier problem's value then its
gradient is the barrier problem's dual, so the *value* reported to the outer
method must include the barrier term. Reporting the unbarriered cost alongside a
barrier-problem dual hands the outer method an inconsistent value/gradient pair
-- the same class of defect as the adaptive-`rho` inconsistency. The
implementation therefore computes `-mu sum log(slack)` over every finite
variable bound explicitly and sets `bound_relax_factor = 0` so that accounting
matches Ipopt's own.

## 6a. The barrier remedy was tested and did NOT work

Corrected run (ieee2522, `T = 12`, `mu = 1e-4`, fixed `rho`): status 8, `du_inf`
oscillating between 15 and 154, `alpha` down to `1e-16`, 17-52 line-search
backtracks per iteration. It never reached the `du_inf = 9.2e-03` that the
UNSMOOTHED run achieved at iteration 14. **Section 6 is not supported by this
test.**

Getting there took two self-inflicted detours worth recording, because both
produced convincing-looking failures that meant nothing:

* `tol` must scale with `mu_target`. Ipopt's NLP error includes complementarity,
  pinned at `~mu`, so `tol = 1e-10` at `mu = 1e-4` is unsatisfiable; it surfaces
  as `SLOW_PROGRESS`, and 5 of 12 stages returned the `1e12` failure sentinel.
  The giveaway was the failure count being non-monotonic in `mu` (9/12 at 1e-4,
  10/12 at 1e-5, 0/12 at 1e-6).
* `bound_relax_factor = 0`, set to make the barrier bookkeeping exact, stalled 5
  of 12 stages on its own.

Also: the "infeasible dispatch" counter is **meaningless under smoothing** (the
barrier holds slacks off zero, so every evaluation trips a `>1e-7` test), and
the barrier term itself was not the problem -- including or excluding it from
the reported value moved value/gradient agreement only from `3.39e-02` to
`3.37e-02`.

## 6b. The gradient identity, verified on the system that fails

Central differences against the balance dual, step-size study to separate
finite-difference noise from genuine error:

| configuration | best relative error | verdict |
|---|---|---|
| cold solves, `rho = 0` | `1.5e-11` | gradient **exact** |
| cold solves, `rho = 1e4` | `1.0e-09` | gradient **exact** |
| warm solves, `rho = 1e4` | `3.0e-04` .. `3.0e-03` | contaminated |

So Section 2's identity holds exactly on ieee2522 `T = 12`, not merely on the
ieee123 `T = 3` case where it was first checked.

The warm-start error follows a clean `1/h` law (`3.020e-05`, `3.020e-04`,
`3.020e-03` at `h = 1e-3, 1e-4, 1e-5`) -- the signature of a constant absolute
error of about `8.7e-6` USD in the **values**. The duals are identical warm and
cold to nine digits, so warm starting corrupts the value, not the gradient. That
matters because FilterDDP's filter accepts or rejects steps by comparing values.

**This is recorded as measured, not as the cause.** Per-step objective decreases
near the breakdown are `~1e-2`, three orders above `8.7e-6`, so the arithmetic
does not obviously support value noise as the blocker. The clean test is a cold
run at `T = 12`; it costs roughly 4x.

## 7. Status

Sections 1-4 are established: the formulation, the gradient identity (now
verified to `1e-11` on the failing system, Section 6b), and the trace. Section 5
remains the best available explanation of the T-dependence but is **inference
from a matching signature, not a verified mechanism** -- no measurement yet
counts active-set changes along the failing trajectory.

Section 6's remedy was implemented and **failed** (Section 6a), so the
piecewise-`C^2` diagnosis has not been confirmed by a successful fix. Either the
diagnosis is incomplete, or barrier smoothing at a fixed `mu` is not the right
form of it -- for example because the outer method runs its own barrier
continuation and the two are not coupled.

The direct test of Section 5 has still not been run: instrument the inner solves
to record which bounds are active, and check whether the active set changes at
exactly the iteration where `du_inf` jumps. That is cheap and would settle the
diagnosis rather than adding another remedy on top of an unverified one.

The broader lesson is structural: a nested value function is generically only
piecewise smooth, so "exact inner solve + second-order outer method" was never a
sound combination. Smoothing the inner problem is not a numerical convenience,
it is what makes the outer model well posed.
