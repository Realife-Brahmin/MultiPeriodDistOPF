# How FilterDDP's stage KKT matrix changes between iterations

Question from the meeting with Anamika (2026-09-26). Reusing a stale stage
factorization diverged in the general case and worked only on the easiest
problems. Which blocks of the KKT matrix move between iterations, and by how
much? Measured 2026-09-26 on the 309 lab PC. Nothing here is timed.

## Setup

- **Runs.** One untimed FilterDDP solve per system at `T = 6`, in the paper's
  current formulation:
  - periodic profile, per-system `C_B`, soft terminal SOC
  - diagonal stage Hessian with the exact assembly rewrites
  - run to FilterDDP's own convergence, with no near-optimality stop:
    ieee123 75 iterations, med2522 81, large10k 128
- **Capture.** The periodic capture hook recorded every stage's
  `K = [H cu'; cu 0]`, `Σ_L + Σ_U` and the incoming `V_xx`. It captured every
  iteration for ieee123 and med2522, and every 10th for large10k (a large10k
  snapshot is about 2.4 GB).
- **Analysis.** `kkt_evolution_analysis.jl` compares each capture with the
  previous one of the same stage. `kkt_evolution_drift.jl` measures drift over
  1-10 iterations from the kept ieee123 captures. `summarize_kkt_evolution.py`
  writes `ddp/results/kkt_evolution/SUMMARY.md` (full per-block tables) and
  `figures/kkt_evolution_<system>.png`. Reproduce with
  `bash ddp/examples/power_system/run_kkt_evolution.sh`, then the summary
  script.
- **Blocks.**
  - **Hessian.** Its diagonal is split into three parts: the barrier terms
    `Σ_L + Σ_U`; the value-function term `dt² diag(V_xx)` on the battery
    powers; and the rest. The rest is the SOCP-multiplier curvature, the
    constant cost terms and the `1e-8` floor. The diagonal is also split by
    variable type.
  - **Jacobian.** It is split by constraint type.
- **Metric.** Changes are `||new − old||_F / ||old||_F`. The Frobenius norm is
  dominated by the largest entries, so SUMMARY.md also gives the largest
  per-entry change factor.

## Result

**1. Most of the matrix cannot change at all.**
- Every linear constraint row is constant in every iteration, on every
  system: power balance P and Q, voltage drop, root voltage and battery
  energy.
- The only Jacobian rows that move are the branch-current (SOCP) rows
  `P² + Q² − v·ℓ + s = 0`, one per line: 23-24% of the constraint rows.
- The constant cost curvature (`C_B`, terminal term) is fixed too.

**2. Median change per iteration by block, split by phase.** "Barrier phase"
means `μ ≥ 1e-6`; "final" means `μ < 1e-6`. large10k compares captures 10
iterations apart and never reached the final phase in a captured iteration.

| block | ieee123 barrier / final | med2522 barrier / final | large10k (per 10 it.) |
|---|---:|---:|---:|
| whole K | 170% / 0.17% | 54% / 7.1% | 74% |
| barrier terms Σ | 190% / 0.17% | 57% / 7.1% | 110% |
| V_xx term (battery powers) | 0.31% / 2e-6 | 1.6% / 3e-6 | 1.6% |
| other Hessian curvature | 4e-5 / 5e-7 | 2e-5 / 6e-6 | 0.14% |
| Jacobian, SOCP rows | 16% / 0.05% | 2.5% / 5e-6 | 54% |

**3. Almost all of the change is in the voltage-limit barrier terms.** The
barrier terms of the squared-voltage variables `v` carry a median of 100%
(99.7% at large10k) of `||ΔK||²`. The other bounded variables change 20-30% per
iteration, but their barrier entries are orders of magnitude smaller.
Examples: `ℓ ≥ 0`, battery power and energy slack bounds, DER reactive limits,
SOCP slacks.

**4. Over several iterations the change compounds.** This is what a frozen
factorization would have to absorb. ieee123, median over stages and
iterations:

| iterations frozen | K (barrier phase) | K (final phase) | non-barrier Hessian (barrier phase) | Jacobian (barrier phase) |
|---:|---:|---:|---:|---:|
| 1 | 1.7x | 0.2% | 0.13% | 13% |
| 2 | 6.9x | 0.2% | 0.26% | 25% |
| 5 | 218x | 1% | 0.9% | 58% |
| 10 | 49,000x | 3% | 1.2% | 83% |

In the barrier phase the barrier entries grow by orders of magnitude within a
few iterations. At large10k the 90th percentile of the 10-iteration change is
`3e4`, concentrated at barrier-parameter updates. Once `μ` reaches its floor
(`1e-8`), the matrix barely moves: ieee123 changes 3% over 10 iterations;
med2522 still changes 7% per iteration below `μ = 1e-6`.

**5. Some tiny entries toggle.** The `P` and `Q` diagonal entries switch
between the `1e-8` floor and `O(1)` curvature from the SOCP multipliers,
changing by up to `5e7x`. They are negligible in norm, but they are pivots,
so a stale factor sees them wrong too.

## Reading

The matrix changes because of the interior-point barrier, not the network
physics. The barrier curvature of a bounded variable is `Σ = z/(distance to
the bound)`.
- **Binding bounds:** as `μ` falls, the distance shrinks and the multiplier
  settles, so `Σ` grows roughly like `1/μ`.
- **Slack bounds:** `Σ` decays.

On these feeders the binding bounds are the voltage limits, which is why the
voltage barrier terms dominate. They jump at every `μ` update and keep
drifting within each barrier subproblem as the iterates approach the bounds.
Everything the network contributes is either exactly constant (the linear
rows) or moves by percent (the SOCP rows, the `V_xx` term).

That explains the stale-factorization result. A factor frozen for even 2-5
iterations in the barrier phase is off by one to two orders of magnitude in
exactly the diagonal pivots that decide the voltage-limited directions. Late
in the solve, or on easy problems where `μ` reaches its floor quickly, the
matrix is nearly fixed, which is where reuse did work.

Two implications, not tested here:
- **Refresh on `μ` updates.** A refresh schedule tied to `μ` updates and
  to the size of the barrier change would be the principled version of
  "reuse when safe".
- **Stale factor as a preconditioner.** Since only the diagonal and the SOCP
  rows move, using a stale factor to precondition an iterative solve could
  work in the final phase. It could not work in the barrier phase, where the
  diagonal changes by orders of magnitude.
