# Two targeted T48 changes after the full-worker comparisons

The exact MA57 full-worker baseline remains T48/rho120000:322 iterations,
334.2166 s core, same profiles/weights/model and current residual/dual target.
T96 remains slow and oscillatory. Its converged trace changed rho46 times,
with116 rounds having a local solve above20 Ipopt iterations (local95th
percentile11). T48 had21 rho changes,53 such rounds (95th percentile13).

Test two arms separately; do not combine changes or weaken accuracy:

1. `warm_retry`: same legacy controller, same primal-dual warm start. Interrupt
   only a warm attempt that reaches20 Ipopt iterations, discard that attempt,
   and restart from the previous fully solved primal with dual warm start off.
   Recovery must reach the original tight local tolerances. Record both attempts'
   solver time and summed iteration counts, interruption and recovery fields.
   The final local status must still be LOCALLY_SOLVED/OPTIMAL. This does not
   accept an inexact subproblem solution. The threshold20 is above the observed
   95th percentiles; it is one tested implementation choice, not a fitted law.
2. `bootstrap_fixed`: retain exact ordinary local solves. Start at120000, divide
   rho by5 every5 iterations, clip at48 and then keep rho fixed; no subsequent
   residual-balancing/watchdog cycle.48 is an exploratory candidate informed by
   the previous T24 fixed-48 result, not a validated horizon-scaling heuristic.
   Rescale scaled duals by old/new rho at every change. This directly tests
   whether repeated penalty adaptation is causing the late oscillations.

Use48 independently pinned workers on one matched Gold node,64 allocated CPUs,
128 GiB, one thread per process, MA57, full warmup/reset and max1000 outer
iterations. Preserve all original milestone definitions, cost and physical
flags, saved dispatches, independent OpenDSS and passive energy checks. Both
arms keep the same initial rho and all model/data settings as the baseline.
Judge total core/loop time to the same goal, including recovery work; retain
failed/budget-limited outcomes. Inspect T48 before extending a candidate to96.

Before submission, force an interruption after one iteration on the nonlinear
unit problem and verify exact recovery, primal/objective agreement, coefficient
updates and complete status logging. Test the deterministic penalty schedule
independently. Existing generic defaults remain unchanged.
