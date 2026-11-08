# Full temporal parallelism and exact local-solve reuse

User correction on November 8: use one worker per period (48 for T48, 96 for
T96), then assess full Ipopt warm starts and objective-update overhead before
any inexact local solves. MA57, physical model, profiles, weights, rho controller
and quality criteria remain those of the preceding frozen sweep.

The completed 24-worker T24/T48 evidence remains intact. The two unfinished
24-worker T96 batches were explicitly superseded and stopped before completion;
partial traces and cancellation/resource records are retained. They are not
converged timing results. Their previous stop-after-first watcher and dependent
audit were cancelled because neither selected arm completed validation.

## Execution

T48 has 48 pinned worker processes plus coordinator on one Gold-6438Y+ node.
T96 has 96 pinned workers (48 per host) plus coordinator on two matched nodes.
Reserve 64 physical CPUs and 128 GiB per node; worker BLAS/Julia/OpenMP threads
are one. Reserve a separate CPU for each node's passive sampler. Reject
duplicate host/CPU assignments or more than one period per worker. Use
SlurmClusterManager 1.1.0: a non-solver probe confirmed that compute-node SSH is
not available, so the worker launch must use srun inside the allocation.

Warm-start variants are explicit: `primal` retains the prior primal-only path;
`primal_dual` caches all primal and constraint/bound dual values before any
model mutation, supplies them to Ipopt and enables warm_start_init_point.
The first local solve uses the original starting values. Centralized reference
measurement keeps its prior full warmup and fresh primal reset protocol; local
warm-start changes do not seed central reference measurements with optimal duals.

`batched` updates the energy linear objective coefficients in a single JuMP
call and changes quadratic coefficients only when rho changes. `legacy`
retains the previous scalar update path for an equivalence gate. Timers now
separately report warm-start transfer and objective updates; both remain inside
the core solve benchmark. No constraint or local solve tolerance is removed.

## Validation and experiment order

1. Exercise equality, inequality, cone, and variable-bound multiplier transfer
   and unchanged/changed rho on a small nonlinear problem. Compare final primal
   values, objective and feasibility for all three implementation variants.
2. A short T96 execution gate uses all 96 workers on two nodes and records
   actual host/CPU/period assignments. Compare legacy and batched primal-only
   numerical paths; exercise the full-dual path independently. These short
   runs are launch/correctness checks, not converged benchmarks.
3. Run T48 at rho120000 with primal-only/batched and primal-dual/batched,
   with full warmup/reset and the existing saved-point/OpenDSS checks. Keep
   both outcomes; select warm-start policy using time to the same quality goal.
4. Run corrected T96 with 96 workers using the selected policy. Retain the
   two originally requested starts (60000 and240000), including budget-limited
   outcomes. Do not change rho tuning in this implementation comparison.

The original centralized reference and each new reference are retained; no
slower central timing is cherry-picked. Record total allocated CPUs, CPU time,
and node-local PSS separately. Node samplers register each actual Julia PID;
their independent samples must not be presented as synchronized global peaks
without an explicit aggregation rule. Validation remains outside core timers.
