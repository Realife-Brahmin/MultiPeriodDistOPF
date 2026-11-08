# July battery weight and slow-subproblem comparison

The user requested restoring July's C_B and investigating the earlier strategy
of cutting off slow local solves. This supersedes the earlier choice of DDP's
system weights for the present tADMM workstream. SPAR-DDP defaults and all
previous frozen results remain unchanged.

## Recovered settings

The July 31 release is commit50b8901c324ba9eefd0d672c6e484e0ededc7a57.
Its config.jl and tadmm_socp.jl set C_B=1e-6*minimum(LoadShapeCost), with
prices constructed on range(0,2pi,length=T). Thus at T24 the exact legacy
coefficient is **8.013987384856764e-8**, the same for all systems. This is
approximately 23709 times smaller than the DDP medium weight and299 times
smaller than the DDP large-system weight. It is horizon-dependent because
the sampled minimum is horizon-dependent. experiment_settings.jl reproduces
that coefficient from the original rule; it does not silently alter today's
periodic input arrays to obtain it.

The July local solver uses max_iter30, tol1e-6, acceptable_tol1e-4,
mu_strategy=adaptive, MUMPS and bound_relax_factor1e-8. It accepts
ITERATION_LIMIT and ALMOST_LOCALLY_SOLVED. The explicit30-step cap was added
in commit7de323d350ac973b5f84baee2b82126f719b5d4e (April5), which reports
a large10k T12 trial. A later split run_tadmm.jl uses a50-step cap. This
is an iteration cap, not evidence of a wall-clock timer killing workers.

The historical reporting criteria also differ: paper-era rows use looser
primal surrogates and a0.5% near-optimal objective rule. Today's criterion is
0.1% objective agreement AND original-model feasibility1e-5 on med/large,
with every local energy copy monitored. July's original r_norm monitors
owned-period copies only. These differences must remain visible in any
historical iteration comparison.

## First controlled experiment

Use current frozen native T24 network exports, periodic load/PV/price curves,
initial SOC, loss-aware power-flow starts, free terminal SOC and MA57.
Retain the exact, proved bound cleanup. Change C_B to the exact July value.
Recompute an exact centralized reference; do not reuse DDP-weight objectives.
Pass the same numeric C_B to the independent OpenDSS checker.

Keep each system's current selected rho setup to isolate the weight change:
med3800 with the existing relative/balance10 policy; large fixed48. This is
explicitly a controlled July-weight experiment, not a replay of July's rho
controller or a claim that these rhos are optimal at the restored weight.

For each system compare two arms on identical data and fresh primal/dual starts:

- exact: current accurate local-solve tolerances, maximum3000 Ipopt iterations;
- cap30: identical tolerances and model, but a30-iteration local budget.

The cap-only arm intentionally does not simultaneously loosen tolerances or
change barrier strategy. This separates the benefit of truncating outliers
from other changes in the historical local-solve recipe. An ITERATION_LIMIT
point is provisional, must have finite primal values, retains its real status,
and records its local primal violation. It is not labeled a solved NLP.
The centralized reference always rejects iteration limits.

Three discarded compile steps precede each assessment, followed by a complete
model/primal/dual reset. Worker-side initial-start SHA256 values must match
between arms. Each arm gets one assessment, up to200 outer steps initially;
NO is scored from its first original-feasible crossing. Actual NO dispatches
are saved, hashed and independently OpenDSS-checked. Strict global criteria,
local truncations, and capped unsuccessful endpoints have separate flags.
These short-warmup observations are not finalized publication timings.

## Current evidence for the outlier hypothesis

In the prior DDP-weight large fixed48 pilot, all2376 local solves after outer
iteration1 used at most24 Ipopt steps. A30-step cap would therefore have no
effect on those later solves. At outer iteration2, median/max solver time
was0.572/0.980 s; at100,0.680/0.818 s. Restoring July's weaker curvature may
change this behavior, so analyze the new matched stage logs rather than assume
the historical outlier problem persists unchanged.

The comparison must report cap-hit counts, affected periods, max/median
stage times, outer iterations to both the100 target and actual NO, objective
gap, original feasibility, and independent OpenDSS checks. No iteration-limit
status or historical surrogate residual substitutes for original feasibility.
