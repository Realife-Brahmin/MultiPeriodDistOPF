# July battery weight and slow-subproblem comparison

## Current stopping-rule decision (November 8)

The user explicitly selected the historical **owned-period RMS consensus**
threshold 1e-4 for the 10k comparison. This is not a relaxation of the maximum
original-model feasibility tolerance. The historical tadmm_socp.jl rule is
owned RMS <=1e-4 AND (unique consensus-change dual RMS <=1e-2 OR objective
range/max <1e-3 over the last five iterations, with positive maximum).
The dual metric counts each battery/time once, unlike today's replicated-copy
dual RMS. Objective agreement with centralized OPF is not part of this exit rule.

historical_stopping.jl reproduces that definition. assess_historical_stop.jl
uses it to terminate one large10k T24 trajectory, reusing the independently
validated July-weight reference and holding fixed rho=48 and exact MA57 locals.
The old stricter NO metric, all-copy residuals, original feasibility and actual
endpoint OpenDSS diagnostics remain separately logged. A historical stop does
not set the strict-convergence or physical-NO success flag.

The previous 200-step assessments remain frozen evidence under their original
stopping rules. They cannot be exactly rescored from the saved CSVs, which lack
the owned-period metric. A fresh prefix is checked against the prior numerical
trace and initial-start hashes. No model or penalty change is mixed into this
stopping-rule assessment. Other existing drivers keep explicit strict defaults.

Source caveat: commit 5e8815d2c changed the monolithic large10k threshold from
1e-4 to 1e-3; the July split run_tadmm.jl defaults to 2e-3 with overrides and
requires the dual condition rather than stagnation. The selected 1e-4 rule is
the user's requested historical variant, not a claim that every July run used
identical thresholds. Local Ipopt acceptable_tol=1e-4 is a different setting.

The user requested restoring July's C_B and investigating the earlier strategy
of cutting off slow local solves. This supersedes the earlier choice of DDP's
system weights for the present tADMM workstream. SPAR-DDP defaults and all
previous frozen results remain unchanged.

## Recovered settings

The July 31 release is commit 50b8901c324ba9eefd0d672c6e484e0ededc7a57.
Its config.jl and tadmm_socp.jl set C_B=1e-6*minimum(LoadShapeCost), with
prices constructed on range(0,2pi,length=T). Thus at T=24 the exact legacy
coefficient is **8.013987384856764e-8**, the same for all systems. This is
approximately 23709 times smaller than the DDP medium weight and299 times
smaller than the DDP large-system weight. It is horizon-dependent because
the sampled minimum is horizon-dependent. experiment_settings.jl reproduces
that coefficient from the original rule; it does not silently alter today's
periodic input arrays to obtain it.

The July local solver uses max_iter30, tol 1e-6, acceptable_tol 1e-4,
mu_strategy=adaptive, MUMPS and bound_relax_factor 1e-8. It accepts
ITERATION_LIMIT and ALMOST_LOCALLY_SOLVED. The explicit30-step cap was added
in commit 7de323d350ac973b5f84baee2b82126f719b5d4e (April 5), which reports
a large 10k T=12 trial. A later split run_tadmm.jl uses a 50-step cap. This
is an iteration cap, not evidence of a wall-clock timer killing workers.

The historical reporting criteria also differ: paper-era rows use looser
primal surrogates and a 0.5% near-optimal objective rule. Today's criterion is
0.1% objective agreement AND original-model feasibility 1e-5 on med/large,
with every local energy copy monitored. July's original r_norm monitors
owned-period copies only. These differences must remain visible in any
historical iteration comparison.

## First controlled experiment

Use current frozen native T=24 network exports, periodic load/PV/price curves,
initial SOC, loss-aware power-flow starts, free terminal SOC and MA57.
Retain the exact, proved bound cleanup. Change C_B to the exact July value.
Recompute an exact centralized reference; do not reuse DDP-weight objectives.
Pass the same numeric C_B to the independent OpenDSS checker.

Keep each system's current selected rho setup to isolate the weight change:
med 3800 with the existing relative/balance10 policy; large fixed 48. This is
explicitly a controlled July-weight experiment, not a replay of July's rho
controller or a claim that these rhos are optimal at the restored weight.

For each system compare two arms on identical data and fresh primal/dual starts:

- exact: current accurate local-solve tolerances, maximum 3000 Ipopt iterations;
- cap30: identical tolerances and model, but a 30-iteration local budget.

The cap-only arm intentionally does not simultaneously loosen tolerances or
change barrier strategy. This separates the benefit of truncating outliers
from other changes in the historical local-solve recipe. An ITERATION_LIMIT
point is provisional, must have finite primal values, retains its real status,
and records its local primal violation. It is not labeled a solved NLP.
The centralized reference always rejects iteration limits.

Three discarded compile steps precede each assessment, followed by a complete
model/primal/dual reset. Worker-side initial-start SHA256 values must match
between arms. Each arm gets one assessment, up to 200 outer steps initially;
NO is scored from its first original-feasible crossing. Actual NO dispatches
are saved, hashed and independently OpenDSS-checked. Strict global criteria,
local truncations, and capped unsuccessful endpoints have separate flags.
These short-warmup observations are not finalized publication timings.

## Current evidence for the outlier hypothesis

In the prior DDP-weight large fixed 48 pilot, all 2376 local solves after outer
iteration 1 used at most 24 Ipopt steps. A30-step cap would therefore have no
effect on those later solves. At outer iteration 2, median/max solver time
was 0.572/0.980 s; at 100,0.680/0.818 s. Restoring July's weaker curvature may
change this behavior, so analyze the new matched stage logs rather than assume
the historical outlier problem persists unchanged.

The comparison must report cap-hit counts, affected periods, max/median
stage times, outer iterations to both the 100 target and actual NO, objective
gap, original feasibility, and independent OpenDSS checks. No iteration-limit
status or historical surrogate residual substitutes for original feasibility.

The first cap arm is instrumented: it evaluates a full JuMP local primal
report at every cutoff. On the medium pair this reduced the median slowest
native Ipopt solve from 0.331 to 0.307 s, but increased critical-worker framework
overhead from 0.024 to 0.066 s. An explicit lightweight timing check disables
that redundant local report with TADMM_LOCAL_PRIMAL_REPORT=0. Finite-value
checks, real cutoff statuses, every original global feasibility check, and
independent OpenDSS at NO remain enabled. Its entire trajectory and initial
starts are compared against the instrumented cap arm, using the same validated
reference. Missing local reports are marked unevaluated/NaN rather than zero.

The completed [assessment and raw evidence](results/july_weight_comparison/README.md)
record medium NO at 126/125 for exact/capped locals and validated OpenDSS at those
actual points. The lightweight cap preserves the 200-step trajectory and gives
a modest observed loop-time saving. The large cap never activates (all solves
use at most 24 steps), and neither arm reaches NO within 200. Its saved endpoint
is limited by accumulated SOC-bound error, while a controlled worst-copy local
probe confirms accurate local solves. These observations do not establish an
optimal penalty for the restored weight or a final benchmark timing.
