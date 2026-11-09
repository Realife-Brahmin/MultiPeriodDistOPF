# Focused large10k T24 early-inexact experiment

User authorized this single focused test. MA57,24 period workers, one Gold6438Y+
node, July C_B, corrected losses, periodic inputs, idle-battery/load-flow starts,
no terminal SOC condition, rho15000 reduced by five every five rounds to48,
primal-dual reuse and batched updates match the completed exact reference31660055.
That reference took49 rounds /50.4500s core /54.7688s loop; central33.7901s;
final maximum reconstructed energy-bound violation0.86962kWh. Strict physical
feasibility was false. Old reference is not relabeled as a fully feasible OPF.

Candidate: first local solve exact; max30 Ipopt iterations during ADMM rounds2–30.
An interrupted primal is accepted only if finite with finite reusable multipliers
and the unscaled JuMP local feasibility report <=1e-6. Rejected attempts continue
from the interrupted point to exact tolerances; attempt, screening and recovery
work are included in core/stage elapsed time. Attempt counters and primal error
are logged. From round31 onward all local solves are exact. Final RMS/dual stop
is prohibited before the exact phase. No local or global final tolerance changes.

Native gates force both acceptance and rejection, verify exact recovery and
phase restoration, and test the observer's exact-final guard. Existing425 checks
also run. Full warmup/reset numerical trajectories must agree. Independently
validate all status exceptions and their phase/feasibility gate, saved milestones,
final SOC consistency and OpenDSS. Use the same centralized reference protocol.
Do not claim a speedup if improved time comes with degraded quality. This is one
exploratory measured trajectory, not a tuning sweep or repeated median benchmark.

Remote: `/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_early_cap`.
Reproduction launcher: repository-root `setup-early-cap-nov09.sh`.
Frozen numerical source db1bb2bfddc2c2fc4386ac1d85a2a82534569c47.
Native gate31663363 passed455 checks (425 existing,26 acceptance/recovery/phase,
4 endpoint-guard). Run31663364 follows the gate; passive energy audit31663365
follows the run. Initial run status RUNNING; no result is asserted by submission.

## Direct finite-iterate follow-up

The screened candidate measured49 rounds /52.3815s core, versus baseline49/
50.4500s, with essentially identical endpoint physics. Seven capped outputs
were accepted; the two expensive period2 attempts at rounds22/23 were rejected
because local feasibility at30 was5.7633e-4/1.3318e-2. They still consumed102/172
total Ipopt iterations after exact recovery. Thus the screen prevented cutting
the target outliers. This negative result is preserved separately.

Test the same30-step early cap with `TADMM_CAP_SCREEN=finite`: allow finite
local primal/multipliers without an early feasibility requirement, record actual
constraint errors, then require all-exact solves from31 and unchanged final
global checks. This is a distinct inexact-ADMM experiment, not loosening final
feasibility or relabeling early dispatches. The candidate selection follows the
observed rejected attempts; do not claim it was an untuned heuristic. All local
error evaluation remains counted conservatively. New native gate explicitly
forces a finite accepted iterate above the screened tolerance and then verifies
restoration of exact feasibility. Remote campaign has suffix `_early_cap_finite`.
