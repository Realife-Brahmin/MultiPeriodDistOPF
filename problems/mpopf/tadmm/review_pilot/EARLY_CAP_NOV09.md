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
