# Persistent inexact local caps, large10k T24

User explicitly removed the exact-local phase and requested the final quality
tradeoff. Test30 and10 Ipopt iterations per local solve from the very first
ADMM round through termination, including discarded warmup. No uncapped first
solve, exact recovery, late exact phase, or final local polishing is permitted.
Ordinary local solves may converge naturally before the maximum. Nonfinite
primal/multipliers abort the candidate rather than trigger exact recovery.

Keep MA57,24 period workers, one Gold6438Y+ node per arm, July C_B, current
periodic inputs, corrected losses, free terminal SOC, idle/load-flow starts,
rho15000 -> fixed48, primal-dual reuse, batched updates and thread policy.
Full centralized solves remain independent quality references, not warm starts
for tADMM. Both candidate arms have full warmup/reset and source/input manifests.

Use the unchanged owned/all-copy RMS1e-5 plus dual1e-2 consensus endpoint or1000
outer-iteration budget. Report objective agreement and original-model violation
separately; reaching consensus is not proof of OPF feasibility. Independently
audit integrated battery SOC, dynamics, power bounds, relaxation gap and OpenDSS
network agreement. Preserve failures, plateaus and invalid early points.

Finite-only persistent acceptance does not compute a local feasibility report
inside each stage: it does not use one to make algorithm decisions. Original
dispatch physical diagnostics remain outside core as before; their wall cost,
all actual core/loop time and warmup/setup are logged. This differs from the
prior early-cap assessment's conservatively timed local reports; report that
boundary change when comparing policies. Persistent30 versus10 uses identical
boundaries. There is no hidden final feasibility gate or exact correction.

Native tests use a one-iteration cap at outer1,2,31,1001, verify bounded work,
finite iterates, no recovery and that the optimizer limit never becomes exact.
Independent collection validation checks every warmup/measured local row,
including the final round, against its cap and allowed statuses. No final exact
guard applies to this policy. Selected30/10 values are experimental settings,
not a validated tuning heuristic or repeated-median paper benchmark.

Exact reference31660055:49 rounds /50.449995s core, central33.790131s; final
energy-bound error0.869624kWh, strict physical flag false. Prior early-cap
follow-up31663479:52 rounds /49.471650s,0.668225kWh, strict physical flag false.

Remote: `/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps`.
Reproduction: repository-root `setup-persistent-caps-nov09.sh`.
