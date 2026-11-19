# Matched MUMPS comparison

After the completed MA57 full-worker comparisons, test MUMPS for both central
and tADMM solves. Start with large10k T48,48 workers, rho120000 and the selected
primal-dual warm-start/batched-update implementation. Keep identical exported
data, July C_B, periodic profiles, free terminal SOC, corrected losses,
screening, initial values, controller, full warmup/reset and quality criteria.
Local and centralized Ipopt tolerances remain tight and local solves exact.

Use the same one-node Gold-6438Y+64-CPU/128-GiB allocation as MA57; one solver
thread per process,49 active compute processes for tADMM, separate sampler CPU.
Record both complete central reference timings and actual tADMM core/loop time.
An OOM, time limit or budget-limited solve is an outcome, not a speedup claim.
Inspect T48 before extending the backend comparison to T96.

The new optional solver argument defaults to MA57, preserving prior launch
behavior. MUMPS omits HSL-specific attributes. The independent validator checks
the actual centralized Ipopt log header against the recorded backend.
Run the repeated-solve/multiplier/objective unit gate natively with MUMPS before
the full T48 comparison. No rho or stopping tolerance tuning is part of this
backend comparison. Large-array/SOC/OpenDSS audits remain outside core timing.
