# Exact DDP transfer and bounded T24 penalty tuning

The two-system model regression passed in job31650190. It compared every
physical metric against frozen C4 and exercised deliberately invalid original
current, import, voltage, battery-power and energy bounds. Local MA57 original
feasibility was 3.85e-16 on med2522 and 4.84e-13 on large10k. All positive
lower-voltage bounds remain. Removed current lower bounds per period:
2521 med2522; 10320 large10k. The shared data-only voltage screen is reused.

Three large10k T24 candidates ran on identical frozen input and numerical
source a07d579d, with free terminal SOC and MA57. Each was capped at 100 steps,
with no full warm-up/repeat and no clean timing claim. Every subproblem status,
iteration/period pair and captured source hash passed validation.

- Job31650191: rho0=384, relative controller balance ratio10. It dropped to
  192 after iteration5 and 96 after iteration10. Final original violation
  1.3958e-4; no NO point by100.
- Job31650218: fixed384. Final original violation 3.977e-4; no NO point by100.
- Job31650238: fixed48. Final original violation 2.327e-5; no NO point by100.

The latter two use the policy-selectable launcher committed at e1d8da71,
archived separately under campaign, while the numerical source remains
unchanged. The actual orchestration SHA is recorded in each result. The
initial gate31650186 failed because its regression modules lacked an include
binding; no optimization ran. Its logs and source bundle are preserved.

These are recorded failed targets, not completed full-convergence results.
No NO/OpenDSS pass is claimed when no qualifying dispatch existed. The
accepted NO rule remains objective error<=0.1% AND original violation<=1e-5,
and the user target is the first qualifying iteration strictly below100.
The C4 reference is a separately fully converged and OpenDSS-validated run;
its first trace crossing was168, with the same frozen inputs and weights.

Full evidence: ddp-transfer-three-calibrations.tgz. It includes source
snapshots, all traces/statuses, hardware/affinity, resource samples, manifests,
gate results, raw failed-harness logs and orchestration. Serialized Julia
checkpoints remain on Kamiak. Archive SHA256:
`336b0a9c08110b1fc01634009ae5b69d8a12cfcc52d3304ce314565deb14a068`.

The separate balance-ratio2 trial uses newly frozen driver source c06bb64b.
Its common model and six shared numerical/validation helpers are byte-identical
to the passed gate by cmp; that model validation is explicitly reused. The
driver exposes the previously hardcoded balance threshold, retaining default10.
Controller changes, first crossing and OpenDSS results must be assessed from
that trial's own completed evidence.
