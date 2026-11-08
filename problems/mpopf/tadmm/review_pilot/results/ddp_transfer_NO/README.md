# Exact DDP transfer and bounded T24 penalty tuning

The two-system model regression passed in job 31650190. It compared every
physical metric against frozen C4 and exercised deliberately invalid original
current, import, voltage, battery-power and energy bounds. Local MA57 original
feasibility was 3.85e-16 on med2522 and 4.84e-13 on large10k. All positive
lower-voltage bounds remain. Removed current lower bounds per period:
2521 med2522; 10320 large10k. The shared data-only voltage screen is reused.

Three large10k T24 candidates ran on identical frozen input and numerical
source a07d579d, with free terminal SOC and MA57. Each was capped at 100 steps,
with no full warm-up/repeat and no clean timing claim. Every subproblem status,
iteration/period pair and captured source hash passed validation.

- Job 31650191: rho=0=384, relative controller balance ratio 10. It dropped to
  192 after iteration 5 and 96 after iteration 10. Final original violation
  1.3958e-4; no NO point by 100.
- Job 31650218: fixed rho=384. Final original violation 3.977e-4; no NO point by 100.
- Job 31650238: fixed rho=48. Final original violation 2.327e-5; no NO point by 100.

The latter two use the policy-selectable launcher committed at e1d8da71,
archived separately under campaign, while the numerical source remains
unchanged. The actual orchestration SHA is recorded in each result. The
initial gate 31650186 failed because its regression modules lacked an include
binding; no optimization ran. Its logs and source bundle are preserved.

These are recorded failed targets, not completed full-convergence results.
No NO/OpenDSS pass is claimed when no qualifying dispatch existed. The
accepted NO rule remains objective error<=0.1% AND original violation<=1e-5,
and the user target is the first qualifying iteration strictly below 100.
The C4 reference is a separately fully converged and OpenDSS-validated run;
its first trace crossing was 168, with the same frozen inputs and weights.

Full evidence: ddp-transfer-three-calibrations.tgz. It includes source
snapshots, all traces/statuses, hardware/affinity, resource samples, manifests,
gate results, raw failed-harness logs and orchestration. Serialized Julia
checkpoints remain on Kamiak. Archive SHA256:
`336b0a9c08110b1fc01634009ae5b69d8a12cfcc52d3304ce314565deb14a068`.

The separate balance-ratio 2 assessment uses frozen driver source c06bb64b.
Its common model and six shared numerical/validation helpers are byte-identical
to the passed gate by cmp; that model validation is explicitly reused. The
driver exposes the previously hardcoded balance threshold, retaining default 10.
Its completed large10k results also miss NO by 100:

- Job 31650253: initial rho=48, relative balance ratio 2; final original violation
  3.11863e-5. Rho changes after iterations 5/10/15/20/30/35 are
  48->96->192->96->48->96->48. The extra adaptation did not improve the tail.
- Job 31650276: initial rho=32, relative balance ratio 2; final original violation
  2.49727e-5. Rho temporarily drops to 32 late and returns to 64, with another
  transient feasibility bump. The plot shows the actual path.

**No tested 10k setting meets the below-100 target.** Fixed48 is the best of
these five capped candidates, not an established optimum. Objective accuracy
passes at each endpoint; original feasibility is the limiting criterion.
These failed targets are retained rather than treated as converged results.

## Current-code medium validation

Job 31650298 is the first medium assessment after the current-bound cleanup.
It uses root+ell screening, rho=3800 with the existing relative balance ratio 10
controller (7600 after iteration 10), MA57, identical C4 input and free terminal
SOC. The system-selectable orchestration is ec5e6c66, archived under campaign;
the numerical source remains c06bb64b and the passed model gate still applies.

First NO is **66**, matching the C4 trace crossing. Original violation there
is 8.119e-6. The actual serialized first NO dispatch passes independent OpenDSS:
maximum voltage difference 3.38e-10 p.u., zero voltage violations, and objective
agreement. It is hashed in no-dispatch.sha256 and retained on Kamiak. All 100
steps have validated stage statuses and source snapshots; no strict-convergence
or warmed timing claim is made for this assessment.

## Saved-point diagnosis and artifacts

Job 31650344 reconstructs all five 10k endpoints from their stage primals and
consensus energy, using the identical common-model variable map. Every
objective and original violation matches its recorded final trace row. This
is read-only evaluation with **zero additional optimization solves**.
The maximum violation equals the assembled battery dynamic-equation error in
all five cases. For best fixed rho=48, original violation is 2.32688e-5,
integrated-power SOC-bound violation3.15845e-6, and cone violation zero.
See [endpoint-components.toml](endpoint-components.toml) for all components.

- [Medium NO validation](med2522-no-validation.png): source cleanup preserves
  the 66-iteration crossing; the saved NO dispatch is OpenDSS-validated.
- [10k initial tuning](large10k-no-comparison.png): best fixed rho=48 and the
  initial rho=32/balance2 trial, alongside the fully converged C4 baseline.
- [10k larger penalties](large10k-high-penalty.png): the rejected larger-rho
  paths. Line styles identify settings; colors keep their established meanings.

The plotter verifies the exact reference-summary SHA and frozen input SHA.
The first two normalized panels jointly determine NO; the dual panel uses
the separate strict threshold. The vertical 100 line marks the requested
iteration target, and the C4 first crossing is shown separately.

Second full archive: ddp-transfer-balance-and-med-no.tgz, including both
trials with balance ratio 2, the medium NO validation, no-solve diagnostics, source
snapshots, orchestration, manifests and raw logs. Serialized Julia files
remain on Kamiak. SHA256:
`0cce4e9ac0143f3b0b31dae6c63a76b0d0475b43cbb2fb28ab2246a6d39436e3`.

The initial assessment and DDP-transfer audit are complete. The 10k below-100
target remains open; these results do not justify a finalized performance
claim or promoting another adaptive policy to default.
