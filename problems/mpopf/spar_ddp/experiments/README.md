# SPAR-DDP studies

These are self-contained research studies, not a second solver distribution.
Start from findings, then consult that study's harness and raw evidence.
Plans and launch scripts describe the recorded experiment; they are not
standing instructions to submit jobs. Production entry points remain in
[examples/power_system](../examples/power_system/README.md).

## Completed studies and current decisions

- [Tree/network-elimination ablation](tree_ablation/README.md): controlled
  network-backend comparison; [collected evidence](tree_ablation/kamiak_oct20/README.md).
  Distinguish full solves from frozen-matrix local linear-algebra replay.
- [Optional-constraint screening](screening_ablation/README.md): eight-mask
  factorial; [final findings](screening_ablation/results_oct20/FINDINGS.md).
  Conditional effects and interactions, not a universal ranking of constraints.
- [Nested battery counts](battery_count/README.md): 250/100 retained physical
  batteries, with matched references and repeated solves;
  [final findings](battery_count/results_oct21/FINDINGS.md). Changing the battery
  population changes the optimization problem; this is not a pure code ablation.
- [Initial value-representation inspection](value_representation/README.md)
  and [short replay](value_representation/replay/README.md): historical
  prerequisites, superseded in coverage by the subsequent audit.
- [Grouped full solves](value_group_solves/README.md): all batteries retained,
  selected cross-group interactions dropped; [findings](value_group_solves/results_oct22/FINDINGS.md).
  Observed iteration penalties motivated retaining coupling instead.
- [Coupling-retaining compact computation](value_coupling_audit/README.md):
  capture, frozen-direction checks, compression, compact composition and direct
  generation. [Final campaign decision](value_coupling_audit/results_campaign_oct23/FINDINGS.md):
  accurate on the tested four-period large10k point, but fresh-input cost gate
  not met. Revise downstream consumption before optimizer integration; no
  compressed optimizer has been tested. Earlier component successes are not
  a whole-solver speedup.

The matched main-paper evidence remains the separate
[October 16 campaign](../results/proper_logging_campaign_oct16/README.md).
Do not replace its rows with these diagnostics.

## Adding or maintaining a study

Use a named directory with a README describing the question, status, input
assumptions, production/experimental boundary and conclusion. Keep source,
validation/collection scripts, manifests and relevant failures together.
Put dated outputs below that directory; pin source/input/environment/hardware
and distinguish local component, Phase B and whole-optimization measurements.
New experiments must be opt-in. A plan is not authorization for a campaign.
Do not reshuffle frozen evidence or rewrite its historical paths; use its
recorded commit when reproducing it. Local scratch does not become a supported
feature merely by being merged.
