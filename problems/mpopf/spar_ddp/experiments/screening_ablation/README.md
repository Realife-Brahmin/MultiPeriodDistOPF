# Optional constraint-screening factorial diagnostic

Authorized 20 October 2026. Final solver configuration, not historical versions.
Keep fixed-substation voltage cleanup enabled: its equality fixes voltage at
the upper bound, leaving no strictly feasible intersection with that bound.
This known formulation/presolve issue is not an optional timing arm.

Independent toggles: bit1 upper-voltage screening, bit2 squared-current
nonnegativity, bit4 substation import nonnegativity. Eight masks form the full
factorial. Non-root lower-voltage screening is never enabled.
Only analytically implied constraints are screened, not arbitrary limits.

All masks on periodic IEEE123 T3/two workers gate the full experiment.
Then three networks T96, two repeats/mask, reversed order on repeat2,
63 workers/coordinator on matching Gold6438Y+64CPU600GiB with occupancy guard.
Sequential warmed solves, max400 iterations and25min/arm;8h bounds campaign.
Preserve partial/failed outcomes; campaign.finished is not all-arms-validated.
No automatic expansion after timeout. Existing MSOPF work is untouched.

Every arm keeps identical physical power-flow initialization, affine line
search, exact curvature, tree backend, condensation, balanced scan, GC and
coefficients. Removing a redundant bound changes barrier trajectories:
identical iteration trajectories/counts across masks are NOT required.
Preserve all objectives, original primal/dual/complementarity, phases, setup
and warmup. An opt-in original-primal audit additionally checks every original
voltage bound, ell>=0 and P_Subs>=0 per iteration irrespective of screening;
retained boxes remain in normal diagnostics. Final original bounds and
independent OpenDSS gate every completed arm. Compare final objectives and
dispatches posthoc; convergence alone is not cross-arm equivalence evidence.

Full-convergence diagnostic first. Near-optimality requires original-feasibility
traces plus a matched centralized reference. Do not mix timer definitions with
paper results. No automatic paper edits or backend-specific memory claims.
The private depot can read the completed tree-campaign depot as secondary
cache; all new environment/cache writes go to this campaign's first depot.

October9 evidence grouped all optional rules, used diagonal curvature and
other horizons, and lacked warmed runtime comparisons. Do not promote its
iteration gains to current isolated-rule speedups. Summarize each rule via
matched on/off pairs conditional on the other two, baseline add-one and
full-screen leave-one-out comparisons. Report interactions and failures.
Two repeats provide initial repeatability, not a broad statistical conclusion.

Initial job31331921 stopped after the mask7 T3 solve passed convergence and
final bound checks but failed the required original-primal log gate. The
process-parallel pass bypassed the sequential diagnostic hook. The v2 source
installs the existing original-model hook on workers, returns one maximum
per time-step with phase-A scalars, and prints the global maximum each pass.
This is opt-in for every arm, changes no solver decisions, and is included
in every arm's measured overhead. Preserve v1 evidence, not a valid timing row.
