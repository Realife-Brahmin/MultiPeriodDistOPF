# October 18, 2026 repository cleanup

This cleanup changes organization and guidance, not solver equations,
benchmark numbers, raw feeder data or published evidence.

## Path map

- `ddp/` → `problems/mpopf/spar_ddp/`: active solver, drivers, notes and evidence.
- `envs/ddp2026/` → `problems/mpopf/spar_ddp/`: active solver environment.
- `envs/ddp/` → `archives/first_order_ddp/original/`: retired first-order code.
- First-order `user_ddp.jl`, `user_ddp_t3.jl` and
  `first_order_dual_average_experiment.jl` → archive `examples/`.
- Root `plots/` → `problems/mpopf/tadmm/plots/`: historical tADMM tools.
- `scripts/` → `scratch/`, except the maintained resource fetcher,
  now `tools/fetch_resources.sh`.
- `assets/networks/` → `figures/networks/`: editable SVGs and PDF exports.
- Tracked editor workspaces and the always-enabled old Cursor sweep rule
  were removed. Current shared instructions live in AGENTS.md.

## Evidence and executable paths

Raw results, logs, hashes, manifests, frozen campaign harnesses and network
figure bytes are preserved. Historical records retain their old paths and
labels. For exact reproduction, use their recorded commit and environment;
do not run an old frozen launcher against renamed folders without adapting
a separate working copy. Historical requests do not authorize new runs.

The current paired and process wrapper copies live under
`problems/mpopf/spar_ddp/examples/power_system/`; the original copies remain in results
for provenance. Active driver paths and the local package source path were
updated. Module/package names, environment variable names and solver modes
called `ddp` are unchanged.

The medium case is displayed as med2522; raw ieee2522C_1ph filenames and
recorded historical labels are retained. The paper repository is separate
and was not edited by this cleanup.

## Maintenance

Review scratch/ at milestones. Promote useful maintained tools, preserve
provenance needed to interpret evidence, and remove abandoned scratch after
checking references. The tADMM plots/results are historical material; fresh
benchmarking is deferred. Keep personal workspace settings out of Git.
