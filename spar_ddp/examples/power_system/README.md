# Power-system model, drivers and diagnostics

The model and shared numerical helpers are colocated here because they are
included directly by the drivers. The many older probes/sweeps are not all
supported launchers. Start with [the workstream guide](../../README.md), not
an arbitrary `run_*.sh` file.

## Current entry points

- `export_ieee123c_data.jl`: OpenDSS instance export (`envs/tadmm`).
- `ieee123c_filterddp.jl`: MPOPF model/generic solver driver (`envs/spar_ddp`).
- `centralized_ipopt_matched.jl`: matched centralized reference.
- `run_paired.sh`, `paired_main.jl`, `procs_main.jl`: explicit sequential,
  paired or one-node process-parallel launch paths.
- `opendss_check.jl`: independent dispatch validation.

Shared model/solver support includes `per_stage_model.jl`, `control_layout.jl`,
`loadflow_start.jl`, `terminal_soc_penalty.jl`, `voltage_screening.jl`,
`tree_kkt*.jl`, `pint_scan.jl`, and the paired/process pass implementations.
The `FILTERDDP_VALUE_GROUPING` and frozen-capture helpers are opt-in diagnostics,
not default compression. No model or algorithm changes are implied by this index.

## Research helpers

`analyze_*`, `benchmark_*`, `probe_*`, `capture_*` and most `run_*` sweeps
belong to particular dated investigations. Read the corresponding
[notes](../../notes/README.md), [results](../../results/README.md) or
[experiments](../../experiments/README.md) before using them. Many encode old
target tolerances, local paths or parameter combinations. Their presence is
for reproducibility and subsequent research, not a request to rerun them.

These files and frozen result harnesses retain their paths to avoid breaking
direct includes, manifests and evidence references. New standalone studies
belong in `experiments/<study>/`; do not add another unattended queue here.
