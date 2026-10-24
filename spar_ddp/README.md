# SPAR-DDP for multi-period distribution OPF

SPAR-DDP combines sparse radial-network elimination with parallel
battery-value recursion. Each time-step keeps network variables local;
battery energy is the state linking time-steps.

[DDP4OPF](DDP4OPF.jl/) is our fork of
[FilterDDP.jl](https://github.com/mingu6/FilterDDP.jl), with the upstream filter
line-search foundation and MPOPF-specific sparse/network and parallel paths.
Retain its [NOTICE](DDP4OPF.jl/NOTICE.md), MIT license and original method
citations (arXiv 2504.08278 and 2606.01487).

## What belongs here

This is the SPAR-DDP **workstream**, not a scratch directory or a Julia
environment. Keep reusable code and the evidence needed to understand and
reproduce this project, including unsuccessful but informative experiments.

- `DDP4OPF.jl/`: maintained solver package and upstream notices.
- [examples/power_system](examples/power_system/README.md): model, drivers,
  validation and shared numerical helpers; older diagnostics are indexed there.
- [experiments](experiments/README.md): bounded opt-in studies, each with
  scope, source/harness, checks, results and a conclusion. Not solver defaults.
- [results](results/README.md): matched benchmark evidence and dated diagnostic
  collections. Logs and manifests are scientific provenance, not loose scratch.
- [notes](notes/README.md): derivations and implementation rationale, including
  dated notes whose settings are no longer current.
- `logs/`, [patches](patches/README.md) and `README_FILTERDDP_EXPERIMENT.md`:
  preserved upstream-evaluation/development history, not a second build recipe.
- `resources/`: reference-source manifest; fetched papers remain local.

Dependencies belong in `../envs/spar_ddp`, common feeder decks in `../rawData`,
and unpromoted one-offs in `../scratch`. Keep binary Newton captures, private
libraries, depots and checkout bundles out of the tracked tree. Do not remove
evidence merely because an approach did not advance; label its status instead.

## Environments and entry points

From the repository root, using Julia 1.12 and Bash:

```bash
julia --project=envs/tadmm -e 'using Pkg; Pkg.instantiate()'
julia --project=envs/spar_ddp -e 'using Pkg; Pkg.instantiate()'
```

The exporter uses OpenDSS through `envs/tadmm`. Solver and centralized
reference drivers use `envs/spar_ddp`. For a small instance export:

```bash
PROFILE_PERIODIC=1 julia --project=envs/tadmm \
  spar_ddp/examples/power_system/export_ieee123c_data.jl ieee123C_1ph 6
```

Exported `.jls` files are local generated inputs under
`results/network_filterddp/`; they are not included in a fresh clone.
For the centralized reference on that same instance:

```bash
REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 \
  julia --project=envs/spar_ddp \
  spar_ddp/examples/power_system/centralized_ipopt_matched.jl ieee123C_1ph 6
```

The current paired/process drivers now live beside the model:

- [run_paired.sh](examples/power_system/run_paired.sh): explicit sequential,
  threaded diagnostic, or process-parallel mode; read its settings first.
- [procs_main.jl](examples/power_system/procs_main.jl): process workers and
  the shared-memory value pass.
- [paired_main.jl](examples/power_system/paired_main.jl): paired-pass wrapper.
- [ieee123c_filterddp.jl](examples/power_system/ieee123c_filterddp.jl): model
  and generic solver entry point.

The launcher usage is:

```bash
bash spar_ddp/examples/power_system/run_paired.sh <raw-system-ID> <T> <count|1|procs> [tag]
```

This runs optimization and a final OpenDSS check. It is not an instruction
to launch a new campaign. Large cases require deliberate worker/memory
configuration; the one-node shared-memory path is not a multi-node backend.
The launcher is an entry point, not a promise to reproduce the paper's
hardware, full campaign settings or timings. Use the frozen campaign
configuration at its evidence commit for that purpose.

## Current benchmark protocol

Main-paper **time to near-optimality** requires both:

- relative objective agreement within **0.1%** of the matched reference;
- maximum **original-model** constraint violation at most **1e-6** for
  ieee123, **1e-5** for med2522 and large10k.

This is a reporting target, not the solver convergence test. The completed
campaign runs continue to declared full convergence; crossing times are
scored posthoc. The reference is the final objective of the matched MA97
first repeat. Use two-repeat medians and the faster centralized MA57/MA97
median per case. Keep internal slack residuals distinct from original-model
violations, and optimizer time distinct from setup, warm-up and replay.

Both paths must use identical exported inputs, periodic profiles, weights,
initial point, screened constraints, hardware/thread policy and timing
boundaries. `REDUCED_CB=system` selects C_B = 1e-3, 1.9e-3 and 2.4e-5
for ieee123, med2522 and large10k respectively. `TERMINAL_SOC_SOFT=1`
selects soft terminal-energy penalties with coefficients 1.75e4, 3.29e4
and 420, fixed per system across horizons. See
[terminal_soc_penalty.jl](examples/power_system/terminal_soc_penalty.jl).

Require full timing/feasibility logs, validated completion/status and final
OpenDSS checks. Historical 0.5% targets and strict-convergence clocks remain
separate. Some old launchers still use those older settings; they do not
define the current paper criterion. No solver defaults changed in the
October 18 folder/documentation cleanup.

## Evidence and organization

- [proper_logging_campaign_oct16](results/proper_logging_campaign_oct16/README.md):
  completed 96-run matched campaign; three systems, T=24/96/144/288, two
  repeats each of MA57, MA97, owner DDP and balanced DDP.
- [ma97_thread_sweep](results/ma97_thread_sweep/): separate large10k T96
  thread diagnostic; one observation per setting, not campaign medians.
- [balanced_stage_b](results/balanced_stage_b/): value-pass assignment
  diagnostics, distinct from full optimizer timings.
- `examples/power_system/`: model, exporter, reference solvers and drivers.
- `DDP4OPF.jl/`: solver package.
- [Experiments index](experiments/README.md): tree/screening ablations, nested
  battery-count study and value-function investigations, with final findings.
- [Compact campaign conclusion](experiments/value_coupling_audit/results_campaign_oct23/FINDINGS.md):
  direct generation and four-stage compact consumption passed accuracy, but
  fresh-input cost did not justify optimizer integration. No compressed
  optimizer or general compact-horizon speedup has been demonstrated.
- `notes/`, `results/`, `logs/`, `patches/`: dated research records.
- [README_FILTERDDP_EXPERIMENT.md](README_FILTERDDP_EXPERIMENT.md): original
  reproducibility study; historical settings are not current benchmark rules.

The original first-order DDP and its restatement are retired in
[archives/first_order_ddp](../archives/first_order_ddp/README.md).
The [path map](../docs/repository-layout-oct18.md) explains old names preserved
in notes, raw evidence and frozen harnesses. Those records must not be
rewritten to make historical commands appear current.
