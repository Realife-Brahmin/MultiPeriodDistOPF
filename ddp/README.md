# DDP: differential dynamic programming for multi-period OPF

This workstream asks whether second-order **differential dynamic programming (DDP)**
can solve full-network multi-period OPF, and how it compares with solving the same
problem centrally with Ipopt.

- **The solver is [DDP4OPF](DDP4OPF.jl/)**, our MIT-licensed fork of
  [FilterDDP.jl](https://github.com/mingu6/FilterDDP.jl) by Mingda Xu: filter
  line-search DDP with a sparse KKT solve at each stage.
  [NOTICE.md](DDP4OPF.jl/NOTICE.md) records what changed. Please cite the original
  method papers (arXiv 2504.08278 and 2606.01487) when you use it.
- **The model** is the branch-flow SOCP of the tADMM work. Each period's network
  variables are DDP controls, and the battery energies are the only states carried
  from one period to the next.
- **Every result is checked against a centralized JuMP/Ipopt solve of the identical
  exported instance** before any timing is compared.

A paper on this work is in preparation.

## Setup

Two Julia environments are involved. `envs/tadmm` exports feeder instances from
OpenDSS; `envs/ddp2026` holds DDP4OPF, JuMP and Ipopt
([its README](../envs/ddp2026/README.md)). From the repository root:

```bash
julia --project=envs/tadmm   -e 'using Pkg; Pkg.instantiate()'
julia --project=envs/ddp2026 -e 'using Pkg; Pkg.instantiate()'
```

## Run one matched case

ieee123 with a 6-period horizon takes about two minutes, including compilation.
Run from the repository root:

```bash
# 1. Export the instance to ddp/results/network_filterddp/ (gitignored)
PROFILE_PERIODIC=1 julia --project=envs/tadmm \
    ddp/examples/power_system/export_ieee123c_data.jl ieee123C_1ph 6

# 2. Settings both solvers must share
export REDUCED_PROFILE=periodic REDUCED_CB=1e-3 TERMINAL_SOC_SOFT=1

# 3. Centralized reference; prints one CENTRAL_IPOPT line
julia --project=envs/ddp2026 \
    ddp/examples/power_system/centralized_ipopt_matched.jl ieee123C_1ph 6

# 4. DDP4OPF with the diagonal stage Hessian and full per-iteration logging
FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8 \
FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 \
julia --project=envs/ddp2026 \
    ddp/examples/power_system/ieee123c_filterddp.jl ieee123C_1ph 6 solve
```

**Expected output.** Ipopt reports `objective=3143.1159...` after 37 iterations.
DDP4OPF prints `FilterDDP objective=3143.1160...` after 75 iterations, a relative
difference of about `4e-8`.

Compare only against the `CENTRAL_IPOPT` line. If tADMM's `run_bf.jl` has been run
on your machine, the DDP driver also prints an `Ipopt objective=... objective_gap=...`
line. That line compares with a tADMM solution of a *different* problem; ignore it.

**Other systems and horizons.** Replace `ieee123C_1ph 6` with `ieee2522C_1ph` or
`large10kC_1ph` and any horizon `T`. Large cases take hours. For large10k at
`T >= 12`, also set `FILTERDDP_FACTOR_BACKED_POLICY=1`. It stores a compact
per-stage policy instead of dense response maps, which gives identical iterates in
far less memory.

**The full protocol** (paired runs, near-optimality detection, background-load
sampling, summaries) is automated in
[`run_matched_ipopt_race.sh`](examples/power_system/run_matched_ipopt_race.sh).
It **commits and pushes** its results after every case, so run it only on a branch
meant to receive them.

## Rules for reportable numbers

A comparison is only reported if all of these hold:

- **One problem.** Both solvers read the same exported instance, with `C_B = 1e-3`
  and the periodic price profile, and their objectives agree before any time is
  quoted.
- **Soft terminal state of charge.** The objective carries
  `gamma * sum_j (B_j^T - B_j^0)^2`, with `gamma` fixed per system and independent
  of `T`: ieee123 1.75e4, ieee2522 3.29e4, large10k 420. The single source is
  [`terminal_soc_penalty.jl`](examples/power_system/terminal_soc_penalty.jl).
- **Time to near-optimality** is the first iteration with objective within 0.5% of
  Ipopt **and** primal infeasibility below the system's threshold: ieee123 `1e-6`,
  ieee2522 `1e-5`, large10k `1e-4` (`NEAR_OPT_PRIMAL_BY_SYSTEM` in the same file).
- **Full logging.** Every timed run sets `FILTERDDP_TIMING_DIAGNOSTIC=1` and
  `FILTERDDP_FEASIBILITY_DIAGNOSTIC=1`. Without the per-iteration trace, the
  near-optimality point cannot be reconstructed afterwards.

The centralized sweep in `results/centralized_ipopt/` predates these rules and uses
a different instance family, so it is not comparable with anything above.

## Layout

```text
ddp/
  DDP4OPF.jl/                 the solver
  examples/power_system/      drivers: instance export, DDP4OPF network model,
                              centralized Ipopt, sweeps
  results/                    one folder per study, each with its raw logs
  notes/                      write-ups of findings and methods
  resources/                  reference papers (MANIFEST.txt; PDFs not committed)
  logs/                       terminal output of the original reproducibility stages
  patches/                    development history of the fork; not a build recipe
  README_FILTERDDP_EXPERIMENT.md  the original reproducibility study
```

The main studies under `results/`:

| Folder | Contents |
|---|---|
| [`matched_ipopt_race/`](results/matched_ipopt_race/) | DDP4OPF against Ipopt on identical instances, all three systems |
| [`centralized_ipopt_matched_knee/`](results/centralized_ipopt_matched_knee/) | Ipopt at growing horizons, up to its memory limit |
| [`hessian_rewrites/`](results/hessian_rewrites/) | exact rewrites of Hessian and KKT assembly, with timings |
| [`diag_hessian_horizon/`](results/diag_hessian_horizon/) | diagonal against exact stage Hessian across horizons (before the soft terminal SOC) |
| [`reduced_space/`](results/reduced_space/) | the network eliminated by an inner OPF |
| [`copper_plate/`](results/copper_plate/), [`official_example/`](results/official_example/) | the original validation stages |

Good entry points in `notes/`:
[`FILTERDDP_DIAGONAL_HESSIAN.md`](notes/FILTERDDP_DIAGONAL_HESSIAN.md),
[`FILTERDDP_HESSIAN_EXACT_REWRITES.md`](notes/FILTERDDP_HESSIAN_EXACT_REWRITES.md),
[`CENTRALIZED_IPOPT_KNEE_SWEEP.md`](notes/CENTRALIZED_IPOPT_KNEE_SWEEP.md),
[`REDUCED_SPACE_INNER_OPF_FEASIBILITY.md`](notes/REDUCED_SPACE_INNER_OPF_FEASIBILITY.md),
and, for how the solver maps to the papers,
[`PAPER_CODE_MAP.md`](notes/PAPER_CODE_MAP.md) and
[`FILTERDDP_API.md`](notes/FILTERDDP_API.md).

## History

Work began in August 2026 as a reproducibility study of the FilterDDP authors'
code on a copper-plate battery problem, validated against three independent
references. That record, stages 1 to 9, is
[`README_FILTERDDP_EXPERIMENT.md`](README_FILTERDDP_EXPERIMENT.md). The network
model, the fork into DDP4OPF, and the matched Ipopt comparisons followed.

Before any of this, a **first-order** DDP scheme for copper-plate MPOPF was written
independently in 2025. It passes the dynamics dual backward one stage per sweep,
with no curvature across stages. It is kept as a read-only reference in
[`envs/ddp/`](../envs/ddp/); a minimal restatement is
[`examples/power_system/user_ddp.jl`](examples/power_system/user_ddp.jl).
