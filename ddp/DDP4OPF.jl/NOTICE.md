# Provenance and attribution

`DDP4OPF` is derived from **FilterDDP.jl** by **Mingda Xu**, released under the
MIT License. The original copyright notice and permission text are retained
verbatim in [LICENSE](LICENSE), as that licence requires.

| | |
|---|---|
| Upstream | https://github.com/mingu6/FilterDDP.jl |
| Forked from commit | `513a104` (2026-06-03) |
| Upstream version | `FilterDDP` v0.6.0, UUID `3fc7acb0-3ec8-433f-906a-2d1e86d217d2` |
| Method papers | Xu et al., arXiv [2504.08278](https://arxiv.org/abs/2504.08278) and [2606.01487](https://arxiv.org/abs/2606.01487) |

The filter line-search differential dynamic programming algorithm is the
authors' work. Please cite their papers when using this package.

## Why this is a separate package

FilterDDP is a general-purpose solver, developed and benchmarked on robotics
trajectory optimisation. Using it for multi-period optimal power flow on
distribution networks required changes deep in its backward pass, KKT assembly
and constraint handling, accumulated across many patches. Keeping those as a
patch series against a gitignored upstream clone stopped being reproducible:
on 2026-09-15 the documented patch recipe was found to fail on a clean checkout
of `513a104` (3 of 9 patches did not apply), and the working solver existed only
as uncommitted modifications on one machine.

This package is that working solver, committed as source.

It has a **new name and a new UUID** (`005a218a-47b4-46c4-a1f7-47229a6d1cc8`)
because FilterDDP is registered in the Julia General registry. Julia identifies
packages by UUID, so a modified copy that kept the upstream UUID would collide
with the registered package in any environment where both could appear.

## What changed relative to upstream

This source is the working tree the solver actually ran from, not a
reconstruction from patches. Net change against `513a104`: 12 files, 581
insertions and 175 deletions, concentrated in `backward_pass.jl` and `solver.jl`.

The changes were developed as the patch series in `ddp/patches/`, each with a
companion note in `ddp/notes/`. That series is kept as **history, not as a build
recipe** -- it does not apply cleanly to `513a104` -- so treat the list below as
the development record of what went in, not as a guarantee that each patch
applies independently:

- **Per-stage data** — stage-varying objectives and constraints
  (`per_stage_data`).
- **Dynamic network scaling** (`dynamic_network_scaling`).
- **Factor-backed bound sensitivities** and **factor-backed policy actions**
  (`factor_bound_sensitivities`, `factor_backed_policy_actions`).
- **Allocation-free update rule** (`no_copy_update_rule`).
- **In-place KKT right-hand sides**, **reusable RHS workspace** and
  **reusable stage-rule buffers** (`in_place_kkt_rhs`, `reuse_kkt_rhs_workspace`,
  `reuse_stage_rule_buffers`).
- **Blocked value right-hand sides** (`blocked_value_rhs`).
- **Typed constraint residual vector** (`type_constraint_residual_vector`).
- **In-place `B` assembly** and **active `B` rows** (`in_place_B_assembly`,
  `active_B_rows`).
- **Periodic capture instrumentation** — opt-in diagnostic dumps of backward-pass
  quantities (`periodic_capture_instrumentation`).

Beyond those, the only change is the rename: `module FilterDDP` became
`module DDP4OPF`. Upstream's `experiments/` (robotics benchmarks), `media/` and
README images were not carried over.

## Deliberately unchanged

- **Environment variables keep the `FILTERDDP_` prefix** (`FILTERDDP_MAX_ITERATIONS`,
  `FILTERDDP_QUIET`, `FILTERDDP_PERIODIC_*`, ...). They are an interface that
  existing sweep and queue scripts set, so renaming them is a separate decision,
  not a side effect of the fork.
- **The ASCII banner** printed at solve start still spells the upstream name.

## Behavioural equivalence

Verified at the fork: re-running the reference cases with this package produced
results identical to the upstream clone it was taken from. See the commit that
introduced this directory for the exact comparison.
