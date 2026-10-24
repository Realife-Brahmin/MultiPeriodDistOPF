# SPAR-DDP result collections

This directory is an evidence library, not generated scratch to clean away.
Start with the current reference and follow a study's own README/notes for
older measurements; not every folder represents comparable or current results.

- [proper_logging_campaign_oct16](proper_logging_campaign_oct16/README.md):
  validated 96-run matched campaign and current near-optimality protocol.
- [ma97_thread_sweep](ma97_thread_sweep/): separate centralized thread diagnostic.
- [balanced_stage_b](balanced_stage_b/), [stage_b_replay](stage_b_replay/):
  scheduling and selected Phase B diagnostics, not whole-solver benchmarks.
- [network_filterddp](network_filterddp/README.md): original network experiments;
  large generated `.jls` inputs/solutions are local and ignored.
- [centralized_ipopt](centralized_ipopt/README.md),
  [ieee123c](ieee123c/README.md), [reduced_space](reduced_space/README.md):
  earlier reference and formulation investigations.

Other named folders retain dated Hessian/KKT, screening, feasibility and
parallelization studies. Consult [the notes index](../notes/README.md) rather
than treating their rows as the main-paper comparison. Newer self-contained
studies keep their own outputs in [experiments](../experiments/README.md).

Preserve raw logs, failure records, CSVs, source/input hashes and frozen
harnesses. Recorded old paths are provenance, not live launch commands.
Do not overwrite a study to make a newer run appear to be the same experiment.
Binary captures, private HSL libraries, package depots and temporary exports
belong outside the tracked collection. Gzipped/archive evidence is intentional.
