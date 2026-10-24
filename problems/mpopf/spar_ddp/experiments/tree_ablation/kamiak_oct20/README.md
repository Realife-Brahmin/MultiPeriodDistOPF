# Kamiak tree-only ablation, 20 October 2026

Completed job 31329692, cn165, Gold 6438Y+, 64 allocated CPUs,
63 Julia workers plus coordinator, 600 GiB request, Julia 1.12.4.
Source: 9b75ae759868e070eb6244641bb4076aa1d3ea02, derived from PR #171.
The node-occupancy guard passed; coordinator affinity covered CPUs 0-63.
Elapsed allocation time: 2:53:51, exit 0. No paper files were changed.

Three alternating-order repeats per backend, T=96, full convergence.
Median optimizer seconds (tree / generic condensed UMFPACK):

- IEEE123: 0.887 / 1.155, a 1.30x ratio.
- med2522: 11.477 / 79.613, a 6.94x ratio.
- large10k: 49.340 / 536.241, a 10.87x ratio.

This isolates network elimination while retaining battery condensation,
exact curvature, balanced temporal scan and matched solver policies.
It is not a comparison against centralized IPOPT, nor a near-optimality
benchmark, nor evidence that this UMFPACK baseline is optimally tuned.
Large10k Phase A median was 4.747 versus 479.179 s; Phase B remained
35.547 versus 33.497 s. The benefit is concentrated in network elimination.

All nine numerical pairs and the periodic T3 correctness pair passed.
Both arms had identical iteration counts: 24, 66 and 16 for the three T96
systems. Maximum relative trajectory disagreements were below 1.5e-13.
All 20 independent OpenDSS replays passed with zero voltage violations.
Final primal/dual/complementarity gates passed at each system's tolerance.
Known OpenDSS precompilation-overwrite diagnostics are preserved; runtime
replay succeeded. They are not being silently treated as solver failures.

`evidence.tgz` contains raw logs, per-run timings, final pair validation,
completion markers, input/source hashes, dependency manifests, hardware
and allocation records. Saved solutions remain on Kamiak in
`/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2`.
`summary.csv` is the generated median comparison.

Slurm batch MaxRSS was 566113772 KiB (~539.9 GiB); this is an allocation-wide
accounting high-water mark, not backend-specific peak memory. Coordinator
RSS alone is not total worker memory. No backend memory-saving claim is made.

Original job 31329138 failed before any benchmark because the snapshot
omitted the shared root-level DSS profiles. `failed_setup.tgz` preserves
that attempt. The new snapshot includes those files, replaces forbidden
`--exclusive` with a full-node CPU request and occupancy guard, and generates
only the missing periodic T3 correctness input privately. The three T96
exports were copied unchanged from the frozen October 16 campaign.
