# Local validation — 20 October 2026

**Kamiak status: not submitted.** This cloud session has no configured Kamiak
SSH connection or Slurm client. No cluster jobs were inspected, cancelled,
changed or launched. Claude's MSOPF simulations were not touched.

Julia 1.12.4 was installed in the transient cloud workspace and dependencies
were installed in a private depot. Package/library manifests are archived.
OpenDSSDirect's build needed `TAR_OPTIONS=--no-same-owner` on this root-owned
cloud filesystem; this is an environment workaround, not a solver change.

## Passed

- Nine synthetic tests: branching network, parent-and-child batteries,
  off-diagonal future curvature, full Newton correction against independent
  UMFPACK, selected battery responses, feedback and value-curvature agreement,
  invalid backend rejection, and capture serialization.
- Frozen synthetic matrix replay: tree, generic condensed UMFPACK, and full
  KKT UMFPACK; three interleaved repetitions each, all numerical gates passed.
  Do not interpret these tiny-fixture times as network performance.
- Full process-parallel IEEE123 smoke: T=3, two workers plus coordinator,
  51 batteries, 791 controls, 562 equalities, balanced scan tasks, exact
  curvature, full warm-up, and full convergence in both arms.
- Both arms converged in 16 iterations, with identical printed objective.
  Maximum relative differences: battery-state trajectory 2.734e-16;
  full time-step optimization-variable trajectory 7.639e-17.
- Both OpenDSS replays converged, with zero lower/upper voltage violations.
- The final end-to-end shell harness produced its `pair.complete` marker;
  neither final arm reported a warm-up/timed-stream mismatch.
- Bash syntax, Python compilation, numerical/log validation, incomplete-run
  rejection, and `git diff --check` passed. OpenDSS gates rejected a
  nonconverged report, voltage violations, and a truncated report.
- The submission script failed safely outside Slurm, before creating a campaign.

This uses the committed **unprofiled T3 IEEE123 export** for local correctness,
not the current T96 periodic campaign data. Its SHA256 is preserved in
`input.sha256`; no original input was modified. Per-system battery cost and
soft terminal penalty were selected by the normal driver.

## Preserved evidence

`local_smoke_evidence.tar.gz` contains the final paired raw solve logs, saved
solutions, independent OpenDSS logs, validation/completion markers, Julia
manifests and synthetic replay CSV/log. `source.sha256` identifies the tested
source files. Source base before this ablation: commit
`789bf3d165745f67d128ff33ab54411b9667aa26`.

The first exploratory correctness pair also agreed numerically, but had a
warm-up output-stream mismatch; none of its times are promoted. The archived
final harness run clears inherited experimental/color settings and rejects
that timing problem. Cloud elapsed times remain smoke evidence only, not a
substitute for the exclusive-node Kamiak comparison.

## Next execution step

Use a separate DDP checkout/worktree of `ddp/oct20-tree-ablation` on the machine
with Kamiak access. Do not switch Claude's MSOPF checkout or environment.
Follow `../README.md` to submit the isolated campaign with the current periodic
exports. The launcher requests one exclusive node and runs all paired arms
sequentially. The original paper results and solver defaults remain unchanged.
