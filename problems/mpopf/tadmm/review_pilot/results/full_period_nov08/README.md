# Full period parallelism, exact warm starts and objective-update reuse

Native MA57 unit tests: 248/248 passed. Full T96 execution gate31653531 passed
with 96 workers on two Gold-6438Y+ hosts, 48 independently pinned workers per
host, one period per worker, plus coordinator. All three variants completed
768 exact local solves each; no iteration-limit exits. Source hashes and
actual host/CPU/period maps were independently rechecked after download.

The eight-iteration primal-only trajectories match exactly for scalar and
batched objective updates (objective, physical error, owned/all-copy residuals
and rho). The full-dual trajectory agrees within the declared numerical audit
tolerances. All three OpenDSS checks converge with zero voltage violations,
maximum voltage disagreement6.30e-11 pu. These are short execution checks;
their final unique dual RMS remains7.10 and they are not converged benchmarks.

The [independent gate audit](gate-summary.json) excludes iteration one from
local timing medians. Legacy/batched/full-dual median Ipopt iterations are
13/13/4; median stage times0.7239/0.6977/0.4138 s. Batched coefficient updates
take0.00543 s median versus0.01971 s for legacy updates. Full-dual transfer
takes0.11693 s versus0.02977 s for primal-only transfer; it is included in
stage/core timing. Modes execute sequentially with different JIT histories;
total gate time ratios are not reportable implementation speedups.

## Complete convergence comparisons now running

T48 jobs31653532 (primal/batched) and31653533 (primal-dual/batched) each use
48 workers on one Gold node. Both use rho120000, full warmup/reset, MA57,
exact local solves and the existing final residual/dual and cost assessments.
No inexact solve policy, stopping tolerance or physical model was changed.
Both fresh centralized references are retained. Node sampling runs on reserved
CPUs; the coordinator is pinned before centralized timing as well as ADMM.

Follow-through job31653565 waits for both T48 jobs. It records both outcomes
and selects the fastest validated implementation at the declared final
owned/all-copy RMS<=1e-5, unique dual<=1e-2 and objective gap<=0.1% goal.
Original-model errors remain separate flags. It then submits the two requested
T96 starts60000/240000 with 96 workers, without changing rho control. If no
validated implementation meets that goal, it records the issue instead of
submitting an unvalidated follow-through trial. Both implementation trials and
their warmups/reference solves count toward exploration cost. Passive energy
audits run after the corresponding measured jobs.

## Provenance

Gate numerical/launcher source: d034d33da55085814995fd3adccad09cbf61dfa2.
Measured-job source:3daf5cf5d226ce8b96e1ce3c966e3cf4b1fce5f0; the post-gate
patch only isolates the coordinator from the sampler before central timing
and adds analysis/docs. Model, numerical driver and actor hashes are unchanged.
Slurm-only dependency is pinned in a separate colocated environment; numerical
Project/Manifest versions remain those of the preceding sweep.

Archive full-period-gate-31653531.tgz SHA256:
a057763ca174af7cba6c1368f63604a378e8bb1a45c76be4dbc192cd7d4b1a35.
It includes successful and failed setup gates, logs, source snapshots, hashes,
actual layouts and scheduler accounting. Serialized large arrays remain on
Kamiak. Setup failures are preserved and excluded from converged solve claims.
The complete submission/selection/passive-audit protocol is in
full-period-submission-protocol.tgz, SHA256
5d900b0fa8600f29b8eb75dd4d851678429c466c6c4fef77be74eb34c1a15b96.
The selector records scheduler CPU/resource accounting alongside both outcomes.
The obsolete 24-worker T96 partial traces remain in the preceding campaign's
horizon-T96-superseded.tgz; no incomplete case is relabeled converged.
