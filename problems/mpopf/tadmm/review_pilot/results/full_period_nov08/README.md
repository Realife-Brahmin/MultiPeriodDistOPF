# Full period parallelism, exact warm starts and objective-update reuse

Native MA57 unit tests: 248/248 passed. Full T96 execution gate31653531 passed
with 96 workers on two Gold-6438Y+ hosts, 48 independently pinned workers per
host, one period per worker, plus coordinator. All three variants completed
768 exact local solves each; no iteration-limit exits. Source hashes and
actual host/CPU/period maps were independently rechecked after download.

The selection job submitted T96 successfully, then its passive-audit stage
failed because a helper had not been copied into the frozen source bundle.
Audit repair31656659 completed separately after all solves; no numerical source
or timed trajectory changed. The exact helper, wrapper, hashes, logs and energy
results are preserved in full-period-audit-repair.tgz, SHA256
564b896887912fc29e710de4a8a4682f7f104b2e6704e3e6aea34d7cdd68b0c0.
Completed-run archive full-period-completed-runs.tgz SHA256:
0a2aa5c5536e748a1dee9f388647aa5fc4b296d7e64b5baa0cfbcc2406a427a3.

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

## Completed comparisons

All measured numerical runs and saved-point network checks completed. Native
validation was independently rerun locally: source hashes, fresh warmup/reset
paths, worker coverage, local statuses, stopping flags and OpenDSS records pass.

- T48 primal-only:322 iterations,393.7843 s core,432.8034 s loop.
- T48 primal-dual:322 iterations,334.2166 s core,373.6063 s loop. Full warm
  starts save15.1% against the48-worker primal-only implementation. Combined
  with full parallelism/batched updates this is1.96x faster than the preceding
  24-worker653.9158-s case. Fresh matched central solves take88.2794/87.0827 s;
  tADMM is still slower. Both final RMS/dual goals are met, objective gaps
  below3.41e-7, and integrated energy overshoot0.35764 kWh /0.01867 SOC points.
- T96 rho60000:1,000-iteration budget reached without the final RMS/dual goal,
  core1586.7922 s. Objective gap0.0644% does not make this a converged solution:
  integrated energy error24.4783 kWh /1.5920 SOC points, dual RMS0.1660.
- T96 rho240000:849 iterations,1073.4769 s core,1228.1118 s loop, final RMS/dual
  goal met. Objective gap7.83e-7; integrated energy overshoot0.68170 kWh /
  0.04666 SOC points. Matched full central solve184.8716 s. This is5.81 times
  the centralized time, not a speedup.

The small engineering errors at the RMS-goal points are explicitly retained;
they do not satisfy the separate original-model1e-5 physical NO threshold.
All39 saved-point/endpoint OpenDSS records across four arms converge with zero
voltage violations and maximum voltage disagreement<=1.50e-10 pu, including
the budget-limited arm. Network agreement does not establish battery consensus.

Slow local solves remain: T48 full-dual median4 inner iterations but max144,
versus primal-only median13/max28. T96 maxima158/112 for the two rho starts.
The total-time comparison includes these outliers and multiplier-transfer cost.
The convergence figure shows persistent oscillation, especially at T96; reducing
local solve time has not fixed the outer controller's horizon dependence.

Measured T48 aggregate PSS is53.79/53.47 GiB, central5.18/5.03 GiB. T96 sums of
independent per-node peak PSS are131.77/133.45 GiB, a conservative upper bound
rather than a synchronized global peak; central9.80 GiB. No aggregate memory
advantage is claimed. Complete workflow times including warmups/build/reference
and validation are1501.3/1373.8 s (T48) and7717.8 s (both T96 starts).

See [convergence plot](analysis-completed/implementation-comparison.png),
[PDF](analysis-completed/implementation-comparison.pdf), [all metrics](analysis-completed/summary.json),
and [T96 SOC plot](T96-SOC.png). Full-worker outcomes are exploratory single
measured trajectories, not repeat medians or a validated universal rho heuristic.

## Execution history

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
