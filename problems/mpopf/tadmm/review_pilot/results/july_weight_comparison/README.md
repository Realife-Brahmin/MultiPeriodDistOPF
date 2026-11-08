# July-weight and local-cutoff assessment

July's coefficient has been restored: C_B=8.013987384856764e-8 at T=24,
recovered from release 50b8901c's original sampled-price rule. The current
periodic exports, corrected losses, MA57, initial SOC, zero terminal penalty,
loss-aware starts and proved bound cleanup are held fixed. The SPAR-DDP
defaults and previous DDP-weight evidence are unchanged.

This isolates the weight change and then the inner-iteration cap. It is not
a literal July replay: July's input sampling, rho controller, MUMPS/inner
tolerances and looser reporting rules differ. Each weight uses its own fully
solved, original-feasible, independently OpenDSS-validated centralized reference.
The current NO rule stays 0.1% objective agreement AND original feasibility 1e-5.

## Matched settings and counts

- med2522: initial rho 3800, relative controller, balance ratio 10; same setup
  as the prior post-screening DDP-weight assessment. The DDP C_B=0.0019 run
  reached NO at 66. With July's weight, exact locals first reach NO at 126;
  capped locals at 125. All three July-weight NO dispatch checks (exact,
  instrumented cap and lightweight cap) pass independent OpenDSS.
- large10k: fixed rho 48, matching the earlier DDP-weight fixed 48 pilot.
  Neither July-weight arm reaches NO within 200. Every local solve finishes
  in at most 24 Ipopt steps, so the 30-step cap never activates; the entire
  exact/capped numerical traces are identical. End original violation is
  1.417868831e-4, despite relative objective gap 1.05872e-7.

The initial rho/policy is deliberately held constant for this ablation;
these values are not established optimal settings for July's weaker curvature.
The older full DDP-weight large run with rho 48->96 reached NO at 168, but
that different policy is not used as a controlled weight-only iteration row.

The weight-comparison plot checks identical input SHA256, screening, initial
rho, rho policy, balance threshold and each reference identity. Objective gaps
are normalized against each problem's own centralized objective. The old
fixed 48 large pilot has only 100 steps; its curve is not extended or extrapolated.

## Cutoff implementation and observed timings

Exact locals use the current accurate settings with a 3000-iteration budget.
The cap arm changes only that budget to 30; other model/tolerance settings are
held fixed. ITERATION_LIMIT remains a provisional status, with finite-value
checks. It is never relabeled LOCALLY_SOLVED. The original global checker
still evaluates every iteration, and accepted NO dispatches get OpenDSS checks.

The medium cap fires 552 times across 4800 local solves,544 after outer
iteration 1. Exact locals need at most 37 steps. The median slowest native
Ipopt solve per round falls from 0.3306 s to 0.3070 s with the instrumented
cap, but full local-feasibility reports increase critical-worker framework
overhead from 0.0244 s to 0.0657 s. They mask the native-solver saving.

The explicit lightweight check disables those additional local reports,
marking them unevaluated/NaN. It retains finite checks, actual cutoff statuses,
all original global checks and independent OpenDSS. Its full 200-step numerical
trajectory matches the instrumented cap arm within the recorded tight tolerance;
its initial-start hashes also match. It again reaches NO at 125.

Observed medium timings, after three discarded compile steps and fresh resets:

- Exact: NO core/loop 49.53/86.47 s; 200-step core/loop 80.15/137.61 s.
- Instrumented cap: NO core/loop 54.43/87.87 s; 200-step86.46/139.66 s.
- Lightweight cap: NO core/loop 48.49/85.39 s; 200-step78.08/135.53 s.

These are single exploratory observations, not full-sequence-warmed repeat
medians. The overall gain is modest; native straggler reductions are clearer
than total-loop savings. The lightweight case's median slowest native solve
is 0.3051 s, with largest post-first-iteration slowest solve 0.3425 s versus
0.5648 s exact. Timing differences alone do not establish a universal speedup.

The large pair takes443.11/438.79 s for 200 steps (exact/capped), core
202.03/201.79 s. Because the cap never fires and traces match, that small
timing difference is not credited to capping. The redundant lightweight large
replay was cancelled during startup; its logs and explicit reason are retained.

NO counts are first crossings. Feasibility is intermittent afterwards: only
31/33 of the 200 medium exact/capped rows qualify, and the observed final
qualifying suffix starts at 200. No arm reaches strict global convergence
within this assessment. A saved valid NO point remains independently checked;
it is not evidence that every later iterate is feasible.

## Endpoint and local-accuracy audit

The large endpoint is reconstructed directly from saved stage primals and
consensus energies. Its objective and original violation match the final CSV.
The maximum is accumulated SOC-bound error 1.41787e-4; battery-dynamics
equation error is 5.78086e-5 and cone-inequality violation is zero.

The largest local-copy discrepancy is at period 6, battery 44082, copied time 7.
Two controlled next-step local probes use the same checkpoint center and
initial primal values: default settings versus explicit unscaled dual tolerance
1e-8. Both solve in 17 Ipopt iterations and return identical energy copies.
Local original constraint error is 2.45e-13, unscaled dual infeasibility
8.22e-14, complementarity 3.25e-10, and maximum raw multiplier 607. This tested
block does not reproduce the earlier inaccurate-local-termination defect.
The probe adds two local solves, not another ADMM trajectory or timing row.

Thus the present large bottleneck is storage consistency under the current
penalty settings, rather than a local solve exceeding30 steps. This does not
isolate the effect of line losses on conditioning; all arms retain corrected
losses. July's penalty schedule/reporting rules remain separate variables.

## Evidence and reproduction

- [Weight comparison](weight-comparison.png), with vector PDF alongside.
- [Medium convergence](med2522-analysis/med2522-july-convergence.png) and
  [medium stragglers](med2522-analysis/med2522-july-stragglers.png).
- [Large convergence](large10k-analysis/large10k-july-convergence.png) and
  [large stragglers](large10k-analysis/large10k-july-stragglers.png).
- Per-round CSVs and source-hashed summary JSON are in each analysis directory.
- [Endpoint/KKT probe](endpoint-audit.toml); raw optimizer logs in its archive.

All figures retain the manuscript's serif fonts, tADMM violet and system
colors. Raw archives preserve source snapshots, current/native Manifest
provenance, input hashes, hardware/affinity, stage statuses, physical traces,
OpenDSS logs, initial-start hashes, source gates and historical July sources.
Serialized checkpoints and actual dispatches remain on Kamiak.

Source c16d43cd: gate 31650587; medium pair 31650588; large pair 31650589.
Source 1d8ba823: gate 31650635; lightweight medium 31650636; redundant large
31650637 cancelled. Audit 31651004 uses the frozen c16d43cd model and the
standalone diagnostic committed in c9ce6e8a. All assessment statuses are
distinct from full-convergence success flags.

Archive SHA256 values:

- med2522-july-pair-31650588.tgz:
  3a1544fd2e3d81e2166e923307a2b0383773f1e42a956beb2e1ff9923a2dd503
- large10k-july-pair-31650589.tgz:
  2a384e92f30f649b944076eeb17702245797fcce2be39c1ea8704830bd26142a
- medium-july-light-31650636.tgz:
  11165f4bb263d23a77a4c3664d2cc90d3e266f38f5850127d29831530c8d9c74
- july-endpoint-audit-31651004.tgz:
  6a53f93ce0ae963c609cbbfc3980fcb5fabfb111c6cc2230124107315916ae25
