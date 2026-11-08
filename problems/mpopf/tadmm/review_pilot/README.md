# Corrected tADMM reviewer-response experiments â€” November 7, 2026

**Current user pivot:** new tADMM comparisons restore July's battery coefficient,
8.013987384856764e-8 at T24, and test exact local solves against a 30-step cap.
See [JULY_COMPARISON.md](JULY_COMPARISON.md). Prior DDP-weight results below
remain historical evidence. Current periodic profiles and corrected losses
are retained; this is a controlled weight/cutoff comparison, not a full replay
of July's inputs and stopping rules. SPAR-DDP defaults are unchanged.

**November 8:** the selected owned-period RMS stop at 1e-4 is now logged
explicitly. The paper's high-rho start (15000 on large T24), primal-first phase
and original 15-step watchdog give a historical exit at 44 iterations / 39.64 s
core, versus 48 / 51.09 s fixed-rho. Endpoint violation improves 9.85x, but the
0.2703% cost gap and 8.67e-4 original violation do not meet the separate physical
NO rule. This is not a matched-quality speedup claim. See the
[rho/diagnostic assessment](results/rho_diagnostics_nov08/README.md) and
[period distributions / overhead](results/kamiak_execution_insights/README.md).
Core solve time excludes independent physical validation; total loop wall time
is reported separately. Typed diagnostics keep every check and are opt-in via
TADMM_DIAGNOSTICS=typed; SPAR-DDP solving/defaults are unchanged.

The first matched July-weight assessment is complete. At T24, med2522 first
passes NO at 126 with exact locals and 125 with the 30-step cap; both actual
NO dispatches pass OpenDSS. The lightweight capped run preserves the full
numerical trajectory and has a modest observed timing gain. Large10k reaches
no NO point by 200; its maximum local iteration count is 24, so the cap never
fires. See [results and plots](results/july_weight_comparison/README.md).

**Numerical-modeling update:** the T=24 plateau audit found redundant fixed-root
voltage bounds creating opposing ~2.66e9 multipliers and inaccurate local
termination under the scaled tolerance. The common model now removes those
redundant bounds and explicitly limits unscaled complementarity. See
`PLATEAU_AUDIT.md` for controlled same-rho evidence and the fresh-run validation
requirement. Prior frozen campaigns are preserved under their original sources.

**DDP transfer and NO target:** see [DDP_TRANSFER_AUDIT.md](DDP_TRANSFER_AUDIT.md)
for the source-mapped review of the Augustâ€“October work. The repaired C4
baseline first reaches the current NO rule at **66 iterations on med 2522** and
**168 on large 10k**. The user's target is below 100, separate from strict
stopping. Redundant current lower bounds are now removed using DDP's proof;
data-only voltage-upper screening follows the October case-specific evidence.
New runs save and OpenDSS-check the actual first NO dispatch.

**Post-transfer assessment completed:** med 2522 job 31650298 again reaches NO
at 66, now with the actual first NO dispatch independently passing OpenDSS.
Five capped 10k penalty candidates all miss NO by 100. Fixed rho=48 is best
among those tested, ending at 2.327e-5 original violation versus the 1e-5 limit.
Saved-point reconstruction identifies battery-dynamics inconsistency as the
maximum violation. See `results/ddp_transfer_NO/README.md` for raw evidence,
the plots in the manuscript's figure style, and the remaining tuning limit.

**Fresh validation passed:** med 2522 T=24 job 31645746 converged in 389 warm-up
and 389 measured iterations with matching traces and independent OpenDSS checks.
The original failed run stopped at 5000. Both use rho=3800, increased to 7600
after iteration 10; no tuning change accounts for this result. The corrected
large 10k T=24 follow-up job 31645790 also passed, with 267 strict iterations.
See results/accuracy_fix_T24_31645746.

Current assessment figures follow the manuscript's actual PGFPlots styles,
`figures/traj_common.tex` and `figures/scaling_inflation_common.tex`, with the
IEEE serif typography from `main.tex`. `plot_style.py` uses Times New Roman
with STIX math, tADMM violet #6A4C9C, feasible centralized-reference green
#2E8B57, and infeasible-reference maroon #8C1D1D. System titles use the paper's
blue/gold/green system colors. Every diagnostic curve in the current NO plots
is tADMM, so it is violet; line styles distinguish settings, and axis labels
identify metrics. Gray lines mark thresholds. PNG and vector PDF versions
are generated together. Original archived figure styles remain preserved.

## Current first-assessment scope (user update)

Use T>=24 to assess performance and penalty policies. T=6 evidence is retained
as correctness/debugging material, not an adequate horizon for those conclusions.
T=12 is a diagnostic fallback only. The user requested one selected backend and
one measured run per case initially, rather than a two-backend/two-repeat sweep.
Compilation warm-up and the measured restart remain distinct; repeated benchmark
trials are deferred until needed for a later finalized comparison.

The user has selected **MA57 only for future tADMM workstream runs**, including
their centralized reference. The first-pass launcher defaults to MA57 and
rejects another backend. Older comparison/diagnostic scripts remain historical
tools and do not authorize further MA97 or MUMPS runs. Available evidence favors
MA57 in these cases; this is a workstream choice, not a universal solver ranking.
The previously active med 2522 T=24 MA57 case in task
31643458_6 is retained. Duplicate MA97 task 31643458_7 was canceled, and pending
tasks 8/9 (large 10k T=6) and 10/11 (old paired T=24 repeats) were canceled.
They are superseded, not numerical failures; their existing evidence is preserved.
Historical job 31643647 instead ran one large 10k T=24 MA57 case via
`kamiak_single_case.sh`, with the same frozen numerical source and CPU policy.

The old active paired launcher would otherwise start MA97 after med 2522 MA57.
A reserved MA97 output directory with `SKIPPED_BY_USER_REQUEST.txt` deliberately
activates its existing refusal-to-overwrite guard. No numerical source is changed.
The enclosing paired job may therefore exit 2 after the MA57 case is validated;
read the MA57 directory's own COMPLETE and validation reports, rather than
classifying this planned orchestration exit as a failed MA57 solve. The skipped
directory has no numerical metadata or COMPLETE marker and is never a result.

Preliminary centralized optimizer timings already available for med 2522 T=24
are 21.322863 s (MA57) and 26.671025 s (MA97), one observation each on the selected
node family and 25-CPU numerical budget. MA97 took about 25% longer here. These
are partial-campaign optimizer timings, not a completed repeated benchmark or
proof that the centralized backend difference is always small. Backend choice
for a later final comparison can be revisited without duplicating every ADMM run.
That earlier option is now superseded by the user's MA57-only direction.

Earlier pre-fix T=24 status: med 2522 task 31643458_6 stopped at its 5000-iteration warm-up
cap with a persistent two-cycle; it has no validated measured result. large 10k
job 31643647 converged in 793 warm-up and measured iterations. Both are
superseded by the validated numerical correction summarized above. The failed
med 2522 plateau was diagnosed rather than rerun with a larger cap.

This is a fresh experimental harness. Historical July drivers and result rows
are preserved. It implements the same overlapping temporal consensus ADMM:
each period owns one branch-flow OPF, local battery-energy copies span that
period and its immediate neighbors, their consensus is averaged, and scaled
duals are updated for every copy.

## Agreed problem

- Use the current shared test decks and SPAR-DDP exporter with periodic profiles.
- The new comparison uses July's C_B rule from `experiment_settings.jl`.
  `TADMM_CB_PROFILE=ddp` explicitly selects the preserved system-specific
  baseline weights: 0.001 IEEE123,0.0019 med 2522,0.000024 large 10k.
- **No terminal-SOC penalty and no terminal equality** (explicit user decision).
  Initial SOC, SOC bounds and lossless battery dynamics remain enforced.
- Include incoming-line real/reactive losses r*ell and x*ell, voltage-drop loss
  terms, nonnegative current square, and the second-order-cone inequality.
- Preserve the current SPAR-DDP benchmark's root-device convention: root balances
  contain only source import and outgoing flows. Root-connected devices are
  excluded from the OpenDSS reinjection check. The med 2522 root resource is
  therefore not an active network injection. This is an explicit matched-model
  convention, not a claim that all raw-deck resources are represented.
- Preserve the shared benchmark's absence of branch ampacity upper bounds.
  Do not describe a successful check as demonstrating thermal-limit compliance.

Existing soft-terminal SPAR-DDP timings are NOT matched comparators for this
free-terminal experiment. The new centralized solve here is the comparator.
Corrected large 10k is also a different physical instance from July's deck.

## Verified July defects and differences

`260731-ias-trans-tadmm:run_tadmm.jl` omits r*ell and x*ell from the nodal
balances. Current master includes those terms, but its stopping residual still
checks only diagonal SOC copies. This harness checks all copies. Its dual
residual repeats consensus changes with the correct copy multiplicity.
The archived `run_bf.jl` centralized baseline omits those loss terms too;
correcting displayed equations alone cannot repair the old numerical evidence.
Both sides of the present comparison use the corrected balances.

The default retains July's phase-aware rho policy: double every five iterations
until primal consensus first passes, then allow residual-based changes; a
15-iteration high-primal watchdog multiplies rho by 2.5. Every change rescales
the scaled duals. rho starts at 3000*sqrt(T/24), bounded between 1 and 1e6.
Phase transitions now use the maximum error across ALL copies. A generic
ten-iteration residual-balancing variant (`TADMM_RHO_POLICY=residual`) is retained
as a diagnostic; it reduced rho prematurely in the first med 2522 warm-up and
that diagnostic was interrupted and preserved, not reported as converged.

Battery-power variables exist only where local dynamics need them; unconstrained
powers for distant periods are removed. Dynamics at t and t+1 constrain the
three-copy local window. At consensus, the network-owned battery power at every
period obeys the original full-horizon energy dynamics. T=1 and T=2 windows are
well-defined. No terminal condition is introduced by an end window.

## Environment and solver layers

tADMM does not call FilterDDP/SPARDDP. It uses JuMP and Ipopt. The small
`Project.toml` and pinned `Manifest.toml` avoid reinstantiating the historical
plotting/interactive environment or requiring Gurobi. The DDP environment and
defaults are unchanged. The sequential pilot uses Julia 1.12.6 on ETRL309-ARYAN;
Kamiak uses 1.12.4 with a separately instantiated Linux environment and recorded
package manifest.

MUMPS is the initial Ipopt linear-system backend. `TADMM_LINEAR_SOLVER=ma57`
or `ma97` and `TADMM_IPOPT_OPTIONS` select alternatives for BOTH the centralized
reference and local OPFs. The local HSL builds are under
`C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build`. Include
`hsllib=<DLL path>;linear_system_scaling=none` when selecting them.
Library existence alone is not evidence that a backend loaded successfully.
`check_backends.jl` subsequently passed an actual Ipopt solve with each of
MUMPS, MA57 and MA97 on this machine; see `backend_checks.toml`.

The large 10k phase-aware MUMPS warm-up is an explicitly interrupted diagnostic.
A separately labeled MA57 pilot uses fixed rho=12, the heuristic energy-cost
curvature scale `2*C_B*kVA_B^2/dt`. That rule is a calibration candidate, not a
claim of an optimal or validated universal tuning rule. Select fixed rho with
`TADMM_RHO_POLICY=fixed` and `TADMM_RHO_INIT=12`. The IEEE123/med 2522 phase-aware
MUMPS cases and this MA57/fixed-rho case are different configurations; do not
present them as a uniform solver/tuning sweep or isolate a MA57 speedup from
this joint change. Every individual comparison uses the same linear backend
for its own centralized reference and ADMM stages.

The fixed-rho large 10k warm-up was also preserved as an interrupted diagnostic:
its dual residual passed while original-model primal feasibility decayed slowly.
The subsequent `_relative_ma57` pilot starts at the same rho=12 and uses relative
residual balancing from Wohlberg, *ADMM Penalty Parameter Selection by Residual
Balancing* (2017), https://arxiv.org/abs/1704.06209, Eqs. 49 and 51. With `E`
replicating the consensus vector, it compares `norm(b-E*z)/max(norm(b),norm(E*z))`
to `norm(rho*E*(z-z_old))/norm(rho*u)`. Denominators are floored at 1e-12.
Every five iterations it doubles/halves rho if the normalized ratio exceeds 10,
with rho bounded in [1e-6,1e6] and scaled duals rescaled. Final ABSOLUTE stopping
and physical validation thresholds are unchanged. This is an explicit controller
variant, not a silently substituted result for the July phase-aware policy.

A native conic backend is a separate future comparison. This harness currently
passes the cone as a smooth quadratic inequality to Ipopt, a general NLP solver.

## Execution and checks

From the repository root in PowerShell:

```powershell
& ./problems/mpopf/tadmm/review_pilot/run_case.ps1 -SystemName ieee123C_1ph -Horizon 6 -RunLabel pilot
```

The machine-local wrapper uses 309's existing Julia executable and depot. Each
new label gets a fresh result directory; existing evidence is never overwritten.
Source snapshots and deck hashes are frozen before execution.

Centralized and tADMM each receive a full warm-up, then restart from their
original battery-idle power-flow starts for the measured solve. Local OPFs reuse
their previous primal iterates within ADMM. No centralized optimum initializes
ADMM. Ipopt must actually terminate successfully; iteration-limit subproblems
are rejected. Warm-up/timed ADMM numerical traces must agree.

Default convergence requires maximum all-copy consensus error <=2.5e-7,
maximum dual residual <=1e-4, and original-model feasibility <=5e-7. The latter
also integrates dispatched battery powers from initial SOC and checks resulting
energy bounds, so small consensus errors cannot hide accumulated SOC drift.
An iteration cap is a failure, not convergence. Full traces allow later scoring
at looser objective/feasibility criteria without stopping the solve early.

`validate.py` requires objective agreement with the centralized reference within
1e-5 relative, independent original-model feasibility, matching warm-up/timed
trajectories, completed OpenDSS solves, no voltage violations at the shared
checker's 1e-4 p.u. limit tolerance, maximum voltage difference <=1e-4 p.u., and
OpenDSS/optimizer objective agreement within 1e-4 relative.

The shared OpenDSS checker reads raw DSS elements independently, reinjects only
the optimized device dispatch, and uses the same 1.05 p.u. source and constant
power device assumptions. Its `OPENDSS_ZERO_LINE_C=1` removes DSS default line
capacitance, because the branch-flow model has no shunt capacitance. This is a
check of the stated benchmark physics, not every native DSS modeling feature.
SOCP equality gaps and fictitious losses are recorded separately; convergence
of the relaxed problem alone is never called AC feasibility.

## Timing interpretation

The historical 309 pilots run sequentially: one Julia thread, one BLAS thread, one
OpenMP thread. No period-worker parallel simulation is performed.

- `central_solver_s`: Ipopt-reported centralized solve time, after warm-up.
- `central_wall_s`: wall time around its optimize call.
- `tadmm_algorithm_s`: measured sequential stage work plus coordination,
  excluding model builds, warm-up, physical diagnostic evaluation and CSV output.
- `tadmm_wall_with_diagnostics_s`: actual measured ADMM loop elapsed time,
  including per-iteration diagnostics/output, excluding build/warm-up/OpenDSS.
- `ideal_parallel_estimate_s`: sum of each iteration's maximum stage elapsed
  time plus measured coordination. This assumes one worker per period and zero
  communications/contention overhead. It is NOT measured parallel time.

Setup, builds, warm-up and validation are reported separately. Process samples
record working set, private bytes and CPU time every 0.5 s. The peak includes
setup/warm-up and resident sequential models; it is not parallel aggregate RAM.
These are exploratory single-repeat pilots, not publication-ready speedup rows.

## Kamiak protocol and frozen campaign

The user requested MA57/MA97 first on Kamiak; no Kamiak MUMPS runs are scheduled.
The numerical source is frozen at commit
`4fbae546bcc7ae6b4d53f6a1ed1d766dc6eb330e` in
`/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_parallel/source`.
Source and native Linux input-export SHA256 manifests are checked before every
campaign task. Orchestration is stored separately in `harness/` and hashed per
task, so later reporting work cannot change a running numerical experiment.

`parallel_workers.jl` uses persistent Julia processes, with local subproblems
assigned to workers and concurrent remote solves followed by ordered consensus
reduction. Workers have one Julia, BLAS and OpenMP thread each. Each worker is
pinned to its own physical CPU; another CPU is reserved for coordination.
The correctness gate uses two workers; system verification and benchmarks use
one worker per period. Workers rebuild models independently for the full warm-up
and the measured restart. Every stage result and actual worker layout is checked.

For timed cases, reserve all 64 physical CPUs on a 64-CPU Gold 6438Y+ node
(`OverSubscribe=NO`); Kamiak rejects the `--exclusive` flag. Numerical work is
restricted to T+1 CPUs: 7 for T=6 and 25 for T=24. The centralized MA97 baseline
may use that same CPU budget through OpenMP; MA57 uses one OpenMP thread.
Julia and BLAS remain single-threaded. This is a comparison under the same CPU
budget, not a claim that both algorithms internally use every available core.
Verification jobs share nodes and must not be reported as benchmark timings.

`optimizer_loop_wall_s` is the primary measured parallel ADMM loop time. It
includes actual worker communication, coordination, diagnostics, checkpointing
and output. `tadmm_algorithm_s` excludes diagnostics and output and is only a
breakdown. The ideal parallel estimate remains hypothetical and must never be
substituted for measured runtime. Warm-up, model construction and independent
OpenDSS validation are separate phases. `phase_events.csv`, per-process CPU/RSS,
aggregate private resident/PSS sampling, and available Slurm cgroup peaks support
resource analysis. Aggregate RSS can double-count shared pages.

Checkpoints every 25 iterations preserve numerical state and source identities.
A restart from a checkpoint is not a fresh contiguous timing result; no automatic
checkpoint resume is used for the reported timed solves.

Completed gates:

- Job 31643403: IEEE123 T=6 serial MA57 and MA97, full convergence and OpenDSS.
- Job 31643439: serial reference plus two-worker MA57/MA97, all 73 iterations;
  MA57 objective identical to serial, MA97 relative difference 6.16e-13 or less.
  All optimization, warm-up/restart, worker-layout, source and OpenDSS checks pass.

Submitted next, with results still pending when this note was written:

- Jobs 31643454 and 31643455: med 2522 and large 10k T=6, six workers, both backends.
- Array 31643458: T=6/T=24, all three systems, two repeats per backend; dependent
  on successful completion of both system gates. At most two task groups run
  concurrently. Backend order reverses on the second repeat. Completion markers
  require both backends' validated runs; partial/preempted attempts are preserved.

Subsequent verification update: MA57 completed both med 2522 (77 iterations) and
large 10k (132 iterations), including matching warm/measured traces and OpenDSS.
MA97 med 2522 has a much longer residual tail. This is not a converged benchmark
result. Saved-checkpoint diagnostic jobs 31643470, 31643476 and 31643484 compare
the same local OPFs without changing or restarting any measured trajectory.
Nominal `tol=1e-9` alone does not enforce equally small unscaled complementarity;
see [Ipopt termination options](https://coin-or.github.io/Ipopt/OPTIONS.html).
The 1e-10 absolute-complementarity diagnostic hit an iteration limit. In the
bounded diagnostic, adaptive barrier with `compl_inf_tol=1e-8` solved all six
subproblems successfully with both backends and much closer objective agreement.
That is evidence for an accuracy-control candidate, not proof of full ADMM
convergence or of the sole cause of the long tail.

Job 31643493 tests that candidate through full warm-up, measured restart and
OpenDSS on med 2522, both backends. Its labels begin `adaptive_verify_`; its case
wrapper and hash are captured separately and all extra Ipopt settings appear in
metadata. The original source, running trajectories and benchmark queue remain
unchanged. A candidate must be validated before any new benchmark protocol is
selected. `diagnose_checkpoint.jl` intentionally records failed local solves in
its diagnostic CSV; those failures are never accepted by the production driver.

Later completion update: both original system gates passed with MA57 and MA97.
The long MA97 tails converged in 1703 iterations on med 2522 and 568 on large 10k;
warm-up/restart and independent OpenDSS checks passed. The adaptive candidate's
MA97 run failed during warm-up with ALMOST_LOCALLY_SOLVED and is rejected.
The original benchmark configuration remains in use. Array 31643458 tasks 0/1
completed the two IEEE123 T=6 repeats; tasks 2/3 advanced to IEEE123 T=24.
See `results/kamiak_t6_update_31643458/README.md` and its archived evidence for
initial timings, timer boundaries and the preserved calibration failure.

The normalized relative-residual controller and curvature-scaled initial rho
were fixed across that original C3 campaign. Current bounded penalty tuning
is documented separately in DDP_TRANSFER_AUDIT.md. These experiments use the corrected
formulation, not a reproduction of July's omitted-loss model. The existing
soft-terminal DDP campaign remains a different optimization problem. A matched
DDP comparison needs the same free-terminal convention and hardware policy.
Keep full-convergence timings separate from any later posthoc 0.1% objective /
original-feasibility scoring. The former two-backend/repeat protocol is
superseded for the present first assessment by the user's MA57-only,
one-observation direction above.
