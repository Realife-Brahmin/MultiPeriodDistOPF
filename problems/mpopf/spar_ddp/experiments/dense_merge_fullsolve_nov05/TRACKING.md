# Tracking

Base a6cf54904fe7276aa6cb377b2cf2c76666647c14 on
ddp/nov04-grouped-phase-b. Remote and workspace were checked before work;
no newer commit was present. The two unrelated untracked battery-count
reference files were preserved.

Active first batch is Kamiak job31574994, snapshot
/home/aryan.jha/kamiak_jobs/dense_merge_fullsolve_nov05.
It requests Gold6438Y+,64 CPUs,300GiB,3h, no requeue, with a node occupancy
guard. It has four full-warm-up T3 correctness arms followed by eight T96
main solves. At the initial check c1 smoke passed numerical, independent
physical and OpenDSS validation; c2 smoke was running.

The snapshot copies battery_count_oct21_campaign/source and links its original
T3/T96/T288 inputs. The original input manifest was verified before preparing
the copy. execution-source.sha256 and harness.sha256 pin the actual patched
files and inputs. Right-kernel algebra comes unchanged from a6cf5490.
Physical start arrays are hashed for per-arm agreement. Numerical controls
remain those of the matched full-solve study, including full warm-up and
pass-level GC.

Collection and analysis remain to be done. Verify full stdout/stderr and all
numerical/physical/OpenDSS markers. Run compare_solutions.jl for each candidate
against c1 from the same repeat, preserve comparison failures as evidence,
and run collect.py strictly only after all required outputs exist. Archive
raw logs, source/input/hardware/affinity checks, Slurm accounting, solution
hashes and analysis. Commit and push the T96 evidence. Report T96 results
before expanding to a valid candidate at T288. If timing ranges overlap, add
a third paired repeat before declaring a winner. No T288 job is queued yet.

Proactive supervision should stay quiet while pending/healthy and notify on
T96 completion/decision, a material failure, and final validated completion.
Send a concise status before8AM PT on November5 even if work remains queued.
Remove the dedicated monitor once the authorized campaign assessment is pushed.

Dedicated heartbeat finish-spar-ddp-dense-merge-full-solves is ACTIVE on this
thread every10 minutes. It carries the explicit T96-before-T288 reporting
gate, scientific validation, source-preservation rules, third-repeat rule and
meeting status deadline. At03:41 PT, c1/c2/c3 T3 smokes had passed all validators
and c4 T3 with two BLAS threads was starting. No main-run result was available.

Nov5 04:12PT post-processing update: c1 and c2 first T96 repeats have independently
passed all validators; c3 is running. The original standalone comparison helper
hit Julia's global soft-scope rule before producing a comparison. Its failed
log comparisons_T96/c2_r1.log is preserved. The committed helper now wraps its
unchanged gates in compare_main(args). The repaired copy is isolated at
SNAPSHOT/posthoc_v2/compare_solutions.jl with posthoc_v2/source.sha256; use this
path for all further comparisons. The active snapshot's original helper,
harness manifest and numerical solver sources were NOT edited.

The repaired c2_r1 comparison passes the unchanged1e-10/same-iteration gate:
both16 iterations, state scaled difference2.78e-15 and control scaled3.74e-16.
Physical x0/ubar hashes are identical for c1/c2. This is an interim comparison,
not the two-repeat campaign decision.

Nov5 T96 completed and strictly validated. Job31574994 COMPLETED0:0 in1h37m36s.
Two-repeat full medians c1/c2/c3/c4 are48.2925/45.9345/43.7375/41.3375s.
c4 wins both pairs with disjoint ranges,14.40% lower full time and34.07%
lower Phase B versus c1. All eight have16 iterations,25 backward passes,
8 mu updates and0 backtracks. All original validators and strict1e-10
solution matches pass, including the reference-repeat comparison.
Physical starting array hashes are identical. The collector now also
verifies every recorded OS-thread mask against disjoint physical-core maps.

The T96 result was reported to the user before extension. Select c4 for
the authorized T288 comparison against c1, two repeats in reverse order.
No third T96 repeat is needed. Commit the results_T96 evidence and comparison
table, write the validation/selection markers, then submit job.sh288.

The authorized T288 extension is submitted as job31575079 in the same immutable
snapshot, Gold6438Y+64CPU300GiB3h/no-requeue. Only c1 and selected c4 run,
ordered c1/c4 then c4/c1. T96.validated contains the strict collector result
and selected_candidate.txt is c4. T96 evidence is pushed at2110ef56.
The numerical sources and original harness were not edited for the extension.
Use posthoc_v2/compare_solutions.jl for comparisons_T288/c4_r1.toml and c4_r2.toml,
then collect.py ROOT OUT --T288. Keep the heartbeat until T288 comparison,
strict validation, final findings and remote push are complete.

FINAL Nov5: T288 job31575079 completed0:0 in1h51m49s. Strict post-processing,
all source/harness hashes, independent physical/OpenDSS validation, runtime
threads/OS masks and same-repeat solution comparisons pass. Full medians
c1=187.5305s and c4=142.4715s (24.03% lower). Phase B150.092112->91.1305385s
(39.28% lower), with A/C+other penalties retained. Both pairs improve and timing
ranges are disjoint, so no third pair is needed. All four runs have18 iterations,
26 backward passes,7 mu reductions and1 backtrack, with strict1e-10 solution
matches. The original comparison-helper scoping failure remains archived.

The authorized campaign is complete. results_T288 contains raw evidence,
strict validation and tables; FINDINGS.md summarizes both horizons, tradeoffs
and limitations. Setup/build/warm-up/process-wall times are kept distinct
from optimizer times. No further optimizer tuning or production/paper changes
were made. Verify the final evidence commit on the remote, then delete the
dedicated finish-spar-ddp-dense-merge-full-solves heartbeat.

NEW REQUEST after completed campaign: user authorized T144 original1020
reference-versus-winner full solves, two repeats each in reverse order.
Snapshot /home/aryan.jha/kamiak_jobs/dense_merge_fullsolve_T144_nov05,
job31575217, Gold6438Y+64CPU300GiB3h/no-requeue/occupancyguard, initially
pending on priority. Source copy inherits the exact validated procs/driver
and right kernel, with the original frozenT144 input linked. Source and
input hashes were rechecked. The standalone root comparison helper in
this NEW snapshot is already the repaired function-scoped version, also
available in posthoc_v2; old snapshots remain untouched.

Collect comparisons_T144/c4_r1.toml and c4_r2.toml against c1 from each
repeat. Run collect.py ROOT OUT --T144 strictly after T144.finished and
all independent validators. The collector's count gate was extended
with verified624/420 merge counts. Preserve any numerical differences,
failure logs and raw solutions/hashes. Archive raw evidence, write the
short matched comparison and push it. The two-horizon findings remain
historical until the new result is available. A new dedicated T144
completion monitor should be deleted after the final evidence push.
