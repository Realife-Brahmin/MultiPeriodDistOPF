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
