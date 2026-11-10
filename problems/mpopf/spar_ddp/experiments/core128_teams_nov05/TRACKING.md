# Tracking

Authorized paired original large10k1020-battery T96 comparison, branch
ddp/nov04-grouped-phase-b. Harness893b7bf5 committed/pushed before submission.
Verified remote/worktree; preserve two unrelated untracked battery reference
files. Prior128-core evidence remains untouched and no other jobs were active.

Kamiak job31580246, snapshot/home/aryan.jha/kamiak_jobs/dense_merge_teams96_nov05.
Xeon6760P128physicalCPU600GiB3h/no-requeue, explicitlycn178 (same node as prior
128-core tests), occupancyguard. Initially pending priority.
Four COMPLETE fullwarmup+measured solves: w48t2_r1,w24t4_r1,w24t4_r2,w48t2_r2.
Both use97physical cores inclcoordinator. No fresh central run, additional
horizons or production-default change.

Exact validated variants/source copied from dense_merge_core128_nov05.
Original input and inherited execution-source manifest verified before copy.
No numerical source modifications; only launcher worker/BLAS selection differs.
Six offline physical-core maps pass, including24x4, no SMT/socket crossing.
Shell/Python syntax passed, UTF8/LF normalized to avoid earlierCR setup issue.
Runtime affinity.py and repaired compare_solutions.jl copied unchanged.
execution-source.sha256 and harness.sha256 pin actual runtime files.
Source-control helpers retain every-merge/A/C thread assertions, before/after
actual mask/backend audits, target _BLAS_DEFAULT and temporary _blas1 restore.

Finish collection/validation rather than stopping at polling:
* Require Slurm0,T96.finished and all4COMPLETE, empty/fully understoodstderr.
* Numerical status0/primaldualcomplementarity<=1e-5/originalbound gates,
  independent original physical/objective and OpenDSS, full finite traces.
* SameDDPx0/ubarhashes across ALLarms/repeats, actual48x2/24x4 worker masks
  disjoint physical groups with each team within a socket, coordinator1.
* Source/input/harness checks before/after, allsha and solution hashes.
* ROOTcompare_solutions.jl withJulia1.12.4 --startup-file=no for w24t4 versus
  w48t2 SAMErepeat; outputscomparisons_T96/w24t4_r1.toml andw24t4_r2.toml.
  Compare w48t2 repeats into reference_repeats.toml. SAMEiterations and
  scaledstate/control<=1e-10 unchanged; preserve false checks/differences.
* Archiveaccounting.txt via sacctJobID State Elapsed ExitCode MaxRSS AllocCPUS-P,
  source-check-final.log from variants/10k_1020/source againstROOTmanifest,
  harness-check-final.log and solution-final.sha256. Binarysolutionsstaycluster.
* Tarallrawroot excludingvariants,julia-depot,solution.jls; collectNEWlocaldir.
  Run THISexperiment collect.py ROOT OUT STRICT (no--mode; --partial isinterim
  only). It adapts priorDDPcollector to w48t2/w24t4 with unchangedgates and
  adds explicitSlurm/source/input/hardware/samesocketaudits. Debugposthoconly
  fromactualoutputs,no loosenings. No automaticpromotionon strictequivalence
  failure. Entire rawtraces should be compared, not onlyfinalsolution.
* Fulloptimizer+A/B/C, acceptedsteps,backwardpasses,mubarriers,backtracks,
  trajectories/GC/RSSand setup/fullwarmup/processwall separately. Reportas
  paired DDPablation, notcentralizedspeedup. Old48x2median26.686s iscontext.
  If desired sharedquality target mayreuse completed SAME128node reference
  fromcore128_t24_nov05/results_T96_core128; NO newcentralizedrun.
* If substantial timing overlap, add one thirdpairedrepeat ONLYwithNEWimmutable
  launcher/r3paths, exactsettings/source; currentjob.shisNOTresume-safe.
  No broaderconfiguration/horizon sweep. Preservefailedattempts and repair
  onlynewimmutablesnapshots,never edit queued/running source.

Commit/pushraw evidence, strict tables and conciseassessmentsamebranch,
verifyremote then remove dedicated completion heartbeat. Stayquiethealthy/
pending; notifyonlymeaningfulfailure,necessarydecision or validatedcompletion.


USER-APPROVED SCHEDULING CHANGE Nov5: exactcn178 restriction removed from PENDING31580246 via scontrol update JobId=31580246 ReqNodeList=. All numerical source/harness,128CPU600GiB3h,no-requeue and6760P-2.20ghz&cores-128 features unchanged. All4runs stillshareONEallocatednode. Before/after scontrol records preserved remotely in scheduling/. Initial ReqNodeList=ALL syntax was rejected as invalid node name and made no change; installed Slurm man page documents blank value to clear. Actual ReqNodeList now(null). Capture31580863 stillafterok31580246; nooverlap or job resubmission. Scheduler estimate remainsNov8 08:03PT immediatelyafterupdate and may change at recomputation. The previous exactcn178 notes are historical; future collection must accept ANYmatching6760Pnode while requiring same-node paired runs and actual audits.

FINAL:31580246 completed0:0 oncn179. Allfour COMPLETE, strict collector andsame-repeat solution/trajectory/scientific/affinity/provenance gates pass. 48x2 median26.7625s versus24x4 29.855s; B improves21.18% but fullsolve worsens11.56%. Raw archive/analysis/FINDINGS committed. Nooverlap-triggeredthirdpair needed. Capture31580863 began afterward; keep combined heartbeat through its full influence audit/report/push.
