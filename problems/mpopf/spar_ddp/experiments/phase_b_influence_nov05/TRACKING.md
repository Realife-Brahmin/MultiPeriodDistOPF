# Tracking

Worktree branchddp/nov04-grouped-phase-b. User's pasted request explicitly
authorizes one targeted baseline capture if consecutive-pair records are missing.
Read prior Oct22/23 final findings; no repeat compression implementation.

Existing capture inventory was read live:
value_coupling_audit_oct22/results/large10kC_1ph_T24_dense_capture/captures:
pass1k0j0,pass5k4j0,pass6k4j1,pass11k7j3,pass23k14j8.
Only3fullS perstate (1/12/24),85.9MB each. Later campaign capture ONLYk14j8
andfullS1:4. No ordinary consecutiveaccepted pairs/all24stage records.
Old baseline full pass log establishes fixed-barrier ordinarypairs0->1j0,
7->8j3,12->13j7 andbarrier4j0->4j1. These are planned keys, not yet
claimed actual new coverage. Eightstates/all24S+leaves~22GB,75GBfree observed.

Created bounded plan and capture helper at7def1879. LocalJulia1.12.7 AND
clusterJulia1.12.4 selector/layout/deepcopy/all24cache/regularizationmetadata
units and Julia source parsing passed. Shellsyntaxpassed,allUTF8/LF.
New runtimecopy only replaces frozen_value_capture.jl and adds one metadata
argument d.reg_last tofc_before call. Original solver arithmetic/decisions
and baseline balanced/currentkernel stay unchanged; no group/right/4BLAS.

New capture job31580863 snapshot
/home/aryan.jha/kamiak_jobs/phase_b_influence_capture_nov05.
Gold6438Y+64CPU300GiB90minno-requeue/occupancyguard,oneT24baselinecoldcapture,
40minprocesscap,24workersBLAS1. PRIVATEfirstdepot+readonlyoldcache.
DEPENDENCYafterok:31580246, so the existing24x4 experiment is untouched and
must finish first. Initially pendingDependency;31580246 pendingPriority.
Diskguard>=30GiB beforepreparation/solve. No expensive overlap.

capture_job validates originalnumerical/bounds, independentphysicalobjective/
primal andOpenDSS, then verify_capture.jl requiresactualeightstates,
192fullSrecords, consecutivepasses, samej/mu/k+1 and ACCEPTED earlierpass
forordinarypairs, samek/j+1/mureductionbarrier, allhashes andrawcapturetags.
If selectionchronology differs, FAILandpreserve;never manufacturecoverage.
All24live current-stage zero/sine/actualSOC fullcontrol/equality/bound probes
are retained. ActualSOC recursively uses liveCURRENTfeedbackdu andS.fx/fu.
PrephaseCfullSdeepcopied; rawnominalprimal/dual/barrier retained. Metadata
records actualreg_lastbeforeandreportedafter (originalprocssetsafter0).
Captured baseline clocks are NOTbenchmarks. Binarycapturessstayimmutablecluster.

Next required work (not finished merely by capture):
1 Finish/report/push24x4pairedstudy withoriginalstrictgates FIRST, preserve
   all results/failures. Do notedititsrunning/queuedsource orchangeaffinity.
2 Verifycompletecapture(stdoutstderr/Slurm/source/runtime/inputhashes/
   physical/OpenDSS/capture_validation.toml+actualrawlogchrono). Archive
   textualmetadata/provenance, not22GBbinaryinGit; hashespin retainedfiles.
3 IMPLEMENT bounded frozen replay perPLAN.md, reuse original productioninplace
   reference and _pp_solve/_tp_post. First reproduceALLsavedlive direction
   probes at unchangedscaled1e-7ANDrelativefloor1e-10 threshold1e-5; dense
   reference sensitivity is its own outcome. Do notattributeit tosubstitution.
4 Three selectedmergelocations fordelta distribution/depth statistics;
   initial24interventions only atinteriormerge depth2,left5,right7:
   H,B_L,B_R,h wholeprevious retention plusonlyHdiagonal stale andonlyB_L
   offdiagonal stale across4pairs. H=x.Mbb+y.Maa,BL=x.Mab,BR=y.Mab,h=x.mb+y.ma.
   These maybeinconsistentNewtonupdates: labelsensitivitydiagnostic.
   Originaluppersemantics/rawasymmetry preserved. Nofixedrankassumption.
5 ReexecuteALLremainingglobalbalancedrounds afterONEsubstitution, track
   dependencyDAG/affectedclosedvalues andALL24stagecurrentS reconstruction+
   fullbatteryforwardpath. Interiorrecord5isreadbylaterprefixrounds; do not
   useonlyrecord1atfinalroundthenmistakeunobservedchangeforharmlessness.
   Truecurrentfuture/full24horizon, no artificialtruncatedterminal.
6 Fullblock/vector norms/actions onactualsavedSOC+additionalcos/topology
   probes,fullcontrol/equality/bound errors withcoordinate+feederoutliers,
   delta singularspectra/ranks,diagonal/offdiag/block/row distributions,
   actualuppercondition/inertia and firstamplification location. Several
   probes/negativecontrol independentconsumer; all gates unchanged.
7 Charge diagnosticchange-detection/SVD/condition/reconstruction/fallback
   costs, not PhaseBbenchmarks. Stop bounded; no approximateoptimizer,
   newcompression/tolerance/worker/horizon/network campaigns.
8 Commitpushcode/RAWmeasurements/provenance/failures/conciseassessment and
   supported next-experiment recommendation. Preserve2untrackedbatteryrefs,
   defaults/paper/mainrows/MSOPF. CompletionheartbeatmustcontinueAFTERfirst
   study and capture untilauditfinalassessment/evidencepushed.

Purely read-only inventory and local/unit preparation while31580246 queued
did not modify or compete with its compute work. Capture is dependencygated.

