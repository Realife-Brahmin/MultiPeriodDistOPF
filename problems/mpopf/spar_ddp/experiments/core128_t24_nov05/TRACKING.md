# Tracking

Source/plan committed dae5cb42 on ddp/nov04-grouped-phase-b. Original T144
predecessor31575217 completed0:0, independently validated and now archived.

Prepared new immutable snapshots with inherited source/input SHA checks and
audit-only centralized runtime logging. Numerical DDP and Ipopt options are
unchanged from dense_merge_fullsolve_nov05. New allocator changes only64 to
the explicit64/128-core budget. Five offline regression maps passed; central
audit driver Julia syntax parsed successfully. Shell syntax passed.

* job31575546: /home/aryan.jha/kamiak_jobs/dense_merge_T24_nov05,
  Gold6438Y+,64CPUs300GiB4h,no-requeue, afterok31575217; now RUNNINGcn154.
  Arms original/latest/ma57/ma97 then reverse; two repeats each.
* job31575547: /home/aryan.jha/kamiak_jobs/dense_merge_core128_nov05,
  Xeon6760P,128CPUs600GiB6h,no-requeue, afterok31575546; pending dependency.
  Armsw63/w96/w48t2/ma57/ma97 then reverse; two repeats each.

All solves have full warmup and inherited50min per-process cap. Occupancy guard
requires the allocated node otherwise unoccupied. No numerical-factor reuse.
Do not edit either queued/running snapshot. source.commit is dae5cb42;
execution-source.sha256 pins all actual numerical/input/audit files;
harness.sha256 pins actual runtime helpers. Raw legacy affinity and checks
are preserved. New snapshots retain repaired compare_solutions.jl in ROOT.

Finish through collection/analysis, not polling only. Implement posthoc strict
collector against real outputs, reusing existing dense_merge_fullsolve and
battery_count extract/quality/validation checks. Existing collector assumes
c1/c4 and64cores so cannot be applied unmodified to these labels/maps.
T24 compare latest_r1/r2 against original SAMErepeat. T96 compare w96 and
w48t2 against w63 SAMErepeat using unchanged1e-10/same-iteration checks.
Central outputs use labels ma57_1/ma97_1 etc, resulting paths
results/10k_1020_T24_ma57_rma57_1, not bare r1.
Preserve comparison failures. Also verify reference-repeat agreement.

All arms require independent numerical/physical/OpenDSS/status/source/input/
actual masks/thread/backend checks. Validate all finite traces, mutual central
objective agreement<=1e-5relative and every final objective gap<=0.1%.
Record all DDP physical start hashes and matched central physical starting
objective. Score feasible<=1e-5 AND objective gap<=0.1% against fresh mutually
agreeing central references, separate from complete convergence. Faster central
backend chosen by two-repeat median, not fastest individual. Retain setup,
warmup/processwall and peak RSS limitations. T24 candidate uses24x2.
No scientific gate loosening to force a speedup or recolor historical rows.

If either batch fails, inspect full logs and preserve failure archive before
new immutable harness repair. Never blind retry/bypass failed dependency.
Record jobs/hardware/source/input/private-library hashes, raw logs/all traces,
solutions hashes, OS masks/accounting, archive binaries on cluster. Archive in
NEW local collection directories, keep final evidence under this experiment.
Add third paired comparisons only if timing overlap requires one. No further
networks/horizons/tuning/production/paper changes authorized by this request.

Completion heartbeat formerly named finish-spar-ddp-t144-full-solves is updated
to supervise BOTH new tests and retain T144 completion. Stay quiet while
healthy; notify meaningful failure/required decision or validated batch/final
results. Commit/push final source/evidence/assessment and verify remote before
deleting the monitor.

