# Completed matched campaign: cloud-workspace evidence

All96 cases completed and passed the frozen campaign's collection validators.
Three systems x four horizons(24,96,144,288) x four configurations
(IpoptMA97,IpoptMA57,DDPowner,DDPbalanced) x two repeats.
Snapshot: /home/aryan.jha/kamiak_jobs/proper_logging_campaign_oct16.
Final job31233788 completed with originalDDP222.012s at large10kT288.

Paper edits remain paused in the309 session; this directory is the data handoff
for the cloud workspace to incorporate into the paper.

- completed_summary.json: all96 full optimizer times, iteration counts,
  configurations, repeats and raw-log provenance directories.
- near_opt_completed.csv: posthoc time to shared objective-gap and primal
  thresholds. Select the intended metric/threshold/gap explicitly; do not mix
  internal and original-physical infeasibility or full and posthoc times.
  Ipopt original residuals retain printed precision, as labelled in the CSV.
- completed_logs_checkpoint.tgz: complete collected textual evidence including
  solve logs, internal Ipopt timing logs, metadata, hashes, logging validation,
  OpenDSS and owner/balanced comparison evidence. Binary solutions are excluded;
  they remain on Kamiak. Unpack to inspect each original log directory.
- source.sha256: original campaign source manifest.
- campaign_plan.csv/submitted_jobs.csv: case plan and job mapping, including
  the earlier first submissions. Raw attempts/preemptions remain on Kamiak;
  final comparison should use validated complete attempts only.

All jobs use the matched periodic soft-terminal formulation with fixed
per-system costs, shared physical initialization and full warm-up. Refer to
each metadata/log for exact CPU, allocation, stopping criteria, preprocessing
and validation. Logged optimizer time excludes warm-up/setup and OpenDSS.
The original collector rechecks logging, completion markers and OpenDSS and
compares owner/balanced original-primal trajectories before exporting.

Illustrative full-time medians: large10kT288 owner222.803s,
balanced198.523s,MA97242.7485s,MA57381.0395s. These are not substituted for
near-optimality times. No MA97 multi-thread result is included here; that is a
separate ongoing diagnostic under ddp/results/ma97_thread_sweep.
