# V3 process-cap interruption: preserve13points, continue ONLY missing2

Final update: continuation31358119 COMPLETED0:0 in10m36. Both missing large
points completed with unchanged code, clean stderr and source/input checks.
Strict final collection validates all15 as13+2 and preserves the original
FAILED attempt, partial files and timeout. See final_analysis/validation.json,
point_sources.csv and FINDINGS.md. This is replay coverage, not validation of
an approximate optimizer or a timing result.

Job31358012 did NOT complete: SlurmFAILED9:0 at20m22s, batchMaxRSS3792744KiB.
The large-case Julia process reached the harness12minute CASE-SERIES limit
and was killed after its grace period. This was not the45minute Slurm batch
limit, preemption or a reported numerical/representation guard failure.

Completed and posthoc-validated13of15 points: all five small T3, all five medium
T24, original1020 early/barrier-before/barrier-after four-period windows. The
interrupted large middle had passed dense saved probes, dense/exact/drop scans
but did not finish all modes or write metadata. Large late had not started.
Do not count either as complete. analyze_replay.py --partial reports13complete,
2missing, complete=false; never substitute a success marker for this attempt.

replay_v3_interrupted_evidence.tgz preserves ALL outputs including the partial
middle files, full timeout stack/Slurm stderr, hardware, source manifests and
the pinned overlay. The original cluster snapshot remains unchanged.

remaining_points_job.sh creates a NEW independent output root for ONLY large
middle and late. It refuses to repeat any point with completed v3 metadata,
requires the inherited small/medium gates and hashes, and calls the EXACT
unchanged v3 replay and numerical source. Separate8minute per-point caps,
18minute batch cap; no optimizer, capture, numerical change or tolerance
relaxation. New private writable depot, old depots/source read-only.

Final coverage must explicitly combine13+2 with per-point source/attempt
provenance. Preserve the failed v3 status/stderr and partial files; require
the continuation's two independent completion records and clean logs. Do not
pretend the old batch succeeded or overwrite its incomplete middle output.
All clocks remain non-benchmark diagnostics, not repeated optimizer timings.
