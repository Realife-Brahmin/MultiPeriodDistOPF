# Corrected dense capture: coverage closed, replay still pending

Job31357503, source e6b8b37b, completed in12m47s on cn162, Gold6438Y+,
64CPU/300GiB. Slurm batch MaxRSS40481996KiB is job accounting, not an
isolated Phase-B peak. All four captures (IEEE123T3 gate, then three original
systems T24) passed numerical stopping, original physical-objective/bounds,
independent OpenDSS and actual sampling/file-hash coverage validators.

Verified middle/late (k,j): small(11,4)/(22,9), medium(26,4)/(52,8),
original1020-battery large(7,3)/(14,8). Each also contains early(k0,j0)
and consecutive first-barrier before/after passes at equal k with mu1->0.2.
Thus coverage is no longer restricted to j0. Exact labels, primal/dual/
complementarity, objectives, raw trajectories and file hashes are archived.

Raw logs retain the known OpenDSS precompilation method-overwrite errors;
these did not prevent its runtime validation. All four actual OpenDSS
checks passed, with maximum voltage discrepancies1.59e-9/1.48e-8/2.40e-6/
8.18e-9 pu (T3/small/medium/large). There were no Slurm stderr messages.

capture_evidence.tgz contains complete textual logs, coverage CSV/JSON,
selected-point metadata, hardware, input/source manifests and completion
markers. Binary leaves/local factors remain at the immutable cluster path
/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22/results; their SHA256
and sizes are in the metadata. They are not copied into Git. Home has7.3GiB
free after capture; preserve existing evidence and avoid another capture.

These are deliberately cold dense trajectories with diagnostic copy/write
overhead, not benchmark timings or evidence of a faster approximate solver.
Production defaults, benchmark result rows and paper are unchanged.

Next: isolated frozen replay, first saved full control/equality probe gate,
then zero-truncation hierarchy, dropped coupling, and adaptive retained
cross-block budgets1e-8/1e-5. Use all five points, medium full T24 suffixes,
and original1020 first-four-period window with its true dense future.
