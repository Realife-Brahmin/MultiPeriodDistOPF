# sf6937 reference configuration

`canonical_configuration.json` records the original Kamiak package snapshots, selected solver backends, hardware allocation, measurements and their known uncertainties. The original environment files and their hashes are unchanged.

`cd_exact_status.csv` is the current per-direction coverage ledger: **all 78 whole-day battery certificates are exact-checked**, with 59 passing in-run checks and 19 passing VD audits. The [six follow-up audits](vd_audits_nov09/README.md) complete directions 2, 11, 14, 16, 17 and 19 at their original CD radii. The new CSVs and Slurm accounting are included; detailed log files remain in the original run archive.

The [historical report](REPORT.md) and original `sf6937_kamiak_reference_2026-11-03.tar.gz` preserve the earlier 72/78 snapshot. They are historical evidence rather than the current coverage ledger. The archive SHA256 in the canonical configuration still identifies those original bytes; the archive was not repacked.

Job environments and submit lines were reconstructed because Slurm 17.11 did not retain them. Unknowns remain labelled. External solver binaries are not redistributed. These follow-up audits were performed separately on Kamiak; integrating their evidence did not launch further jobs.
