# sf6937 reference configuration

`canonical_configuration.json` records the original Kamiak package snapshots, selected solver backends, hardware allocation, measurements and their known uncertainties. The environment files are under `problems/msopf/environments/kamiak/`; their SHA256 hashes match the author-supplied reference archive. No new optimization or exact-audit jobs were run to prepare this release metadata.

`cd_exact_status.csv` is the original per-direction coverage ledger from the read-only campaign audit. An independent scan of the archived logs confirmed 59 passing in-run exact checks plus 13 passing VD audits, leaving directions 2, 11, 14, 16, 17 and 19 without an exact recheck. Preserve these original numerical certificates and this coverage distinction. The release notes summarize why those six were outside the audit batch.

The configuration is based on text launch records, logs, Slurm accounting, node descriptions, environment files and checkout history. Job environments and submit lines were reconstructed because Slurm 17.11 did not retain them. Unknowns are labelled rather than inferred as observations. External solver binaries are not redistributed.
