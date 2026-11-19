# Saved med2522 T96 chronological controls

Replay job31667879 completed with lossless batteries, no leakage and energy
carried through96 quarter-hour intervals. Native cap/storage gates passed.
Original saved optimization job31660059 used two nodes; replaying its endpoints
does not supply a new single-node timing benchmark.

Centralized passes chronological execution and the rating-normalized0.05
engineering target. Maximum SOC-state difference is1.64e-8 percentage points;
maximum dispatch error is5.55e-8 percent of rated power.

The saved tADMM endpoint is a failed trajectory:1000 iterations, final RMS
not met, objective disagreement38.97%. Chronological replay confirms failure:
96.61055kWh maximum requested bound violation, up to6.385kW dispatch refusal,
and maximum voltage difference0.05615pu. The independently normalized maxima
are32.69 SOC percentage points of state disagreement and100% of rated power
for dispatch tracking error. These maxima occur at different locations/times.

Independent maximum absolute SOCP equality gaps are7.79199e-7pu squared for
centralized and8.19689e-7 for the failed endpoint. P95 relative gaps are0.23945%
and0.38098%; maximum relative gaps approach100% at tiny flows. The denominator
floor is1e-12pu squared. Report absolute and relative gaps together. Tight
absolute cone slack does not establish battery or original-model feasibility.

`med2522-saved-sequential-nov09.tgz` contains complete battery/bus/branch/period
tables, native summaries and gate logs. SHA256 is in `provenance.json`.
`frozen-Storage.dss` is the exact ratings input used by the frozen campaign;
the replay overrides its95% efficiency values to100%, preserving zero idle
loss. Each saved-case directory's `socp-storage-audit.json` contains independent
table-hash checks, normalized errors, cone distributions and worst locations.
The derived full branch-gap CSV is reproducible with repository-root
`audit-socp-storage-nov09.py REPLAY_FOLDER frozen-Storage.dss` after extracting
the archive. Original tables and logs are never rewritten by this audit.

Fresh single-node cap10/30 runs31667880/31667882 remain queued at the last
check. Explicit SOCP/storage/memory follow-up jobs31667909/31667910 are attached
to their completion. Centralized remains fully solved with unchanged settings.
Optimization sampling retains RSS/PSS/private-memory/CPU records and phase
events; follow-ups summarize sampled peaks and Slurm memory accounting.
No fresh capped2522 results or scalability claims are made here.
