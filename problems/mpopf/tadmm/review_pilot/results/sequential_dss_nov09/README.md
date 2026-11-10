# Chronological OpenDSS verification — November 9, 2026

These are passive replays of saved large10k, T=24 endpoints. No optimization
was rerun, and replay time is excluded from all solve-time benchmarks.
The primary battery model has 100% charging and discharging efficiency,
zero idling/leakage loss, and no inverter efficiency curve.

Centralized passes the full chronological schedule check. The cap10 and cap30
tADMM endpoints fail it, despite excellent network-snapshot agreement.

- Centralized: maximum energy-accounting error2.974701e-7kWh; maximum
  engine-energy versus OPF-state difference4.003886e-7kWh; no SOC violation.
- Cap10: maximum energy-accounting error0.798033kWh; maximum unclipped
  requested-schedule SOC-bound violation0.818639kWh; maximum engine-state
  versus OPF-state difference0.737307kWh;1,276 clamped battery-period records.
- Cap30: maximum energy-accounting error0.644060kWh; maximum unclipped
  requested-schedule SOC-bound violation0.668225kWh; maximum engine-state
  versus OPF-state difference0.436059kWh;841 clamped battery-period records.

All reported errors are maxima over individual batteries/intervals, not
averages or fleet totals. The energy tolerance is0.01kWh (1e-5 times the
1000kVA base). Raw tables retain numerical precision for other tolerances.

## Why the old snapshot pass was insufficient

The shared snapshot checker substitutes constant-power injections for Storage.
It verifies network agreement but cannot certify sequential energy execution.
The new replay keeps real Storage enabled, initializes energy once, commands
signed power in `DispMode=External`, and carries energy through exactly one
`FinishTimeStep()` per interval. Actual terminal power, stored energy, OPF
states and independently integrated requested energy are all compared.

The installed engine can deliver requested AC power over a coarse interval and
then clip energy at the reserve. Bounded engine SOC therefore does not suffice:
the energy-accounting check detects the inconsistency. The capped cases fail
this check even with perfect efficiency and no leakage.

The historical Julia implementation did use live Storage. The switch recalled
by the author was Follow to **External**, documented in commits
[92a6cf8](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/commit/92a6cf804b829f451394ea2215b06d5792aa6212)
and [a974207](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/commit/a974207318ab23e363d36c00cd5de4ea0e97cf57).
The historical detailed Julia validator is retained at
[28ed0c2](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/blob/28ed0c233841872455237687a2872054d8596a63/src/openDSSValidator.jl).

## Artifacts and reproduction

`chronological-discrepancies.png` and `.pdf` show the primary comparisons.
The three `*-model-worst-battery-discrepancies.csv` files list the100 worst
battery records per case. Complete primary tables are stored in
`sequential-dss-{central,cap30,cap10}-model-nov09.tgz`; extract each archive
with `tar -xzf ARCHIVE` into its own directory. Each contains compressed
battery, bus, branch and period CSV tables plus a TOML summary. The compressed
CSV files can be read directly by pandas. Battery energy is in kWh, active
power in kW, reactive power in kvar; voltages and squared-current variables
are in per unit as identified by their column names.

`sequential-dss-metadata-nov09.tgz` preserves frozen source, Storage input,
logs and combined summaries. `sequential-dss-archives-nov09.json` records
SHA256 for all seven campaign archives. The three additional `*-dss-*`
archives retain a separate95%-efficiency mismatch diagnostic locally/on
Kamiak; they are not the requested lossless model and are not bundled here.

Successful main replay: Kamiak31667060, source
ab3713e33022f3f78c5bffa202d5f8872536e20d, schema1.
Latest-source gate:31667199, source97f68c75cfcc282d787a06584f7a2cdb3664a000,
schema2; logs and the repeated centralized summary are in
`metadata/current-gate/`. Schema2 adds explicit metadata for efficiencies,
base, network size and idle loss; primary physics is unchanged.

The latest native gate passed73 assertions, including idle no-leakage,
reserve/capacity clamp detection, signs, clock continuity, quarter-hour
durations, efficiency diagnostics and prevention of double-applied profiles.
`independent-table-verification.json` records independent re-scoring of
3,119,328 rows across six cases. Passing that integrity check does not mean
every schedule passed physical validation. `report-scope-tests.json` confirms
snapshot-only legacy reports cannot imply full-MPOPF validation.

Reproduce native setup with repository-root
`setup-sequential-dss-nov09.sh`; test latest-source integration with
`setup-sequential-current-gate-nov09.sh`. Independent checks and plots use
`verify-sequential-tables-nov09.py`, `test-sequential-report-scopes-nov09.py`
and `plot-sequential-discrepancies-nov09.py`. Original failed setup logs are
preserved; the corrected unity loadshape is explicit in frozen source.
