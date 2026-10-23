# Public distribution test systems (OpenDSS), as published

Downloaded 2026-10-20 by `bash tools/fetch_test_systems.sh` (re-runnable on any machine) and
solved once by `julia --project=envs/msopf tools/summarize_test_systems.jl`, which writes
`summary.csv` here and in each folder. No electrical data is changed. Single-phase conversion and
our own PV and batteries are a separate step (as done for the 9500-node system:
`envs/multi_poi/root_level/convert_to_1ph.jl`, decks `rawData/<system>_1ph`, checked by `check_1ph.jl`).

The raw files (30 MB in all, none over 50 MB) are committed as the source of truth (user,
2026-10-20); the fetch script reproduces them.

| Folder | Systems | Source |
|---|---|---|
| `epri_test_circuits/` | EPRI Ckt5, Ckt7, Ckt24; J1, K1, M1 (in `epri_dpv/`) | OpenDSS SVN trunk `Version8/Distrib/EPRITestCircuits`, via its GitHub mirror `dss-extensions/electricdss-tst` |
| `ieee8500/` | IEEE 8500-node | `GRIDAPPSD/Powergrid-Models`, branch `cimhub_refactor`, `models/feeders/OpenDSS/IEEE/IEEE8500` |
| `smartds_sfo_p1u/` | NREL SMART-DS v1.0, San Francisco sub-region P1U: 4 substations on a 69 kV sub-transmission network | OEDI bucket `s3://oedi-data-lake/SMART-DS/v1.0/2018/SFO/P1U/scenarios/base_timeseries/opendss/` over HTTPS |

Not downloaded: IEEE8500_3subs (published only as CIM XML; the user works in OpenDSS only, and the
9500-node system, already in `rawData/ieee9500_source`, is its successor); SMART-DS annual load
profiles (`P1U/profiles/`, more than 1000 files and at least 0.65 GB; the snapshot solve does not
need them); ANEEL BDGD (later, if at all).

## Results of the snapshot solve (nominal load)

| System | Buses | Loads | Regulators | Capacitors | Switches (open) | Load (kW) | Losses (kW) | Converged |
|---|---|---|---|---|---|---|---|---|
| EPRI Ckt5 | 2,998 | 1,379 | 0 | 4 | 73 (0) | 7,133 | 266 | yes |
| EPRI Ckt7 | 1,255 | 906 | 0 | 2 | 11 (0) | 76,608 | 230 | yes |
| EPRI Ckt24 | 6,058 | 3,891 | 1 | 3 | 0 (0) | 50,942 | 1,007 | yes |
| EPRI J1 (no PV) | 3,434 | 1,385 | 9 | 5 | 18 (0) | 10,950 | 469 | yes |
| EPRI K1 (no PV) | 1,282 | 321 | 1 | 1 | 2 (0) | 12,455 | 199 | yes |
| EPRI M1 (no PV) | 2,596 | 1,470 | 1 | 3 | 88 (0) | 15,467 | 191 | yes |
| IEEE 8500 (balanced) | 4,876 | 1,177 | 12 | 10 | 43 (0) | 10,773 | 1,210 | yes |
| IEEE 8500 (unbalanced) | 4,876 | 2,354 | 12 | 10 | 43 (0) | 10,773 | 1,275 | yes |
| SMART-DS P1U, whole region | 18,585 | 17,207 | 6 | 16 | 2,910 (0) | 121,926 | 7,769 | yes |
| P1U substation hs0 | 5,270 | 4,373 | 2 | 8 | 970 (0) | 49,882 | 1,821 | yes |
| P1U substation hs1 | 9,025 | 9,272 | 2 | 3 | 1,153 (0) | 37,883 | 2,033 | yes |
| P1U substation hs2 | 3,489 | 3,200 | 1 | 2 | 553 (0) | 17,503 | 617 | yes |
| P1U substation hs3 | 791 | 362 | 1 | 3 | 224 (0) | 16,659 | 366 | yes |

Buses count every bus at every voltage level (customer secondaries included). Load is the sum of
the loads' nominal kW; "regulators" counts RegControls. None of these systems has a normally-open
switch in its published configuration: the substations of SMART-DS P1U meet only through the
69 kV network, and the EPRI and 8500 circuits each have one source. Ckt7's nominal load
(76.6 MW, against 230 kW of losses) is as published and has not been looked into.

## How each master is run (in memory only)

Each master's commands run one by one from its folder, skipping those that run a study or open a
window (`Solve`, `Plot`, `Show`, `Export`, `Buscoords`, `Set mode`); then `Calcvoltagebases`,
snapshot mode, `maxiterations=100` (the 8500 needs more than OpenDSS's default 15; its own
`Run_8500Node.dss` sets 20), `maxcontroliter=200`, and one solve. For SMART-DS, `LoadShapes.dss`
reads profiles that are not downloaded, so each load shape it names is defined as a flat 1.0
instead. File-name capitalization (Ckt5, Ckt7, Ckt24) needs no fix on Windows; on Linux check the
`Redirect` targets against the file names.

## Citations

- EPRI Feeder J1 / K1 / M1, DPV Monitoring and Feeder Analysis, dpv.epri.com, 2013.
  EPRI Ckt5, Ckt7, Ckt24: EPRI test circuits distributed with OpenDSS
  (https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Version8/Distrib/EPRITestCircuits/).
- R. F. Arritt and R. C. Dugan, "The IEEE 8500-node test feeder," IEEE PES Transmission and
  Distribution Conference and Exposition, 2010. Files: https://github.com/GRIDAPPSD/Powergrid-Models.
- NREL, SMART-DS: Synthetic Models for Advanced, Realistic Testing: Distribution Systems and
  Scenarios, v1.0, https://data.openei.org/submissions/2981.
