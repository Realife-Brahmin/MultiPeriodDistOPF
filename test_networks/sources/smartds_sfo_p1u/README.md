# NREL SMART-DS v1.0, San Francisco sub-region P1U (OpenDSS)

- Source: public OEDI bucket, `s3://oedi-data-lake/SMART-DS/v1.0/2018/SFO/P1U/scenarios/base_timeseries/opendss/`
  (293 files, 15 MB), fetched over HTTPS without sign-in on 2026-10-20, unmodified.
  Dataset page: https://data.openei.org/submissions/2981.
- Contents: the whole region (`Master.dss`: a 230 kV source, the 69 kV sub-transmission network in
  `subtransmission/`, and four 12.47 kV substations `p1uhs0_1247` ... `p1uhs3_1247`), plus a master
  for each substation and each feeder.
- Not downloaded: the annual load profiles the load shapes read (`P1U/profiles/`, more than 1000
  files and at least 0.65 GB). For the snapshot solve each load shape is defined in memory as a
  flat 1.0; the masters' own `Solve mode=yearly`, `Export` and `Plot` lines are skipped.
- File fixes: none.
- Citation: NREL, "SMART-DS: Synthetic Models for Advanced, Realistic Testing: Distribution
  Systems and Scenarios," v1.0, https://data.openei.org/submissions/2981.
- Snapshot results: `summary.csv`.
