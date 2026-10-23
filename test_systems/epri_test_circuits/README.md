# EPRI test circuits (OpenDSS)

- Source: OpenDSS SVN trunk, `Version8/Distrib/EPRITestCircuits`
  (https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Version8/Distrib/EPRITestCircuits/),
  downloaded 2026-10-20 through its GitHub mirror https://github.com/dss-extensions/electricdss-tst
  (branch master), unmodified.
- Contents: `ckt5/`, `ckt7/`, `ckt24/`, and `epri_dpv/J1`, `K1`, `M1`.
- Masters used: `ckt5/Master_ckt5.dss`, `ckt7/Master_ckt7.dss`, `ckt24/master_ckt24.dss`,
  `epri_dpv/J1/Master_noPV.dss` (this mirror has a no-PV master for J1 as well as
  `Master_withPV.dss`), `epri_dpv/K1/Master_NoPV.dss`, `epri_dpv/M1/Master_NoPV.dss`.
- File fixes: none. Capitalization of Ckt5/Ckt7/Ckt24 file names needs no fix on Windows; on Linux
  check `Redirect` targets against the file names.
- Citation: "EPRI Feeder X, DPV Monitoring and Feeder Analysis, dpv.epri.com, 2013" (J1, K1, M1);
  EPRI test circuits distributed with OpenDSS (Ckt5, Ckt7, Ckt24).
- Snapshot results: `summary.csv`.
