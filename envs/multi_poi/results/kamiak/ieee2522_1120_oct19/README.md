# ieee2522, second substation at bus 1120: T = 12, 24, 96 (Kamiak, 2026-10-19/20)

Stage `ieee2522` with `LARGE_AT=1120`, PV reactive power dispatched, MA57, one job per horizon
(vcea partition, 2 cores, 16 GB). Each folder: `metadata.txt` (commit, node, resources, wall
time), `run.log`, and that job's rows; merged into `../../multiperiod/large_windows.csv`.
T = 12 reproduces the committed Surface result exactly.

| T | idle batteries (deg) | cyclic batteries (deg) | cost-optimal day | cyclic solve | peak memory |
|---|---|---|---|---|---|
| 12 | [-9.590, +8.013] | [-12.736, +10.691] | $8544.57 | 283 s | 2.1 GB |
| 24 | [-9.430, +7.851] | [-12.508, +10.473] | $8630.92 | 2980 s | 2.7 GB |
| 96 | [-9.162, +7.620] | [-12.344, +10.322] | $8717.31 | 30364 s | 6.9 GB |

- The window narrows slightly with resolution: from T = 24 to T = 96 by 0.27 / 0.23 deg (idle) and
  0.16 / 0.15 deg (cyclic), 1.3-3%. Finer steps resolve the peak hour's import better.
- Every edge sits at zero import (no-backflow); with cyclic batteries both edges also touch the
  0.95 pu lower limit. Predicted from the import model: idle [-7.59, +7.28], cyclic [-10.33, +9.92]
  at T = 96, conservative by 4-20%.
- OpenDSS replay of each cost-optimal day: |dP| at most 1.4e-4 kW, |dV| 7.9e-9 pu.
- Cost: the cyclic window at T = 96 took 8.4 h (no `MULTIPERIOD_LEAN` memory trouble: 6.9 GB).
