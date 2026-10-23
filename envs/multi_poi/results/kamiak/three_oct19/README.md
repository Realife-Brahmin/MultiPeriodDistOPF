# Three-substation windows on ieee2522 and large10kD, Kamiak, 2026-10-19/20

Stages `three_ieee2522` (second and third substations at buses 947 and 1120) and
`three_large10kd` (buses 1120 and 1060), T = 24, PV reactive power dispatched, batteries idle
(`nobatt`) and cyclic, 16 directions each. One folder per Slurm job: `metadata.txt` (commit,
node, resources, wall time), `run.log`, and that job's rows.

- `*_g0`...`*_g3` (jobs 31315117-32): meant to run four directions each, but `THREE_DIRS=0,4,8,12`
  was split by `sbatch --export` at the commas, so each ran its first direction only (0, 90, 180,
  270 deg). Valid rows for those directions.
- `*_d1_2_3` ... `*_d13_14_15` (31316439-56): the other twelve directions, colon-separated.
- `*_bracket` (31316447/56, 31319022/23, 31321117): ieee2522 directions where the maximize-angle
  solve failed its edge check, redone by feasibility tests bracketing the predicted reach
  (idle 67.5, 90, 112.5; cyclic 67.5, 90, 112.5, 135 deg). large10kD uses the bracket throughout.

Merged result (bracket rows replace failed ray rows; none left with failed checks):
`../../multiperiod/three_large_windows.csv`.
