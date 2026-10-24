# MSOPF on Kamiak

How MSOPF studies run on WSU's Kamiak cluster (set up 2026-10-19).

## Layout on Kamiak

- `~/msopf/repo`: a clone of `MultiPeriodDistOPF` on the MSOPF branch, kept at the remote's
  tip (`git fetch && git reset --hard origin/<branch>`). Separate from `~/MultiPeriodDistOPF`,
  which belongs to the DDP work; never run MSOPF jobs from that one.
- `~/msopf/jobs/<study>_<date>/`: one folder per study, holding a copy of `msopf_job.sh`,
  `logs/` (Slurm's own output; create it before `sbatch`, or the job fails at once), and one
  `out_<label>_<jobid>/` per job.
- Julia 1.12.7: `~/opt/julia-1.12.7/bin/julia`. HSL MA57 (Linux build of the user's licensed
  sources): `~/opt/hsl-private/install/ma57-3.11.3/lib/libma57.so`, loaded with
  `LD_PRELOAD` of its `libfakemetis.so` and `ma57_pivot_order=2` (AMD ordering; see
  `~/opt/hsl-private/README.md`). Keep HSL out of Git.
- Compute nodes reach the internet: `Pkg.instantiate()` of `problems/msopf` ran in a job.

## Submitting

[`msopf_job.sh`](msopf_job.sh) `<label> <julia args ...>` runs Julia in `problems/msopf` from the
clone and writes `out_<label>_<jobid>/`: `metadata.txt` (commit, node, CPU, solver, Slurm
resources, exit code, wall time, peak memory from `sstat`), `run.log`, `run.err`. It sets
`MULTIPERIOD_OUT_DIR` to that folder, so concurrent jobs never share a CSV, and
`SOLVER=ma57` selects MA57. Example, from the study folder:

```bash
mkdir -p logs
sbatch --partition=kamiak --account=dubey --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G \
  --time=03:00:00 --job-name=msopf_2522_T24 --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --export=ALL,SOLVER=ma57,LARGE_AT=1120,LARGE_T=24,MULTIPERIOD_LEAN=1 \
  msopf_job.sh ieee2522_1120_T24 problems/msopf/root_level/multiperiod_study.jl ieee2522
```

`kamiak` is the backfill partition: jobs can be pre-empted and requeued, which these short
single-process studies tolerate. `/usr/bin/time` does not exist on the compute nodes.

## Portability (2026-10-19)

`portability_check.jl` recomputes committed IEEE123 results (T = 24, PV reactive power free;
`main_horizon_windows.csv`, `main_days.csv`). Evidence in
`problems/msopf/results/kamiak/portability_oct19/`, commit `5f91b4cc`:

| run | windows vs committed | imports vs committed | day cost | wall |
|---|---|---|---|---|
| 309 PC, MUMPS | 2.5e-6 deg | 5.0e-5 kW | $2433.8267 | 141 s |
| Kamiak cn185, MUMPS (job 31310694) | 2.5e-6 deg | 5.0e-5 kW | $2433.8267 | 156 s |
| Kamiak cn185, MA57 (job 31310695) | 2.5e-6 deg | 5.0e-5 kW | $2433.8267 | 106 s |

No failed edge check in any run; peak memory 1.4 GB.
