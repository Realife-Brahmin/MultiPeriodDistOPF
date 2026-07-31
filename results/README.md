# Curated results

A small, version-controlled snapshot of the **final (winner) results** behind the
IAS-Trans tADMM paper. Full per-run output (every `rho` in each sweep, solution `.jls`
vectors, model dumps) is *not* tracked — it is regenerated locally under
`envs/tadmm/processedData/` by `run_tadmm.jl` / `run_rho_sweep.jl`.

## Contents

- **`summary.csv`** — one row per reported `(system, T)` cell: the winning initial
  penalty `rho0`, the monolithic Branch-Flow (`BF_time_s`) time, the tADMM near-optimal
  effective time (`tADMM_NO_time_s`), and the resulting `speedup`. For `large10k, T=48`
  the BF time is the wall-clock at MUMPS memory-exhaustion failure (BF never converges).

- **`trajectory/`** — the raw convergence data for the two BF-vs-tADMM trajectory figures
  (Figs. 2–3 of the paper): `med2522_T144/` (winner `rho0 = 50000`) and
  `large10k_T48/` (winner `rho0 = 30000`). Each contains the winning tADMM run's
  `convergence_data.csv`, `near_optimal_summary.csv`, `results_socp_tadmm.txt` plus the
  BF `results_socp_bf.txt` and `ipopt_bf.log`. These feed
  `scripts/_export_traj_csv.py`, which produces the PGFPlots CSVs in the paper repo.

## Reproducing the full results

Winning `rho0` values per cell are in `summary.csv` (and the paper's Appendix A). To
regenerate any cell's full sweep:

```bash
SYSTEM_OVERRIDE=ieee2552C_1ph T_OVERRIDE=144 julia run_rho_sweep.jl
```
