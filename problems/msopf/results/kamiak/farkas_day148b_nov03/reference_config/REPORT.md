# sf6937 Kamiak reference configuration (paper runs of 2026-11-03/04)

Read-only audit, 2026-11-09. Nothing was submitted, rerun, installed or modified on Kamiak. Sources: the job
folder `~/msopf/jobs/farkas_day148b_nov03`, the per-task `.meta` files written by the job script, Slurm
accounting (`sacct`), `scontrol show node`, the checkout's reflog, the environment files, and the Claude Code
session transcript for the submit lines. Kamiak runs Slurm 17.11, which stores neither the submit line nor
the job's environment, so the `FB_*` settings below come from the recovered `sbatch` commands
(`launch/submit_commands_recovered.txt`). `sbatch` exports the caller's environment by default.

## 1. Paper data → phase → result files

| Paper file (`figure-data/farkas_rated_nov03/`) | Phase | Kamiak result (`sfp1u_four/`) | Local copy (`problems/msopf/results/kamiak/farkas_day148b_nov03/sfp1u_four/`) |
|---|---|---|---|
| `sf6937_hour.csv` | A: idle edge | `fb_A_<d>.csv` | `fb_A_<d>.csv` |
| `sf6937_bracket.csv` | C with derived bounds (FB_CTAG=d) | `fb_Cd_<d>.csv` | `fb_Cderived_<d>.csv` (renamed: Windows case collision with `fb_CD`) |
| `sf6937_day.csv` | AD: battery reach | `fb_AD_<d>.csv` | `fb_AD_<d>.csv` |
| `sf6937_daycert.csv` | CD: whole-day certificate | `fb_CD_<d>.csv` | `fb_CD_<d>.csv` |
| none (audits) | V (ratings-only C), VD (CD) | `fb_V_*`, `fb_VD_*` | same |

`fb_C_*` (ratings-only certificates), `logs_C_rated_ratingsonly/` and `C_unverified/` are superseded and not in
the paper.

## 2. Configuration per phase

Shared by every phase:
- **Launcher:** `farkas_job.sh` (sha in `launch/`), from the checkout `~/msopf/repo_rated`.
- **Julia:** `~/opt/julia-1.12.7/bin/julia --startup-file=no --project=$HOME/msopf/envs/${FB_ENV:-sdp_hsl}`.
- **Environment variables:**
  - `FB_RATINGS=1 FB_LOADDAY=148 FB_SYSTEM=sfp1u_four`
  - `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 MKL_NUM_THREADS=${NTH:-1}`
  - `MULTIPERIOD_SOLVE_LOG=0`
- **Ipopt:** `linear_solver=ma57; hsllib=~/opt/hsl-private/install/ma57-3.11.3/lib/libma57.so; ma57_pivot_order=2; ma57_pre_alloc=3.0; linear_system_scaling=none`, plus the code's `tol = constr_viol_tol = 1e-9` and `max_iter = 3000`. `LD_PRELOAD=.../libfakemetis.so`.
- **Clarabel** (certificate phases): `max_iter 400, static_regularization_constant 1e-7, tol_gap_abs 1e-6, tol_gap_rel 1e-6, tol_feas 1e-6`, unchanged across all four commits. `FB_EXACT` default 1, `FB_EQTOL` 1e-5, `FB_FLOOR` 1e-4.
- **Workers:** one Julia process per Slurm array task, one task per direction. No Distributed workers, and Julia threads = 1.

| Phase | Script | Env | Conic LDL (FB_LDL) | Other settings | Slurm request | Array (job ID) | Commit(s) at task start |
|---|---|---|---|---|---|---|---|
| A | farkas_bracketing.jl | sdp_hsl | (Ipopt only) | | kamiak, 1 CPU, 6G, 6 h, `%40` | 31541120 (68 tasks) + 309 PC (dirs 69–78) | 335a30e3 |
| AD | farkas_dispatched.jl | sdp_hsl | (Ipopt only) | | kamiak, 1 CPU, 12G, 24 h, `%40` | 31541132 | 335a30e3 |
| C (paper: Cd) | farkas_bracketing.jl | sdp_hsl | qdldl (default) | FB_DERIVED=1, FB_CUTS=kcl0,jbound,lsec, FB_CTAG=d | kamiak, 1 CPU, 12G, 8 h, `%20` | 31552027 | 865e26d1 (all 78) |
| **CD** | farkas_dispatched.jl | **sdp_pardiso** | **mkl** (Clarabel + MKL Pardiso) | NTH=8 (MKL_NUM_THREADS=8), FB_CUTS=kcl0,jbound,lsec, FB_CD_BISECT=0 (grid 1.02/1.05/1.10), no derived bounds | kamiak, 8 CPU, 80G, 24 h, `%8` | 31541133 | 335a30e3 (19), 4f29d098 (22), 865e26d1 (37) |
| V (audit of ratings-only C) | farkas_bracketing.jl | sdp_hsl | qdldl | | kamiak, 1 CPU, 12G, 3 h | 31544190, 31544191 | 5712c808 (74), 4f29d098 (4) |
| VD (audit of CD) | farkas_dispatched.jl | sdp_pardiso | mkl | NTH=8, FB_CUTS=kcl0,jbound,lsec | kamiak, 8 CPU, 80G, 12 h | 31545814 (13 tasks) | 4f29d098 |

**Commits.** All four are on `msopf-nov02-current-limits`; the checkout was moved between them while the arrays
ran (reflog in `repo/`). Between 335a30e3 and 865e26d1 the only changes are:
- `farkas_bracketing.jl`, `farkas_dispatched.jl` and `current_bounds.jl`: the exact gate, phases V/VD, and the derived bounds;
- two new scripts.

The model code (`multiperiod_full_angle.jl`, `feasibility_certificate.jl`, `full_angle_pf.jl`) and the deck
`test_networks/smartds_p1u_xf_1ph` are unchanged across them. The checkout is clean today. The reflog shows only
checkouts, but uncommitted edits during the runs cannot be ruled out: **unknown, none recorded**.

**309 PC (directions 69–78 of phase A).** These ran on the 309 PC, `etrl309-aryan`:
- Hardware: Intel i9-10900X, 32 GB RAM, Windows 11.
- Software: Julia 1.12.7 (juliaup), with `C:/Users/.../msopf_envs/sdp_hsl` byte-identical to Kamiak's `sdp_hsl` (same SHA256), and the Windows MA57 build `Documents/hsl/build/libma57.dll` (SHA256 5b4a5bc3…).
- Code: the commit recorded is 335a30e3, the session worktree's HEAD. That worktree also held unrelated uncommitted files, so a clean tree is **unknown**.
- Launcher: `launch/local_worker.sh`, with `launch/sync_kamiak.sh` copying its rows to Kamiak.

## 3. Environments, packages and libraries

| File | SHA256 |
|---|---|
| `sdp_pardiso/Project.toml` (2026-11-01) | f44fd8dd9b8e88c752f0786d26d8bb2ad0dbca6044a4e99b7324205697b2724d |
| `sdp_pardiso/Manifest.toml` | b2a86d4c478ae4dd29d5365f8413a670727accba81715b56aff1245a7d10cc9e |
| `sdp_hsl/Project.toml` (2026-10-26) | 08d3012c07c72b74a345f2d83108b7f5463dbca70ec3b3ff7e55a343fdfd4390 |
| `sdp_hsl/Manifest.toml` | 1510a10c16e835dba8aca50b841e2c2fce42404fd5ecf5999ac4d37a2c0b5c2b |

Both files predate the runs and were not modified afterwards. `sdp_pardiso` is `sdp_hsl` plus Pardiso.jl; every
version that both environments share is the same.

**Julia 1.12.7**, official binary, commit 6d172b025e4befc4d274d9fbc9339917a8a86b65, Manifest format 2.0.

| Package | Version | Used by the paper phases? |
|---|---|---|
| JuMP | 1.30.0 | yes |
| MathOptInterface | 1.50.0 | yes |
| Ipopt / Ipopt_jll | 1.14.1 / 300.1400.1901+0 | yes (A, AD; MUMPS_seq_jll 500.800.200 installed, not selected) |
| Clarabel | 0.11.1 | yes (C, CD, V, VD) |
| Pardiso.jl / MKL_jll / IntelOpenMP_jll | 1.1.2 / 2025.2.0+0 / 2025.2.0+0 | CD and VD only (sdp_pardiso) |
| QDLDL | 0.4.1 | C and V (the conic LDL), sdp_hsl |
| OpenDSSDirect | 0.9.9 | yes (A, AD OpenDSS replays) |
| DSS C-API (OpenDSS engine bundled with ODD) | 0.14.3 (`deps/linux/libdss_capi.so`) | yes |
| HSL.jl / HSL_jll | 0.5.6 / 4.0.8+0 | installed; **not** the MA57 Ipopt loaded |
| Gurobi / Gurobi_jll | 1.9.2 / 13.0.0 | installed, not used |
| libblastrampoline / OpenBLAS_jll / OpenBLAS32_jll | 5.15.0 / 0.3.29 / 0.3.30 | Julia BLAS |

**HSL actually loaded.** The private build `~/opt/hsl-private/install/ma57-3.11.3/lib/libma57.so` (MA57 3.11.3,
SHA256 6125be8d…e2dee5), set through Ipopt's `hsllib` in every phase's `.meta`:
- **Build:** 2026-10-14, gfortran `-O2 -fPIC`, `--with-blas=-lopenblas`, AMD ordering.
- **BLAS it links:** the system OpenBLAS, `/lib64/libopenblas.so.0` → `openblas-serial-0.3.26-2.el9` on the login node. The compute nodes' copy was not inspected: **unknown**, presumably the same Rocky 9.6 image.
- **METIS:** a dummy `libfakemetis.so` (SHA256 6d1db857…) preloaded.

MA97 was built (`hsl_ma97-2.8.1*`) but not used. Build scripts and `config.log` are in `hsl_build_text/`; the
binaries are not included.

**Backends actually used, and how we know:**
- **CD/VD:** every CD and VD log names the solver as `(…, mkl)`, which is Clarabel's `direct_solve_method = :mkl` (MKL Pardiso via Pardiso.jl's extension).
- **Peak memory:** the CD logs record a max RSS of 43.5 GB median and 46.4 GB maximum over 107 relaxations.
- **Threads:**
  - MKL_NUM_THREADS was 8. Clarabel calls `Pardiso.set_nprocs!(ps, max_threads)` with `max_threads = 0` (its default, "automatic"), so the effective Pardiso thread count was not logged: **unknown**.
  - Measured CPU use (TotalCPU over elapsed) was 2.2 busy cores median and 2.85 maximum, out of the 8 allocated.
- **Julia's BLAS:** OpenBLAS via libblastrampoline with OPENBLAS_NUM_THREADS=1. Nothing in the code forwards LBT to MKL (no `using MKL`, no `lbt_forward`). The BLAS backend in use at run time was not printed: **inferred, not recorded**.

## 4. Compute resources (from the job records)

All Kamiak tasks ran on partition `kamiak` (the backfill queue), account dubey.

| Phase | Tasks | CPUs | Requested memory | Peak RSS (sacct, batch step), median / max | Wall time, median / max | Nodes |
|---|---|---|---|---|---|---|
| A | 68 (+10 on the 309 PC) | 1 | 6 GB | 1.6 / 1.6 GB | 0.32 / 0.52 h | Granite Rapids 6760P 2.2 GHz (cn176–187) |
| C (Cd) | 78 | 1 | 12 GB | 7.1 / 8.3 GB | 0.53 / 1.00 h | Granite Rapids |
| AD | 78 | 1 | 12 GB | 4.6 / 5.2 GB | 1.21 / 2.75 h | Granite Rapids |
| **CD** | 78 | **8** | **80 GB** | **43.7 / 45.7 GB** | **2.71 / 9.03 h** | 74 on Granite Rapids 6760P (128 cores, 1–1.5 TB); 4 on Sapphire Rapids Gold 6438Y+ 2.0 GHz (64 cores, 512 GB) |
| V | 78 | 1 | 12 GB | 3.0 / 3.3 GB | 0.19 / 0.26 h | Granite Rapids |
| VD | 13 | 8 | 80 GB | 43.6 / 44.4 GB | 1.52 / 3.58 h | 12 on Sapphire Rapids, 1 on Granite Rapids |

- Nodes are shared and not exclusive, and hyperthreading is off. Per-task data: `slurm/per_task.csv`.
- Node specs come from `scontrol show node`: `slurm/compute_nodes.txt`.
- Job IDs of every task: `slurm/sacct_all_tasks.txt`.

## 5. Finding: six whole-day certificates were never checked in exact arithmetic

The exact gate for CD went in at 5712c808, checked out at 11:30 PST on 2026-11-03. CD tasks that started earlier
ran 335a30e3, which has no exact check (19 directions). VD audited the 13 of them that had finished by then, and
all 13 were certified. The other six finished later without an exact check or an audit:

| Direction | Solver status | Float bound |
|---|---|---|
| 2 | INFEASIBLE (Farkas ray) | 5.29e4 |
| 14 | INFEASIBLE (Farkas ray) | 5.13e4 |
| 16 | INFEASIBLE (Farkas ray) | 5.20e4 |
| 19 | INFEASIBLE (Farkas ray) | 5.10e4 |
| 11 | SLOW_PROGRESS | U −331.8 kW |
| 17 | ITERATION_LIMIT | U −3.5e36 kW |

The paper's sentence "each one with the batteries in use, was also recomputed in exact rational arithmetic" is
therefore not yet true for sf6937: 72 of 78 are exact-checked, 59 in-run plus 13 by VD.

Closing it needs VD on directions 2, 11, 14, 16, 17 and 19: six 8-CPU, 80 GB tasks of about 1.5–3.6 h each. That
was **not run**, because this audit is read-only. Per direction: `cd_exact_status.csv`.

## 6. Unknown or not recoverable

- **Environment variables:** the exact environment of each job. Slurm 17.11 does not store it; reconstructed from the submit lines.
- **Working-tree state:** whether `repo_rated` had uncommitted edits during the runs (now clean), and the 309 PC worktree's state.
- **Threads and BLAS at run time:** the effective MKL Pardiso thread count, and the BLAS library actually bound at run time. Neither was logged.
- **Compute-node system libraries:** the system OpenBLAS version on the compute nodes.
- **Earlier runs:** only the last run of each (phase, direction) is in the `.meta`/`.log` files, which are overwritten on rerun. Earlier attempts are visible only through `sacct` and the submit lines. `logs_C_rated_ratingsonly/` preserves the overwritten ratings-only C logs.

## Archive contents (`sf6937_kamiak_reference_2026-11-03.tar.gz`)

- `envs/`: both Project/Manifest pairs and SHA256SUMS.
- `launch/`:
  - `farkas_job.sh`;
  - `submit_commands_recovered.txt`;
  - the 309 PC's `local_worker.sh` and `sync_kamiak.sh`.
- `hsl_build_text/`: README, build scripts and MA57 `config.log`.
- `libs/`: the Julia commit, the DSS C-API version and `ldd` of libma57.
- `slurm/`: sacct records, per-task CSV and JSON, and node specs.
- `repo/`: reflog and commit list.
- `logs/`:
  - every `.meta`, `.log` and `.err` of the campaign;
  - the superseded ratings-only C logs;
  - `C_unverified/`.
- `cd_exact_status.csv`.

No binaries, no result archives and no certificate files are included.
