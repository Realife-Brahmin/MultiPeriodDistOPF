# Release history and archived notes

This file preserves publication history separately from the study overview. The Releases tab presents one current MSOPF release. Removing an older release entry does not remove its Git tag, source tree, commit history or archived results. No release entries listed here have attached binary assets.

## Maintaining links

- Use `https://github.com/Realife-Brahmin/MultiPeriodDistOPF/releases/latest` for navigation to the current release.
- Cite a fixed Git tag or commit for archived code and results. The PESGM 2027 study uses `msopf-pesgm2027-v1.0.1`; its source URL is `https://github.com/Realife-Brahmin/MultiPeriodDistOPF/tree/msopf-pesgm2027-v1.0.1`.
- Preserve existing tags when consolidating the Releases listing. Do not move a cited tag to newer code.
- If a future study uses a new snapshot, record the new tag here and update that study's citation explicitly. Existing citations continue to identify their original snapshot.

## Notes before consolidation (10 November 2026)

### MSOPF PESGM 2027 — complete 78/78 exact checks (v1.0.1)

Git tag: `msopf-pesgm2027-v1.0.1`. [Preserved source tree](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/tree/msopf-pesgm2027-v1.0.1). Published: 2026-11-10T07:56:51Z.

All **78/78 sf6937 whole-day battery certificates now have passing exact-arithmetic checks: 59 in-run and 19 by VD audit**. This release includes the six final Kamiak audits for directions 2, 11, 14, 16, 17 and 19, integrated through #189. It supersedes v1.0.0’s 72/78 coverage record.

The six CSV hashes match the documented SHA256 values for their LF Git bytes, all six audited angles exactly match the original published CD points, and Slurm records show successful completion. Direction 11 is certified by a negative exact objective upper bound despite a slow-progress floating-point solver status; the other five final audits use Farkas certificates. The result CSVs and accounting are included; detailed logs remain in the original run archive.

Updated: per-direction coverage ledger, canonical metadata, result index, viewer text and release documentation. **The numerical points, map geometry, model, network inputs and recorded solver environments are unchanged.** The original report and reference tarball retain their earlier audit snapshot and are clearly identified as historical records. Packaging did not launch further jobs or rerun the rational verifier.

Start with the [interactive maps](https://realife-brahmin.github.io/MultiPeriodDistOPF/farkas_maps.html), [paper-result index](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/blob/msopf-pesgm2027-v1.0.1/problems/msopf/RESULTS.md), [six-audit evidence](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/blob/msopf-pesgm2027-v1.0.1/problems/msopf/results/kamiak/farkas_day148b_nov03/reference_config/vd_audits_nov09/README.md), and [hardware/software reference](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/blob/msopf-pesgm2027-v1.0.1/problems/msopf/environments/kamiak/README.md).

Download **Source code (zip)** for the focused MSOPF snapshot. The released HTML includes its numerical data and loads Plotly.js from a CDN. Connected display surfaces remain interpolations between sampled points, not certificates of continuous regions or evidence of superiority to radial reconfiguration.

### MSOPF PESGM 2027 reference snapshot — v1.0.0

Git tag: `msopf-pesgm2027-v1.0.0`. [Preserved source tree](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/tree/msopf-pesgm2027-v1.0.0). Published: 2026-11-10T07:36:12Z.

This release packages the rated-current, multi-source, multi-period distribution OPF study with prescribed substation voltage phasors, PV reactive dispatch, and batteries returning to their initial charge. The dedicated `msopf-pesgm2027-release` branch presents the MSOPF code, shared network inputs, and archived results without the separate DDP workstreams. The full research history remains on `master`; preparation was merged in #188.

**Start here:** [interactive result maps](https://realife-brahmin.github.io/MultiPeriodDistOPF/farkas_maps.html), [paper-result index](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/blob/msopf-pesgm2027-v1.0.0/problems/msopf/RESULTS.md), and [release notes](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/blob/msopf-pesgm2027-v1.0.0/problems/msopf/RELEASE.md). A versioned HTML snapshot is included under `docs/`; its plotting library loads from a CDN.

**Recorded reference:** sf6937 on Kamiak, Julia 1.12.7; JuMP 1.30.0, MathOptInterface 1.50.0, Ipopt 1.14.1 with an external MA57 3.11.3 build, Clarabel 0.11.1, and OpenDSSDirect 0.9.9 / DSS C-API 0.14.3. Whole-day certificates used Pardiso.jl 1.1.2 with MKL 2025.2.0. Both original Julia environments, package credits, environment hashes, and external-library requirements are included. Proprietary solver binaries are not distributed.

**Hardware:** the 78 sf6937 whole-day tasks requested 8 CPUs and 80 GB RAM each, on Intel Xeon 6760P and Xeon Gold 6438Y+ nodes. Recorded median wall time was 2.71 h per direction and median Slurm peak memory was 43.7 GB. Allocation is not a measured parallel-speedup claim; further details and uncertainties are documented with the environments.

**Validation coverage:** sf6937 retains 78 numerical whole-day certificates. Exact-arithmetic checks cover 72/78, all passing. Asynchronous completion left six directions (2, 11, 14, 16, 17, 19) outside the earlier audit batch; their original numerical results are preserved without a subsequent exact recheck. Per-direction coverage is included. No new solves or audits were run for this release. The displayed maps interpolate sampled results and do not certify every point between samples or establish superiority to radial reconfiguration.

Download **Source code (zip)** below for the focused snapshot. Browsing results requires no Julia installation; recomputation requires the documented numerical environment and external libraries. See the included software and network credits before reuse.

### tADMM — IAS Transactions reproduction

Git tag: `260731-ias-trans-tadmm`. [Preserved source tree](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/tree/260731-ias-trans-tadmm). Published: 2026-07-31T19:09:54Z.

Repository state behind the IAS Transactions paper *"Scalable Multi-Period
Optimal Power Flow in Distribution Systems Using Temporal Decomposition"*.

**tADMM** decomposes the multi-period OPF along time into per-period branch-flow
OPFs coordinated by an augmented-Lagrangian consensus, solving horizons the
monolithic formulation cannot.

## Quickstart

```bash
git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
julia --project=envs/tadmm -e 'using Pkg; Pkg.instantiate()'

# cheapest cell showing a real speedup (~30 s)
SYSTEM_OVERRIDE=ieee2552C_1ph T_OVERRIDE=6 RHO_OVERRIDE=4000 \
  JULIA_NUM_THREADS=16 julia run_tadmm.jl
```

## Reproducing a specific result

Every reported cell — winning penalty `rho0`, BF time, tADMM time, speedup — is
in [`results/summary.csv`](results/summary.csv). Take a row's `rho0_winner` and
`T`, pass them as `RHO_OVERRIDE` and `T_OVERRIDE`.

| Paper name | `SYSTEM_OVERRIDE` | Buses |
|---|---|---|
| ieee123  | `ieee123C_1ph`  | 128 |
| med2522  | `ieee2552C_1ph` | 2,522 |
| large10k | `large10kC_1ph` | 10,321 |

| Cell | tADMM | Centralized baseline |
|---|---|---|
| `ieee123C_1ph`, any `T` | seconds | seconds — solution-quality check only |
| `ieee2552C_1ph`, `T=6..48` | 0.5–9 min | 0.6–12 min |
| `ieee2552C_1ph`, `T=144` | ~23 min | ~3.5 h — the 9.31× headline |
| `large10kC_1ph`, `T=48` | ~1.6 h | **fails after ~4.9 h** |

> ⚠️ **The `large10k, T=48` baseline is *supposed* to fail.** That is the paper's
> headline result: the monolithic solve exhausts memory during MUMPS
> factorisation (~10.5 GB) and never converges, while tADMM completes. If you run
> `run_bf.jl` on that cell expecting a solution, nothing is misconfigured — the
> reported "BF time" is the wall-clock at failure. Budget ~16 GB RAM.

## Notes

- Package versions are pinned by the tracked `envs/tadmm/Manifest.toml`.
- Times assume **16 threads**; tADMM parallelises over periods, so fewer threads
  scale proportionally.
- This repository is broader than tADMM — `envs/ddp` and `envs/multi_poi` are
  other methods, not exercised by these results.
- Prefer **Source code (zip)** below (~33 MB) over a full clone; the history is
  large for reasons unrelated to reproduction.


### PESGM2026-MultiSourceMultiPeriodOPF

Git tag: `251117-pesgm-multi-source-mpopf`. [Preserved source tree](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/tree/251117-pesgm-multi-source-mpopf). Published: 2025-11-18T00:52:35Z.

<img width="1358" height="307" alt="image" src="https://github.com/user-attachments/assets/57255ed2-1d4c-46a5-90a3-1a31ee579971" />

Go to `multi_poi_test.jl` and try it out for yourself!

This repo has two other algorithms currently in development for solving the MPOPF problem itself, but they're all segregated with their own environments.

### v1.0.0

Git tag: `v1.0.0`. [Preserved source tree](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/tree/v1.0.0). Published: 2024-12-01T22:58:12Z.

The point between completion and testing of the first temporally decoupled MPOPF algorithm and a stable temporally brute forced MPOPF algorithm.

