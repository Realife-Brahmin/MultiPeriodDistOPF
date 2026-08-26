# Archived feeder data

Feeder models and data files that no code in this repository currently
references. They are kept, not deleted: several are intended for future work
(the 730-bus systems in particular).

Nothing here is needed to reproduce the IAS-Trans tADMM results. The systems
those results use — `ieee123C_1ph`, `ieee2552C_1ph`, `large10kC_1ph` — remain in
`rawData/`, alongside the other systems still referenced by `envs/ddp` and
`envs/multi_poi`.

## Contents

| Directory | Notes |
|---|---|
| `ieee730_1ph/`, `ieee730_1ph_rahul/`, `ieee730_1ph_basecase_OpenDSS_rahul/` | 730-bus feeder variants, earmarked for future work |
| `ieee729_1ph/` | 729-bus variant |
| `ieee123_1ph/` | earlier IEEE 123 variant, superseded by `ieee123C_1ph` |
| `8500-Node/` | IEEE 8500-node test feeder |
| `ads3_1ph/`, `ads10_1ph/` | small ADS test cases |
| `case533mt_hi/` | 533-bus case |
| `Data_Anup.xlsx` | source spreadsheet |
| `LoadShapePSubsCostDefault0.dss` | superseded load-shape variant |
| `ieee123A_1ph/` | earlier IEEE 123 variant, superseded by `ieee123C_1ph` |
| `ieee2552_1ph/` | earlier IEEE 2552 variant, superseded by `ieee2552C_1ph` |
| `ieee2552_rahul/` | original MATLAB/Excel source data (base case + 10/30/50% PV power flows) that `ieee2552_1ph` was converted from |
| `large10k_1ph/`, `large10kB_1ph/` | earlier 10k-node variants, superseded by `large10kC_1ph` |
| `convert_all_areas_10k.jl`, `convert_flexparams_10k_to_opendss.jl`, `convert_ieee2552_to_opendss.jl`, `convert_simple_10k.jl` | conversion scripts for the superseded variants above; no longer needed since the `C` variants are checked in directly |

## Restoring one

```bash
git mv rawData/archived/<name> rawData/<name>
```
