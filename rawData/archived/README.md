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

## Restoring one

```bash
git mv rawData/archived/<name> rawData/<name>
```
