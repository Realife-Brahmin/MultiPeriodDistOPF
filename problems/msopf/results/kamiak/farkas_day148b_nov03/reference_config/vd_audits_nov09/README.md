# sf6937 VD audits of the last six CD directions (2026-11-09)

These audits close the gap described in `../REPORT.md` ("72 of 78 are exact-checked"). With them,
**all 78 sf6937 rated whole-day battery certificates (CD) are exact-checked: 59 in-run and 19 by VD**
(13 in array 31545814 and 6 here). Every VD verdict is CERTIFIED.

## Results

| dir | r_cert (deg) | float status | exact check | verdict | VD wall |
|---|---|---|---|---|---|
| 2  | 5.259779 | INFEASIBLE / INFEASIBILITY_CERTIFICATE | Farkas 52936.9 | CERTIFIED | 1:19:58 |
| 11 | 2.148050 | SLOW_PROGRESS / OTHER_RESULT_STATUS | objective -355.402 kW | CERTIFIED | 2:07:53 |
| 14 | 1.520360 | INFEASIBLE / INFEASIBILITY_CERTIFICATE | Farkas 51147.5 | CERTIFIED | 1:49:29 |
| 16 | 4.810860 | INFEASIBLE / INFEASIBILITY_CERTIFICATE | Farkas 51949.6 | CERTIFIED | 1:27:20 |
| 17 | 4.454110 | INFEASIBLE / INFEASIBILITY_CERTIFICATE | Farkas 51426.1 | CERTIFIED | 1:57:53 |
| 19 | 1.299720 | INFEASIBLE / INFEASIBILITY_CERTIFICATE | Farkas 51130.8 | CERTIFIED | 1:52:42 |

Direction 11 has no Farkas ray (the float solve stopped on slow progress, as its CD run did); it is
certified by the exact objective bound being negative. Direction 17's CD run hit the iteration
limit; its VD solve returned a Farkas ray. The radii are the CD radii; VD re-solves at them.

## Provenance

- Slurm array **31667126** (tasks 2, 11, 14, 16, 17, 19), partition `kamiak`, account `dubey`,
  8 CPUs, 80 GB, 12 h limit; submitted and started 2026-11-09 20:06:13 PST, all exit 0:0.
  Per-task step IDs in the `.meta` files (31667138-31667142, 31667126). Nodes cn149, cn153 (x3),
  cn169, cn178. Peak RSS 43.6-44.9 GiB (`sacct_31667126.txt`, MaxRSS in KiB).
- Commit **865e26d19a8a1d35f0b72881ce63e83ff61ed280** (checkout `~/msopf/repo_rated`). The earlier
  VD array 31545814 ran 4f29d098; the diff between them only touches `FB_DERIVED`, which is off
  for VD, so the VD code path is identical.
- Environment `sdp_pardiso` (same Project/Manifest hashes as `../REPORT.md`), Clarabel with MKL
  Pardiso (`FB_LDL=mkl`), exact rational check of the returned certificate.
- Submit command (from `~/msopf/jobs/farkas_day148b_nov03`):

  ```
  FB_SCRIPT=farkas_dispatched.jl FB_ENV=sdp_pardiso FB_LDL=mkl NTH=8 FB_CUTS=kcl0,jbound,lsec \
    sbatch --parsable -p kamiak -A dubey -c 8 -t 12:00:00 --mem=80G --array=2,11,14,16,17,19 \
    -J kVD_sfp1u_four -o logs_kVD_sfp1u_four_%a.out farkas_job.sh sfp1u_four VD
  ```

## Files

| file | where |
|---|---|
| `fb_VD_{2,11,14,16,17,19}.csv` | `../../sfp1u_four/` beside the 13 earlier `fb_VD_*.csv` (new, uncommitted) |
| `logs/VD_*.log`, `.meta`, `.err` | here (copied from Kamiak `~/msopf/jobs/farkas_day148b_nov03/sfp1u_four/logs/`) |
| `sacct_31667126.txt` | here |

The Slurm stdout files `logs_kVD_sfp1u_four_*.out` on Kamiak are empty (output goes to the logs above).
The `.err` files hold only the known OpenDSSDirect 0.9.9 precompile warnings/error on Julia 1.12;
the package loads and runs.

SHA256 of the result CSVs (identical on Kamiak and here):

```
8beb460a9e9fff4b30fce0683d8486de750545ca9a90603688e6878878210a85  fb_VD_2.csv
06fe7f313d2ddf82e8380eb061ef14ea43e017f5b18bf9e46e135ccb4e53aeae  fb_VD_11.csv
0388f92eec3a1df7d6fb658ad6d3cdd2fd9bd3f6ec4cd0e127a2120792cd8cdd  fb_VD_14.csv
55b9fab9cb935f200a75f63ed8a74299d807034caf532561ba45171d03ac9d5a  fb_VD_16.csv
a266848cf88d32f30cfe77ff060aff4af21b13cb65686660e257103c69e0e064  fb_VD_17.csv
8be84824bfd6c766b7cc54d9a1c7ad78b7f3f58a977cb51ec17728a06de90439  fb_VD_19.csv
```

## What still says 72 of 78

Not changed by these audits: `../REPORT.md` (lines near "72 of 78" and the VD resource row),
`cd_exact_status.csv` inside `../sf6937_kamiak_reference_2026-11-03.tar.gz`, the paper, and the
maps page. They need updating to 78 of 78 (59 in-run, 19 by VD) if the release records these runs.
