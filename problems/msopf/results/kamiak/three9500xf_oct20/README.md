# Three-substation window on ieee9500xf_1ph (service transformers modeled), Kamiak, 2026-10-20

The sensitivity check on secondaries: `ieee9500xf_1ph` is `ieee9500_1ph` with each of its 1,275
service transformers modeled (a line with its leakage impedance to its secondary bus, where its
loads sit); PV and batteries are identical. Stage `three_ieee9500xf`, the five ties between
feeders closed, T = 24, PV reactive power dispatched, 16 directions, bracket search (jobs
31322625-32, commit b4743437). One folder per job: `metadata.txt`, `run.log`, rows; merged into
`../../multiperiod/three_large_windows.csv`.

Against `ieee9500_1ph` (same directions, `three9500_oct20b`), every boundary point moves outward:

| case | mean | range |
|---|---|---|
| batteries idle | +1.71% | +1.59% to +1.84% |
| cyclic batteries | +1.27% | +1.15% to +1.48% |

The binding limit is unchanged (no-backflow at every edge). The transformers' losses add about
1% to each substation's import, so it takes slightly more angle to drive one to zero. Voltages at
the secondary buses fall to 0.950 pu at the cyclic edges (0.954 idle) without becoming the
binding limit. So the primary-only deck understates the window by 1.2-1.8%, never overstates it.
Cyclic directions took 47-178 min each (3,980 buses), idle 4-12 min.
