# Three-substation window on ieee9500_1ph, Kamiak, 2026-10-20

Stage `three_ieee9500`: the deck's own three substations, the five normally-open switches
between feeders closed, T = 24, PV reactive power dispatched, batteries idle (`nobatt`) and
cyclic, 16 directions, edges by feasibility tests bracketing the predicted reach. Jobs
31321352-60 at commit e0eaaa1d (switches as lines, substation impedance inside the Vsource;
electrically the same circuit as the later explicit-transformer layout). One folder per job:
`metadata.txt`, `run.log`, rows. Merged into `../../multiperiod/three_large_windows.csv`.

All 32 edges are no-backflow (smallest import 0.00 kW) with voltages 0.953-1.045 pu, and lie
3-11% beyond the region predicted from the angle-to-import sensitivities (median +8%), idle
and cyclic alike: the prediction is conservative on this system throughout. Cyclic batteries
widen the reach by 1.29-1.36x. Idle directions took 4-8 min each, cyclic 46-125 min.
