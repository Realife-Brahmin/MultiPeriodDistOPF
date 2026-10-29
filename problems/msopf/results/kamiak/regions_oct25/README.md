# Regions on the decks with service transformers, Kamiak, 2026-10-25/26

The paper's region method on the two decks with service transformers, the 0.95-1.05 band at customer buses only
(2026-10-24 rule). This replaces the matched campaign of 2026-10-25 (`matched_secondary_oct25`), whose fixed
360 s cap left about half the directions unresolved.

* **Stages** (`root_level/multiperiod_study.jl`):
  * `four_smartds_p1u_typeiii`: San Francisco P1U, its four substations, Type III ties closed.
  * `four_ieee9500xf`: the 2705-bus feeder with its 1,275 service transformers, three substations plus one
    added at l2804256.
  * `three_ieee9500xf`: the same feeder, its own three substations.
* **Settings:**
  * T = 24, PV reactive power dispatched, batteries idle and cyclic.
  * Directions: 78 (fib78) for four sources, 16 for three.
  * Edge search: import-margin secant with `MULTI_VMARGIN=20000`, falling back to bracketing.
  * Time caps: the reference solve uncapped, later solves capped at 8x its time, each retried once at 4x.
  * MA57, one thread per task.
* **Run:** commit `4a5856e6`, job arrays 31420559 (28 idle tasks) and 31420560 (86 cyclic tasks), all 114
  COMPLETED. Job script [`../../../kamiak/regions_job.sh`](../../../kamiak/regions_job.sh); task lists
  `tasks_fast.txt` and `tasks_cyclic.txt`. One folder per task (`metadata.txt`, `run.log`, `run.err`, rows).
* **Merged** into `../../multiperiod/multi_windows_fib78.csv` (merge_multi.jl) and
  `../../multiperiod/three_large_windows.csv`. The latter replaces the 2026-10-20 ieee9500xf rows, which banded
  every bus.

## Results

No direction is unresolved. "At 20 deg" counts directions still feasible at the search limit; those are kept out
of the ratios.

| System | Batteries | Directions | At 20 deg | Measured / predicted: median (range) | Below prediction | Set by |
|---|---|---|---|---|---|---|
| San Francisco P1U, 4 sources | idle | 78 | 11 | 1.20 (1.03-1.30) | 0 | no-backflow 53, voltage 10, none 4 |
| San Francisco P1U, 4 sources | cyclic | 78 | 25 | 1.19 (1.13-1.27) | 0 | no-backflow 38, voltage 13, none 2 |
| 2705-bus with secondaries, 4 sources | idle | 78 | 0 | 1.05 (1.03-1.38) | 0 | no-backflow 78 |
| 2705-bus with secondaries, 4 sources | cyclic | 78 | 0 | 1.05 (0.95-1.40) | 6 | no-backflow 78 |
| 2705-bus with secondaries, 3 sources | idle | 16 | 0 | 1.09 (1.06-1.10) | 0 | no-backflow 16 |
| 2705-bus with secondaries, 3 sources | cyclic | 16 | 0 | 1.08 (1.04-1.11) | 0 | no-backflow 16 |

* **Band at customers only.** With it, voltage no longer cuts inside the prediction on the 2705-bus feeder:
  the primary-only deck had 12 idle points below it, down to x0.35. Only 6 cyclic points fall below, by at most 5%.
* **Three-source edges** move by at most 0.01 deg against the all-bus band of 2026-10-20: all are set by
  no-backflow.
* **San Francisco** is never inside its prediction (x1.03 at worst). Its substation-4 face is barely sampled:
  substation 4 imports the most (11.6 MW mean), and most directions toward it are still feasible at 20 deg
  (user, 2026-10-26: 20 deg is already far beyond what is needed; no rerun).
* **Voltages.** Customer voltages reach 0.950 pu at many edges, together with no-backflow. Where voltage alone
  binds, the edge is still beyond the prediction.
* **Time per direction:**

  | System | Idle: median (max) | Cyclic: median (max) |
  |---|---|---|
  | San Francisco | 13 min (36 min) | 44 min (4.9 h) |
  | 2705-bus feeder | 3 min (4 min) | 14 min (56 min) |

* **Currents** are not recorded here: the current tabs (`eac3c458`) came after this campaign.
