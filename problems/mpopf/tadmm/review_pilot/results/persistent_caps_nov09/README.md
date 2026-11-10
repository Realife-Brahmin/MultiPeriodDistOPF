# Persistent local caps: completed large10k T24 comparisons

Frozen numerical source8d688457cbccd78d1af95818a095503988f0a5af. Gate31665844
passed501 assertions, including first-round and late-round caps with no exact
recovery or phase switch. Cap30 job31665887/audit31665888 and cap10 job31665889/
audit31665890 completed. Both native collection validators and independent
local revalidation pass: all source/input hashes, initial values, warmup/reset
numerical paths, worker coverage, capped local rows and allowed final statuses.

No exact cleanup was used. Every tADMM local solve, including first and final,
has its declared maximum30 or10; finite-only acceptance never triggers exact
recovery. Ordinary convergence can occur before the cap. At the final round
both arms happen to have zero truncated local solves: all finish naturally
within their limit. Fully solved centralized references are separate and never
initialize tADMM.

Both arms use MA57,24 single-thread workers, one Gold6438Y+ node with128GiB
requested memory, identical current periodic exports/corrected losses, July C_B,
free terminal SOC, shared load-flow/idle-battery starts, full warmup/reset and
rho15000 decreasing to48. Core includes starts/updates/communication/coordination.
All original-model diagnostics, passive SOC audits and OpenDSS remain outside
core. Persistent caps have no local feasibility decision/report inside each
stage. Prior early-cap experiments counted those local reports conservatively;
do not attribute all cross-policy timing differences to the cap alone. The
exact reference has no such timed local report, so this boundary matches it.

## Measured outcomes

Exact baseline31660055:49 outer iterations /50.449995s core /54.768778s loop;
central33.790131s; objective disagreement1.5439e-6; energy-bound error0.869624kWh
/0.065228 SOC percentage points; dynamics mismatch0.612496kWh. Separate strict
original-model physical feasibility was false.

Persistent30:52 rounds /48.454048s core /53.691627s loop, central32.485390s.
This is3.96% less core time than the single exact reference. Objective agreement
2.1488e-6 (0.000215%); max integrated energy-bound violation0.668225kWh /
0.036451 SOC percentage points; max dynamics mismatch0.671864kWh. There are14
truncations among1248 local solves. Local iteration median/95th/max4/23/30;
local solver time median/95th/max0.1897/0.9375/1.2290s.

Persistent10:53 rounds /37.262515s core /41.571119s loop, central33.793986s.
This is26.14% less core time than the exact reference and23.10% less than cap30,
but remains10.26% slower than its matched full centralized solve. Objective
agreement3.4308e-6 (0.000343%); max integrated energy-bound violation0.818639kWh
/0.050850 SOC percentage points; max dynamics mismatch0.341810kWh. There are244
truncations among1272 local solves. Local iteration median/95th/max4/10/10;
local solver time median/95th/max0.1940/0.4407/0.7234s. Power-bound violation0.

Both arms reach owned/all-copy RMS<=1e-5 and unique dual RMS<=1e-2. Neither
reaches physical NO (original maximum violation<=1e-5 plus0.1% cost agreement)
or the stricter driver flag. Thus the measured solver-work reduction is not a
fully feasible centralized speedup. These are single measured trajectories,
not repeat medians or a proven optimal cap/heuristic. Objective values are
slightly below the central reference; assembled dispatches remain inexact,
so objective agreement alone is not an optimality certificate.

All17 saved network records including central/endpoints converge with zero
voltage violations, maximum optimizer/OpenDSS voltage difference<=1.48e-10pu.
Normalized cone gaps are4.67e-7/4.68e-7. Network agreement is separate from
battery consistency. No thermal upper constraints are modeled by this pilot.

The exact baseline has9/1152 warmed local solves over30 Ipopt iterations
(0.78%), versus201/1152 over10 (17.45%). Thus30 only cuts rare outliers;10
changes a broader portion of local work and performs better here with similar
endpoint engineering errors. No generalization to other feeders/horizons or
smaller caps is established by these two trials.

## Evidence and preservation

`summary.json` includes exact endpoint flags, local distributions, truncation
counts and independent energy/network audits. `cap30-convergence` and
`cap10-convergence` PNG/PDF/CSV/JSON restore the July objective/r/s/rho view,
explicitly labeling persistent caps. Raw numeric curves are not smoothed.
Compact archive `persistent-caps-completed-nov09.tgz` SHA256:
`b5b9a65e4b01e1abf3a85a512a32ad3bce26f094b64f65f9456b87e2f1379cf4`.
Full serialized instances/dispatch checkpoints remain on Kamiak and are
excluded from this archive; source snapshots and logs remain included.

Pre-solve jobs31665845/31665847 failed only because generated Python bytecode
changed after manifest creation. All scientific source/input hashes matched;
only `__pycache__` entries were excluded in the repaired manifest. Original
manifest, failed logs and repair explanation are preserved. No numerical
outcomes were removed or relabeled. No more persistent-cap arms were launched.
