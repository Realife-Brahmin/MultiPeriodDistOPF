# Focused large10k T24 cutoff experiments

Both variants preserve July C_B, current periodic inputs/corrected loss model,
rho15000 -> fixed48,24 workers, MA57 single-thread, Gold6438Y+ hardware,
idle-battery/load-flow starts, free terminal SOC and full warmup/reset.
Exact reference31660055:49 iterations,50.449995s core,54.768778s loop;
central33.790131s. Maximum integrated energy-bound violation0.869624kWh;
strict original-model physical feasibility false. These are single trajectories.

## Screened cap completed

Source db1bb2bfddc2c2fc4386ac1d85a2a82534569c47; gate31663363 passed455 checks;
run31663364 and energy audit31663365 completed. Native collection and independent
local validation pass, including unchanged initial hashes, warmup/reset path,
phase/status exceptions, final exact round and saved network/energy audits.

Max30 local iterations in rounds2–30, accepted only if unscaled local primal
<=1e-6; exact from31. Result49 rounds /52.381487s core /57.306291s loop,
central34.354851s. Cost gap1.5439e-6; max integrated energy-bound error0.869624kWh
and0.065228 SOC percentage points, essentially unchanged. Strict physical flag
false. The measured core time is3.83% above the historical exact trajectory;
no repeat-median or statistical significance claim.

Seven capped solves were accepted. The two expensive period2 attempts at22/23
were rejected: primal violations5.76328e-4/1.33182e-2 at30 local iterations.
Exact recovery still consumed102/172 total inner iterations. Thus this screened
policy did not actually cut the key outliers. Local screening cost (roughly
0.17–0.26s per cap attempt) and all recovery work were counted.

Compact raw/source/log archive `early-cap-screened-completed.tgz` SHA256:
`d6a7c655ca170c88d3ea54251390a6221938a1495b58e07d4c47024c89dd1a86`.
Full exported instances/checkpoint arrays remain on Kamiak; `.jls` arrays are
excluded from the compact archive. `screened-summary.json` records individual
attempts and all quality flags. The July-style convergence plot explicitly marks
the exact-phase boundary; it does not turn the endpoint into a feasibility claim.

## Direct finite-iterate follow-up

Source4b5d07e5af2305912c17ca0b3e8e1b64028fd1b3. Gate31663478 passed472 tests,
including finite acceptance above the screened tolerance and exact restoration.
Run31663479 with passive audit31663480 is pending/running at submission.
This follow-up accepts finite capped primal/multipliers in rounds2–30 and logs
actual constraint errors instead of applying the1e-6 early screen. Exact solves
from31 and final quality criteria are unchanged. It is selected after inspecting
the rejected attempts, not an untuned heuristic. No completed result is claimed
until all numerical and independent validation records are available.
