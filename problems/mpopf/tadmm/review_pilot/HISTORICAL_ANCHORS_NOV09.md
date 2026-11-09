# Historical anchor cases, November 9

User selected large10k T24 and med2522 T96 to focus on cases with successful
July references. The archived summary records MUMPS centralized/tADMM-NO times
of 4332.0/2264.6 seconds and 1622.20/738.30 seconds respectively, with winning
initial rho 15000 and 20000. That summary does not contain outer iteration
counts. Do not infer iteration counts from historical times or substitute the
large10k T48 centralized OOM time for a completed reference. There is no July
large10k T96 result in that summary.

Fresh exploratory comparisons use MA57 and MUMPS for both algorithms at each
case, one worker per period (24 or 96), same current modified feeders and
periodic profiles, July C_B, no terminal SOC requirement, corrected nodal line
losses, exact local solves, full warmup/reset and shared load-flow/idle-battery
initialization. No optimized initial battery schedule is used. Both backends
use identical data, screening, starts, rho policy and quality definitions.

Use one initial rho per case, from the July summary: 15000/20000. The current
improved bootstrap-fixed policy divides rho by five every five iterations
to floor 48, then holds it fixed. Transfer of that floor to these cases is an
experiment, not an established tuning rule or a literal July-controller replay.
No tuning sweep is hidden in this comparison. All 1000-iteration-budget or
solver failures must be preserved. Cases are submitted sequentially to limit
resource demand; native gates run first.

Continue to the owned AND all-copy RMS <=1e-5 and unique dual <=1e-2 endpoint.
Also report first historical residual stops at 1e-3 and 1e-4, reconstructing
the five-objective stagnation/dual alternative from the unsmoothed traces.
Historical-style residual plus cost milestones are distinct from that stopping
rule. Report physical errors and independent OpenDSS validation separately.
Core, loop, build/initialization and warmup costs remain separately logged.
These are current-model comparisons at historical horizons, not identical-input
replays of July and not repeated-median paper benchmarks.

Remote campaign:
`/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_historical_anchors`.
Reproduction launcher: repository-root `setup-historical-anchors-nov09.sh`.

The existing T48 MUMPS job 31660026 is retained. Pending large10k T96 job
31660028 and dependent audit 31660029 were cancelled before running at this
change of focus; they are not failed or completed numerical comparisons.
