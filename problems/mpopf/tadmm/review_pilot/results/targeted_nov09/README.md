# Targeted T48 experiments

Source4a074c167403ba9decb31a578a8df4de1e47f0c0, gate31657404 passed425 tests:
174 forced-interruption/exact-recovery checks,3 fixed-penalty schedule checks,
and248 existing repeated-solve/objective/multiplier checks.

## Completed T48 comparisons

Both native collection validators and independent local revalidation pass:
source hashes, warmup/reset numerical paths, worker coverage, exact final local
statuses, milestone definitions and OpenDSS checks. Passive energy audits also
completed. These remain single measured trajectories, not repeat medians.

- Exact recovery at 20 iterations: 322 outer iterations, 344.1063 s core,
  379.5507 s loop. It triggered 208 cold retries and discarded 4,160 warm-attempt
  iterations; all work is timed. Maximum total local iteration count fell from
  144 to 37, but total time did not beat the 334.2166-s baseline. No timing gain
  is claimed for this candidate. Final integrated energy error is 0.35764 kWh.
- Fixed-48 schedule: 129 iterations, 160.1589 s core, 176.7187 s loop, versus
  322 / 334.2166 s for the legacy controller. Same final owned/all-copy RMS
  <=1e-5 and unique dual <=1e-2 goal; final cost gap 3.84e-6. The matched full
  centralized solve is 90.6251 s, so the tighter goal still takes longer than
  centralized. Energy overshoot is 0.57835 kWh / 0.04338 SOC percentage points.

The looser owned-RMS <=1e-3 plus 0.5% cost milestone occurs at iteration 42 /
59.3809 s, below centralized time, but energy error is 25.5262 kWh / 1.56149 SOC
points. This is a relaxed historical-style milestone, not an equal-physical-
quality speedup. Every-copy RMS <=1e-3 plus 0.1% cost occurs at 66 / 102.3090 s,
with 18.2865 kWh / 1.09710 SOC points. At <=1e-4 plus 0.1% cost, it takes 94 /
128.6625 s. All physical flags remain separate; even the final small engineering
error does not meet the separate original-model 1e-5 strict physical threshold.

See [comparison figure](analysis-T48/implementation-comparison.png),
[PDF](analysis-T48/implementation-comparison.pdf), [all metrics](analysis-T48/summary.json)
and [final SOC trajectory](T48-fixed-SOC.png).

## T96 transfer check running

Job 31657825 extends only the fixed-48 schedule to T96 with 96 workers on two
matched Gold nodes, initial rho 240000 and exact ordinary local solves. Same
source, model, weights, targets and full warmup/reset. No retry policy is mixed
in. Audit 31657826 follows. The candidate and initial rho are recorded before
T96 results; this is an exploratory horizon-transfer check, not a universal
rho heuristic or an untouched-data validation claim.

Completed T48 evidence and the predeclared T96 extension are preserved in
targeted-T48-completed.tgz, SHA256:
2bc1a09ea36e20dcef5627ae2e418fa59ee1bcea9fc33b17635f7c3ad46cbc61.

## Execution history

Two separate full-warm/reset MA57 comparisons were submitted:

-31657405: legacy rho controller with exact primal-only recovery after a warm
 attempt reaches20 Ipopt iterations. Both attempts' work is timed and counted.
-31657407: ordinary exact local solves with high initial rho120000 ramped
 down by5 every5 iterations to48, then frozen; no later penalty cycling.

Both use48 workers, same modified large10k T48 data/weights/physical model,
Gold node64-CPU/128-GiB allocation, same stopping/quality checks and max1000
outer iterations. Passive energy audits31657406/31657408 follow the jobs.
Compare each to the334.2166-s baseline before combining changes or extending
to T96. No performance gain is claimed from the gate or submission.

See [predeclared protocol](../../TARGETED_NOV09.md). Exact setup and native gate
output are preserved in targeted-initial-protocol.tgz, SHA256:
a3407f265ea5b16fbd095d1318d48259a71cf8fcd6bdfa090e92ea72b1bb6b0c.
