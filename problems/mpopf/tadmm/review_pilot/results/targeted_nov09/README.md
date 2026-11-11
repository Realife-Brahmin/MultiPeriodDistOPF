# Targeted T48 experiments

## November 10 persistent cap10 follow-up

Kamiak job 31679860 started on cn168 with all 64 physical cores allocated
and 48 time-step workers. Campaign path:
`/home/aryan.jha/kamiak_jobs/tadmm_review_nov10_large10k_t48_cap10`.
The frozen solver source is b1ba396977702b264eef679ce981f4881c15f787,
with the existing cap and chronological replay gates and verified input and
execution-source hashes. The run uses MA57, the fixed-48 bootstrap schedule
starting at rho 120000, a persistent 10-iteration cap on tADMM subproblems,
and a fresh fully converged centralized reference on the same allocation.
Jobs 31679861 and 31679862 audit energy feasibility and explicit cone,
chronological storage and RAM evidence after the solve. No results are
claimed at submission.

### Exact comparison before selecting another cap

Job 31679947 repeats T48 with exact local solves and the same frozen source,
hardware class, penalty schedule, centralized reference procedure, and
validation as job 31679860. The campaign directory ends in
`tadmm_review_nov10_large10k_t48_exact`. Audit jobs 31679948 and 31679949
follow it. The comparison changes the local solve policy only.

The prior exact fixed-48 trajectory (31657407) contains 6,192 measured
subproblem solves over 129 rounds. Ipopt iterations were 4 at the median,
17 at p95, and 137 at the maximum. A cap10 would have cut 762 solves and
the longest elapsed subproblem in 106 rounds. A cap20 would have cut 163
solves and the longest subproblem in 39 rounds. The final round had a
maximum of 11 iterations. These counts describe the exact trajectory and
cannot predict the runtime or convergence of a capped trajectory.

The fresh cap10 measured trace reached the consensus target after 176
rounds, 154.70 s core and 175.98 s loop, versus 87.81 s centralized.
Three last-round subproblems hit the cap. Its original-model violation
was 8.12e-4 and its maximum raw cone slack was 1.52e-4 per-unit squared.
At the last iteration before half the centralized core time, the all-copy
RMS residual was 2.85e-3 and the original-model violation 1.83e-2.
The earlier first crossing of 1e-3 had an objective about 6.8% above
centralized, so the residual alone cannot select the stopping time.
These are measured solve results pending completion of independent audits.

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

## Completed T96 transfer check

Job 31657825 and passive audit 31657826 completed; native and independent local
validators pass, including source hashes, warmup/reset trajectory agreement,
96 unique period workers on two hosts, exact final local statuses and all 12
saved-point/endpoint OpenDSS records (zero voltage violations, maximum optimizer
voltage difference <=1.50e-10 pu).

The fixed-48 schedule reaches the tight RMS/dual goal in 379 iterations /
474.1936 s core (555.3876 s loop), versus 849 / 1073.4769 s with the legacy
controller at the same initial rho240000. This is a 2.26x core-time improvement.
Final cost gap is7.41e-6, integrated energy overshoot0.80747 kWh /0.04037 SOC
percentage points. The separate strict original-model feasibility flag remains
false. The matched fully converged centralized reference is188.5880 s, so the
tighter-goal tADMM result is still2.51 times slower than central.

Relaxed milestones retain all quality flags:

- Owned RMS<=1e-3 plus0.5% cost:64 iterations /98.6301 s, but25.9641 kWh SOC
  error /1.29821 SOC points; all-copy RMS2.03e-3 and dual RMS0.275.
- Every-copy RMS<=1e-3 plus0.1% cost:92 /144.7510 s, below the full central
  time, but38.8749 kWh SOC error /2.91591 SOC points and dual RMS0.0987.
- Adding dual<=1e-2 to the1e-3 RMS test:152 /216.8124 s, cost gap0.0738%,
  energy error5.0857 kWh /0.30512 SOC points. Still slower than central.
- Every-copy RMS<=1e-4 plus0.1% cost:188 /258.2386 s, energy error5.5845 kWh.

Thus the earlier time advantages are qualified relaxed-target results, not
equal-physical-quality speedups. Per-iteration core time at the final goal is
about1.25 s for T96 versus1.24 s for T48; most remaining horizon penalty is in
outer iteration count rather than the cost of one parallel iteration.

See [T96 comparison](analysis-T96/implementation-comparison.png),
[PDF](analysis-T96/implementation-comparison.pdf), [metrics](analysis-T96/summary.json)
and [SOC trajectory](T96-fixed-SOC.png). Complete evidence is in
targeted-T96-completed.tgz, SHA256:
f91cccfaa061f4f9a27245dbfb140adcbbe0bdd08a93b27d1ade87d3a04acd8f.

### Transfer-test history

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
