# Kamiak period distributions and execution overhead

The [period boxplots](period-distributions.png) ([PDF](period-distributions.pdf))
use the corrected-model July-weight T24 exact-local assessments, jobs 31650588
and 31650589. Each plot uses outer iterations 2-200; the aggregate accounting
below uses all 200 iterations. Raw source archives are indexed in the
[July-weight assessment](../july_weight_comparison/README.md).

Medium local Ipopt counts have median 26, p95 33, maximum 37. Periods 18-20
have medians 34, 32 and 33 and native median times about 0.31-0.32 s. Period 12
has median 19 steps and 0.196 s. Large counts have median 11, p95 16, maximum
24; periods 5-6 have median 16/15 steps and 0.702/0.657 s. Easier large periods
are around 0.46-0.50 s. These are descriptive observations, not a causal
attribution to a specific load/PV event.

For medium, summed per-round maximum native Ipopt time is 66.56 s; summed
maximum local-stage time is 74.96 s; core algorithm time is 80.15 s. The core
remainder outside the maximum local-stage durations is 5.19 s, or 6.47%.
For large the corresponding values are 155.47, 184.48 and 202.03 s, leaving
17.54 s, or 8.68%. Median per-round remainders after the first round are
0.0218 s medium and 0.0823 s large.

**That remainder is not pure communication.** It combines process dispatch,
serialization, worker start skew, scheduling and consensus/dual calculations.
Communication that overlaps worker execution is not isolated by subtraction.
There are no network-only or inter-node MPI measurements in these results.

All periods run concurrently in 24 Julia processes on one Gold 6438Y+ node,
with a separately pinned coordinator. Each process has one Julia/BLAS/OpenMP
thread. The allocation reserves 64 physical CPUs, but the active configured
compute budget is 25; it is not a 64-worker or inter-node experiment.

The native workload balance ratio, sum of all native solve times divided by
24 times the sum of the per-round slowest native times, is 0.800 medium and
0.705 large. This quantifies load imbalance, not measured CPU utilization or
speedup over an independently timed sequential implementation. The recurring
slow periods make straggler handling relevant, especially medium periods 18-20.

Physical validation and full dispatch retrieval are outside the core solve
timer. Native Ipopt timers are separate again. Total loop wall time includes
those diagnostics, logging and checkpoint work and is not substituted for
core solve time. The new diagnostic implementation changes neither the
frequency nor the mathematical checks; its speedup is reported separately.

[Summary and per-period statistics](summary.json) and source hashes accompany
the figures. Repeat medians and fully warmed final benchmark claims are not
made for these exploratory assessments.
