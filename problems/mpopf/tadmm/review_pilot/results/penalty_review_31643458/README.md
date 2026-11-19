# Penalty-policy diagnosis for med2522 and large10k

This report preserves completed validated runs under the frozen relative-residual
controller. It does not change rho or solver options in an ongoing trajectory.
The data archive contains the underlying traces, captured sources, validation
reports and a reproducible extractor. Serialized Julia instances and dispatches
remain on Kamiak. The PNG is a scientific diagnostic of the shared-node T=6
verification trajectories; its elapsed times are not benchmark timing claims.

## What the observed traces show

med2522 T=6 starts rho=950. With MA57, updates after iterations 5 and 10 raise
rho to 1900 and 3800; it then stays fixed and full convergence takes 77 iterations.
With MA97, it follows those same early updates, then increases after iterations
55, 60, 90 and 125, reaching 60800. It remains there until iteration 1425, halves
to 30400, and returns to 60800 after iteration 1450. Full convergence takes 1703
iterations. This trajectory reproduces in both completed benchmark repeats.
It is an observed poor-performing configuration, not an intrinsic iteration
requirement of tADMM or evidence of fabricated iteration counts.

For the repository's posthoc optimization criterion (0.1% absolute relative
objective agreement and original-model feasibility <=1e-5), MA57 first qualifies
at iteration 28 and stays qualified. MA97 first qualifies at 31, subsequently
leaves that tolerance, and stays qualified from iteration 92 through the end.
The original stricter stopping rules still require all-copy consensus <=2.5e-7,
dual residual <=1e-4 and original feasibility <=5e-7. These are different metrics;
none should silently replace the others. OpenDSS checks validate the final saved
dispatch, not every intermediate posthoc criterion-crossing point.

large10k T=6 starts rho=12. MA57 reaches rho=384 in five doublings and converges
in 132 iterations. MA97 reaches rho=192 in four doublings and converges in 568.
The posthoc criterion first qualifies at 59 / 116 respectively, but remains
qualified through the end only from 111 / 540. Thus a first crossing alone can
give a misleading impression of sustained convergence, especially on large10k.

The plots demonstrate a long residual plateau in the medium MA97 case and
different penalty sensitivity on large10k. They do not isolate the controller
from inner-Ipopt accuracy and backend effects. Saved-subproblem diagnostics in
the earlier evidence archive show that nominal local solver success can have
different unscaled complementarity accuracy. Do not conclude from these plots
alone that raising rho caused the entire plateau.

## Current rule and useful next comparisons

The initial rho heuristic is 2*C_B*base^2/dt: 950 / 12 for med2522 / large10k
at T=6, and 3800 / 48 at T=24. This is an energy-cost curvature heuristic, not a
proved optimum. Every five iterations, normalized residual balancing doubles
rho when relative primal exceeds ten times relative dual, or halves rho for the
reverse imbalance; scaled duals are rescaled. The current permitted range is
[1e-6,1e6]. Exact observed decisions are in penalty_changes.csv.

For the next explicitly selected calibration, useful controlled comparisons are:

- Hold the initial rho fixed, to separate initialization from update behavior.
- Compare the current adaptive policy with a capped or frozen policy after the
  early adjustment phase, keeping the starting rho and all other settings fixed.
- Only then vary the initial rho around its current scale. Keep inner accuracy,
  profiles, weights, initialization, CPU policy and validation requirements fixed.

No such new sweep has been launched by this report. Report near-optimal first
crossing, sustained qualification, full convergence and failed local solves
separately, with the complete rho trajectory and wall time for every candidate.

## Snapshot status

Both med2522 T=6 benchmark repeats are included and validated with both backends.
Both large10k T=6 shared-node verification runs are included and validated;
their benchmark repeats are queued. med2522 T=24 runs are active. The first
IEEE123 T=24 MA97 attempts failed with LOCALLY_INFEASIBLE on a local subproblem
and are not represented as completed results. This is a solver termination
status, not a proof that the full benchmark optimization problem is infeasible.

The benchmark task concurrency cap was increased from two to three to allow
10k coverage alongside med2522 T=24, with unchanged numerical settings and
allocation policy. Scheduler priority still controls actual start times.

Archive SHA256:
`0da1df92b13cdc6e8f9257553a2542915f9bff814d47e4818114b30bb1517473`
