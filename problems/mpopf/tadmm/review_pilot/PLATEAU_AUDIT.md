# T=24 plateau: local accuracy and redundant source-voltage bounds

The failed med2522 run is preserved at checkpoint iteration 5000. Its penalty
was fixed at 7600 after iteration 10. Do not attribute its persistent two-cycle
to ADMM tuning before accounting for the numerical-modeling defect below.

## Controlled evidence

MA57-only saved-state replay 31645686 holds rho=7600 and the network, objective,
initial state, consensus and dual update unchanged. Eight default steps reproduce
the odd/even objective and original-feasibility oscillation. The quadratic and
linear augmented-objective coefficients match their independently calculated
values exactly. The interior consensus dual-sum identity holds to about 1.4e-17.
Rewriting the objective as a complete square gives identical replay trajectories.

Default local solves report unscaled complementarity as large as 7.10e-6 even
with nominal tol=1e-9. Ipopt's tol is a scaled criterion; absolute dual and
complementarity limits are separate ([official options](https://coin-or.github.io/Ipopt/OPTIONS.html)).
Setting only compl_inf_tol=1e-8, with the same rho and algorithm, breaks the
cycle: after eight steps original-model violation falls to 6.11e-8 and all-copy
residual to 3.34e-8. Dual residual falls to 3.79e-4 and is still above its final
1e-4 threshold, so this is diagnostic improvement, not full convergence.
An adaptive-barrier accuracy variant gives similar damping and is not required
for the production fix.

## Source-voltage redundancy

The original model imposes v[root]=1.05^2 and also gives that variable an upper
bound 1.05^2. In period 7 at the preserved state, those two constraints have
opposing multipliers of magnitude about 2.658e9. This degeneracy distorts the
scaled local termination test. The rest of the sampled largest multipliers are
about 228 in magnitude.

Replay 31645728 removes only the redundant root voltage bounds, retaining its
equality and checking that 1.05 lies inside the original voltage limits. The
feasible set is unchanged. The huge multiplier pair disappears; the largest
sampled multiplier becomes about 228. Period-7 unscaled dual error falls to
5.75e-11 and complementarity to 2.27e-10. The plateau begins damping without
changing rho. Raw logs and multiplier lists are authoritative.

## Production change and verification scope

The common centralized/local model now excludes root voltage bounds already
implied by its fixed voltage equality, and explicitly sets compl_inf_tol=1e-8.
Original physical metrics still check root voltage and all original voltage
limits. No line-loss terms, profiles, costs, terminal convention, penalty rule
or ADMM update is changed. The centralized reference and local subproblems use
the same revised physical model. A fresh med2522 T=24 run must pass full warm-up,
measured restart, objective/physical checks and independent OpenDSS validation
before claiming an end-to-end fix. Checkpoint replay is never a fresh timing.

Fresh validation outcome: job 31645746 passed full warm-up and measured restart
in 389 iterations each, matching numerical traces, all original feasibility /
consensus / dual criteria and independent OpenDSS checks. Initial rho=3800 and
its single increase to 7600 after iteration 10 match the failed run. The
controlled tests and the fresh run therefore support the numerical-modeling /
local-accuracy explanation without a penalty-policy change. Full measured loop
time is 304.555781 s; centralized wall time 14.971893 s. Convergence is repaired;
this is not a claim that the retained penalty is optimal or performance is yet
competitive. See results/accuracy_fix_T24_31645746 for frozen evidence and plots.
Job 31645790 is the one corrected-source large10k T=24 follow-up, gated on that
successful validation. The earlier backend iteration counts are confounded by
the local-accuracy defect and should not be treated as a universal backend ranking.

The earlier 10k T=24 MA57 run has now completed all validation in 793 iterations;
its original-source evidence is retained. It is not silently relabeled as a run
of this revised source.

Diagnostic development attempts are preserved: 31645650 logged one valid step
then invalidated cached JuMP results by changing its log destination; corrected
instrumentation is in 31645686. Attempt 31645674 failed environment loading from
an incorrect depot path before any solve. Attempt 31645711 failed multiplier
inspection because a local array shadowed the dual function; 31645728 qualifies
JuMP.dual explicitly. These are audit-harness failures, not production ADMM
termination results. Their files are not overwritten or counted as converged.
