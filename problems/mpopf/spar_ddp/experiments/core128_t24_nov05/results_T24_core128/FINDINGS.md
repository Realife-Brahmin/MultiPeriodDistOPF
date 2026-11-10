# Final matched128-core T24 result

All six original large10k1020-battery solves completed and passed strict
validation on the SAME Xeon6760P128-core nodecn178,600GiB allocation.
Forward/reversed two-repeat medians:

| Method | Full optimizer,s | Feasible0.1% target,s | Full range,s |
|---|---:|---:|---:|
| Groupedb2/right24x2 DDP |13.6840|7.2454094|13.684-13.684|
| MA97 |17.7835|9.6831455|17.679-17.888|
| MA57 |18.4775|9.8584685|18.464-18.491|

Against the faster central BACKEND MEDIAN (MA97 for both clocks), DDP is
**1.2996x at full convergence and1.3365x at feasible near-optimality**.
The target is objective gap<=0.1% AND original primal<=1e-5, not objective
alone. Observed target-time ranges are7.2321-7.2587s for DDP versus
9.6028-9.7635s for MA97, so no third paired repeat is needed. No interpolation
or fastest-single-run selection is used. Full-time reduction versus MA97
is23.05%; do not confuse that with the1.3365x target-speedup ratio.

DDP median A1.678795s,B9.4873795s,C1.493787s; other optimizer work1.0240385s.
Both DDP runs have14iterations/accepted steps,23backward passes,
8barrier reductions,0backtracks and status0. Final state/control arrays
match bit-identically across repeats. Each also passes the unchanged
same-iteration/1e-10 scaled state+dispatch comparison against the original
usual-node T24 reference: state2.221e-15, control3.361e-16 maximum scaled
differences. All physical batteries/constraints/coupling remain intact.

Every run passes original-model physical residual/bounds and independent
OpenDSS validation with no voltage violations over247680bus-periods under
the existing1e-4 tolerance. DDP independent original primal max6.822e-13;
central maxima4.867e-9. DDP independent objective2682722.226822585,
central reference2682722.010619225: relative agreement~8.06e-8, well below
0.1%. The central references mutually agree to floating-point precision.
Original numerical dual/complementarity checks, finite full objective/
feasibility trajectories and matched physical starting objectives pass.

Actual DDP settings are24workers,2BLAS threads per worker,1coordinator
thread/core, with49physical cores in use within the128-core allocation.
Every recorded OS-thread mask matches disjoint same-socket core groups;
phase/merge assertions and before/after audits confirm runtime settings.
Central BLAS/OpenMP stays1 under the inherited matched protocol. This
is not a centralized threading sweep or a global best-hardware claim.

Full times cover the warmed optimizer, including its ordinary overhead,
not model construction/full warmup/validation. DDP build54.855-55.406s,
fullwarmup76.402-78.464s, total process235-264s remain separately logged.
Worker historical peak RSS1.65GiB-scale includes shared mappings; summed
worker peaks are not simultaneous private memory. Batch MaxRSS38689616KiB
is the whole-job historical high-water mark. A shortT24 solve cannot exploit
all128cores through one worker per stage, and no cross-generation speed
difference is attributed solely to core count.

Job31576169 COMPLETED0:0 in15m29s. Source/harness/input hashes pass;
stderr is empty. The additional independent posthoc original-input recheck
is explicitly timestamped in posthoc/input-provenance.txt; it did not alter
running source or any scientific threshold. The LF-only preparation repair
and earlier setup failure remain preserved separately. No failed optimizer
or scientific gate relaxation occurred.

Raw logs, full trajectories, settings, audits/masks, hardware, HSL hashes,
source/input manifests, accounting, strict solution comparisons and solution
hashes are in dense_merge_core128_T24_final_evidence.tgz. Binary solutions
remain immutable on Kamiak. Reproduce with collect.py ROOT OUTPUT
--mode core128t24. No active source/default, published benchmark row, paper,
MSOPF or dispatch-aggregation change was made. This bounded campaign is complete.
