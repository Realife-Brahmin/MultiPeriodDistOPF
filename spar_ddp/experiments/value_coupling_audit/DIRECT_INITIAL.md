# Bounded direct initial-factor investigation

Final: results_direct_initial_oct23/FINDINGS.md. Direct generation passes for
the available period1 record; required chain2:4 elimination records are missing.
Their discovery and true-future discovery remain charged. Validated zero-JIT
three-pair partial-pipeline total is2.59433s vs1.33727s dense, not a full direct
pipeline test or optimizer result. No further capture/campaign was launched.

Base cfff0291. Same immutable original1020 large-late four-period window;
no optimizer/new trajectory/capture/medium/longer-horizon extension.

## Availability is an explicit limitation

Capture code fc_cache_stage/fc_write_stage saved full S only at1,12,24.
All24dense quadratic leaves exist, but they cannot uniquely recover omitted
network factors or local primal/dual Newton states. Local nominal in stage1
packet is NOT a full24stage trajectory. Required chain leaves2,3,4 therefore
cannot use a verified direct generator from these immutable records. Do not
substitute stage1 numerical factors into other stages or assume iterate reuse.
The incoming true future at5 is a composition of the remaining20periods, whose
elimination factors are likewise unavailable. Its dense-input discovery remains
explicitly charged; no zero/guessed/free boundary or new capture is authorized.

Smallest justified path: generate period1's leaf directly, within this SAME
window, and feed its actions into the initial-state/battery consumer and closed
period1 value probes. Leaves2:4 retain charged existing discovery, then use the
already validated compact chain. This is a PARTIAL direct pipeline, not removal
of the1.42s discovery obstacle for the three unavailable leaves plus future.

## Direct algebra (no dense boundary discovery)

Recreate only local Dg/energy-elimination blocks from S.K/nzc, Bd, PN/NP/Ninv.
The raw canonical battery Schur is block-local minus U Mroot^-1 U'. Scatter
using b0.bat, preserving separate raw left/right factors, not symmetrization.
Apply the dynamics scaling p=1/f. Build sparse local pieces and separator
factors for each boundary block:

    Mbb = A + Uc Vc'
    Mab = -A - R - Uc Vc'
    Maa = A + R + R' + D + Uc Vc'

R=Rhat' Diagonal(p); vectors use the original rhat/d formula. The generator
never reads existing dense S.b0.Sh0 or constructs global dense S0/Maa/Mab/Mbb.
Local blocks and sparse cross interfaces are retained; QR/small-core SVD may
recompress their KNOWN factors at1e-12. This is not an SVD-free algorithm:
no full dense input cross-block SVD is needed for the directly generated leaf.
Ranks/interfaces are measured, not fixed. Dense cross-block fallbacks remain
explicit and charged; existing original-coordinate upper completion is used
only when production uses it. No changes to physical variables/constraints.

## Gates and accounting

Hash-check actual source/input/capture/stage data. Verify actual saved stage
inventory and all existing ordinary/production, saved full-probe and independent
consumer gates. Compare every generated block/vector to immutable raw leaf1;
record raw skew separately (not an arbitrary symmetrization). Poison the saved
dense Sh0 with NaN and require generated output unchanged: a read-dependency
negative test. Compare chained two-merge/closure values and composed period1
actions/vectors to production. Recover battery/full-control/equality/bound
directions with unchanged1e-7 scaled AND1e-5 relative/floor1e-10 gates. Corrupting
the composed gradient must fail full direction checks.

Three alternating warmed SAME-input dense/partial-direct pairs in one process.
Exact helper/type paths and condition probes are warmed explicitly; every
phase records @timed compile_time/recompile_time, both required ZERO by strict
collection. Failed logging attempt and mixed-compilation pilot are preserved
and excluded, not repaired in place or promoted as warmed timings.
Charge leaf1generation, unavailable leaf2:4 discovery, incoming future discovery,
QR/recompression, chained updates, condition probes, period1value actions,
consumers and cleanup/in-phase GC; record allocations. Dense leaf1generation
uses _pp_elem with its existing saved dense Sh0; explicitly report this given
input distinction. Saved stage/network factors are common phaseA INPUTS here,
not free new network factorizations or a whole-PhaseA benchmark. Existing
topology is structural setup, cold clocks are not warmed benchmarks.
No cached numerical-factor assumption across Newton iterates.

Generated matrix reconstruction for validation stays outside execution timers
and does not feed execution. Coefficient/object sizes are not RSS savings;
both paths and reference/capture objects coexist. A partial pipeline win/loss
cannot establish a full direct pipeline, PhaseB wall or optimizer performance.
Report missing record/boundary costs before considering any further campaign.

New immutable small overlay, reused read-only numerical/capture roots, private
first depot, Gold6438Y+64CPU300GiB,Julia/BLAS/OMP1,occupancyguard,no-requeue,
35min batch/25min process cap. No expensive overlap; preserve failures and only
repair safe harness issues in new snapshots. Commit/push code/raw evidence and
assessment; production defaults, benchmark rows and paper untouched.
