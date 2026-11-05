# Completed T24 and 128-core T96 results (November5)

The renewed exact Phase B operations improve original large10k T24 enough
to exceed the existing centralized reference at the shared feasible0.1%
objective target. On128-core hardware, one worker per timestep is NOT the
fastest tested allocation:48x2 beats96x1 on the same node with the same97
physical cores in use including the coordinator.

## Usual hardware: T24,1020 batteries

| Method | Full optimizer median,s | Feasible0.1% target median,s |
|---|---:|---:|
| Original balanced/current24x1 |22.6705|11.9664841|
| Groupedb2/right24x2 |14.4070|7.6207093|
| Existing MA97 reference |18.6050|10.2142905|
| Existing MA57 reference |20.6110|10.9933650|

The candidate reduces full optimizer time36.45% and target time36.32%
against the fresh original-DDP control. Against the REUSED October16 MA97
two-repeat median, full speedup1.2914x and shared-quality speedup1.3403x,
versus0.8207x and0.8536x for the fresh original-DDP control.
Thus the targeted near-optimality comparison exceeds1 in this follow-up.
This is a new separately identified result, not an overwritten paper row.

Fresh DDP reference range22.595-22.746s; candidate14.277-14.537s. Both pairs
improve with disjoint ranges; no third pair is needed. Phase A1.709481->
1.712798s, B17.3450555->10.095429s, C2.5218035->1.5989255s.
AtT24 both arms have24workers; the candidate uses2BLAS threads each.
Every run has14accepted steps/iterations,23passes,8barrier reductions and
0backtracks. Same physical start hashes and strict same-iteration/state/
control1e-10 comparisons pass, including reference-repeat equivalence.
Max scaled state difference2.221e-15, control3.361e-16.

The existing central rows are from the completed validated October16
Gold6438Y+ campaign, same original frozen instance and matched settings.
Their source paths, exact rows and whole-CSV SHA are retained. No new
centralized run was necessary for this comparison.

Scope correction: the initially submitted batch DID run four extra fresh
centralized solves on the usual node. By the time the approved replacement
was applied, both original batches had already completed; cancellation could
no longer remove that work. Those outputs are preserved in the raw archive
but are EXCLUDED from the T24 comparison above. They are not represented as
historical references or as newly requested work.

## Same128-core node: T96,1020 batteries

| Configuration | Full,s | A,s | B,s | C,s | Feasible0.1% target,s |
|---|---:|---:|---:|---:|---:|
| Grouped/right63x1 |39.8380|4.1603|27.7296|4.8578|21.2098237|
| Grouped/right96x1 |37.3270|2.6235|29.1440|2.8274|19.8339246|
| Grouped/right48x2 |26.6860|4.0084|16.5776|3.3579|14.3644872|
| Fresh MA97 |88.4595||||62.3987935|
| Fresh MA57 |93.9290||||64.9217950|

Job31575547 ran all ten solves on Xeon6760P cn178,128physicalcores600GiB.
The same scan/kernel/source/model are held fixed across the DDP arms.
96x1 improves A/C relative to63x1 but not B. The48x2 arm is28.51% faster in
full time than96x1 and33.01% faster than63x1.48x2 and96x1 each use97physical
cores including the coordinator, while63x1 uses64. This supports a tradeoff
between stage-level and matrix-level parallelism, not a general ideal of
one worker per timestep.

The48x2 range26.657-26.715s is separated from96x1's35.897-38.757 and63x1's
39.789-39.887; no third pair is needed to identify the best tested arm.
Against same-node fresh MA97,48x2 is3.3148x at full convergence and4.3440x
at the shared target. Centralized BLAS/OpenMP remains1 under the inherited
reference protocol; this is NOT a centralized thread sweep or a claim of
globally optimal central configuration. All six DDP solves have16iterations/
accepted steps,25passes,8barrier reductions and0backtracks.96x1 solutions
match63x1 bit-identically;48x2 passes the unchanged strict gate with maximum
scaled state3.765e-15 and control3.338e-16. Reference repeats also match.

## Validation and timing scope

Strict collectors passed actual completed logs, original numerical primal/
dual/complementarity1e-5 and physical/bound checks, independent OpenDSS with
zero voltage violations, finite full objective/feasibility traces, source/
harness/input checks before and after, matched initial objectives, identical
DDP start hashes, runtime BLAS/backend and every OS-thread mask against disjoint
physical-core groups. Fresh central objectives mutually agree within1e-5
relative; all final objectives agree within0.1%. Scientific gates are unchanged.

Reported full time is the warmed optimizer loop, not total process wall.
Per-run files retain model build/full warmup/whole process and lifetime RSS.
Summed worker RSS includes shared mappings and is not private resident RAM.
The128-core batch MaxRSS180187264KiB is a whole-job high-water mark, not an
arm's peak; cross-generation speed differences cannot be assigned to core count.
Job31575546 completed0:0 in27m15s;31575547 completed0:0 in1h20m12s.
Raw solutions stay immutable on Kamiak, with hashes and raw evidence archived.

The adapted DDP collector initially reused a name for the reference arm and
reference row, skipping intended checks during development. Renaming the arm
variable restored all unchanged start/objective/trajectory checks. Only the
corrected strict outputs are reported. No solver source was changed.

## Remaining authorized follow-up

T24 on128-core hardware is newly submitted as job31576169 in
/home/aryan.jha/kamiak_jobs/dense_merge_core128_T24_nov05_v2.
It runs latest24x2 plus MA57 and MA97, two each forward/reversed, SAME node
within one allocation, preserving all numerical/settings/physical gates.
The inherited source and CPU maps are unchanged. The first preparation
attempt stopped BEFORE submission because a trailing blank line contained
CR; its source/failure archive is preserved. A new immutable LF-only
preparation succeeded. No numerical gate or arithmetic was changed.
No more centralized T24 work is queued on the usual hardware.
