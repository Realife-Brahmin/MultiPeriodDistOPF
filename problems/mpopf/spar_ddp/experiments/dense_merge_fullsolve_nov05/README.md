# Dense merge full-solve ablation

The first batch compares original large10k at T96 with four configurations.

| Arm | Scan | Merge kernel | Workers | BLAS threads per worker |
|---|---|---|---:|---:|
| c1 | original balanced | current | 63 | 1 |
| c2 | grouped b=2 | current | 63 | 1 |
| c3 | grouped b=2 | right-sided | 63 | 1 |
| c4 | grouped b=2 | right-sided | 31 | 2 |

All configurations share one copied numerical snapshot, the original frozen
physical model and the same settings inherited from the completed battery-count
campaign. Full warm-up is followed by a fresh measured solve from the same
physical start. Initial energy and control arrays are hashed before warm-up.
The eight main runs execute c1/c2/c3/c4 in repeat1 and c4/c3/c2/c1 in repeat2.
Four T3 numerical, physical and OpenDSS checks must pass first. Smoke runs use
min(T,requested_workers), with the same kernel and BLAS setting as the arm.

The copied driver normally forces one BLAS thread during worker setup. The
patch replaces that reset with the explicit experimental target and updates
the existing _BLAS_DEFAULT reference. Existing temporary one-thread regions
inside tree elimination retain their normal restoration policy. Each worker
asserts the selected thread count and default at A/C boundaries and before
each merge. Runtime audits after setup and after the solve record thread
counts, BLAS backend, OS process/thread affinity and RSS. The coordinator uses
one BLAS thread and its own core. Worker core groups are disjoint, inside one
socket, with at most64 physical cores including the coordinator.

Only pw_combine selects the opt-in right-sided merge. Other consumers retain
their original algebra. The grouped schedule is the previously verified b=2
schedule, including the conceptual terminal record and within-round safety
checks. The right-sided kernel is copied from a6cf5490 without algebra changes.
Original upper-triangle semantics and failed-Cholesky fallback are retained.

Every run uses the original numerical, bound, independent physical-residual
and OpenDSS validators. Full optimizer and phase clocks, iteration/barrier/
backtracking records, objective and original feasibility traces are retained.
Warm-up, setup and validation are outside the reported optimizer clock. Raw
solutions stay on Kamiak and their hashes must accompany the archived logs.
RSS is historical peak with shared mappings, rather than simultaneous private
memory. Cgroup/Slurm accounting should be archived separately.

compare_solutions.jl records the existing same-iteration and1e-10 scaled
dispatch/state match gate against the matched c1 reference. A mismatch is
reported without deleting the solutions or forcing counts to match. Each arm
must independently pass the original1e-5 numerical/physical gates. The existing
0.1% objective-agreement threshold is recorded as an additional comparison,
without claiming a centralized speedup or a new optimality certificate.

collect.py ROOT OUTPUT validates completed main runs and compares their traces
at common iteration/barrier keys. Before strict collection, run the standalone
Julia comparison once for every candidate/repeat against c1 from that repeat,
writing comparisons_T96/c2_r1.toml and analogous files. Use --partial only for
an interim view. It never signifies campaign completion.

Report T96 results before any T288 extension. If a valid candidate beats c1 in
both repeat pairs and its timing range is clearly separated, test that candidate
versus c1 at T288 with two repeats in reverse order. Keep grouped b=2 fixed in
that extension. If the full-solve timing ranges overlap or the candidates are
too close to choose, add a third paired repeat before naming a winner. Do not
weaken the numerical or equivalence gates to advance. Preserve trajectory and
convergence differences, even when objective/feasibility checks pass.

prepare.sh and job.sh identify the immutable historical dependencies and the
new snapshot. Never edit a queued/running snapshot. Harness repairs require a
new immutable attempt and preservation of the old logs. Expensive jobs run
sequentially on an otherwise unoccupied Gold6438Y+ node with64 CPUs and300GiB.

Production defaults and published benchmark rows are outside this experiment.
