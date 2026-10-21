# One-node balanced combination scheduling

Uses the immutable first-pass large10k T96 capture from stage_b_replay_oct15.
Compares production left-stage ownership with task-index round-robin scheduling
across 32,48,63 workers on one matched Gold6438Y+ node. All stages remain in
shared memory; stage input ownership and reset behavior are unchanged.

Each combination writes its left stage only. Every round asserts distinct
left/right operand sets, preserving the existing round barrier. Every worker
has private scratch and one Julia/BLAS computation thread. Same-node shared
memory only; this is not a cross-node scheduling change.

One warm-up per scheduler then five timed repeats, alternating scheduler order
between repeats. Both schedules replay exactly the same captured matrices,
using the same combination arithmetic and all 96 stage results checked against
ordinary value recursion to relative infinity-norm error below 1e-11.
Reset, GC and numerical verification are outside the Phase B timer; collection
is held off during each pass. Logs report busy workers, maximum tasks assigned
per worker, each round's wall time, compute sums, memory and CPU/affinity.
These are Phase B diagnostics, not full optimizer solve times.
