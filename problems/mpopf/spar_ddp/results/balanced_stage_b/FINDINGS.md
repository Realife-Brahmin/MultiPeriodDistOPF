# Balanced shared-memory Stage B scheduling

Completed 2026-10-15. Job31189926, matched Gold6438Y+ node cn165, 64 CPUs
allocated, 63 persistent single-threaded workers plus coordinator. No other
jobs on node at startup. Job elapsed4m11s, status COMPLETED/exit0. No runtime
errors in stdout or stderr. All64 registered processes reported Gold6438Y+;
CPU masks and memory high-water marks are preserved in raw stdout.

At32 workers, existing ownership took2.879132s versus2.089921s balanced:
27.41% reduction. At48:1.989472s versus1.338464s,32.72% reduction. At63:
1.947521s versus1.325754s,31.93% reduction. Medians use five timed repeats,
after one warm-up per scheduling method. Methods alternate first/second
within each repeat; same job, same inputs, same shared storage and arithmetic.

All36 replay passes checked every stage against ordinary sequential value
recursion; worst relative infinity-norm error2.977617591605e-14, against
required1e-11. Source SHA checks passed. The launcher checks the original
capture SHA before starting. Input is the first real large10k T96 backward
pass, so these are Stage B replay times, not optimizer solve times and not
an assertion about every iteration or horizon.

## What changed

Stage ownership and matrix storage remain fixed. Existing scheduling sends a
combination to owner(left stage). Balanced scheduling sends combination q to
worker (q-1)modW. Shared-memory mapping permits that worker to read operands
and overwrite the left stage. Its private scratch is unchanged. Each round
asserts disjoint left/right operand sets and uses the same barrier.

At48 workers the first five rounds now use48 workers once each, rather than
only24 workers twice each in the first four rounds. For measured repeat5,
balanced first-five round times were0.1813,0.1831,0.1827,0.1823,0.1829s.
This confirms ownership-induced imbalance was real and correctable. The last
two rounds have64 tasks, and still require at most2 tasks per worker at48.

Balanced63 is only0.95% faster than balanced48 in this input. Thus once the
assignment is improved,48 workers are already sufficient for most of this
particular phase's benefit, but this is an empirical result for this schedule
and input, not a universal T/2 worker limit.

## Scope and memory

The tested change is implemented in the diagnostic replay harness. Production
procs_pass.jl and complete optimizer benchmarks have not yet been changed or
rerun. No two-node claim follows: moving computation there requires explicit
data transfers and locality-aware placement.

REPLAY_MEMORY gives sum of active workers' historical maxRSS, not simultaneous
node peak, and does not separate owner/balanced memory: reported sums at32,
48,63 were46201311232,77117014016,96337387520 bytes. Repeated configurations
retain processes and shared-memory mappings, so these are harness high-water
marks, not solver RAM savings. Reset, validation and per-pass collection are
outside the timer; GC is held off inside Phase B. Full raw evidence, original
source hashes and capture hashes are in evidence/.
