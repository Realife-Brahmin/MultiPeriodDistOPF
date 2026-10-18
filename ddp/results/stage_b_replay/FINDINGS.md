# Stage B replay: ownership and transfer costs

Completed 2026-10-15. Diagnostic only: replay of the first real large10k T=96
backward pass's Stage A outputs, not full optimization or new OpenDSS results.
The original full solutions were separately validated. Each configuration has
one warm-up and three measured replays, reset to identical inputs. Every
stage's resulting Vx/Vxx is compared with the ordinary sequential recursion
outside the timer. Final two-node worst relative infinity-norm discrepancy:
2.977617591605e-14, against a required 1e-7.

## Measured result

One-node shared memory: 7.841, 5.092, 2.916, 1.945, 1.973 seconds at
7,16,32,48,63 workers. Same-node explicit fetch: 3.717, 2.584, 2.649 seconds
at 32,48,63 workers. Two-node explicit fetch: 5.714, 4.959, 3.958, 4.695
seconds at 32,48,63,96 workers. Complete repeated values are in comparison.csv.

At 48 workers, explicit fetch on the same node adds 32.8% relative to shared
memory; spreading those workers over two nodes adds another 91.9% relative to
same-node fetch. These are changes in complete replay wall time, not a pure
measurement of network latency: placement, serialization, memory access and
source-worker availability all contribute. No claim of a universal one-node
advantage or optimal worker count follows from one input and one allocation.

## Why the 48-to-63 plateau is not ideal T/2 saturation

The fixed ownership map is owner(t)=(t-1) mod W. Each combination executes at
the left stage's owner. For W=48, rounds 1--4 each have 48 combinations but
only 24 busy workers, each doing two combinations; round 5 uses all 48 once.
At W=63, busy worker counts in the seven rounds are 48,40,36,34,33,33,63.
Thus ownership causes round-dependent collisions and idle cores even when the
number of workers exceeds the number of simultaneously useful tasks.

For the third measured shared-memory replay, W48 round walls were
0.337,0.335,0.348,0.346,0.181,0.255,0.157 seconds; W63 walls were
0.181,0.341,0.344,0.341,0.341,0.254,0.166 seconds. W63 speeds round 1 but
slows round 5, almost cancelling the gain. Merely counting 48 combinations
and assuming 48 distinct workers is incorrect. The final rounds have 64
combinations. A potential next improvement is scheduling those combinations
more evenly while accounting for operand locality; no implementation change
or full-solve speedup from that idea is claimed here.

## Provenance and limitations

Capture job31187159 took 6m54s; one-node31187170 took5m19s. Original two-node
31187171 launched surplus Slurm tasks and is excluded from comparison.csv.
The v2 launcher fixed the task count, but job31187707 printed unhandled lazy
peer-connection errors during warm-up despite ultimately completing numerical
checks; its results are also excluded. Their full logs remain preserved.

Final two-node job31187814 completed in4m27s, exit0, with no unhandled task
errors or surplus-task exits. Its stderr does contain a SlurmClusterManager
warning about the requested lazy=false keyword. Therefore do not claim the
experiment proves eager connections fixed the earlier failure. The explicit
srun task limit is effective; the source records the attempted keyword too.
All 97 registered processes (96 workers and coordinator) reported Gold6438Y+.
Per-process CPU masks and RSS are in the final stdout. Worker maxima are
historical per-process maxima and their sum is not a simultaneous node peak;
unused persistent workers and previous configurations remain part of job RAM.
The original one-node harness did not record individual CPU masks.

Captures and sources are hashed. Every replay launcher checked the capture
checksums before starting; final source hashes were independently rechecked.
The snapshots remain on Kamiak; local evidence includes raw logs and hashes.
No expensive runs overlapped. Same-node explicit fetch and cross-node fetch
use the same matrix combination routine; remote fetch also includes waiting
and serialization. Summed task/fetch times overlap and are not wall time.
The harness's scheduling_wait_est_s is only an imperfect residual estimate,
not a direct per-worker barrier-wait measurement. Cold warm-ups are excluded
from medians. These first-pass replays do not establish variation throughout
the optimization trajectory.
