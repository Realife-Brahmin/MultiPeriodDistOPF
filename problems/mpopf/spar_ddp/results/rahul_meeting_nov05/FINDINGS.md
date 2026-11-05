# Rahul meeting campaign, November 5, 2026

The complete solve queue and the original four-case profile are finished.
All 22 full solves have completion markers, numerical checks, independent
original physical residual checks below 1e-5, and OpenDSS voltage checks with
zero violations under the existing tolerance. Four are separate sampled
profiling solves. Two are a single balanced/grouped pair. The other 16 form
eight two-repeat timing cases. Raw solutions remain on Kamiak.

Full-convergence optimizer medians for original 1,020-battery large10k are
18.518s at T12, 28.1705s at T48, 42.027s at T72, 50.1125s at T96,
111.496s at T144 and 199.038s at T288. T24 remains the previously committed
October 16 point. Workers are min(T,63), BLAS/Julia/OMP threads are one,
and all new jobs use Gold6438Y+ nodes. These DDP-only runs do not establish
new centralized speedups or shared feasible-near-optimal crossing times.

The newly paired complete T96 solves take 50.028s balanced and 46.731s
grouped, with Phase B 35.673633s and 32.828737s respectively. This is one
pair, with identical independently checked final objectives and no
production-default change. The 6.59% full-time improvement is preliminary.

The four-case comparison gives the following two-repeat medians.
100 batteries, T96: full5.224s, A2.519008s, B0.2218825s, C0.543095s.
1020 batteries, T96: full50.1125s, A5.078579s, B35.852710s, C5.8581665s.
100 batteries, T144: full8.0935s, A3.969331s, B0.420307s, C0.839960s.
1020 batteries, T144: full111.496s, A8.8851615s, B87.0040525s, C9.786764s.
Changing physical battery count also changes installed storage and the
optimization problem. It motivates computational investigation but is
not a pure computational ablation.

## Fixed-input Phase B worker study

The original T96 first-pass capture is reused. Both schedules use the
same dense kernels and shared slots. There are five alternating warmed
timing pairs per worker count, plus warmups, and every resulting gradient
and Hessian passes the existing 1e-11 scaled infinity-norm check.
Sequential control: median7.347037s, including the control routine's
deepcopy and stored suffix outputs.
32 workers: balanced2.066641s, grouped1.455091s.
48 workers: balanced1.355613s, grouped1.236975s.
63 workers: balanced1.319444s, grouped1.225511s.

At 63 workers, summing per-round merge time on the busiest worker explains
99.5-99.7% of round wall time. This quantity includes memory stalls inside
the dense kernels. It is not a bandwidth or communication measurement.
The average round's completion spread is about42ms balanced and21ms
grouped. Completion spread includes differences in task count and arrival
order; it is not an independent cost to add to the critical path.

The replay timer field 'compute' includes minor slot lookup and return
overhead around the merge. GC is disabled inside these replay passes and
runs between passes, matching the earlier replay protocol. This differs
from arbitrary optimizer iterations and cannot establish a universal
bottleneck fraction.

## Operation timer correction

Job31562015 completed and passed numerical checks, but its source
instrumentation missed Cholesky because a keyword argument contains a
semicolon. Its zero Cholesky column is invalid. Preserve all its raw data,
but use the corrected replay for the complete operation breakdown.
The correction explicitly times each complete Cholesky call, retaining
the same arithmetic and accuracy gates. No solver tolerance changed.

## Evidence

phaseb_full_profile_nov04_text_evidence.tgz includes the complete four-case
logs, stack profiles, hardware, input/source checks and solutions hashes.
rahul_nov05_text_evidence.tgz includes the horizon runs, paired solve,
worker study and first operation-timer attempt. SHA256 values are recorded
in evidence.sha256. The final corrected profiler has its own archive.

analysis/full_runs.csv contains all22 solve records.
analysis/fresh_full_solve_profile.png shows the phase comparison and the
fresh original-battery horizon series.
harness/ contains the launchers and analysis scripts.

To reproduce the analysis, extract the archives into a temporary directory
and run summarize_full.py from that directory. For the per-merge analysis,
run analyze.py beside the extracted operation replay, adjusting the input
directory to the corrected attempt as documented in its source.

Historical tree/screening ablations were reused. The interface-dimension
ablation was deferred by the user. No new centralized runs or longer
frozen T288 replay were conducted, and no T288 frozen first-pass capture
was found in the inspected replay locations.

Kamiak snapshots are immutable under /home/aryan.jha/kamiak_jobs/.
Main job IDs are31554990-31554993 (four-case campaign),
31559808/31559809/31559810/31559811/31559842/31559879
(horizon/pair/worker/freshT288),31562015 (first timer attempt), and
31574751 (corrected timer attempt).

## Corrected operation breakdown

Corrected job31574751 COMPLETED0:0 in2m36s. All12 warm/timed passes
passed the same gradient/Hessian gate, max2.98e-14 versus1e-11.
There are1840 timed balanced merges and1280 timed grouped merges,
with zero fallback in every merge. No operand arithmetic was changed.

Summed worker operation-time shares for balanced are39.64% matrix triangular
solves,20.83% symmetric rank-k updates,16.77% cross-block matrix products,
9.06% Cholesky,6.51% transpose,2.27% matrix assembly,2.16% triangle mirror,
1.83% copies,0.69% vector matrix products and0.23% vector triangular solves.
The unaccounted timer/dispatch-within-kernel residue is about0.04%.
These are elapsed times inside operations and include cache/memory stalls.

For the busiest worker in each round, summed across rounds then averaged
over five repeats, balanced operation times are0.5230s matrix triangular
solves,0.2775s symmetric updates,0.2074s cross-block product,0.1224s Cholesky,
0.1030s transpose,0.0311s assembly,0.0207s mirroring,0.0166s copies,
0.0087s vector products and0.0032s vector solves. This is an attribution
of the observed critical worker workload, not an independently removable
serial budget or a guaranteed speedup estimate.

Instrumented medians are1.320941s balanced and1.219395s grouped. The
corresponding earlier clean63-worker replay medians are1.319444/1.225511s.
Timing drift across jobs is larger than any measurable instrumentation
penalty here, so do not subtract them to infer negative timer overhead.

Two-sided merges average173.62ms balanced/167.57ms grouped, while closing
into a terminal one-boundary value averages84.09/79.01ms. This difference
is consistent with two-sided merges requiring two matrix solves, two
symmetric Hessian updates and a cross-block product, compared with one
matrix solve and one symmetric update when closing against a suffix.
These are concurrent-worker latencies on this fixed first-pass input.

The original profile's missing Cholesky timer was repaired before drawing
the complete conclusion. Its archive remains intact. Cholesky acceleration
alone can only address a modest part of measured merge work. The largest
target is applying dense factors to whole matrix right-hand sides and
forming the resulting dense Hessian blocks. Further approximation or
structured/matrix-free consumption must include setup and recovered
directions. Hardware counters or isolated locality experiments are still
required to distinguish arithmetic throughput from memory bandwidth/stalls.

All queued meeting jobs and the corrected profiler are complete. No
dedicated campaign heartbeat was created in this session, so it did not
provide proactive completion updates while the chat was idle.
