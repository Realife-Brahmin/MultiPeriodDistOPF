# Dense merge and execution configuration study

This opt-in experiment replays the original 96-stage, 1,020-battery Stage A
capture using the grouped b=2 schedule. Production solver code is read from
the pinned historical snapshot. The candidate kernel is defined separately.

The current kernel stores Y1 = U' \ Mab' and Y2 = U' \ Mab_right.
The candidate stores Y1 = Mab / U through BLAS.trsm! with side=R,
uplo=U, trans=N, diag=N. It changes the left Hessian update to
Y1*Y1', the cross block to -Y1*Y2 and the gradient update to -Y1*g.
The right-side Hessian update remains Y2'*Y2. Upper-triangle updates and
mirroring match production. Both failed-Cholesky branches call the original
pint_combine and copy back exactly the same fields.

PairWork remains three nB-by-nB matrices and one nB vector. At nB=1020 this
is24,977,760 coefficient bytes per worker before object metadata. The shared
96-stage storage is2,398,648,320 bytes in every execution configuration.

run.sh launches the synthetic algebra/fallback test, the original saved late
policy check, then three persistent-worker processes in sequence. kernel63
alternates current/right kernels; threads31 alternates current_1/current_2;
threads15 alternates current_1/current_4. Every mode has a warm pass and
five clean timed passes. Each also has an allocation-diagnostic warm pass
and one allocation-diagnostic measured pass outside clean timings.
Slot reset, between-pass GC, CPU pinning, source checks, logging and reference
validation are outside clean Phase B timers. Matrix assembly/copies,
triangular solves and Hessian/gradient updates are inside every merge timer.

CPU maps choose64 physical CPUs from the allocated cpuset. One is reserved
for the coordinator. Each active worker receives a disjoint group of1,2 or4
cores inside a socket. taskset applies the group to all existing OS threads.
The thread count is explicitly set on every worker before each pass, then
its affinity/backend/thread count are logged. The initial worker-to-socket
assignment balances the available core counts. Memory placement and NUMA
first touch are not separately controlled in this pilot.

The exact-dispatch log and the pinned Julia source excerpt document that
ldiv! reaches LAPACK.trtrs! for these strided Float64 matrix right-hand
sides. Native stack samples in the preceding committed profile identify
OpenBLAS dtrsm inside that LAPACK path. Explicit right-sided TRSM also
changes which matrix operands are transposed in subsequent BLAS calls,
so performance differences cannot be attributed solely to the removed
transpose operation.

Every warm, timed and diagnostic replay checks all96 gradients and Hessians
against the same ordinary dense reference using scaled infinity-norm error
below1e-11. The saved late T24 trajectory uses all24 leaves and the true
future suffix, with four full stage packets and stored probes. Values use
the unchanged1e-7 gate; battery/control/equality/bound directions use both
scaled1e-7 and relative1e-5 with1e-10 denominator floor. The synthetic test
also exercises asymmetric raw input blocks and the unchanged fallback.

choose.py can trigger one combined kernel/configuration comparison when
both initial studies find a lower median. It uses the same five paired
repeats and reruns the late check at that selected BLAS width. All candidate
results must then be judged for repeatability from the raw pairs.
No optimizer invocation is part of this study.

To reproduce, create a fresh Kamiak snapshot containing these files, copy
the existing capture.sha256, record source.sha256, and submit run.sh with
64 physical cores,300GiB, Gold6438Y+ and a75-minute cap. Frozen paths in
run.sh and the source manifests identify the exact historical dependencies.
Do not modify an existing queued or running snapshot.
