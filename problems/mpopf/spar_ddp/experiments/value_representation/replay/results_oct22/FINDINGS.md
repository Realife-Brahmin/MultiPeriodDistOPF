# Eight-period original1020 representation probe

Job31356047 COMPLETED/exit0 in21seconds, cn164 Gold6438Y+,64CPU allocation,
32GiB request,Julia1.12.7,one computing process/BLAS1. No other jobs at startup;
full stdout/stderr and hardware/affinity archived. Sourcee592655d. No full
optimizer, new capture, benchmark rerun, production change or paper edit.

Input: eight existing first-pass original1020 T96 stage records. Every capture
matched its original SHA. Original input SHA:
02d3b377324a7c02314799292fe4393f8ef119c2fe05071acbfc95da7c2616f1.
The parent graph yields102 groups of10 individual batteries. No state was
removed or merged; approximate path sets cross-group quadratic terms to zero.
This is NOT an exact separator solver or a feasibility guarantee.

## Numerical result

All seven approximate compositions had maximum relative infinity error across
the three quadratic blocks and two vectors <=9.04e-15 versus dense composition.
Balanced versus left-fold DENSE association error1.16e-14 passed the1e-7 gate.
Closing the eight-period quadratic with zero terminal value gave approximate
Hessian error2.29e-14 and gradient error1.64e-14. No Cholesky fallback occurred.
This is a synthetic closing of the first8periods, NOT the actual96-period future.

All22950 singular-value rows (15intervals x3blocks x510values) are retained in
the archive. Across the selected whole-group half partition, the largest cross
singular value is only~3.3e-21 of full-matrix Frobenius norm in the first leaf
and~0.9--2.7e-20 after8periods. At relative-to-cross-block1e-8, apparent ranks
grow from1 to11/12/11 for Maa/Mab/Mbb. These tiny blocks can be dominated by
rounding/cancellation; do NOT interpret the apparent ranks as a demonstrated
fixed physical rank or a universally justified truncation tolerance.

## Kernel cost and storage limitations

Warmed allocating dense production pint_combine took median0.158--0.162s per
combination; groupwise uses the same routine on each block and took0.00274--
0.00284s (~56--58x kernel difference). Each has three repetitions, alternating
order, with GC before each pair; reported within-call GC waszero. The ordering
is mildly asymmetric across an odd number of repeats. These are single-process
allocation-heavy diagnostic kernels, NOT the production in-place parallel scan,
NOT a complete Phase B replay, and absolutely NOT a whole-solve speedup.

The hypothetical packed group-block payload is261120bytes versus24985920bytes
for dense3blocks+2vectors (95.7x smaller). The test actually embeds group
results in full dense zero arrays for comparison: that storage saving is NOT
implemented or measured. Observed allocations are26.29MB/groupwise call versus
49.99MB/dense call; neither is the footprint of an optimized compact backend.
Julia process maxRSS1,904,365,568bytes (~1.77GiB) includes inputs, spectra,
intermediate copies and analysis. Slurm sampled MaxRSS2796KiB misses the
short-lived Julia peak; do not use that number as the diagnostic's true memory.

## Decision

Promising on the first-pass slice, not sufficient for production. Next evidence
needed is separately captured early/middle/late passes under current exact
curvature, including barrier transitions, with scaled residual/inertia/direction
checks and rank growth. Inspect existing saved states first; reconstructing
matrices at them is setup/diagnostic work, never part of benchmark solve time.
Only then decide between exact separator representation and an explicitly
controlled approximation. Do not promote the old diagonal-Hessian structured
path or these near-machine-error first-pass results into an exactness proof.

The v1 source package had CRLF shell lines and failed bash-n before any Slurm
submission. It remains at value_representation_oct22; correctedv2 pinned LF and
passed syntax validation. No numerical failure was hidden and no gate relaxed.

Evidence archive SHA256:
f55682c5d4e018051e505af0cd7fc7fbc6dda3a65cf8beb80f8305225911ff3c.
Source snapshot /home/aryan.jha/kamiak_jobs/value_representation_oct22_v2.
Logs, complete spectra, rounds.csv, hardware and source hashes are preserved.
