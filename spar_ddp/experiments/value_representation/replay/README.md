# Authorized short follow-up, October22

First eight existing original1020 first-pass captures; no optimizer, physical
battery removal, new capture, model changes, or paper changes. One bounded
single-node diagnostic allocation: same Gold6438Y+,64CPUs,32GiB,15min limit,
Julia1.12.7/BLAS1. Only one Julia process computes; unused cores are requested
to enforce uncontended matched-node diagnostics, not to imply parallel speedup.

The test measures complete cross-partition singular spectra at composition
depths0--3, storage counts, and a deliberately APPROXIMATE blockwise path.
Groups come from the network parent graph, never hardcoded feeder counts.
All1020 state rows and both gradient vectors remain; cross-group interactions
are dropped, not batteries aggregated. This is NOT the proposed exact
separator backend. Its errors and cost establish a diagnostic reference for
whether developing that backend is worthwhile. A dense balanced/left-fold
association comparison is the exact arithmetic regression gate (<1e-7).
No small approximation-error gate is imposed: a large error is a scientific
outcome, not reason to relax the validator or silently change the algorithm.

Dense production pint_combine and groupwise production pint_combine are warmed,
then each timed three times with alternating order. Allocations/GC recorded;
these allocation-heavy kernels are not optimized production scan benchmarks.
No claim about shared-memory scheduling, full-solve speed, feasibility or
later-iteration structure can follow. A zero-terminal closing checks value
Hessian and gradient errors but is not the actual full96-period future.

The older network structured path and separator-size algebra remain candidates
for a subsequent exact implementation, not something already validated here.
Ranks across multiple network groups can grow through elimination. Inspect
singular values relative to both the cross-block norm and whole matrix; tiny
first-pass couplings can misleadingly appear compressible. Masking out diagonal
blocks can itself inflate rank, so spectra use intact cross-partition blocks.

No new full optimizer/capture campaign or automatic extension is authorized.
Preserve source/input SHA, all stdout/stderr and adverse outcomes. If queued,
monitor quietly; report completion/failure and only expand after assessment.
