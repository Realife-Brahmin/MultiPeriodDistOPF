# Full-solve value-group diagnostic, authorized October22

Not a mock battery study: all physical batteries and their individual limits,
states, locations and exact network constraints remain. Opt-in research flag
FILTERDDP_VALUE_GROUPING selects dense (unchanged), root-branch grouping,
or chunks16 (deterministic sorted-bus chunks of16, a stress partition that can
cut strongly coupled branches). No root count is hardcoded. This tests the
concern that only large10k's particular root-branch structure makes dropping
cross-group terms benign. Existing small/medium networks may have only one
root branch, so they also need the stress partition to exercise approximation.

One serial Gold6438Y+64CPU300GiB allocation,Julia1.12.4/BLAS1. Unit tests plus
all3 IEEE123T3/two-worker full-solve smokes gate3systemsT24 x3modes=9 full
diagnostic solves, one repetition only. Full warm-up per arm;150iterations,
15min per optimizer process including build/warm-up. Batch2hcap. Scientific
approximation failures are preserved and later systems still tested. Dense
baseline failure stops the job. No default production flags are changed.

The optional helper acts ONLY in scan combinations: subblocks use production
combination math, cross-group quadratic blocks become zero, vectors remain
individual. Phase C still solves the original exact stage block with the
approximate future value. Shared storage remains dense: no compact-memory
implementation or claim. Diagnostic reference computations at iterations
0/5/10/15 (barrier index0) perform the exact ordinary value recursion at the
SAME iterates and report suffix Vx/Vxx errors at every period. Dense mode must
agree<1e-7. Approximate errors are scientific observations, not masked by a
dense fallback. Selected exact leaf matrices are saved at first/middle/last
period in the timed diagnostic solve for early/middle/late structure analysis.

Audits/captures add overhead INSIDE the printed solve clock. These are
convergence/accuracy/robustness diagnostics, not publishable solve-time
comparisons. Preserve pass/iteration/barrier counts, raw errors/GC/timing and
all unsuccessful outcomes; do not claim whole-solve speedup from this batch.
Require status, original primal/dual/complementarity, original bounds,
independent physical objective/residual re-evaluation and OpenDSS. Compare
each approximation to its own same-system dense reference; no assumption of
identical trajectories or automatic equivalence. No tolerance relaxation.

No T96/full-horizon escalation, centralized retuning, new physical network,
or paper edits. Snapshot is immutable after submission. Setup-only repairs
need a new snapshot. Full-convergence failures can be interesting outcomes;
campaign.finished does not mean all arms validated. Do not relaunch arms into
existing output paths after preemption: preserve and recover in a new version.

Completed job31356155: all9full cases and3smokes validated. See
results_oct22/FINDINGS.md. Root partition is trivial (onegroup) on small/medium;
stress partition converged but increased iterations41--50%. Diagnostic clocks
are not benchmark timings. The j==0 sampling condition prevented intended late
barrier audits: actual coverage is k0 small/large and k0/5 medium. Later-state
structure remains unverified and requires a separate corrected diagnostic.
