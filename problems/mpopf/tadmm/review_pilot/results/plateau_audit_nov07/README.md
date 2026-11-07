# Causal audit of the med2522 T=24 plateau

The archive preserves raw solver/KKT logs, replay traces, source hashes and all
diagnostic source versions, including unsuccessful audit-harness development
attempts. Checkpoint/instance JLS files remain on Kamiak. See ../../PLATEAU_AUDIT.md
for the findings, assumptions, failed-attempt classification and model fix.

The controlled-replay PNG compares the same iteration-5000 state, fixed rho=7600
and MA57. Metric colors follow the repository convention. Eight original-accuracy
steps reproduce the two-cycle; eight tighter-complementarity steps restore
damping; four steps with only redundant source-voltage bounds removed produce
similar damping. These short replays do not establish full convergence or timing.

Archive SHA256:
`0fe1e965d25b631e0d0b230c90344be2e889eeb93766a46afb36ee2d0f90ccf9`

A fresh corrected-source med2522 T=24 case is job 31645746, in the isolated
`tadmm_review_nov07_accuracy_fix` campaign. Numerical source commit is d4491c1c;
the native Linux Manifest and frozen input exports are inherited explicitly
from the original campaign, with execution-source/input SHA256 manifests.
The earlier input/algorithm verification gates are referenced by hash, not
relabeled as validation of the new source. Full validation is required on the
new result itself. No rho initialization or update policy was changed.
