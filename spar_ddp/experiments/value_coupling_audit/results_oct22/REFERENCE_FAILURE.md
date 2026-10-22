# Preserved v1 replay failure and independently confirmed reference repair

Replay31357649 (source96154b8a) FAILED after12m56s, not a successful campaign.
All five T3 gates and four medium T24 points completed; medium late failed
BEFORE any approximation at stage1's sine equality-direction probe. Required
unit-floor Frobenius error<1e-7; observed2.1699319339532596e-6. The guard blocked
all large-case work. No tolerance, physical constraint or gate was relaxed.
All v1 text/CSV outputs, manifests and logs are in replay_v1_failed_evidence.tgz
(SHA2560bc4d04a0b4d114faf0e73de6a15f0eda8b2fba70a6092a6eddb23e4b489ff68).

An independent ONE-POINT frozen diagnostic31357837 (sourcedcfbab49,
/home/aryan.jha/kamiak_jobs/value_reference_diagnosis_oct22) completed in21s.
It uses the same saved medium-late leaves, factors, nominal point and zero/sine
probes, and compares three numerical paths without optimization:

- Ordinary dense recurrence: worst saved-probe error2.16993e-6.
- Allocating pairwise scan: worst error3.55322e-8.
- Production in-place scan, actual round pattern: ALL12 saved control/equality
  probe errors exactly0 across periods1/12/24.

Thus the v1 reference did not reproduce the production floating-point path.
Small value-recurrence differences can be amplified in a late equality-
multiplier feedback direction. This is not evidence that dropping or retaining
cross-group coupling caused that failed gate: no approximation had run there.
It is also a useful warning that a Hessian/value norm alone cannot certify
direction accuracy near convergence.

Repair ONLY in a new immutable replay snapshot: use production in-place scan
for saved-policy checks and matched dense/candidate combinations, keep ordinary
recursion as an independently reported all-stage value comparison, and retain
the original1e-7 policy gate. Re-run the small gate before remaining cases.
No benchmark or paper row is changed; cold diagnostic clocks are not speedups.
