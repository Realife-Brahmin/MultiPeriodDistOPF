# T144 complete-solve result

The groupedb2/right31x2 candidate reduces the complete warmed optimizer
median from107.5640 to78.2205s: **27.28% lower**, saving29.3435s.
The unchanged original large10k model retains all1020 physical batteries.

| Configuration | Full, s | A, s | B, s | C, s | Full range, s |
|---|---:|---:|---:|---:|---:|
| Balanced/current63x1 |107.5640|7.1122|85.8404|9.4035|107.300-107.828|
| Groupedb2/right31x2 |78.2205|10.2414|52.0445|9.2596|78.093-78.348|

Two forward/reversed repeats on the same otherwise unoccupied Gold6438Y+
nodecn154,64CPUs300GiB. Phase B saves33.7958s (39.37%).
Phase A increases3.1292s, C decreases0.1439s, and other optimizer work
increases about1.4670s. Both paired comparisons improve with disjoint ranges;
no third repeat is needed. Balanced624merges/8rounds versus grouped420/9.

All four runs pass original primal/dual/complementarity/bound gates,
independent physical/objective and OpenDSS validation. Each has18accepted
steps/iterations,27backward passes,8barrier reductions and0backtracks.
Independent original primal maxima<=1.024e-12. Final physical objective
difference is about4.66e-10absolute (1.72e-16relative).
The unchanged1e-10 state/control and same-iteration comparison passes in
both pairs: maximum scaled state5.165e-15, control7.637e-16.
Reference repeats match exactly. All27 objective/feasibility pass keys match;
maximum recorded objective relative difference0 and original-primal
difference3.411e-13. Physical starting-array hashes are identical.

Actual thread/backend/OS-mask audits confirm63x1 versus31x2, coordinator1,
disjoint physical cores and same-socket BLAS teams. All source/harness/input
hashes pass. Job31575217 completed0:0 in1h03m47. Batch maxRSS185027220KiB
is the whole-job historical high-water mark, not per-configuration private RAM.

Full optimizer excludes model construction/full warmup/independent validation.
Those costs remain in per_run.csv: baseline process wall980-1019s,
candidate913s; build~424-429s; warmup~162-203s. These are not clean
compilation-free preprocessing benchmarks. No whole-process speedup claim,
centralized comparison or changed main-paper row follows from this ablation.

Raw logs, audits, accounting, immutable provenance and solution hashes are
in dense_fullsolve_T144_evidence.tgz (binary solutions remain on Kamiak).
Strict collect.py ROOT OUTPUT --T 144 reproduces the tables.
All scientific gates are unchanged. Production defaults/paper are unchanged.
