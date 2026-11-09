# Corrected med2522 T96 hardware allocation

User requires one128-core physical node for med2522 T96,96 single-thread period
workers and distinct coordinator/sampler CPUs. The former two64-core Gold-node
measurements remain historical and are not relabeled. New matched central and
tADMM timings use the same Granite Rapids6760P node type, one128-core allocation,
128GiB requested memory, same numerical model/data/rho policy and warmup/reset.
MA57 then MUMPS run sequentially. Source/input hashes, actual CPU model, node
count, unique worker pins, solver statuses and OpenDSS are validated/logged.
No claim that hardware alone fixes the previous1000-iteration nonconvergence.
Unit gates from the anchor campaign apply to the unchanged numerical kernels;
the only execution-source change is the launcher CPU-capacity/model check.
Reproduction: repository-root `setup-med2522-single-node-nov09.sh`.
Remote: `/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_med2522_single_node`.

## Existing iteration-cost evidence

`results/historical_anchors_nov09/iteration-costs.json` and repository-root
`collect-anchor-iteration-cost-nov09.py` preserve observed core iteration times,
local solver distributions (excluding first iteration), exact40-iteration
prefix times and explicitly hypothetical count-only break-even estimates.
All estimates exclude build/warmup/validation as in the established core timer.
An observed40-iteration prefix is not a converged solution. Reducing outer
iterations can change local solve costs, so projected times are conditional.

Large10k T24 MA57: mean1.030s/round, local solve median0.194s, median slowest
stage0.703s. First40 core44.807s versus central33.790s; observed last prefix
below central is29. MUMPS: mean2.044s, local median0.687s, first40 core89.131s
versus69.721s; last prefix below central30. Mean-based hypothetical break-even
counts32/34 differ from observed prefix break-even because costs vary by iteration.

Med2522 T96, superseded two-node setup: MA57 mean0.301s/round, local median0.117s,
median slowest stage0.209s; first40 core12.777s versus76.446s central. MUMPS
mean0.438s, local median0.237s; first40 core20.165s versus125.121s central.
Both prefixes remain unconverged; these demonstrate timing headroom only.

Large10k T48 MUMPS: mean2.282s/round, local median0.700s; first40 core90.810s
versus180.519s central. Reaching matched quality in40 at that cost would beat
central; the actual tested trajectory required129 iterations for the RMS/dual
endpoint and retained a separate strict physical-feasibility failure.
