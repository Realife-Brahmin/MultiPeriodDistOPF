# Planned large10k horizon/rho sensitivity campaign

User authorization: two starting rhos at T24/48/96; retain July C_B, corrected
model/current periodic profiles, free terminal SOC, MA57 and the restored
primal-first controller/watchdog. Run past early residual crossings to 1e-5.

The fixed plan is T24 [15000,60000], T48 [30000,120000], T96 [60000,240000].
The first two anchors are July winners; no July T96 large-system result is
recorded. Its 60000 anchor is the explicit extrapolation rho0=625*T. The second
candidate is four times that anchor. These are empirical candidates, not an
oracle-derived or theoretically optimal rule. Later fitting must disclose
which instances were used; independent validation must be distinguished.

Keep controller primal gate 1e-4, mu=5, tau_incr=2, tau_decr=5, interval 5,
and the original watchdog (15 high-residual counts, multiplier 2.5). Terminate
only when BOTH owned/all-copy primal RMS <=1e-5 and unique dual RMS <=1e-2,
or at the 1000-iteration budget. The objective-stagnation shortcut does not
terminate this campaign. Report cap/budget exits as such, not convergence.

Record first raw owned/all-copy RMS crossings at 1e-3/1e-4/1e-5, their
dual-qualified crossings, historical 0.5%-cost and current 0.1%-cost crossings,
and the separate original-physical NO rule. Raw crossings can be transient.
Save each unique first-crossing dispatch and evaluate it independently with
OpenDSS after the timed trajectory. Report cone equality gap/fictitious losses,
constraint errors and voltage/import/loss comparisons; relaxation exactness
or full operational feasibility is not assumed from RMS or voltage alone.

Use 24 pinned one-thread period workers plus one pinned coordinator for every
horizon. Each worker owns 1/2/4 periods for T24/48/96; groups solve sequentially
within a worker. Core time is the actual process-parallel stage round plus
coordination, not a hypothetical one-worker-per-period maximum. The legacy
ideal-time diagnostic now accounts for sums within each assigned worker group.
Monitoring uses a separate reserved CPU. Reserve the same full Gold 6438Y+
node CPU allocation (64 CPUs) and 128GB memory limit for all three horizon jobs.
The central MA57 solve remains one-threaded, with that same available allocation;
do not describe allocated CPU capacity as actual centralized CPU utilization.

Fully warm both centralized and each tADMM configuration, then reset original
primal/dual starts. Use one measured trajectory per starting rho, preserving
all outcomes. The coordinator/worker pool is reused between configurations.
Physical checks, full dispatch retrieval, snapshots and OpenDSS are outside
core solve timers. Record total loop, construction/warmup, full campaign and
tuning cost separately. Report single-run exploratory observations, not repeat
medians. Sampling records aggregate process-tree RSS/private/PSS and CPU time;
phase timestamps distinguish centralized and distributed memory footprints.

Central tolerances are not relaxed. Preserve full Ipopt time and iteration
trace, including its original infeasibility table. Posthoc early-quality
estimates are distinct from independently validated saved central dispatches.
DDP's October campaign did credit both methods' first near-optimal iterate,
without weakening Ipopt tolerances; MA57 large10k medians were about 11.0 vs
20.6 s full at T24 and 72.4 vs 105.3 s full at T96. These are different weighted,
soft-terminal instances and cannot be reused as this campaign's baselines.

Reviewer 3 explicitly challenges best-of-sweep times, excluded tuning cost,
centralized solver choice, concurrency and aggregate memory. The editor repeats
these requests, and reviewers 3-5 flag owned-only consensus residuals. This
campaign supplies measured sensitivities and implementation evidence. A dedicated
conic-solver comparison, battery-model extension and literature/novelty response
remain separate reviewer tasks; this sweep does not claim to close them.
