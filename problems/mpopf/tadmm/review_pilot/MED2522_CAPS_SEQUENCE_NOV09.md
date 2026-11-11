# med2522 capped-local-solve campaign — November 9, 2026

Authorized follow-up to the large10k T24 persistent-cap experiments. Frozen
launcher/source commit: b1ba396977702b264eef679ce981f4881c15f787.
Remote campaign: `/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_med2522_caps_sequence`.

Submitted jobs:

- 31667879: native persistent-cap and chronological-storage gates, then
  lossless chronological replays of saved T96 centralized/MA57 endpoints.
- 31667880: fresh T96 persistent cap10 experiment; audit31667881.
- 31667882: fresh T96 persistent cap30 experiment; audit31667883.

Last observed state: gate waiting for scheduler priority; cap10 depends on
the gate, cap30 depends on cap10. These submissions are not completed results.
Existing exact single-node jobs31660388/31660390 remain unchanged and queued.

Both new optimization jobs require one physical128-core Granite Rapids6760P
node,128GiB requested memory,96 single-thread workers and MA57. They use
current periodic profiles, July C_B, no terminal SOC penalty, shared
load-flow/idle-battery initialization, warmup/reset, rho20000 with the existing
bootstrap policy/floor48, and unchanged residual stopping rules. T96 covers
24 hours at15-minute resolution. No numerical model or defaults were changed.

The10/30 caps apply only to tADMM time-step subproblems, throughout all rounds,
with no exact recovery phase. Centralized uses the existing exact policy,
max_iter3000, tight tolerances and mandatory solved status. A truncated
centralized solve is an error, not an accepted benchmark reference.

Fresh runs include centralized and endpoint chronological OpenDSS replay,
100% charge/discharge efficiency, zero idle loss and carried battery state.
All such validation and passive energy checks are outside the solve timers.
Detailed compressed battery/bus/branch/period tables are retained on Kamiak.
Saved T96 endpoints from31660059 are replayed only for early physical insight;
their former two-node timings are excluded from the new single-node benchmark.

For context, the large10k T24 exploratory runs had tADMM core/loop times:
exact50.449995/54.768778s, cap3048.454048/53.691627s and
cap1037.262515/41.571119s. Independently rerun full centralized references were
33.790131,32.485390 and33.793986s, respectively. This is run-to-run timing
variation under identical settings; these are single trajectories, not repeat
medians. Centralized was never subjected to the local caps.
