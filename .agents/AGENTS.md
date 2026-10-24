# Repository instructions

## Scope and workflow

This repository contains SPAR-DDP, tADMM and multi-source OPF (MSOPF) workstreams; MSOPF notes are in envs/multi_poi/CLAUDE.md.
Check the remote, branch and working tree before changing files. Code changes
use a workstream branch and PR; make separate focused commits. Do not launch
simulations merely to initialize a session or clean up the repository.

The SPAR-DDP paper is maintained separately at
https://github.com/Realife-Brahmin/TPWRS-2027-Sparse-Parallel-DDP-for-MPOPF.
Paper edits go directly to that repository's main branch when requested.
Deliver compiled paper PDFs as main.pdf unless the user requests a new name.

## Current SPAR-DDP evidence

DDP means Differential Dynamic Programming. The current solver is DDP4OPF,
our fork of FilterDDP.jl; retain its upstream license, NOTICE and citations.
The medium display name is med2522; the raw model ID remains ieee2522C_1ph.

Use envs/spar_ddp for solver work and envs/tadmm for OpenDSS instance export.
The original first-order DDP is retired in archives/first_order_ddp and is
not an alternative active solver or environment.

The completed October 16 campaign in spar_ddp/results/proper_logging_campaign_oct16
contains 96 validated runs. Main-paper reporting uses objective agreement
within 0.1% AND original-model feasibility (1e-6 ieee123; 1e-5 med2522 and
large10k), scored posthoc from runs continued to full convergence. Use two
repeat medians, selecting the faster MA57/MA97 median per matched case.
Keep full-convergence times, historical 0.5% targets and standalone diagnostics
separate. Do not promote exploratory timings into matched benchmark results.
The separate MA97 thread sweep is at spar_ddp/results/ma97_thread_sweep.

Comparisons require identical exported data, periodic profile, initialization,
weights, hardware/thread policy and timer boundaries. C_B is system-specific:
1e-3 ieee123, 1.9e-3 med2522, 2.4e-5 large10k. Terminal SOC is a soft penalty
with per-system coefficient fixed across horizons (1.75e4, 3.29e4, 420).
See spar_ddp/examples/power_system/terminal_soc_penalty.jl for implementation.
Do not alter solver defaults or result rows during documentation cleanup.
Require full timing/feasibility logging, explicit completion/status validation
and independent OpenDSS checks for newly authorized reportable runs.

## Preservation and historical material

Preserve raw logs, CSVs, source manifests, hashes and pinned campaign harnesses.
Historical evidence retains original paths and names. Consult dated context
in docs/history/ only when relevant: its old branch names, pending sweeps and
launch/publish instructions are historical records, not current authorization.
The obsolete always-enabled Cursor sweep rule has been retired.

Respect all workstreams. tADMM is retained for reference; fresh testing and
research may be needed, but neither is authorized by repository cleanup.
Use Bash for Julia one-liners; PowerShell quoting previously corrupted Julia
commands. Pkg strips Project.toml comments, so put rationale in a README.
