# Validated original-source large10k T=24 MA57 case

Job 31643647 completed its warm-up and measured restart in 793 iterations each,
with matching numerical traces and independent OpenDSS validation of centralized
and tADMM dispatches. This run uses the original frozen source 4fbae546, including
the redundant source-voltage bounds identified by the later plateau audit. It is
preserved as observed evidence, not relabeled as the revised model's performance.

Actual measured full ADMM loop time is 3489.930362 s (58.17 minutes), including
communication, coordination, diagnostics, checkpointing and output, after warm-up.
Core algorithm time is 2446.574890 s. Centralized optimize wall time is 104.605683 s.
Both use the same 25-CPU budget on a fully allocated 64-CPU Gold node; ADMM has
24 single-threaded workers plus the coordinator. These are single observations,
not repeat medians, near-optimal times or a matched SPAR-DDP comparison.

Final relative objective gap is 9.6442e-9; original violation 4.0088e-7; all-copy
residual 2.4945e-7; dual residual 1.0873e-5. Maximum OpenDSS voltage difference
is 1.61e-10 p.u., with no checker voltage violations. SOCP equality slack and
fictitious loss remain separate reported diagnostics. The archive preserves
source snapshots, input identities, traces, solver/validation logs, worker layout
and resource samples; large JLS instances/dispatches remain on Kamiak.

Archive SHA256:
`5b061e78b7b39b569c6285408778e1636e5fa73e9ff8828f88985f02a20388fc`
