# Validated med2522 T=24 numerical-model correction

Job 31645746 completed full warm-up and measured restart in 389 iterations each.
Optimization, matching warm/measured traces, worker layout, source snapshot and
independent OpenDSS validation all passed. Numerical source commit is d4491c1c.
The native Linux Manifest, identical frozen network exports and referenced
earlier input/algorithm gates are explicitly hashed in the archive.

The common centralized/local model removes voltage bounds already implied by
the source-voltage equality and sets an explicit absolute complementarity cap
1e-8. The physical feasible set is unchanged. Original physical checks still
enforce all voltage limits, including the fixed source voltage. Neither rho
initialization/update nor profiles, costs, initial state, terminal convention or
ADMM consensus/dual update was changed.

Initial rho=3800 increases after iteration 10 to 7600 and remains fixed, exactly
as in the failed original-source run. The comparison figure verifies equality
of input SHA256, C_B, terminal-SOC coefficient, backend, rho policy and initial
rho before plotting. The preserved original warm-up failed at 5000; it is not
a successful timing baseline. No speedup is claimed relative to that failure.

Measured full ADMM loop: 304.555781 s, including actual communication,
coordination, physical diagnostics, output and checkpointing after warm-up.
Core algorithm: 195.444191 s. Centralized optimize wall: 14.971893 s. This is one
observation under the requested 25-CPU budget (24 pinned single-threaded period
workers and one coordinator) on a fully allocated Gold node. The remaining
runtime gap is a tuning/scaling question, not resolved by repairing convergence.

Final original violation: 1.8421e-8; all-copy residual: 7.3974e-9; dual residual:
9.7967e-5; relative objective gap: 2.2387e-8. Maximum OpenDSS voltage difference:
3.41e-10 p.u., with zero voltage violations at the shared checker's 1e-4 tolerance.
The SOCP equality-slack and fictitious-loss diagnostics remain separate; the
independent validation uses the stated constant-power/no-line-shunt benchmark.

The archive contains source snapshots, complete CSV traces, status/solver and
OpenDSS logs, metadata, worker/CPU details, resource samples and manifests.
Large serialized Julia instances and dispatches remain on Kamiak.

Archive SHA256:
`c6cb3ca38068bee150784c2daf03d7469debe1db2b7cae64cc6d11a8a427adab`
