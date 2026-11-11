# Validated large10k T=24 numerical-model correction

Job 31645790 completed full warm-up and measured restart in 267 iterations
with matching trajectories. Worker layout, optimization, source snapshots and
independent OpenDSS validation all passed. Numerical source is d4491c1c,
matching the med2522 corrected C4 campaign. Inputs and the native Linux
Manifest are hashed in the archive.

Initial rho=48 increased to 96 after iteration 10 and remained there. First NO
was iteration 168, sustained through the end, using the repository rule of
0.1% objective agreement and maximum original feasibility <=1e-5. Its measured
loop clock was 538.362682 s; this is a posthoc trace crossing, with OpenDSS
validation of the final dispatch rather than a separately saved NO dispatch.
At iteration 100 objective accuracy already passed, but original feasibility
was 1.49484e-4. This motivates the subsequent bounded rho calibration.

Full measured loop was 853.424708 s, including diagnostics, communication,
coordination, output and checkpoints. Core algorithm was 499.109143 s;
centralized optimize wall was 67.990670 s. These are one warmed observation
on the fully allocated Gold node with 24 single-threaded period workers plus
one coordinator, not a repeat-median final comparison against SPAR-DDP.

Final original violation was 3.90483e-7, all-copy residual 2.19064e-7,
dual residual 9.80293e-5 and relative objective gap 4.86832e-9.
Maximum OpenDSS voltage difference was 1.41e-10 p.u.; zero voltage violations
at the shared checker tolerance. Constant-power/no-line-shunt convention and
free terminal SOC are retained. SOCP equality slack and fictitious loss
remain separately reported in summary.toml.

The archive preserves all raw traces, source snapshots, logs, metadata,
resource samples and manifests. Serialized instances/dispatches remain on
Kamiak. no_baseline contains a newly scored, source-hashed NO summary.

Archive SHA256:
`bcefbef8471703cf1c5318d914be748dd400f24bd2d07857a56d7d20d13f958d`
