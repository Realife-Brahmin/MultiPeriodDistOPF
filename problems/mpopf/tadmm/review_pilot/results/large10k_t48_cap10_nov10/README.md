# Large10k T48 persistent cap10

Kamiak solve 31679860 completed on cn168 in 18m15s. Energy audit 31679861
and explicit cone, storage and memory audit 31679862 completed. The fresh
exact comparison 31679947 is running on the same node. Expanded physical
and timing audit 31679974 depends on both original explicit audits.

Measured cap10 tADMM took 176 rounds, 154.6995 s core and 175.9808 s loop.
The fully converged centralized solve took 87.8136 s. Of 8,448 subproblem
solves, 948 reached the cap, including three in the final round. Maximum
sampled process-tree PSS was 53.94 GiB. Validation is excluded from solve
timing.

The all-copy RMS endpoint was 9.15e-6, and objective disagreement was
1.25e-7 relative. The original-model maximum violation was 8.12e-4 and
the strict feasibility and chronological schedule flags were false.
Maximum SOC-state disagreement was 0.5852 kWh and 0.02926 percentage
points of rated capacity. Requested SOC-bound violation reached 0.8123 kWh
and 0.04062 points. Dispatch tracking error reached 0.03480 kW and 0.008843%
of rated power. The maxima in physical and normalized units are evaluated
independently.

The independent raw maximum cone gap was 1.51805e-4 pu squared for
tADMM and 1.55854e-5 pu squared for centralized. The exploratory
absolute target of 1e-5 was exceeded by both. The separate normalized
statistics retain the scaling denominator.

`round-timing.csv` derives the local time distributions from the raw
per-period trace. The compact evidence archive retains traces, source and
input manifests, validation summaries, logs and memory samples. SHA256:
`67a44e392f4a2f05b8252ff0ddde8bf4e51cff316413aa71f161dca6624e2e6f`.
Full dispatches and replay tables remain in the frozen Kamiak run directory.
