# Derivations and implementation notes

These documents explain design decisions; they are not uniformly current
configuration specifications. Main-paper measurement rules live in
[AGENTS.md](../../AGENTS.md) and the [workstream guide](../README.md).

Useful starting points:

- [PAPER_CODE_MAP.md](PAPER_CODE_MAP.md): connections between theory and code.
- [BATTERY_BLOCK_REDUCTION_EXPLAINED.md](BATTERY_BLOCK_REDUCTION_EXPLAINED.md):
  network/battery elimination rationale.
- [PARALLEL_IN_TIME.md](PARALLEL_IN_TIME.md): temporal decomposition history.
- [REF21_USER_INDEXING_DERIVATION.md](REF21_USER_INDEXING_DERIVATION.md):
  two-sided value indexing and composition.
- [CONSTRAINT_SCREENING.md](CONSTRAINT_SCREENING.md): screening discussion.
- [FILTERDDP_API.md](FILTERDDP_API.md): solver interface background.
- [VALUE_REPRESENTATION_ASSESSMENT_OCT21.md](VALUE_REPRESENTATION_ASSESSMENT_OCT21.md):
  original hypotheses; later [experiment findings](../experiments/README.md)
  supply the actual subsequent evidence and limitations.

The `FILTERDDP_*`, `KKT_*` and centralized-Ipopt notes record individual
optimization attempts, including rejected approaches. Old branches, paths,
patch instructions and intended next tests remain historical, not automatic
authorization. Keep reusable insights here; live queues/session handoffs do
not belong in maintained documentation.
