# Feeder locality of FilterDDP sensitivity maps

## Question

Could FilterDDP reduce its sensitivity work by treating a battery-state
perturbation as a local feeder event and discarding responses beyond a small
topological radius?

One exact large10k `T = 3`, stage-1 KKT solution was grouped by the graph
distance from each perturbed battery bus.  Each control and equality-multiplier
family was assessed separately so quantities with different units were not
mixed.  For every state column, squared sensitivity magnitudes were normalized
within the family; the reported radius encloses 90% or 99% of that energy.

## Result

Only the direct battery-power response is genuinely local: both its median
90% and 99% radii are zero.  The network-mediated responses are broad:

- real branch power `beta`: median 90% radius 18 edges;
- reactive branch power: 27 edges;
- voltage: 31 edges;
- branch current: 34 edges;
- battery energy slack: 47 edges;
- real-power-balance multiplier `omega`: 53.5 edges;
- reactive-power-balance multiplier: 36 edges.

At radius 8, a nominally local neighborhood captures median energy fractions
of only 44.0% for real branch power, 5.2% for reactive power, 42.6% for
voltage, 41.1% for current, and 0.18% for the real-power-balance multiplier.
Even radius 16 captures only 0.42% of the median real-power-balance multiplier
energy.  Radius 32 is needed for roughly 88--96% of most network-control
families, yet still captures only 20.8% of the median real-power-balance
multiplier energy.

The tails also vary strongly by battery.  The 90% radius reaches maxima of 82
edges for real power, 77 for voltage/current, 117 for energy slack, and 125
for real-power-balance multipliers.  A fixed small-radius truncation would
therefore damage some state columns much more severely than the median column.

## Decision

A naive spatial cutoff or independent local-feeder KKT solve is rejected.  It
would preserve the direct battery-power response while discarding substantial
network and multiplier information that FilterDDP uses in its value recursion
and rollout.  No approximate solver change is made.

This does not rule out spatial decomposition with an explicit boundary or
interface message.  A credible future method would need to communicate the
long-range balance and voltage effects rather than assume they vanish.  The
raw summaries are `large10k_stage1_locality_distance_energy.csv` and
`large10k_stage1_locality_locality_radius.csv`; the reproducible analysis is
`ddp/examples/power_system/analyze_network_sensitivity_locality.jl`.

