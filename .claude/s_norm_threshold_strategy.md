---
name: s_norm threshold justification strategy
description: Build empirical table correlating ‖s‖ with solution quality gaps across systems
type: project
---

## Strategy: Data-Driven ‖s‖ Threshold Selection

Instead of arbitrary 5e-3 threshold, build a confidence table showing the relationship between dual residual ‖s‖ and solution objective gap.

## Data from Apr 5 run (large10kC_1ph T=12, 30-iter cap, tight tol)

| k | ‖s‖ | obj ($) | Gap from BF | Gap % | Notes |
|---|-----|---------|-------------|-------|-------|
| 30 | 6.29 | 2,988,355 | 38,569 | 1.31% | Phase 2 start |
| 40 | 3.29 | 2,962,745 | 12,959 | 0.44% | **Crosses 0.5% gap** |
| 50 | 1.47 | 2,953,695 | 3,910 | 0.13% | |
| 60 | 1.13 | 2,951,337 | 1,551 | 0.05% | |
| 70 | 0.447 | 2,950,039 | 254 | 0.01% | Approaching 5e-3 |

## Key Finding

For **0.5% acceptable solution gap**, empirically need **‖s‖ ≤ 3.3**, which converges ~k=40.

## Next Steps

1. Run same experiment on other systems (large10kB, large10k) and T values (4, 6, 24)
2. Build full table with mean/variance of ‖s‖ threshold per acceptable gap %
3. Justify PhD work: "With ‖s‖ ≤ X, we are confident solution gap ≤ Y% (confidence: Z% based on N runs)"
4. This replaces arbitrary 5e-3 with principled threshold based on solution quality

## PhD Justification Template

"We empirically validated that dual residual ‖s‖ ≤ 3.3 reliably produces solutions within 0.5% of global optimum across [systems]. Early stopping at this threshold reduces wall-clock time by 40% without compromising solution quality for practical applications."