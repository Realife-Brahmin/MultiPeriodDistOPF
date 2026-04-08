---
name: Convergence Analysis
description: Detailed analysis of tADMM convergence behavior and adaptive rho issues at T=12
type: project
---

## T=12 Adaptive ρ Convergence (large10kB_1ph, 1020 batteries)

```
k    objective       ‖r‖       ‖s‖        ρ       note
1    $2,925,260     6.3e-02   2.15e+01   2121
5    $3,126,776     8.3e-08   8.39e+00   1061    ρ halved, primal converged
10   $3,116,730     6.3e-04   7.06e+00   1061
15   $3,111,163     3.9e-05   4.66e+00    530    ρ halved
20   $3,108,513     5.6e-06   1.79e+00    265    ρ halved
25   $3,107,739     4.3e-04   4.31e-01    265    obj ~converged here
30   $3,107,950     7.6e-05   2.65e-02    133    ρ halved
35   $3,107,919     8.9e-06   9.83e-04     66    ρ halved
37   CONVERGED                             66
```

- Final objective: $3,107,919 (matches BF exactly)
- BF time: 1,071s; tADMM: 1,330s (0.81x)
- Iterations k=25-37 are practically meaningless for objective quality

## Per-Iteration Ipopt Counts (confirming uniform difficulty)
- T=4: k=1 avg 164 iters, k=2-9 avg 10-14 iters (stable)
- T=12: k=1 avg 106 iters, k=2-37 avg 7-15 iters (stable)
- Conclusion: each ADMM iteration is ~equally expensive after warmup

## Adaptive ρ Code Location
- Config: `tadmm_socp.jl` lines 103-114
- ρ update logic: around line 1688 (`if adaptive_rho && k % update_interval == 0`)
- Residual computation: lines 1570-1587
- Convergence check: line 1822 (`if r_norm ≤ eps_pri && s_norm ≤ eps_dual`)
