---
name: Next tADMM Run Plan
description: Settings and validation checklist for next experimental run
type: project
---

## Next Run Configuration

### Ipopt Settings (Subproblem Solver)
- **max_iter**: 30 (keep this)
- **tol**: 1e-2 (CHANGE: loosen from 1e-6 to match pre-solve tolerances)
- **acceptable_tol**: 1e-3 (CHANGE: loosen from 1e-4)
- **Rationale**: Tight tolerance + 30-iter cap causes hitting iter limit. Loose tolerance will converge faster within 30 iters.

### tADMM Settings
- **max_iter_tadmm**: 100 (CHANGE: from 500, as you said "either converges fast or not at all")
- **Early stop criterion**: (‖r‖ ≤ 1e-4) AND (objective_improvement_over_10_iters < 0.1%)
- **Rationale**: System/T-agnostic stopping rule. ‖r‖ ensures feasibility. Objective stagnation (<0.1% change in 10 iters) catches diminishing returns automatically without needing hardcoded ‖s‖ threshold. Works across all systems and time horizons.

### Output Settings
- **Print frequency**: EVERY iteration (currently every 5-10, need granular data)
- **Print format**: k, obj, ‖r‖, ‖s‖, ρ, wall-clock per iteration
- **Rationale**: Build complete Phase 2 trajectory to identify exact early-stop point and diminishing returns curve

### Post-Convergence Validation
- Run feasibility checker on early-stopped solution
- Report: obj gap %, max voltage violation, max power imbalance, max SOC violation
- **Rationale**: Ensure 3x speedup comes with acceptable constraint violations

## Expected Results (Target)
- **Convergence time**: 12-18 min (vs 24.8 min BF, aiming for 2-2.5x speedup improvement over this run's 3x)
- **Solution quality**: Within 0.5% of BF objective (0.44% gap achieved at k=40 last run)
- **Feasibility**: Max violations < engineering tolerance (TBD after validation)

## Systems to Test
1. large10kC_1ph T=12 (baseline, already have Apr 5 data)
2. large10kC_1ph T=4, T=6, T=24 (explore T sensitivity)
3. large10kB_1ph T=12 (smaller battery set)
4. large10k_1ph T=12 (original system)

## Data Capture for PhD
Build confidence table across systems:
| System | T | ‖r‖_threshold | k_converge | wall-clock (s) | obj_gap % | max_feasibility_violation |
|--------|---|---|---|---|---|---|

This justifies: "Early stopping at ‖r‖ ≤ 1e-4 achieves 2-3x speedup with solution quality within 0.5% and feasibility violations < X%"
