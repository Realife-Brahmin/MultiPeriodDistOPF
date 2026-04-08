# Context for Apr 7 Session

## Project State

**Current Branch:** `hoping-for-the-best-apr05`
**Last Commit:** f484c5a3 (Add session notes: Apr 7 T=24 ADMM tuning and results)

## Key Changes This Session

### Code Changes (COMMITTED)

1. **eps_pri_tadmm**: 1e-4 → **1e-3**
   - Relaxed primal residual threshold to enable stagnation-based early stopping
   - Was too strict, blocking convergence detection

2. **stagnation_window**: 5 → **20**
   - Check improvement over 20 iterations (was 5)
   - More patient detection of when algorithm stops improving meaningfully

3. **progress_interval**: 5 → **1**
   - Print every ADMM iteration (was every 5th)
   - Helps monitor convergence in real-time

4. **max_iter_tadmm**: 500 → **100**
   - Cap iterations at 100 (stagnation should trigger before this)

5. **.gitignore**: Added `*.jls` and `**/model_summary*.txt`
   - Large serialized Julia objects not needed for reproducibility
   - Model summary files are verbose debugging output

### T=24 Results (large10kC_1ph)

**Configuration:**
- System: 10321 buses, 1020 batteries, 24 time periods
- Solver: Ipopt (no Gurobi)
- Threads: 16

**Results:**
```
tADMM:
  Iterations: 100 (hit max_iter cap)
  Objective: $2,976,617
  Time: 176.7 min wall-clock
  Gap to BF: 0.22%

BF:
  Objective: $2,970,084
  Time: 77.7 min

Speedup: 0.44x (tADMM slower)
```

**CRITICAL INSIGHT:** 0.5% margin was hit at **k=51 (50.9 min) = 1.53x speedup over BF**
- Algorithm didn't stop because stagnation criterion didn't trigger
- Consider: why didn't 20-iter stagnation window detect plateau?

### Files Saved Locally (NOT in git)

```
envs/tadmm/processedData/large10kC_1ph_T24/
├── convergence_data.csv                    # Iteration history (objfun, residuals, rho)
├── results_socp_tadmm.txt                  # Summary with timing
├── results_socp_bf.txt                     # BF comparison
├── subproblem_timing_details.csv           # Per-iteration timing breakdown
├── convergence/
│   ├── tadmm_convergence_socp.png          # Convergence plot
│   └── tadmm_timing_socp.png               # Timing plot
└── sol_socp_tadmm.jls, sol_socp_bf.jls     # Serialized solutions (62 MB, in .gitignore)
```

## Comparison Table

| System | T | tADMM Time | BF Time | Speedup | Gap |
|--------|---|-----------|---------|---------|-----|
| large10kC | 6 | 11.2 min | 11.1 min | 0.99x | 0.0% |
| large10kC | 24 | 176.7 min | 77.7 min | 0.44x | 0.22% |

## Known Issues

1. **stagnation_window=20 didn't trigger at k=100**
   - Objective was essentially flat from k=30 onwards
   - Why? Check if threshold calculation is correct
   - Possible: stagnation threshold (0.1%) still too strict for final polishing

2. **Archives directory (9700+ files)**
   - Not committed (untracked)
   - Old experimental results
   - Has long filenames causing Windows git issues
   - Safe to delete: `rm -rf archives`

3. **Git can't push due to file size limits**
   - 194 MB model_summary file was blocking
   - Now in .gitignore, but was in history from earlier commits
   - Code itself is fine and pushable

## For Next Session

1. Investigate stagnation threshold:
   - Is 0.1% improvement over 20 iters reasonable?
   - Test: loosen to 0.05% or 0.01%?
   - Or: use max(relative, absolute) improvement?

2. Early stopping target:
   - Goal: get tADMM to 20-30 min at 0.5% above BF
   - Current: 50.9 min at 0.5% (1.53x speedup)
   - Need: 1.5-2.6x speedup

3. Possible fixes:
   - Adaptive stagnation threshold (tight early, loose late)
   - Hybrid: stop if BOTH objfun stagnates AND residuals near threshold
   - Increase ρ more aggressively when stagnant

## Git Status

✅ **Ready to Push:**
- All code changes committed (5 commits ahead)
- Clean working directory
- .gitignore configured properly

⚠️ **Why can't push yet:**
- archives/ has Windows long-filename issues (not blocking, just untracked)
- Solution: `rm -rf archives` then push

