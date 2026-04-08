# Apr 7 Session Summary

## T=24 Run Results

**Configuration Changed:**
- `eps_pri_tadmm`: 1e-4 → 1e-3 (relaxed primal residual gate)
- `stagnation_window`: 5 → 20 (patient convergence detection)
- `progress_interval`: 5 → 1 (every iteration)
- `max_iter_tadmm`: 500 → 100 (cap iterations)

**Results:**
- Converged in 100 iterations (hit max_iter)
- Objective: $2,976,617 (0.22% above BF's $2,970,084)
- Time: 176.7 min (BF: 77.7 min = 0.44x)

**KEY INSIGHT:**
- 0.5% target hit at k=51 (50.9 min) = **1.53x speedup over BF**
- Algorithm didn't stop because stagnation didn't trigger at k=100

## Files Committed

✅ tadmm_socp.jl (all parameter changes)
✅ .gitignore (exclude *.jls, model_summary)
✅ CLAUDE_SESSION_APR07.md (high-level notes)
✅ CONTEXT_APR07.md (detailed technical context)

## Results Saved Locally

envs/tadmm/processedData/large10kC_1ph_T24/:
- convergence_data.csv (iteration history)
- results_socp_tadmm.txt (summary)
- convergence/*.png (plots)
- subproblem_timing_details.csv

NOT in git (by design): *.jls files, model_summary

## Next Session Focus

1. Why didn't stagnation trigger at k=100?
   - Objective was flat from k=30+ 
   - 0.1% threshold over 20 iters should have caught it
   - Investigate threshold calculation

2. Target: 20-30 min at 0.5% above BF
   - Currently: 50.9 min at 0.5%
   - Need: 1.7x-2.5x speedup

3. Ideas to test:
   - Loosen stagnation threshold (0.05%?)
   - Adaptive threshold (tight early, loose late)
   - Hybrid with residual gates

## Git Status

Branch: hoping-for-the-best-apr05
Commits: 5 ahead of origin
Status: Clean, ready to use

