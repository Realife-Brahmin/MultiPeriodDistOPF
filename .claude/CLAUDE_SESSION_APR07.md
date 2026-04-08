# Apr 7 Session: T=24 ADMM Tuning & Results

## Current Status

**T=24 Results (large10kC_1ph):**
- tADMM converged in 100 iterations
- Objective: $2,976,617 (0.22% above BF)
- Time: 176.7 min (vs BF: 77.7 min)
- **Key insight:** 0.5% target hit at k=51 (50.9 min) = 1.53x speedup

**Code Changes:**
- `eps_pri_tadmm`: 1e-4 → 1e-3
- `stagnation_window`: 5 → 20
- `progress_interval`: 5 → 1 (every iteration)
- `max_iter_tadmm`: 500 → 100

**Results Files (saved locally, not in git):**
```
envs/tadmm/processedData/large10kC_1ph_T24/
├── convergence_data.csv              # Iteration history
├── results_socp_tadmm.txt            # Summary
├── convergence/tadmm_convergence_socp.png
└── subproblem_timing_details.csv
```

## Git Status

✅ **Committed & Clean:**
- All code changes (.gitignore, tadmm_socp.jl tuning)
- Ready to push (except blocked by archive long filenames)

❌ **Not Committed:**
- T=24 result files (194 MB model_summary + 62 MB .jls files)
- archives/ directory (Windows git issue with long filenames)

## Next Session

1. T=24 stagnation didn't trigger at k=100 → investigate why
2. Test looser stagnation thresholds
3. Consider extracting/archiving old results separately
4. Focus on achieving 20-30 min runtime at 0.5% above BF

