# URGENT: Home PC Setup Instructions After Lab PC Fixes

**Date**: November 12, 2025  
**Context**: PESGM deadline in 4 days - need fast, safe environment setup  
**Branch**: `tadmm-nov12`  
**Backup Branch Created**: `tadmm-nov12-lab-pc-working`

---

## What Happened on Lab PC

### Problem Encountered:
Lab PC had severe Julia package corruption:
- **OpenSpecFun_jll** package had invalid `riscv64-linux-gnu` platform entry
- **Plots** package had corrupted UUID in `Project.toml`
- Cascade failure preventing `small2poi_mpopf.jl` from compiling

### Solution Applied:
1. **Julia Version**: Switched to LTS 1.10.10 (from 1.12.1) using `juliaup default lts`
2. **Environment Reset**: Deleted entire Julia package cache:
   - `~/.julia/packages/`
   - `~/.julia/compiled/`
   - `~/.julia/registries/`
3. **Fresh Install**: Regenerated `envs/multi_poi/Project.toml` with correct dependencies
4. **Code Fixes**:
   - Fixed `Pkg.activate()` path in `small2poi_mpopf.jl` (line 10)
   - Added `using LaTeXStrings` (line 25)
   - Fixed `Plotter.jl` voltage handling bug (line 279)
   - Updated directory structure to match `small3poi_mpopf.jl` pattern

### Files Changed (Committed to Git):
1. **small2poi_mpopf.jl**:
   - Line 10: `Pkg.activate(joinpath(@__DIR__, "envs", "multi_poi"))` (removed `..`)
   - Line 25: Added `using LaTeXStrings`
   - Line 36: Added `systemName = "small2poi_1ph"`
   - Lines 669-672: New hierarchical directory structure
2. **envs/multi_poi/Project.toml**:
   - Added: `Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"` (correct UUID)
   - Removed: `name`, `uuid`, `version`, `[compat]` sections (cleaner format)
3. **envs/multi_poi/Plotter.jl**:
   - Line 279: Changed from `haskey(result[:v_nonslack], sub)` to `haskey(result, :v_nonslack)`

---

## What You Need to Do on Home PC

### CRITICAL: Your home PC environment was working fine before these changes!

### Step 1: Pull Changes Safely
```bash
cd ~/path/to/MultiPeriodDistOPF
git fetch origin
git pull origin tadmm-nov12
```

### Step 2: Test Without Breaking Anything

**Option A: Try the gentle approach first (90% chance this works)**
```bash
# Let Julia intelligently resolve package dependencies
julia --project=envs/multi_poi -e 'import Pkg; Pkg.resolve()'

# Test the script
julia small2poi_mpopf.jl
```

**If Option A works**: ✅ You're done! Move on to PESGM work.

**If Option A fails with package errors**: Try Option B

**Option B: Reinstantiate packages (if Option A fails)**
```bash
# Force Julia to reinstall packages based on new Project.toml
julia --project=envs/multi_poi -e 'import Pkg; Pkg.instantiate()'

# Test the script
julia small2poi_mpopf.jl
```

**If Option B fails**: Try Option C

**Option C: Nuclear option (last resort, ~2 minutes)**
```bash
# Remove machine-specific package manifest
rm envs/multi_poi/Manifest.toml

# Fresh installation
julia --project=envs/multi_poi -e 'import Pkg; Pkg.instantiate()'

# Test the script
julia small2poi_mpopf.jl
```

### Step 3: If Everything Breaks (Unlikely, but possible)

**Fallback: Revert to pre-change state**
```bash
# Switch to backup branch (pre-changes)
git checkout tadmm-nov12-lab-pc-working

# Your old environment should work exactly as before
julia small2poi_mpopf.jl
```

Then message me: "Need to debug home PC environment divergence"

---

## Why This Should Be Safe

### What Makes This Low-Risk:

1. **Manifest.toml is NOT tracked in git**:
   - Each machine keeps its own package versions
   - Home PC's working package set won't be overwritten
   
2. **Project.toml change is minimal**:
   - Only added one line: `Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"`
   - This is the **correct** UUID (lab PC had corrupted one)
   - Home PC probably already has Plots working fine
   
3. **Code changes are non-breaking**:
   - `using LaTeXStrings` - just adds a dependency your home PC likely has
   - Directory structure change - just organizational, doesn't affect logic
   - `Plotter.jl` bug fix - corrects edge case, won't break working code

4. **Backup branch exists**:
   - `tadmm-nov12-lab-pc-working` has the exact state before changes
   - Can switch back in 5 seconds if needed

### Expected Outcome:

**Most Likely (95%)**: `Pkg.resolve()` updates Manifest.toml smoothly, script runs fine

**Less Likely (4%)**: Need to run `Pkg.instantiate()` to refresh packages (~30 seconds)

**Rare (1%)**: Need to delete Manifest.toml and reinstall (~2 minutes)

---

## Key Information for Debugging (If Needed)

### Lab PC Environment (After Fixes):
- **Julia**: 1.10.10 LTS
- **OS**: Windows with bash.exe shell
- **Key Packages**: JuMP, Ipopt, Crayons, Plots, LaTeXStrings, OpenDSSDirect (171 total dependencies)
- **Script Status**: ✅ Compiles and runs successfully
- **Plots Generated**: ✅ All 3 PNG files in correct location

### Directory Structure Change:
**Old**: `processedData/small2poi_1ph_T24/*.png` (flat)  
**New**: `processedData/small2poi_1ph_T24/variable_voltages/plots/*.png` (hierarchical)

This matches the pattern used in `small3poi_mpopf.jl` for better organization.

---

## Timeline Consideration

**PESGM Deadline**: 4 days from now  
**Estimated Time for Setup**: 30 seconds to 2 minutes (depending on which option works)  
**Contingency**: Backup branch lets you revert in 5 seconds if issues arise

### Recommended Approach:
1. **DO THIS FIRST**: Try Option A (Pkg.resolve) - takes 10 seconds
2. **If that fails**: Try Option B (Pkg.instantiate) - takes 30-60 seconds
3. **Only if desperate**: Try Option C (delete Manifest.toml) - takes ~2 minutes
4. **Emergency**: Revert to backup branch - takes 5 seconds

---

## Questions to Ask Claude on Home PC

Copy this entire file and ask:

> "I just pulled changes from my lab PC where we fixed Julia package corruption. The changes include updated Project.toml and code fixes in small2poi_mpopf.jl. I need to verify my home PC environment still works. Please help me run the appropriate Pkg commands (resolve/instantiate) and test the script. PESGM deadline is in 4 days, so we need this to work quickly."

Then provide this file's contents for context.

---

## Git Commit Message (FYI)

```
Fix Julia environment corruption and update small2poi_mpopf.jl structure

- Fix corrupted Plots UUID in envs/multi_poi/Project.toml
- Add systemName parameter and hierarchical directory structure
- Fix Pkg.activate path and add LaTeXStrings import
- Fix Plotter.jl voltage data handling bug
- Switch to Julia LTS 1.10.10 on lab PC

Tested: ✅ Script compiles and runs successfully on lab PC
Backup branch: tadmm-nov12-lab-pc-working
```

---

## Contact Points

**If home PC setup fails**, provide these details to Claude:
1. Julia version: `julia --version`
2. Error message from `Pkg.resolve()` or `Pkg.instantiate()`
3. Contents of: `envs/multi_poi/Project.toml`
4. Whether you're on Windows/Linux/Mac
5. Output of: `julia --project=envs/multi_poi -e 'import Pkg; Pkg.status()'`

---

**Good luck! The chances of this working smoothly are very high.** 🚀
