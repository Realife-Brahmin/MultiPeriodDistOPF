# Instructions for Home PC

## Changes Made on Lab PC (Nov 12, 2025)

### What Changed:
1. **Project.toml**: Added Plots package with correct UUID
2. **small2poi_mpopf.jl**: 
   - Added `using LaTeXStrings`
   - Updated directory structure for plots
3. **Julia Version**: Switched to LTS 1.10.10 on lab PC

### When You Pull These Changes on Home PC:

**If everything works fine:**
- Just run the script normally
- Julia will update packages automatically

**If you get package errors:**

```bash
# Option 1: Let Julia resolve it automatically
julia --project=envs/multi_poi -e 'import Pkg; Pkg.resolve()'
julia small2poi_mpopf.jl

# Option 2: Fresh instantiate (if Option 1 fails)
julia --project=envs/multi_poi -e 'import Pkg; Pkg.instantiate()'
julia small2poi_mpopf.jl

# Option 3: Nuclear option (if both above fail)
rm envs/multi_poi/Manifest.toml
julia --project=envs/multi_poi -e 'import Pkg; Pkg.instantiate()'
```

### Backup:
If things break, you can always revert to the pre-change state:
```bash
git checkout tadmm-nov12-lab-pc-working
```

### Lab PC Was Fixed By:
- Corruption in OpenSpecFun_jll (riscv64 platform issue)
- Reset entire Julia environment
- Switched to Julia LTS 1.10.10
- Fresh package installation
