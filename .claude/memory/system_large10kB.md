---
name: large10kB System
description: 10321-bus system variant with 1020 batteries (10 per area) for high temporal coupling tests
type: project
---

Created 2026-03-31 to test hypothesis that more batteries increase tADMM advantage.

- **Buses**: 10,321 (same as large10k_1ph)
- **Batteries**: 1,020 (10 per area × 102 areas, vs 102 in original)
- **PVs**: 1,022 (same as large10k_1ph)
- **Files**: `rawData/large10kB_1ph/` (copied BranchData/Loads/PVSystem from large10k_1ph, new Storage.dss)
- **Storage.dss**: Generated with Julia, seed=42, kVA cycles [80,90,100,110,120]
- **Results**: `envs/tadmm/processedData/large10kB_1ph_T{4,6,12}/`

**Why:** Original 102 batteries showed tADMM advantage only at T=4. Need more temporal coupling to test if tADMM scales.
**How to apply:** Use `systemName = "large10kB_1ph"` in tadmm_socp.jl.
