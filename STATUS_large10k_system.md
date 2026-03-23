# Status: Large 10k Bus System Conversion (2026-03-22)

## Goal
Test tADMM on a truly large system (10k buses) to see if "going bigger" makes tADMM win over brute force.

## SUCCESS: 10,340 Bus System Generated ✓

### Final System Statistics
- **Buses**: 10,340
- **Lines**: 10,319
- **Loads**: 10,339 (0.1 MW, 0.01 MVAr each)
- **PV Systems**: 1,022 (0.07 MW, ~10% penetration)
- **Base**: 1000 kVA, 12.47 kV L-L, 7.1996 kV L-N

### Source Data Structure
- **Location**: `../FlexibleParams-DistOPF/`
- **Main backbone**: `linedata.txt` (100 lines, 101 buses)
- **Area data**: `Area Data/Area1-Area101/linedata.txt` (101 areas, 10,119 total lines)
- **Cross-boundary**: `CB_full.txt` (100 connections)

**Key Discovery**: CB_full.txt area IDs (101-120) are decomposition labels for distributed OPF, NOT directory indices. The actual network is built from all 101 Area directories.

## Conversion History

### First Attempt (FAILED)
- **File**: `rawData/convert_flexparams_10k_to_opendss.jl`
- **Problem**: Only loaded areas mentioned in CB_full.txt (Area101-Area120)
- **Result**: Only Area101 existed → 221 buses instead of 10k
- **Error**: KeyError on bus 1001 (radial tree parser couldn't find parent)

### Second Attempt (TOO SMALL)
- **File**: `rawData/convert_simple_10k.jl`
- **Approach**: Just use main backbone linedata.txt, ignore areas
- **Result**: Only 101 buses (not 10k!)

### Third Attempt (SUCCESS) ✓
- **File**: `rawData/convert_all_areas_10k.jl`
- **Approach**: Load ALL 101 area directories + main backbone + CB connections
- **Bus numbering**: Area X, local bus Y → global bus X*1000+Y
- **Result**: 10,340 buses as expected!

## Generated Files

```
rawData/large10k_1ph/
├── Master.dss          - Circuit definition
├── BranchData.dss      - 10,319 lines (100 backbone + 100 CB + 10,119 area)
├── Loads.dss           - 10,339 loads
├── PVSystem.dss        - 1,022 PV systems (every 10th bus)
└── Storage.dss         - Empty template (can add batteries later)
```

## Current Configuration

### tadmm_socp.jl
```julia
systemName = "large10k_1ph"
T = 6  # Number of time steps
```

## Current Status: Testing BF T=6

**Running**: Brute force MPOPF with T=6 to verify system is well-formed
- If this works: we have a valid 10k system
- Next step: Add batteries and test tADMM to see if large system makes tADMM win

## Files in Repo

### Converter Scripts
- `rawData/convert_flexparams_10k_to_opendss.jl` - First attempt (wrong)
- `rawData/convert_simple_10k.jl` - Second attempt (too small)
- `rawData/convert_all_areas_10k.jl` - Final version (SUCCESS) ✓

### System Files
- `rawData/large10k_1ph/` - Generated OpenDSS system (10,340 buses)
- `STATUS_large10k_system.md` - This status file

### Modified Files
- `tadmm_socp.jl` - Set to use large10k_1ph, T=6

## Key Insights

1. **FlexibleParams structure**: CB_full.txt uses area_id 101-120 as decomposition labels, not directory references
2. **Actual network**: Must load all 101 Area directories to get the full 10k bus system
3. **Notes.txt was misleading**: Said "CB is 9 areas, ~1000 buses" but CB_full has 20 area labels and full system is 10k
4. **User was right**: "You overestimated the systematic-ness of my assortment of files" 😅

## References
- Source: `../FlexibleParams-DistOPF/`
- Notes: `../FlexibleParams-DistOPF/Notes.txt`
