# Status: Large 10k Bus System Conversion (2026-03-20)

## Goal
Test tADMM on a truly large system (10k buses) to see if "going bigger" makes tADMM win over brute force.

## What We Did

### 1. Identified Source Data
- **Location**: `FlexibleParams-DistOPF/` (separate repo, sibling directory)
- **Structure**: 101 areas with ~100 buses each
  - Main backbone: `linedata.txt` (100 buses)
  - Area data: `Area Data/Area1-Area101/linedata.txt` (~100 buses per area)
  - Cross-boundary: `CB_full.txt` (connections from backbone to areas)
- **Total size**: ~10,340 buses, 10,319 lines

### 2. Created Converter Script
- **File**: `rawData/convert_flexparams_10k_to_opendss.jl`
- **What it does**:
  - Reads main feeder + 101 area linedata files
  - Reads CB_full.txt for cross-boundary connections
  - Generates OpenDSS format: Master.dss, BranchData.dss, Loads.dss, PVSystem.dss
  - Bus numbering: Main feeder (1-100), Area X bus Y → X*1000+Y (e.g., Area 101 bus 5 → 101005)

### 3. Generated System
- **Output**: `rawData/large10k_1ph/`
- **Statistics**:
  - 10,340 buses
  - 10,319 lines
  - 10,339 loads (0.1 MW, 0.01 MVAr each)
  - 1,020 PV systems (0.07 MW, 50% penetration)
  - Base: 1000 kVA, 12.47 kV L-L, 7.1996 kV L-N

## Current Problem

### Error When Running
```
ERROR: LoadError: KeyError: key 1001 not found
  @ tadmm_socp.jl:281
```

### Root Cause
The network topology is **non-radial**:
- CB_full.txt creates direct connections from backbone bus 1 to many area buses
- Example: Bus 1 connects to 101002, 102007, 103012, etc.
- The OpenDSS parser (`parse_opendss.jl`) expects a **radial tree** structure
- Parser builds `parent[]` dictionary assuming each bus has exactly one parent
- Bus 1001 exists in BranchData.dss but isn't in the parent dictionary

### Topology Understanding (Unclear)
The original MATLAB system structure from `FlexibleParams-DistOPF` isn't well-documented:
- CB_full.txt format: `[backbone_bus, area_local_bus, area_id]`
- Not clear if this represents actual physical topology or just spatial decomposition for distributed OPF
- Area internal structure: radial trees within each area
- Cross-boundary connections: might be creating multiple roots or mesh topology

## Configuration Changes Made

### tadmm_socp.jl
```julia
# Line 89-93
systemName = "large10k_1ph"
T = 6  # Number of time steps

# Line 1970 (disabled tADMM for initial test)
if false && !isempty(data[:Bset])  # Only run tADMM if there are batteries
```

## Next Steps (When Resuming)

### Option 1: Fix the Topology
- Investigate actual FlexibleParams-DistOPF network structure
- May need to simplify cross-boundary connections to maintain radial tree
- Or: modify parser to handle non-radial (mesh) networks

### Option 2: Use a Different Large System
- IEEE 8500-node system already in `rawData/8500-Node/` (OpenDSS format, native)
- But: need to add batteries and configure for MPOPF
- Advantage: Already in working OpenDSS format

### Option 3: Synthetic Radial System
- Generate a large synthetic radial feeder (10k buses)
- Simple repeating pattern with loads, PV, and batteries
- Guaranteed to work with existing parser

## Files Created/Modified

### New Files
- `rawData/convert_flexparams_10k_to_opendss.jl` - Converter script
- `rawData/large10k_1ph/Master.dss` - OpenDSS master file
- `rawData/large10k_1ph/BranchData.dss` - Network topology (623 KB)
- `rawData/large10k_1ph/Loads.dss` - Load definitions (1.3 MB)
- `rawData/large10k_1ph/PVSystem.dss` - PV systems (134 KB)
- `rawData/large10k_1ph/Storage.dss` - Battery template
- `STATUS_large10k_system.md` - This file

### Modified Files
- `tadmm_socp.jl`:
  - systemName = "large10k_1ph"
  - T = 6
  - Disabled tADMM (line 1972)

## Key Insight from User
> "Yeah it's possible you overestimated the systematic-ness of my assortment of files lol"

The FlexibleParams-DistOPF data may not have clean documentation about actual network topology vs. decomposition structure for distributed optimization algorithms.

## References
- Source data: `../FlexibleParams-DistOPF/`
- Original MATLAB code: `../FlexibleParams-DistOPF/NL_OPF_dist.m`
- Network generation: `../FlexibleParams-DistOPF/Large_net_form.m`
- Notes file: `../FlexibleParams-DistOPF/Notes.txt` ("CB is really just a collection of 9 Areas... like a 1000 bus network")
