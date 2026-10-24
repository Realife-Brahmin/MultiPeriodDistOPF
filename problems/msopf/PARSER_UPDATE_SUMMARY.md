# Parser Update Summary: Network Sets for Multi-POI MPOPF

## Date
November 13, 2025

## Objective
Update the OpenDSS parser to correctly construct network sets for Multi-Period Optimal Power Flow (MPOPF) formulation with multiple Points of Interconnection (POI).

## Problem Statement
The previous parser created network sets that only contained distribution buses, excluding substation buses. For MPOPF formulation, we need:
- **Substation buses** (slack bus candidates) included in network sets
- **POI connection lines** included in line sets
- **Proper separation** between different types of buses and lines
- **Mixed data types** (String for substations, Int for distribution buses)

## Solution Implemented

### New Network Sets Structure

#### Bus Sets
1. **`Sset`** - Substation buses only
   - Type: `Vector{String}`
   - Example: `["1s", "2s", "3s", "4s", "5s"]`
   - Used for: Slack bus selection, angle reference

2. **`Nset`** - ALL buses (substations + distribution)
   - Type: `Vector{Any}` (mixed String and Int)
   - Example: `["1s", "2s", "3s", "4s", "5s", 1, 2, 3, ..., 128]`
   - Used for: Voltage variables, power balance for all buses

3. **`Nm1set`** - Distribution buses only (excluding substations)
   - Type: `Vector{Int}`
   - Example: `[1, 2, 3, ..., 128]`
   - Used for: Array sizing for loads/PV/batteries (indexed by Int only)

#### Line Sets
1. **`Lset`** - ALL lines (POI + distribution)
   - Type: `Vector{Tuple{Any, Any}}`
   - Example: `[(1s, 1), (2s, 33), ..., (1, 2), (2, 3), ...]`
   - Used for: All power flow equations, unified impedance access

2. **`L1set`** - POI connection lines only
   - Type: `Vector{Tuple{Any, Any}}`
   - Example: `[(1s, 1), (2s, 33), (3s, 120), (4s, 116), (5s, 85)]`
   - Used for: POI power exchange, angle coordination
   - Note: First element is substation (String), second is distribution bus (Int)

3. **`Lm1set`** - Distribution lines only (excluding POI)
   - Type: `Vector{Tuple{Int, Int}}`
   - Example: `[(1, 2), (2, 3), ..., (127, 128)]`
   - Used for: Standard branch flow model constraints

### Modified Functions

#### 1. `parse_network_topology!`
**Changes:**
- Creates `Sset` for substation buses
- Creates `Nset` as `Vector{Any}` combining substations and distribution buses
- Creates `Nm1set` for distribution buses only
- Parses ALL lines (POI + distribution) into appropriate sets
- Handles mixed bus name types (String for substations, Int for distribution)

**Key Code:**
```julia
# Create Sset: substation bus names
Sset = substation_buses  # ["1s", "2s", "3s", "4s", "5s"]

# Create Nset: ALL buses (substations + distribution)
Nset = Vector{Any}(vcat(substation_buses, dist_bus_numbers))

# Create Nm1set: only distribution buses
Nm1set = dist_bus_numbers  # [1, 2, 3, ..., 128]

# Parse lines and classify as POI (L1set) or distribution (Lm1set)
for line_name in line_names
    # Determine bus types and create line tuple
    is_poi_line = startswith(lowercase(line_name), "poi")
    if is_poi_line
        push!(L1set, line_tuple)
    else
        push!(Lm1set, line_tuple)
    end
end
```

#### 2. `parse_line_impedances!`
**Changes:**
- Uses `Dict{Any, Float64}` for unified impedance dictionaries
- Stores impedances for ALL lines (POI + distribution) in `rdict_pu` and `xdict_pu`
- Maintains separate POI impedance dicts for backward compatibility

**Key Code:**
```julia
# Unified dictionaries for ALL lines
rdict_pu = Dict{Any, Float64}()
xdict_pu = Dict{Any, Float64}()

for line_name in line_names
    # Determine bus types
    bus1 = endswith(bus1_base, "s") ? bus1_base : parse(Int, bus1_base)
    bus2 = endswith(bus2_base, "s") ? bus2_base : parse(Int, bus2_base)
    
    # Store in unified dictionary
    line_tuple = (bus1, bus2)
    rdict_pu[line_tuple] = r_pu
    xdict_pu[line_tuple] = x_pu
end
```

#### 3. `parse_loads!` and `parse_pv_generators!`
**Changes:**
- Use `Nm1set` instead of `Nset` for array sizing
- Arrays remain indexed by Int only (distribution buses)

**Key Code:**
```julia
# Use Nm1set (distribution buses only) to determine array size
Nm1set = data[:Nm1set]
max_bus = maximum(Nm1set)
p_L_pu = zeros(max_bus, T)
```

#### 4. `parse_voltage_limits!`
**Changes:**
- Uses `Dict{Any, Float64}` for voltage limit dictionaries
- Creates separate substation voltage limit dicts for clarity

**Key Code:**
```julia
# Use Dict{Any, Float64} to handle both String and Int bus names
Vminpu = Dict{Any, Float64}()
Vmaxpu = Dict{Any, Float64}()

# Default voltage limits for all buses
for bus in data[:Nset]
    Vminpu[bus] = 0.95
    Vmaxpu[bus] = 1.05
end

# Also store limits specifically for substation buses
Vminpu_sub = Dict{String, Float64}()
Vmaxpu_sub = Dict{String, Float64}()
```

## Verification

### Test Results
```
--- Network Topology ---
  Substation buses (Sset): 5
  Distribution buses (Nm1set): 128 (1 to 128)
  Total buses (Nset): 133 (substations + distribution)
  POI lines (L1set): 5
  Distribution lines (Lm1set): 127
  Total lines (Lset): 132

--- Network Sets for MPOPF ---
  Sset (substations): 5 buses
  Nset (all buses): 133 buses (substations + distribution)
  Nm1set (distribution only): 128 buses
  Lset (all lines): 132 lines
  L1set (POI lines): 5 lines
  Lm1set (distribution lines): 127 lines
```

### Relationships Verified
- `|Nset| = |Sset| + |Nm1set| = 5 + 128 = 133` ✓
- `|Lset| = |L1set| + |Lm1set| = 5 + 127 = 132` ✓
- POI connections: `1s→1, 2s→33, 3s→120, 4s→116, 5s→85` ✓
- All impedances accessible via unified dictionaries ✓

## Files Modified

1. **`parse_opendss_multi_poi.jl`**
   - `parse_network_topology!()` - Complete rewrite for proper set construction
   - `parse_line_impedances!()` - Updated to handle mixed bus types
   - `parse_loads!()` - Changed to use `Nm1set` for array sizing
   - `parse_pv_generators!()` - Changed to use `Nm1set` for array sizing
   - `parse_voltage_limits!()` - Updated to use `Dict{Any, Float64}`

2. **`test_parser_multi_poi.jl`**
   - Updated to display new network sets
   - Fixed voltage limit printing
   - Added verification of all key data structures

3. **`multi_poi_mpopf.jl`**
   - No changes needed (parser provides all required data)

## New Data Dictionary Keys

- `:Sset` - Substation bus names
- `:Nset` - All buses (substations + distribution)
- `:Nm1set` - Distribution buses only
- `:Lset` - All lines
- `:L1set` - POI lines only
- `:Lm1set` - Distribution lines only
- `:m` - Total number of lines
- `:m_poi` - Number of POI lines
- `:m_dist` - Number of distribution lines
- `:num_dist_buses` - Count of distribution buses
- `:Vminpu_sub` - Voltage limits for substations (separate dict)
- `:Vmaxpu_sub` - Voltage limits for substations (separate dict)

## Usage in MPOPF Formulation

### Example: Define Voltage Variables for All Buses
```julia
# All buses including substations
for i in data[:Nset]
    @variable(model, data[:Vminpu][i] <= V[i] <= data[:Vmaxpu][i])
end
```

### Example: POI Power Flow
```julia
# POI lines connect substations to distribution buses
for (i, j) in data[:L1set]
    # i is String (substation), j is Int (distribution bus)
    r_ij = data[:rdict_pu][(i, j)]
    x_ij = data[:xdict_pu][(i, j)]
    @constraint(model, P_ij == ...)
end
```

### Example: Distribution Network Power Flow
```julia
# Distribution lines connect distribution buses
for (i, j) in data[:Lm1set]
    # Both i and j are Int
    @constraint(model, P_ij + r_ij * l_ij == P_j + P_D_j - p_L_j)
end
```

### Example: Load/PV Array Access
```julia
# Arrays indexed by distribution bus number (Int only)
for i in data[:Nm1set]
    for t in 1:T
        P_load = data[:p_L_pu][i, t]  # Correct indexing
    end
end
```

## Benefits

1. **Correct MPOPF Formulation**: Network sets now match mathematical formulation requirements
2. **Type Safety**: Clear separation between String (substation) and Int (distribution) bus names
3. **Unified Access**: Single impedance dictionary works for all lines
4. **Flexibility**: Easy to loop over specific bus/line types as needed
5. **Clarity**: Explicit naming (Sset, Nm1set, L1set, Lm1set) makes code intent clear

## Next Steps

1. Implement multi-POI MPOPF optimization model using these network sets
2. Handle slack bus selection from `Sset`
3. Implement angle coordination across 5 substations
4. Verify power flow results against OpenDSS baseline

## Documentation
- Network sets documentation: `NETWORK_SETS_DOCUMENTATION.md`
- This summary: `PARSER_UPDATE_SUMMARY.md`
