# DSS-Based Parsing Functions

This document describes the new DSS-based parsing approach for extracting system configuration data directly from OpenDSS.

## Overview

Previously, the parser required dedicated `.dss` files (like `BranchData.dss`, `Loads.dss`, `PVsystem.dss`, `Storage.dss`) with specific formatting. The new approach:

1. **Loads the system directly into OpenDSS** using the `Master.dss` file
2. **Extracts data via the OpenDSS API** rather than parsing text files
3. **Maintains compatibility** with existing MPOPF model builder

## Key Advantages

- ✅ **Single source of truth**: Only need to maintain one `Master.dss` file per system configuration
- ✅ **Flexibility**: Easy to change PV percentage, battery configurations without multiple file versions
- ✅ **Robustness**: Uses OpenDSS's internal data structures rather than parsing text
- ✅ **Validation**: Can run power flow to validate the circuit before optimization

## Functions

### 1. `load_system_in_dss(systemName; kwargs...)`

Loads a system into the OpenDSS engine.

**Purpose**: Simple but fundamental function that clears OpenDSS and loads the specified system's `Master.dss` file.

**Arguments**:
- `systemName::String`: Name of the system (e.g., `"ads10A_1ph"`, `"ieee123_1ph"`)
- `rawDataFolderPath::Union{Nothing, String}`: Path to rawData folder (optional)
- `runPowerFlow::Bool`: Whether to run power flow after loading (default: `false`)

**Returns**: 
- `dss_file::String`: Full path to the loaded Master.dss file

**Example**:
```julia
using MultiPeriodDistOPF

# Load system into OpenDSS
dss_file = load_system_in_dss("ads10A_1ph")
# Output: "Loaded system: ads10A_1ph from c:/Users/.../rawData/ads10A_1ph/Master.dss"
```

### 2. `get_system_config_from_dss(systemName, T; kwargs...)`

Extracts all system configuration data required for MPOPF optimization.

**Purpose**: The main workhorse function that loads the system and extracts all necessary data structures.

**Arguments**:
- `systemName::String`: Name of the system
- `T::Int`: Number of time steps for simulation

**Optional Keyword Arguments**:
- `numAreas::Int`: Number of areas (default: 1)
- `objfun0::String`: Primary objective ("subsPowerCostMin", "lineLossMin", "subsPowerMin")
- `objfun2::String`: Secondary objective (default: "scd")
- `temporal_decmp::Bool`: Enable temporal decomposition (default: false)
- `algo_temporal_decmp::String`: Decomposition algorithm ("DDP", "tENApp")
- `solver::String`: Solver name ("Ipopt", "Gurobi", "Juniper")
- `linearizedModel::Bool`: Use linearized model (default: false)
- `tSOC_hard::Bool`: Hard terminal SOC constraints (default: false)
- `relax_terminal_soc_constraint::Bool`: Relax terminal SOC (default: false)
- And more...

**Returns**: 
- `data::Dict`: Dictionary with all MPOPF parameters

**Example**:
```julia
# Extract data for 24-hour optimization
data = get_system_config_from_dss("ads10A_1ph", 24,
    linearizedModel=false,
    temporal_decmp=false,
    solver="Ipopt",
    tSOC_hard=false,
    relax_terminal_soc_constraint=true)

# Access extracted data
@unpack Nset, Lset, NLset, Dset, Bset = data
println("Buses: $(length(Nset)), Branches: $(length(Lset))")
println("Loads: $(length(NLset)), PVs: $(length(Dset)), Batteries: $(length(Bset))")
```

## Extracted Data Structures

The `get_system_config_from_dss` function extracts the following data:

### Network Topology
- `Nset`: Vector of all bus numbers
- `Lset`: Vector of all branches (tuples `(from_bus, to_bus)`)
- `LTset`: Transformer branches
- `LnotTset`: Non-transformer branches
- `N`, `m`: Counts of buses and branches
- `parent`, `children`: Tree structure dictionaries

### Substation Sets
- `N1set`: Substation bus (typically `[1]`)
- `Nm1set`: Non-substation buses
- `Nc1set`: Buses directly connected to substation
- `Nnc1set`: Buses not directly connected to substation
- `L1set`: Branches connected to substation
- `Lm1set`: Branches not connected to substation

### Branch Impedances
- `rdict`, `xdict`: Resistance and reactance in Ohms
- `rdict_pu`, `xdict_pu`: Per-unit impedances
- `Z_B_dict`: Base impedance per branch

### Loads
- `NLset`: Buses with loads
- `N_L`: Number of loads
- `p_L_R`, `q_L_R`: Rated active/reactive power (kW, kVAr)
- `p_L_R_pu`, `q_L_R_pu`: Per-unit rated powers
- `p_L[(j,t)]`, `q_L[(j,t)]`: Time-varying load profiles
- `p_L_pu[(j,t)]`, `q_L_pu[(j,t)]`: Per-unit time-varying profiles
- `Vminpu_L`, `Vmaxpu_L`: Voltage limits per load bus

### PV Systems (DERs)
- `Dset`: Buses with PV systems
- `n_D`, `DER_percent`: Count and percentage of PVs
- `p_D_R`, `S_D_R`: Rated active/apparent power
- `p_D_R_pu`, `S_D_R_pu`: Per-unit rated powers
- `p_D[(j,t)]`, `p_D_pu[(j,t)]`: Time-varying PV generation profiles
- `irrad`: Irradiance values
- `Vminpu_D`, `Vmaxpu_D`: Voltage limits per PV bus

### Batteries (Storage)
- `Bset`: Buses with batteries
- `n_B`, `Batt_percent`: Count and percentage of batteries
- `B0`, `B0_pu`: Initial state of charge (kWh, pu)
- `B_R`, `B_R_pu`: Rated energy capacity (kWh, pu)
- `P_B_R`, `P_B_R_pu`: Rated power (kW, pu)
- `S_B_R`, `S_B_R_pu`: Rated apparent power (kVA, pu)
- `eta_C`, `eta_D`: Charging/discharging efficiencies
- `soc_min`, `soc_max`, `soc_0`: SOC parameters (%)
- `Bref`, `Bref_pu`, `Bref_percent`: Reference SOC values
- `Vminpu_B`, `Vmaxpu_B`: Voltage limits per battery bus

### Consolidated Voltage Limits
- `Vminpu[bus]`, `Vmaxpu[bus]`: Per-bus voltage limits (maximum of component minimums, minimum of component maximums)

### Base Values
- `kVA_B`, `kV_B`, `MVA_B`, `Z_B`: System-wide base values
- `kVA_B_dict[bus]`, `kV_B_dict[bus]`: Per-bus base values
- `MVA_B_dict`, `Z_B_dict`: Dictionaries for all buses/branches

### Time Parameters
- `T`: Number of time steps
- `Tset`: Vector `[1, 2, ..., T]`
- `delta_t`: Time step size (hours)

### Solver & Objective Configuration
- `solver`: Solver name
- `objfun0`, `objfun2`: Objective function types
- `objfunString`, `objfunUnit`: Display strings
- `linearizedModel`: Boolean flag
- And many more metadata parameters...

## Usage Examples

### Example 1: Basic Usage
```julia
include("./src/setupMultiPeriodDistOPF.jl")

# Parse system data
systemName = "ads10A_1ph"
T = 24
data = Parser.get_system_config_from_dss(systemName, T)

# Build and solve MPOPF model
modelDict = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data)
optimize!(modelDict[:model])
```

### Example 2: With Temporal Decomposition
```julia
# Use DDP with linearized model
data = Parser.get_system_config_from_dss("ads10A_1ph", 24,
    linearizedModel=true,
    temporal_decmp=true,
    algo_temporal_decmp="DDP",
    solver="Gurobi")

# Solve with DDP
modelDict = DDPLinear.optimize_MPOPF_1ph_L_DDP(data, maxiter=100)
```

### Example 3: Custom Configuration
```julia
# Custom DER/battery configuration
gedDict = Dict(
    :systemName0 => "ads10A_1ph",
    :DER_Percent_ud => 30,
    :DER_Rating_factor_ud => 1.0,
    :Batt_Percent_ud => 25,
    :Batt_Rating_factor_ud => 0.5
)

data = Parser.get_system_config_from_dss("ads10A_1ph", 24,
    gedDict_ud=gedDict,
    tSOC_hard=true,
    relax_terminal_soc_constraint=false)
```

## Migration from Old Parser

**Old approach** (parsing dedicated `.dss` files):
```julia
data = Parser.parse_all_data(systemName, T,
    numAreas=1,
    solver="Ipopt",
    ...)
```

**New approach** (loading from OpenDSS engine):
```julia
data = Parser.get_system_config_from_dss(systemName, T,
    numAreas=1,
    solver="Ipopt",
    ...)
```

Both return the same `data` dictionary structure, so existing code using `data` remains unchanged!

## Testing

Run the test script to verify functionality:
```bash
julia test_dss_parser.jl
```

This will:
1. Load the system into OpenDSS
2. Extract all configuration data
3. Verify data structures
4. Test compatibility with model builder
5. Optionally build an MPOPF model

## Known Limitations

1. **Voltage levels**: Currently uses system-wide base voltage. Proper transformer voltage tracking is TODO.
2. **DER/Battery percentages**: Extracted as 100% by default; needs enhancement for automatic calculation.
3. **Load shapes**: Uses default load shape files; custom shapes need to be specified.
4. **Multi-phase**: Only single-phase systems are currently supported.

## Future Enhancements

- [ ] Implement proper voltage level tracking through transformers
- [ ] Auto-calculate DER and battery percentages from circuit
- [ ] Support for multi-phase systems
- [ ] Extract custom load shapes from OpenDSS
- [ ] Parallel system loading for multi-area problems
- [ ] Support for other OpenDSS elements (capacitors, regulators, etc.)

## File Location

- **New module**: `src/Parser/parseFromDSS.jl`
- **Test script**: `test_dss_parser.jl`
- **Exports**: Added to `src/Parser/parseOpenDSSFiles.jl` and `src/MultiPeriodDistOPF.jl`

---

**Note**: The old file-based parsing functions (`parse_all_data`, `parse_branch_data`, etc.) are still available and functional. Both approaches can coexist, giving you flexibility in how you set up your system data.
