module parseFromDSS

export 
    load_system_in_dss,
    get_system_config_from_dss

using OpenDSSDirect
using Parameters: @unpack, @pack!

include("../helperFunctions.jl")
import .helperFunctions as HF

#region load_system_in_dss
"""
    load_system_in_dss(systemName::String; rawDataFolderPath=nothing, runPowerFlow=false)

Load a system into the OpenDSS engine by redirecting to its Master.dss file.

This function clears any existing circuit in OpenDSS and loads the specified system
from its Master.dss file. Optionally, it can run a power flow simulation after loading.

# Arguments
- `systemName::String`: The name of the system to load (e.g., "ads10A_1ph", "ieee123_1ph").
- `rawDataFolderPath::Union{Nothing, String}`: Path to the rawData folder. If `nothing`, 
    it will be constructed relative to this file's location (default: nothing).
- `runPowerFlow::Bool`: Whether to run a power flow after loading the system (default: false).

# Returns
- `dss_file::String`: The full path to the loaded Master.dss file.

# Example
```julia
dss_file = load_system_in_dss("ads10A_1ph")
```
"""
function load_system_in_dss(systemName::String; 
    rawDataFolderPath=nothing, 
    runPowerFlow=false)
    
    # Construct path to Master.dss
    if isnothing(rawDataFolderPath)
        wd = @__DIR__
        rawDataFolderPath = joinpath(wd, "..", "..", "rawData")
    end
    
    dss_dir = joinpath(rawDataFolderPath, systemName)
    dss_file = joinpath(dss_dir, "Master.dss")
    
    # Check if file exists
    if !isfile(dss_file)
        error("Master.dss file not found at: $dss_file")
    end
    
    # Clear any existing circuit and load the new one
    OpenDSSDirect.Text.Command("Clear")
    OpenDSSDirect.Text.Command("Redirect \"$dss_file\"")
    
    # Optionally run power flow
    if runPowerFlow
        OpenDSSDirect.Text.Command("Solve")
    end
    
    println("Loaded system: $systemName from $dss_file")
    
    return dss_file
end
#endregion

#region get_system_config_from_dss
"""
    get_system_config_from_dss(systemName::String, T::Int; kwargs...)

Extract all system configuration data from OpenDSS for MPOPF optimization.

This function loads the system into OpenDSS and extracts all necessary data including
buses (Nset), branches (Lset), batteries (Bset), PVs (Dset), loads (NLset), and various
parameters required for the MPOPF solver.

# Arguments
- `systemName::String`: The name of the system to parse.
- `T::Int`: The number of time steps for the simulation.

## Optional Keyword Arguments (for metadata)
- `numAreas::Int`: Number of areas in the system (default: 1).
- `objfun0::String`: Primary objective function type (default: "subsPowerCostMin").
- `objfun2::String`: Secondary objective function type (default: "scd").
- `temporal_decmp::Bool`: Enable temporal decomposition (default: false).
- `algo_temporal_decmp::String`: Temporal decomposition algorithm (default: "DDP").
- `PSubsMax_kW::Float64`: Maximum substation power in kW (default: Inf).
- `inputForecastDescription::String`: Input forecast description (default: "nonspecific").
- `solver::String`: Solver to use (default: "Ipopt").
- `tSOC_hard::Bool`: Hard terminal SOC constraints (default: false).
- `relax_terminal_soc_constraint::Bool`: Relax terminal SOC constraint (default: false).
- `linearizedModel::Bool`: Use linearized model (default: false).
- `gedDict_ud::Union{Nothing, Dict}`: User-defined DER/battery configuration (default: nothing).
- `alpha_fpi::Float64`: Alpha parameter for FPI (default: 0.43).
- `gamma_fpi::Float64`: Gamma parameter for FPI (default: 0.8).
- `warmStart_mu::Union{Nothing, String}`: Warm start strategy (default: nothing).
- `threshold_conv_iters::Int`: Convergence iteration threshold (default: 3).
- `kVA_B::Float64`: Base power in kVA (default: 1000).
- `kV_B::Float64`: Base voltage in kV (default: 2.4018).
- `rawDataFolderPath::Union{Nothing, String}`: Path to rawData folder (default: nothing).

# Returns
- `data::Dict`: A dictionary containing all parsed system configuration data including:
    - Network topology: `Nset`, `Lset`, `N`, `m`, `parent`, `children`
    - Substation sets: `N1set`, `Nm1set`, `Nc1set`, `Nnc1set`
    - Branch sets: `L1set`, `Lm1set`, `LTset`, `LnotTset`
    - Branch impedances: `rdict`, `xdict`, `rdict_pu`, `xdict_pu`
    - Loads: `NLset`, `N_L`, `p_L_R`, `p_L_R_pu`, `q_L_R`, `q_L_R_pu`, `p_L`, `q_L`
    - PVs: `Dset`, `n_D`, `DER_percent`, `p_D_R`, `p_D_R_pu`, `S_D_R`, `S_D_R_pu`, `p_D`, `p_D_pu`
    - Batteries: `Bset`, `n_B`, `Batt_percent`, `B0`, `B_R`, `P_B_R`, `S_B_R`, `eta_C`, `eta_D`
    - Voltage limits: `Vminpu_*`, `Vmaxpu_*` for loads, PVs, and batteries
    - Base values: `kVA_B`, `kV_B`, `MVA_B`, `Z_B`, `kVA_B_dict`, `kV_B_dict`, `Z_B_dict`
    - Time parameters: `T`, `Tset`, `delta_t`
    - Solver configuration and objective function parameters

# Example
```julia
data = get_system_config_from_dss("ads10A_1ph", 24)
```
"""
function get_system_config_from_dss(systemName::String, T::Int;
    numAreas=1,
    objfun0="subsPowerCostMin",
    objfun2="scd",
    temporal_decmp=false,
    algo_temporal_decmp="DDP",
    PSubsMax_kW=Inf,
    inputForecastDescription="nonspecific",
    solver="Ipopt",
    tSOC_hard=false,
    relax_terminal_soc_constraint=false,
    linearizedModel=false,
    gedDict_ud=nothing,
    alpha_fpi=0.43,
    gamma_fpi=0.8,
    warmStart_mu=nothing,
    threshold_conv_iters=3,
    kVA_B=1000,
    kV_B=2.4018,
    rawDataFolderPath=nothing)
    
    # Step 1: Load system into OpenDSS
    dss_file = load_system_in_dss(systemName; rawDataFolderPath=rawDataFolderPath, runPowerFlow=false)
    
    # Step 2: Extract data from loaded DSS circuit
    data = Dict{Symbol, Any}()
    
    # Extract base values and system parameters
    MVA_B = kVA_B / 1000
    Z_B = (kV_B^2) / MVA_B
    Tset = collect(1:T)
    
    @pack! data = systemName, T, Tset, kVA_B, kV_B, MVA_B, Z_B
    
    # Extract network topology (buses and branches)
    topology_data = extract_network_topology_from_dss()
    data = merge(data, topology_data)
    
    # Extract base values per bus/branch (handles transformers)
    @unpack Nset, Lset = data
    baseValuesDict = compute_base_values_from_dss(Nset, Lset, kVA_B, kV_B)
    data = merge(data, baseValuesDict)
    
    # Extract loads
    load_data = extract_loads_from_dss(T, baseValuesDict)
    data = merge(data, load_data)
    
    # Extract PV systems
    pv_data = extract_pv_from_dss(T, baseValuesDict)
    data = merge(data, pv_data)
    
    # Extract storage (batteries)
    battery_data = extract_batteries_from_dss(baseValuesDict)
    data = merge(data, battery_data)
    
    # Evaluate voltage limits across all components
    voltage_limits_data = evaluate_voltage_limits_from_components(data)
    data = merge(data, voltage_limits_data)
    
    # Extract simulation parameters
    sim_params = extract_simulation_parameters_from_dss()
    data = merge(data, sim_params)
    
    # Add metadata parameters
    metadata = create_metadata_dict(
        systemName, T, numAreas, objfun0, objfun2, temporal_decmp, algo_temporal_decmp,
        PSubsMax_kW, inputForecastDescription, solver, tSOC_hard, 
        relax_terminal_soc_constraint, linearizedModel, gedDict_ud,
        alpha_fpi, gamma_fpi, warmStart_mu, threshold_conv_iters
    )
    data = merge(data, metadata)
    
    # Add cost data
    cost_data = HF.generateBinaryLoadShape(T)
    data = merge(data, cost_data)
    
    # Add folder paths
    if isnothing(rawDataFolderPath)
        rootFolderPath = dirname(dirname(@__DIR__))
        rawDataFolderPath = joinpath(rootFolderPath, "rawData")
    else
        rootFolderPath = dirname(dirname(rawDataFolderPath))
    end
    processedDataFolderPath = joinpath(rootFolderPath, "processedData")
    srcFolderPath = joinpath(rootFolderPath, "src")
    
    @pack! data = processedDataFolderPath, rawDataFolderPath, rootFolderPath, srcFolderPath
    
    # Post-process data (compute derived quantities)
    data = post_process_dss_data(data)
    
    return data
end
#endregion

#region extract_network_topology_from_dss
"""
Extract network topology (buses, branches, parent-child relationships) from OpenDSS circuit.
"""
function extract_network_topology_from_dss()
    Nset = Set{Int}()
    Lset = Set{Tuple{Int,Int}}()
    LTset = Set{Tuple{Int,Int}}()  # Transformer lines
    LnotTset = Set{Tuple{Int,Int}}()  # Non-transformer lines
    
    rdict = Dict{Tuple{Int,Int}, Float64}()
    xdict = Dict{Tuple{Int,Int}, Float64}()
    
    parent = Dict{Int, Union{Int, Nothing}}()
    children = Dict{Int, Vector{Int}}()
    
    # Extract lines (non-transformer branches)
    line_id = OpenDSSDirect.Lines.First()
    while line_id > 0
        bus1_name = OpenDSSDirect.Lines.Bus1()
        bus2_name = OpenDSSDirect.Lines.Bus2()
        
        # Extract bus numbers (before any phase designation)
        bus1 = parse(Int, split(bus1_name, ".")[1])
        bus2 = parse(Int, split(bus2_name, ".")[1])
        
        # Add to sets
        push!(Nset, bus1, bus2)
        push!(Lset, (bus1, bus2))
        push!(LnotTset, (bus1, bus2))
        
        # Get impedances (in Ohms)
        r_matrix = OpenDSSDirect.Lines.RMatrix()
        x_matrix = OpenDSSDirect.Lines.XMatrix()
        
        # For single-phase, take first element
        rdict[(bus1, bus2)] = r_matrix[1]
        xdict[(bus1, bus2)] = x_matrix[1]
        
        # Update parent-child relationships
        parent[bus2] = bus1
        if !haskey(children, bus1)
            children[bus1] = Int[]
        end
        push!(children[bus1], bus2)
        
        get!(children, bus2, Int[])
        get!(parent, bus1, nothing)
        
        line_id = OpenDSSDirect.Lines.Next()
    end
    
    # Extract transformers
    xfmr_id = OpenDSSDirect.Transformers.First()
    while xfmr_id > 0
        # Get winding buses
        wdg_voltages = OpenDSSDirect.Transformers.WdgVoltages()
        
        # For 2-winding transformers
        if length(wdg_voltages) == 2
            buses = OpenDSSDirect.CktElement.BusNames()
            bus1 = parse(Int, split(buses[1], ".")[1])
            bus2 = parse(Int, split(buses[2], ".")[1])
            
            push!(Nset, bus1, bus2)
            push!(Lset, (bus1, bus2))
            push!(LTset, (bus1, bus2))
            
            # Get transformer parameters
            kva = OpenDSSDirect.Transformers.kVA()
            xhl = OpenDSSDirect.Transformers.Xhl()
            pct_r = OpenDSSDirect.Transformers.R()  # %R
            
            kv1 = wdg_voltages[1]
            kv2 = wdg_voltages[2]
            
            # Calculate resistances and reactances (referred to secondary)
            r_sec = (pct_r / 100) * kv2^2 / (kva / 1000)
            x_sec = (xhl / 100) * kv2^2 / (kva / 1000)
            
            rdict[(bus1, bus2)] = r_sec
            xdict[(bus1, bus2)] = x_sec
            
            # Update parent-child relationships
            parent[bus2] = bus1
            if !haskey(children, bus1)
                children[bus1] = Int[]
            end
            push!(children[bus1], bus2)
            
            get!(children, bus2, Int[])
            get!(parent, bus1, nothing)
        end
        
        xfmr_id = OpenDSSDirect.Transformers.Next()
    end
    
    # Create substation-related sets (assuming bus 1 is substation)
    N1set = Set{Int}([1])  # Substation bus
    Nm1set = setdiff(Nset, N1set)  # All buses except substation
    Nc1set = Set{Int}()  # Buses connected directly to substation
    L1set = Set{Tuple{Int,Int}}()  # Branches with substation
    Lm1set = Set{Tuple{Int,Int}}()  # Branches without substation
    
    for (i, j) in Lset
        if i == 1 || j == 1
            push!(L1set, (i, j))
            if i == 1
                push!(Nc1set, j)
            else
                push!(Nc1set, i)
            end
        else
            push!(Lm1set, (i, j))
        end
    end
    
    Nnc1set = setdiff(Nm1set, Nc1set)  # Non-substation buses not directly connected
    
    # Compute counts
    N = length(Nset)
    m = length(Lset)
    N1 = length(N1set)
    Nm1 = length(Nm1set)
    Nc1 = length(Nc1set)
    Nnc1 = length(Nnc1set)
    m1 = length(L1set)
    mm1 = length(Lm1set)
    
    # Sort sets
    Nset = sort(collect(Nset))
    Lset = sort(collect(Lset))
    LTset = sort(collect(LTset))
    LnotTset = sort(collect(LnotTset))
    N1set = sort(collect(N1set))
    Nm1set = sort(collect(Nm1set))
    Nc1set = sort(collect(Nc1set))
    Nnc1set = sort(collect(Nnc1set))
    L1set = sort(collect(L1set))
    Lm1set = sort(collect(Lm1set))
    
    return Dict(
        :Nset => Nset,
        :Lset => Lset,
        :LTset => LTset,
        :LnotTset => LnotTset,
        :rdict => rdict,
        :xdict => xdict,
        :parent => parent,
        :children => children,
        :N1set => N1set,
        :Nm1set => Nm1set,
        :Nc1set => Nc1set,
        :Nnc1set => Nnc1set,
        :L1set => L1set,
        :Lm1set => Lm1set,
        :N => N,
        :m => m,
        :N1 => N1,
        :Nm1 => Nm1,
        :Nc1 => Nc1,
        :Nnc1 => Nnc1,
        :m1 => m1,
        :mm1 => mm1
    )
end
#endregion

#region compute_base_values_from_dss
"""
Compute per-unit base values for each bus and branch, handling transformers.
"""
function compute_base_values_from_dss(Nset, Lset, kVA_B_sys, kV_B_sys)
    kVA_B_dict = Dict{Any, Float64}()
    kV_B_dict = Dict{Any, Float64}()
    MVA_B_dict = Dict{Any, Float64}()
    Z_B_dict = Dict{Tuple{Int,Int}, Float64}()
    
    rdict_pu = Dict{Tuple{Int,Int}, Float64}()
    xdict_pu = Dict{Tuple{Int,Int}, Float64}()
    
    # For now, use system-wide base values
    # TODO: Implement voltage-level tracking through transformers
    for bus in Nset
        kVA_B_dict[bus] = kVA_B_sys
        kV_B_dict[bus] = kV_B_sys
        MVA_B_dict[bus] = kVA_B_sys / 1000
    end
    
    for (i, j) in Lset
        kVA_B_dict[(i, j)] = kVA_B_sys
        kV_B_dict[(i, j)] = kV_B_sys
        MVA_B_dict[(i, j)] = kVA_B_sys / 1000
        Z_B_dict[(i, j)] = kV_B_sys^2 / (kVA_B_sys / 1000)
    end
    
    return Dict(
        :kVA_B_dict => kVA_B_dict,
        :kV_B_dict => kV_B_dict,
        :MVA_B_dict => MVA_B_dict,
        :Z_B_dict => Z_B_dict
    )
end
#endregion

#region extract_loads_from_dss
"""
Extract load data from OpenDSS circuit.
"""
function extract_loads_from_dss(T, baseValuesDict)
    @unpack kVA_B_dict = baseValuesDict
    
    NLset = Set{Int}()
    p_L_R = Dict{Int, Float64}()
    p_L_R_pu = Dict{Int, Float64}()
    q_L_R = Dict{Int, Float64}()
    q_L_R_pu = Dict{Int, Float64}()
    Vminpu_L = Dict{Int, Float64}()
    Vmaxpu_L = Dict{Int, Float64}()
    
    p_L = Dict{Tuple{Int,Int}, Float64}()
    p_L_pu = Dict{Tuple{Int,Int}, Float64}()
    q_L = Dict{Tuple{Int,Int}, Float64}()
    q_L_pu = Dict{Tuple{Int,Int}, Float64}()
    
    # Generate default load shape
    LoadShapeLoad = HF.generateLoadShape(T, filenameLoadShape="LoadShapeDefault.dss")
    
    load_id = OpenDSSDirect.Loads.First()
    while load_id > 0
        bus_names = OpenDSSDirect.CktElement.BusNames()
        bus = parse(Int, split(bus_names[1], ".")[1])
        push!(NLset, bus)
        
        # Get rated powers
        p_L_R[bus] = OpenDSSDirect.Loads.kW()
        q_L_R[bus] = OpenDSSDirect.Loads.kvar()
        
        p_L_R_pu[bus] = p_L_R[bus] / kVA_B_dict[bus]
        q_L_R_pu[bus] = q_L_R[bus] / kVA_B_dict[bus]
        
        # Voltage limits (using defaults)
        Vminpu_L[bus] = 0.90
        Vmaxpu_L[bus] = 1.10
        
        # Create time-varying profiles
        for t in 1:T
            p_L[(bus, t)] = p_L_R[bus] * LoadShapeLoad[t]
            p_L_pu[(bus, t)] = p_L_R_pu[bus] * LoadShapeLoad[t]
            q_L[(bus, t)] = q_L_R[bus] * LoadShapeLoad[t]
            q_L_pu[(bus, t)] = q_L_R_pu[bus] * LoadShapeLoad[t]
        end
        
        load_id = OpenDSSDirect.Loads.Next()
    end
    
    NLset = sort(collect(NLset))
    N_L = length(NLset)
    
    return Dict(
        :NLset => NLset,
        :N_L => N_L,
        :p_L_R => p_L_R,
        :p_L_R_pu => p_L_R_pu,
        :q_L_R => q_L_R,
        :q_L_R_pu => q_L_R_pu,
        :Vminpu_L => Vminpu_L,
        :Vmaxpu_L => Vmaxpu_L,
        :p_L => p_L,
        :p_L_pu => p_L_pu,
        :q_L => q_L,
        :q_L_pu => q_L_pu,
        :LoadShapeLoad => LoadShapeLoad
    )
end
#endregion

#region extract_pv_from_dss
"""
Extract PV system data from OpenDSS circuit.
"""
function extract_pv_from_dss(T, baseValuesDict)
    @unpack kVA_B_dict = baseValuesDict
    
    Dset = Set{Int}()
    p_D_R = Dict{Int, Float64}()
    p_D_R_pu = Dict{Int, Float64}()
    S_D_R = Dict{Int, Float64}()
    S_D_R_pu = Dict{Int, Float64}()
    irrad = Dict{Int, Float64}()
    Vminpu_D = Dict{Int, Float64}()
    Vmaxpu_D = Dict{Int, Float64}()
    
    p_D = Dict{Tuple{Int,Int}, Float64}()
    p_D_pu = Dict{Tuple{Int,Int}, Float64}()
    
    # Generate default PV shape
    LoadShapePV = HF.generateLoadShape(T, filenameLoadShape="LoadShapePVDefault.dss")
    
    pv_id = OpenDSSDirect.PVsystems.First()
    while pv_id > 0
        bus_names = OpenDSSDirect.CktElement.BusNames()
        bus = parse(Int, split(bus_names[1], ".")[1])
        push!(Dset, bus)
        
        # Get rated powers
        p_D_R[bus] = OpenDSSDirect.PVsystems.Pmpp()
        S_D_R[bus] = OpenDSSDirect.PVsystems.kVARated()
        
        p_D_R_pu[bus] = p_D_R[bus] / kVA_B_dict[bus]
        S_D_R_pu[bus] = S_D_R[bus] / kVA_B_dict[bus]
        
        irrad[bus] = OpenDSSDirect.PVsystems.Irradiance()
        
        # Voltage limits (using defaults)
        Vminpu_D[bus] = 0.90
        Vmaxpu_D[bus] = 1.10
        
        # Create time-varying profiles
        for t in 1:T
            p_D[(bus, t)] = p_D_R[bus] * LoadShapePV[t]
            p_D_pu[(bus, t)] = p_D_R_pu[bus] * LoadShapePV[t]
        end
        
        pv_id = OpenDSSDirect.PVsystems.Next()
    end
    
    Dset = sort(collect(Dset))
    n_D = length(Dset)
    DER_percent = 100  # Placeholder
    
    return Dict(
        :Dset => Dset,
        :n_D => n_D,
        :DER_percent => DER_percent,
        :p_D_R => p_D_R,
        :p_D_R_pu => p_D_R_pu,
        :S_D_R => S_D_R,
        :S_D_R_pu => S_D_R_pu,
        :irrad => irrad,
        :p_D => p_D,
        :p_D_pu => p_D_pu,
        :Vminpu_D => Vminpu_D,
        :Vmaxpu_D => Vmaxpu_D,
        :LoadShapePV => LoadShapePV
    )
end
#endregion

#region extract_batteries_from_dss
"""
Extract battery/storage data from OpenDSS circuit.
"""
function extract_batteries_from_dss(baseValuesDict)
    @unpack kVA_B_dict = baseValuesDict
    
    Bset = Set{Int}()
    B0 = Dict{Int, Float64}()
    B0_pu = Dict{Int, Float64}()
    B_R = Dict{Int, Float64}()
    B_R_pu = Dict{Int, Float64}()
    P_B_R = Dict{Int, Float64}()
    P_B_R_pu = Dict{Int, Float64}()
    S_B_R = Dict{Int, Float64}()
    S_B_R_pu = Dict{Int, Float64}()
    eta_C = Dict{Int, Float64}()
    eta_D = Dict{Int, Float64}()
    soc_min = Dict{Int, Float64}()
    soc_max = Dict{Int, Float64}()
    soc_0 = Dict{Int, Float64}()
    Vminpu_B = Dict{Int, Float64}()
    Vmaxpu_B = Dict{Int, Float64}()

    # Get all battery names
    all_names = OpenDSSDirect.Storages.AllNames()
    for batt_name in all_names
        # Select the storage element
        OpenDSSDirect.Storages.Name(batt_name)
        bus_names = OpenDSSDirect.CktElement.BusNames()
        bus = parse(Int, split(bus_names[1], ".")[1])
        push!(Bset, bus)

        # Retrieve parameters using Text.Command
        # kWh rated
        kWhRated = parse(Float64, OpenDSSDirect.Text.Command("? Storage.$batt_name.kWhrated"))
        B_R[bus] = kWhRated
        B_R_pu[bus] = kWhRated / kVA_B_dict[bus]

        # kW rated
        kWRated = parse(Float64, OpenDSSDirect.Text.Command("? Storage.$batt_name.kWrated"))
        P_B_R[bus] = kWRated
        P_B_R_pu[bus] = kWRated / kVA_B_dict[bus]

        # kVA (if available, else approximate)
        kVA_str = OpenDSSDirect.Text.Command("? Storage.$batt_name.kva")
        S_B_R[bus] = isnothing(tryparse(Float64, kVA_str)) ? kWRated * 1.1 : parse(Float64, kVA_str)
        S_B_R_pu[bus] = S_B_R[bus] / kVA_B_dict[bus]

        # Efficiencies
        eff_charge_str = OpenDSSDirect.Text.Command("? Storage.$batt_name.%effcharge")
        eta_C[bus] = isnothing(tryparse(Float64, eff_charge_str)) ? 0.95 : parse(Float64, eff_charge_str) / 100
        eff_discharge_str = OpenDSSDirect.Text.Command("? Storage.$batt_name.%effdischarge")
        eta_D[bus] = isnothing(tryparse(Float64, eff_discharge_str)) ? 0.95 : parse(Float64, eff_discharge_str) / 100

        # SOC parameters
        soc_0_str = OpenDSSDirect.Text.Command("? Storage.$batt_name.%stored")
        soc_0[bus] = isnothing(tryparse(Float64, soc_0_str)) ? 0.625 : parse(Float64, soc_0_str) / 100
        soc_min_str = OpenDSSDirect.Text.Command("? Storage.$batt_name.%reserve")
        soc_min[bus] = isnothing(tryparse(Float64, soc_min_str)) ? 0.3 : parse(Float64, soc_min_str) / 100
        soc_max[bus] = 0.95  # Default max

        B0[bus] = soc_0[bus] * B_R[bus]
        B0_pu[bus] = B0[bus] / kVA_B_dict[bus]

        # Voltage limits
        vminpu_str = OpenDSSDirect.Text.Command("? Storage.$batt_name.vminpu")
        Vminpu_B[bus] = isnothing(tryparse(Float64, vminpu_str)) ? 0.90 : parse(Float64, vminpu_str)
        vmaxpu_str = OpenDSSDirect.Text.Command("? Storage.$batt_name.vmaxpu")
        Vmaxpu_B[bus] = isnothing(tryparse(Float64, vmaxpu_str)) ? 1.10 : parse(Float64, vmaxpu_str)
    end

    Bset = sort(collect(Bset))
    n_B = length(Bset)
    Batt_percent = 100  # Placeholder

    # Reference SOC (default to initial)
    Bref = B0
    Bref_pu = B0_pu
    Bref_percent = Dict(j => Bref_pu[j] / B_R_pu[j] for j in Bset)

    return Dict(
        :Bset => Bset,
        :n_B => n_B,
        :Batt_percent => Batt_percent,
        :B0 => B0,
        :B0_pu => B0_pu,
        :Bref => Bref,
        :Bref_pu => Bref_pu,
        :Bref_percent => Bref_percent,
        :B_R => B_R,
        :B_R_pu => B_R_pu,
        :P_B_R => P_B_R,
        :P_B_R_pu => P_B_R_pu,
        :S_B_R => S_B_R,
        :S_B_R_pu => S_B_R_pu,
        :eta_C => eta_C,
        :eta_D => eta_D,
        :soc_min => soc_min,
        :soc_max => soc_max,
        :soc_0 => soc_0,
        :Vminpu_B => Vminpu_B,
        :Vmaxpu_B => Vmaxpu_B
    )
end
#endregion

#region evaluate_voltage_limits_from_components
"""
Consolidate voltage limits from all components at each bus.
"""
function evaluate_voltage_limits_from_components(data)
    @unpack Nset, NLset, Dset, Bset = data
    @unpack Vminpu_L, Vmaxpu_L = data
    @unpack Vminpu_D, Vmaxpu_D = data
    @unpack Vminpu_B, Vmaxpu_B = data
    
    Vminpu = Dict{Int, Float64}()
    Vmaxpu = Dict{Int, Float64}()
    
    for bus in Nset
        vmin_list = Float64[]
        vmax_list = Float64[]
        
        if bus in NLset
            push!(vmin_list, Vminpu_L[bus])
            push!(vmax_list, Vmaxpu_L[bus])
        end
        
        if bus in Dset
            push!(vmin_list, Vminpu_D[bus])
            push!(vmax_list, Vmaxpu_D[bus])
        end
        
        if bus in Bset
            push!(vmin_list, Vminpu_B[bus])
            push!(vmax_list, Vmaxpu_B[bus])
        end
        
        if isempty(vmin_list)
            Vminpu[bus] = 0.90
            Vmaxpu[bus] = 1.10
        else
            Vminpu[bus] = maximum(vmin_list)
            Vmaxpu[bus] = minimum(vmax_list)
        end
    end
    
    return Dict(
        :Vminpu => Vminpu,
        :Vmaxpu => Vmaxpu
    )
end
#endregion

#region extract_simulation_parameters_from_dss
"""
Extract simulation parameters like voltage source settings and time step.
"""
function extract_simulation_parameters_from_dss()
    # Get substation bus and voltage
    substationBus = 1  # Default assumption
    V_Subs_pu = 1.0
    delta_t = 1.0  # hours
    
    # Try to get vsource information
    vsource_id = OpenDSSDirect.Vsources.First()
    if vsource_id > 0
        V_Subs_pu = OpenDSSDirect.Vsources.PU()
    end
    
    return Dict(
        :substationBus => substationBus,
        :V_Subs_pu => V_Subs_pu,
        :delta_t => delta_t
    )
end
#endregion

#region create_metadata_dict
"""
Create metadata dictionary with solver configuration and objective function parameters.
"""
function create_metadata_dict(
    systemName, T, numAreas, objfun0, objfun2, temporal_decmp, algo_temporal_decmp,
    PSubsMax_kW, inputForecastDescription, solver, tSOC_hard, 
    relax_terminal_soc_constraint, linearizedModel, gedDict_ud,
    alpha_fpi, gamma_fpi, warmStart_mu, threshold_conv_iters)
    
    # Objective function strings
    if objfun0 == "subsPowerCostMin"
        objfunString = "Cost of Substation Power"
        objfunSense = "Min"
        objfunPrefix = "subsPowerCost_min"
        objfunUnit = "\$"
    elseif objfun0 == "lineLossMin"
        objfunString = "Line Losses"
        objfunSense = "Min"
        objfunPrefix = "lineLoss_min"
        objfunUnit = "kW"
    elseif objfun0 == "subsPowerMin"
        objfunString = "Substation Power"
        objfunSense = "Min"
        objfunPrefix = "subsPower_min"
        objfunUnit = "kW"
    else
        objfunString = "unknown objective"
        objfunSense = "Min"
        objfunPrefix = "unknown_obj"
        objfunUnit = ""
    end
    
    objfunAppendix = (objfun2 == "scd") ? "with_scd" : ""
    objfunConciseDescription = objfunPrefix * "_" * objfunAppendix
    
    # Temporal decomposition strings
    if temporal_decmp
        if algo_temporal_decmp == "DDP"
            temporalDecmpString = "Temporally Decomposed via DDP"
            temporalDecmpAppendix = "tmprl_dcmpsd"
        elseif algo_temporal_decmp == "tENApp"
            temporalDecmpString = "Temporally Decomposed via tENApp"
            temporalDecmpAppendix = "tmprl_dcmpsd_tENApp"
        else
            temporalDecmpString = "Temporally Decomposed"
            temporalDecmpAppendix = "tmprl_dcmpsd"
        end
    else
        temporalDecmpString = "Temporally Brute-forced"
        temporalDecmpAppendix = "tmprl_bruteforced"
    end
    
    # Linearized model strings
    if linearizedModel
        linearizedModelString = "LinDistFlow 1ph"
        linearizedModelAppendix = "ldf_1ph"
    else
        linearizedModelString = "BranchFlowModel 1ph"
        linearizedModelAppendix = "bfm_NL_1ph"
    end
    
    # Spatial decomposition strings
    if numAreas > 1
        spatialDecString = "Spatially Decomposed into $(numAreas) areas"
        spatialDecAppendix = "spat_dcmpsd_$(numAreas)_areas"
    else
        spatialDecString = "Spatially Centralized"
        spatialDecAppendix = "spat_centr_system"
    end
    
    machine_ID = gethostname()
    simNatureString = temporalDecmpString * ", " * spatialDecString
    simNatureAppendix = temporalDecmpAppendix * "_" * spatialDecAppendix
    
    macroItrsCompleted = 0
    solution_time = -1
    
    return Dict(
        :alpha_fpi => alpha_fpi,
        :inputForecastDescription => inputForecastDescription,
        :gedDict_ud => gedDict_ud,
        :machine_ID => machine_ID,
        :macroItrsCompleted => macroItrsCompleted,
        :numAreas => numAreas,
        :solution_time => solution_time,
        :gamma_fpi => gamma_fpi,
        :linearizedModel => linearizedModel,
        :linearizedModelAppendix => linearizedModelAppendix,
        :linearizedModelString => linearizedModelString,
        :objfun0 => objfun0,
        :objfun2 => objfun2,
        :objfunString => objfunString,
        :objfunSense => objfunSense,
        :objfunPrefix => objfunPrefix,
        :objfunAppendix => objfunAppendix,
        :objfunConciseDescription => objfunConciseDescription,
        :objfunUnit => objfunUnit,
        :PSubsMax_kW => PSubsMax_kW,
        :relax_terminal_soc_constraint => relax_terminal_soc_constraint,
        :simNatureAppendix => simNatureAppendix,
        :simNatureString => simNatureString,
        :solver => solver,
        :spatialDecAppendix => spatialDecAppendix,
        :spatialDecString => spatialDecString,
        :tSOC_hard => tSOC_hard,
        :temporal_decmp => temporal_decmp,
        :algo_temporal_decmp => algo_temporal_decmp,
        :temporalDecmpString => temporalDecmpString,
        :temporalDecmpAppendix => temporalDecmpAppendix,
        :warmStart_mu => warmStart_mu,
        :threshold_conv_iters => threshold_conv_iters
    )
end
#endregion

#region post_process_dss_data
"""
Post-process extracted data to compute derived quantities and add final metadata.
"""
function post_process_dss_data(data)
    @unpack rdict, xdict, Z_B_dict = data
    
    # Compute per-unit impedances
    rdict_pu = Dict{Tuple{Int,Int}, Float64}()
    xdict_pu = Dict{Tuple{Int,Int}, Float64}()
    
    for (branch, r_val) in rdict
        rdict_pu[branch] = r_val / Z_B_dict[branch]
        xdict_pu[branch] = xdict[branch] / Z_B_dict[branch]
    end
    
    @pack! data = rdict_pu, xdict_pu
    
    # Compute DER and battery percentages (if N_L exists)
    if haskey(data, :N_L) && haskey(data, :n_D)
        @unpack N_L, n_D = data
        DER_percent = Int(ceil(n_D / N_L * 100))
        @pack! data = DER_percent
    end
    
    if haskey(data, :N_L) && haskey(data, :n_B)
        @unpack N_L, n_B = data
        Batt_percent = Int(ceil(n_B / N_L * 100))
        @pack! data = Batt_percent
    end
    
    # Add GED string and appendix
    if haskey(data, :DER_percent) && haskey(data, :Batt_percent)
        @unpack DER_percent, Batt_percent = data
        gedString = "$(DER_percent)% PVs and $(Batt_percent)% Batteries"
        gedAppendix = "pv_$(DER_percent)_batt_$(Batt_percent)"
        @pack! data = gedAppendix, gedString
    end
    
    return data
end
#endregion

end # module parseFromDSS
