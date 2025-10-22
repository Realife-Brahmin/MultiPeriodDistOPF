"""
Standalone OpenDSS Parser for tADMM
Minimal dependencies, extracted from MultiPeriodDistOPF Parser module
"""

using OpenDSSDirect
using Parameters: @unpack

"""
    load_system_in_dss(systemName::String; rawDataFolderPath=nothing, runPowerFlow=false)

Load a system into OpenDSS engine.
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


"""
    parse_system_from_dss(systemName::String, T::Int; kwargs...)

Parse system data from OpenDSS for MPOPF.
Simplified version with essential data only.
"""
function parse_system_from_dss(systemName::String, T::Int;
    kVA_B=1000.0,
    kV_B=2.4018,
    LoadShapeLoad=nothing,
    LoadShapePV=nothing,
    LoadShapeCost=nothing,
    C_B=1e-6 * 0.08,
    delta_t_h=1.0,
    rawDataFolderPath=nothing,
    kwargs...)
    
    # Load system into OpenDSS
    load_system_in_dss(systemName; rawDataFolderPath=rawDataFolderPath, runPowerFlow=true)
    
    # Initialize data dictionary
    data = Dict{Symbol, Any}()
    
    # Store simulation parameters
    data[:systemName] = systemName
    data[:T] = T
    data[:Tset] = 1:T
    data[:delta_t_h] = delta_t_h
    data[:C_B] = C_B
    data[:kVA_B] = kVA_B
    data[:kV_B] = kV_B
    
    # Default load shapes if not provided
    if isnothing(LoadShapeLoad)
        LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2
    end
    if isnothing(LoadShapePV)
        LoadShapePV = zeros(T)  # No PV by default
    end
    if isnothing(LoadShapeCost)
        LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2
    end
    
    data[:LoadShapeLoad] = LoadShapeLoad
    data[:LoadShapePV] = LoadShapePV
    data[:LoadShapeCost] = LoadShapeCost
    
    # Parse network topology
    parse_network_topology!(data)
    
    # Parse line impedances
    parse_line_impedances!(data, kVA_B, kV_B)
    
    # Parse loads
    parse_loads!(data, T, LoadShapeLoad, kVA_B)
    
    # Parse PV/DERs
    parse_pv_generators!(data, T, LoadShapePV, kVA_B)
    
    # Parse batteries
    parse_batteries!(data, T, kVA_B)
    
    # Parse voltage limits
    parse_voltage_limits!(data)
    
    return data
end


"""
    parse_network_topology!(data::Dict)

Parse network topology: buses, branches, parent-child relationships.
"""
function parse_network_topology!(data::Dict)
    # Get circuit name and substation bus
    circuit_name = OpenDSSDirect.Circuit.Name()
    data[:substationBus] = 1  # Assume bus 1 is substation
    
    # Get all bus names and extract bus numbers
    bus_names = OpenDSSDirect.Circuit.AllBusNames()
    bus_numbers = [parse(Int, split(name, ".")[1]) for name in bus_names]
    N = maximum(bus_numbers)
    
    # Node sets
    data[:N] = N
    data[:Nset] = sort(unique(bus_numbers))
    
    # Get all line/branch names
    line_names = OpenDSSDirect.Lines.AllNames()
    data[:m] = length(line_names)
    
    # Parse branches and build topology
    Lset = Tuple{Int,Int}[]
    parent = Dict{Int, Union{Nothing, Int}}()
    children = Dict{Int, Vector{Int}}()
    
    # Initialize children dict
    for i in 1:N
        children[i] = Int[]
    end
    
    # Substation has no parent
    parent[data[:substationBus]] = nothing
    
    # Parse each line
    for line_name in line_names
        OpenDSSDirect.Circuit.SetActiveElement("Line.$line_name")
        
        # Get bus names for this line
        bus1_full = OpenDSSDirect.CktElement.BusNames()[1]
        bus2_full = OpenDSSDirect.CktElement.BusNames()[2]
        
        # Extract bus numbers (e.g., "7.1" -> 7)
        i = parse(Int, split(bus1_full, ".")[1])
        j = parse(Int, split(bus2_full, ".")[1])
        
        # Add to branch set
        push!(Lset, (i, j))
        
        # Build parent-child relationships
        # Assume branches point away from substation
        if i == data[:substationBus] || (haskey(parent, i) && !haskey(parent, j))
            parent[j] = i
            push!(children[i], j)
        elseif j == data[:substationBus] || (haskey(parent, j) && !haskey(parent, i))
            parent[i] = j
            push!(children[j], i)
        else
            # Use distance from substation or assume i -> j
            parent[j] = i
            push!(children[i], j)
        end
    end
    
    data[:Lset] = Lset
    data[:parent] = parent
    data[:children] = children
    
    # Classify branches
    substation_bus = data[:substationBus]
    L1set = Tuple{Int,Int}[]  # Branches directly from substation
    Lm1set = Tuple{Int,Int}[]  # All other branches
    
    for (i, j) in Lset
        if i == substation_bus
            push!(L1set, (i, j))
        else
            push!(Lm1set, (i, j))
        end
    end
    
    data[:L1set] = L1set
    data[:Lm1set] = Lm1set
    
    # Non-substation nodes
    Nm1set = filter(n -> n != substation_bus, data[:Nset])
    data[:Nm1set] = Nm1set
end


"""
    parse_line_impedances!(data::Dict, kVA_B::Float64, kV_B::Float64)

Parse line resistance and reactance.
"""
function parse_line_impedances!(data::Dict, kVA_B::Float64, kV_B::Float64)
    Z_B = kV_B^2 / (kVA_B / 1000)  # Base impedance in Ohms
    
    rdict = Dict{Tuple{Int,Int}, Float64}()
    xdict = Dict{Tuple{Int,Int}, Float64}()
    rdict_pu = Dict{Tuple{Int,Int}, Float64}()
    xdict_pu = Dict{Tuple{Int,Int}, Float64}()
    
    line_names = OpenDSSDirect.Lines.AllNames()
    
    for line_name in line_names
        OpenDSSDirect.Lines.Name(line_name)
        
        # Get bus numbers from bus names (e.g., "7.1" -> 7)
        buses = OpenDSSDirect.CktElement.BusNames()
        i = parse(Int, split(buses[1], ".")[1])
        j = parse(Int, split(buses[2], ".")[1])
        
        # Get R and X (Ohms)
        R_ohm = OpenDSSDirect.Lines.RMatrix()[1]  # First element of R matrix
        X_ohm = OpenDSSDirect.Lines.XMatrix()[1]  # First element of X matrix
        
        rdict[(i,j)] = R_ohm
        xdict[(i,j)] = X_ohm
        rdict_pu[(i,j)] = R_ohm / Z_B
        xdict_pu[(i,j)] = X_ohm / Z_B
    end
    
    data[:Z_B] = Z_B
    data[:rdict] = rdict
    data[:xdict] = xdict
    data[:rdict_pu] = rdict_pu
    data[:xdict_pu] = xdict_pu
end


"""
    parse_loads!(data::Dict, T::Int, LoadShapeLoad::Vector, kVA_B::Float64)

Parse load data.
"""
function parse_loads!(data::Dict, T::Int, LoadShapeLoad::Vector, kVA_B::Float64)
    load_names = OpenDSSDirect.Loads.AllNames()
    
    NLset = Int[]  # Set of buses with loads
    p_L_R_pu = Dict{Int, Float64}()  # Base load per bus
    q_L_R_pu = Dict{Int, Float64}()
    
    for load_name in load_names
        OpenDSSDirect.Loads.Name(load_name)
        
        # Get bus number (e.g., "7.1" -> 7)
        bus_full = OpenDSSDirect.CktElement.BusNames()[1]
        bus_num = parse(Int, split(bus_full, ".")[1])
        
        # Get load power (kW, kVAR)
        kW = OpenDSSDirect.Loads.kW()
        kVAR = OpenDSSDirect.Loads.kvar()
        
        # Convert to per-unit
        p_base_pu = kW / kVA_B
        q_base_pu = kVAR / kVA_B
        
        # Store base values (if multiple loads on same bus, add them)
        if haskey(p_L_R_pu, bus_num)
            p_L_R_pu[bus_num] += p_base_pu
            q_L_R_pu[bus_num] += q_base_pu
        else
            p_L_R_pu[bus_num] = p_base_pu
            q_L_R_pu[bus_num] = q_base_pu
            push!(NLset, bus_num)
        end
    end
    
    # Create 2D arrays for time-varying load (indexed by [bus, time])
    # Get max bus number from already-parsed data
    Nset = data[:Nset]
    max_bus = maximum(Nset)
    p_L_pu = zeros(max_bus, T)
    q_L_pu = zeros(max_bus, T)
    
    # Apply time-varying load shape
    for bus in NLset
        for t in 1:T
            p_L_pu[bus, t] = p_L_R_pu[bus] * LoadShapeLoad[t]
            q_L_pu[bus, t] = q_L_R_pu[bus] * LoadShapeLoad[t]
        end
    end
    
    data[:NLset] = sort(NLset)
    data[:N_L] = length(NLset)
    data[:p_L_pu] = p_L_pu
    data[:q_L_pu] = q_L_pu
    data[:p_L_R_pu] = p_L_R_pu
    data[:q_L_R_pu] = q_L_R_pu
end


"""
    parse_pv_generators!(data::Dict, T::Int, LoadShapePV::Vector, kVA_B::Float64)

Parse PV/generator data.
"""
function parse_pv_generators!(data::Dict, T::Int, LoadShapePV::Vector, kVA_B::Float64)
    gen_names = OpenDSSDirect.PVsystems.AllNames()
    
    Dset = Int[]  # Set of buses with PV
    p_D_R_pu = Dict{Int, Float64}()
    S_D_R = Dict{Int, Float64}()  # Apparent power rating (pu)
    
    # Get max bus number for array sizing
    Nset = data[:Nset]
    max_bus = maximum(Nset)
    
    # Initialize 2D array for time-varying PV (indexed by [bus, time])
    p_D_pu = zeros(max_bus, T)
    
    # Only parse if generators exist (check for "NONE" which indicates no generators)
    if !isempty(gen_names) && gen_names[1] != "NONE"
        for gen_name in gen_names
            OpenDSSDirect.PVsystems.Name(gen_name)

            # Get bus number (e.g., "7.1" -> 7)
            bus_full = OpenDSSDirect.CktElement.BusNames()[1]
            bus_num = parse(Int, split(bus_full, ".")[1])
            
            # Get generator rating
            kW_rated = OpenDSSDirect.PVsystems.kW()
            kVA_rated = OpenDSSDirect.PVsystems.kVARated()

            # Convert to per-unit
            p_rated_pu = kW_rated / kVA_B
            S_rated_pu = kVA_rated / kVA_B
            
            p_D_R_pu[bus_num] = p_rated_pu
            S_D_R[bus_num] = S_rated_pu
            push!(Dset, bus_num)
        end
        
        # Apply time-varying PV shape
        for bus in Dset
            for t in 1:T
                p_D_pu[bus, t] = p_D_R_pu[bus] * LoadShapePV[t]
            end
        end
    end
    
    data[:Dset] = sort(unique(Dset))
    data[:n_D] = length(data[:Dset])
    data[:p_D_pu] = p_D_pu
    data[:p_D_R_pu] = p_D_R_pu
    data[:S_D_R] = S_D_R
end


"""
    parse_batteries!(data::Dict, T::Int, kVA_B::Float64)

Parse battery/storage data.
"""
function parse_batteries!(data::Dict, T::Int, kVA_B::Float64)
    storage_names = OpenDSSDirect.Storages.AllNames()
    
    Bset = Int[]
    B0_pu = Dict{Int, Float64}()  # Initial SOC (pu, in kWh/kVA_B)
    B_R = Dict{Int, Float64}()  # Energy capacity (kWh)
    B_R_pu = Dict{Int, Float64}()  # Energy capacity (pu)
    P_B_R = Dict{Int, Float64}()  # Power rating (kW)
    P_B_R_pu = Dict{Int, Float64}()  # Power rating (pu)
    S_B_R_pu = Dict{Int, Float64}()  # Apparent power rating (pu)
    soc_min = Dict{Int, Float64}()
    soc_max = Dict{Int, Float64}()
    
    E_BASE = kVA_B * 1.0  # Energy base: kVA_B * 1 hour = kWh
    
    # Only parse if storage exists (check for "NONE" which indicates no storage)
    if !isempty(storage_names) && storage_names[1] != "NONE"
        for storage_name in storage_names
            OpenDSSDirect.Storages.Name(storage_name)
            
            # Get bus number (e.g., "7.1" -> 7)
            bus_full = OpenDSSDirect.CktElement.BusNames()[1]
            bus_num = parse(Int, split(bus_full, ".")[1])
            
            # Get storage parameters using Text.Command
            kWh_rated = parse(Float64, OpenDSSDirect.Text.Command("? Storage.$storage_name.kWhrated"))
            kW_rated = parse(Float64, OpenDSSDirect.Text.Command("? Storage.$storage_name.kWrated"))
            kVA_rated = parse(Float64, OpenDSSDirect.Text.Command("? Storage.$storage_name.kva"))
            pct_stored = parse(Float64, OpenDSSDirect.Text.Command("? Storage.$storage_name.%stored"))
            pct_reserve = parse(Float64, OpenDSSDirect.Text.Command("? Storage.$storage_name.%reserve"))
            
            # Convert to per-unit
            B_R[bus_num] = kWh_rated
            B_R_pu[bus_num] = kWh_rated / E_BASE
            P_B_R[bus_num] = kW_rated
            P_B_R_pu[bus_num] = kW_rated / kVA_B
            S_B_R_pu[bus_num] = kVA_rated / kVA_B
            
            # Initial SOC
            B0_kWh = kWh_rated * (pct_stored / 100.0)
            B0_pu[bus_num] = B0_kWh / E_BASE
            
            # SOC limits
            soc_min[bus_num] = pct_reserve / 100.0
            soc_max[bus_num] = 0.95  # Default max SOC
            
            push!(Bset, bus_num)
        end
    end
    
    data[:Bset] = sort(unique(Bset))
    data[:n_B] = length(data[:Bset])
    data[:B0_pu] = B0_pu
    data[:B_R] = B_R
    data[:B_R_pu] = B_R_pu
    data[:P_B_R] = P_B_R
    data[:P_B_R_pu] = P_B_R_pu
    data[:S_B_R_pu] = S_B_R_pu
    data[:soc_min] = soc_min
    data[:soc_max] = soc_max
end


"""
    parse_voltage_limits!(data::Dict)

Parse voltage limits for all buses.
"""
function parse_voltage_limits!(data::Dict)
    N = data[:N]
    Vminpu = Dict{Int, Float64}()
    Vmaxpu = Dict{Int, Float64}()
    
    # Default voltage limits
    for i in 1:N
        Vminpu[i] = 0.95
        Vmaxpu[i] = 1.05
    end
    
    data[:Vminpu] = Vminpu
    data[:Vmaxpu] = Vmaxpu
end
