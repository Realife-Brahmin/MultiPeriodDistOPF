"""
OpenDSS Parser for Multi-POI Systems
Adapted from tADMM parser to handle multiple substation buses
"""

using OpenDSSDirect

"""
    load_system_in_dss(systemName::String; rawDataFolderPath=nothing, runPowerFlow=false)

Load a multi-POI system into OpenDSS engine.
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
    
    println("✓ Loaded system: $systemName from $dss_file")
    
    return dss_file
end


"""
    parse_system_from_dss_multi_poi(systemName::String, T::Int; kwargs...)

Parse multi-POI system data from OpenDSS for MPOPF.
This version handles multiple substation buses (POIs).

Key differences from single-POI parser:
- Identifies all substation buses (ending in 's': 1s, 2s, 3s, etc.)
- Identifies POI connection lines (Line.poi1, Line.poi2, etc.)
- Stores multi-substation information for later use
"""
function parse_system_from_dss_multi_poi(systemName::String, T::Int;
    kVA_B=1000.0,
    kV_B=11.547,  # 20kV line-to-line -> 11.547kV line-to-neutral
    LoadShapeLoad=nothing,
    LoadShapePV=nothing,
    LoadShapeCost_dict=nothing,  # Dict mapping substation number to cost vector
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
    if isnothing(LoadShapeCost_dict)
        # Default: single cost profile for all substations
        default_cost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2
        LoadShapeCost_dict = Dict(1 => default_cost)
    end
    
    data[:LoadShapeLoad] = LoadShapeLoad
    data[:LoadShapePV] = LoadShapePV
    data[:LoadShapeCost_dict] = LoadShapeCost_dict
    
    println("\n" * "="^80)
    println("PARSING MULTI-POI SYSTEM FROM OPENDSS")
    println("="^80)
    
    # Parse multi-POI specific elements (substations and POI lines)
    parse_multi_poi_elements!(data)
    
    # Parse network topology (branches, buses, parent-child relationships)
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
    
    println("="^80)
    println("✓ PARSING COMPLETE")
    println("="^80)
    
    return data
end


"""
    parse_multi_poi_elements!(data::Dict)

Parse multi-POI specific elements: substation buses and POI connection lines.
Identifies buses ending in 's' as substation buses (1s, 2s, 3s, etc.).
"""
function parse_multi_poi_elements!(data::Dict)
    println("\n--- Parsing Multi-POI Elements ---")
    
    # Get all bus names
    all_bus_names = OpenDSSDirect.Circuit.AllBusNames()
    
    # Identify substation buses (ending in 's')
    substation_buses = String[]
    substation_bus_numbers = Int[]
    
    for bus_name in all_bus_names
        # Extract base name (without phase suffix like ".1")
        base_name = split(bus_name, ".")[1]
        
        # Check if it ends with 's' (substation bus)
        if endswith(base_name, "s")
            push!(substation_buses, base_name)
            # Extract substation number (e.g., "1s" -> 1, "2s" -> 2)
            sub_num_str = base_name[1:end-1]  # Remove 's'
            sub_num = parse(Int, sub_num_str)
            push!(substation_bus_numbers, sub_num)
        end
    end
    
    # Sort by substation number
    sort_indices = sortperm(substation_bus_numbers)
    substation_buses = substation_buses[sort_indices]
    substation_bus_numbers = substation_bus_numbers[sort_indices]
    
    num_subs = length(substation_buses)
    
    println("  Number of substations: $num_subs")
    for (i, (bus_name, bus_num)) in enumerate(zip(substation_buses, substation_bus_numbers))
        println("    Substation $i: Bus $bus_name (number $bus_num)")
    end
    
    # Identify POI connection lines
    all_line_names = OpenDSSDirect.Lines.AllNames()
    poi_lines = String[]
    poi_connections = Tuple{String, Int}[]  # (substation_bus, distribution_bus)
    
    for line_name in all_line_names
        if startswith(lowercase(line_name), "poi")
            push!(poi_lines, line_name)
            
            # Get bus connections for this POI line
            OpenDSSDirect.Lines.Name(line_name)
            buses = OpenDSSDirect.CktElement.BusNames()
            
            # buses[1] should be substation (e.g., "1s.1"), buses[2] should be distribution bus
            bus1_base = split(buses[1], ".")[1]
            bus2_base = split(buses[2], ".")[1]
            
            # Determine which is substation and which is distribution
            if endswith(bus1_base, "s")
                sub_bus = bus1_base
                dist_bus = parse(Int, bus2_base)
            else
                sub_bus = bus2_base
                dist_bus = parse(Int, bus1_base)
            end
            
            push!(poi_connections, (sub_bus, dist_bus))
        end
    end
    
    println("\n  POI Connection Lines: $num_subs")
    for (line_name, (sub_bus, dist_bus)) in zip(poi_lines, poi_connections)
        println("    $line_name: $sub_bus → $dist_bus")
    end
    
    # Store in data dictionary
    data[:num_substations] = num_subs
    data[:substation_buses] = substation_buses  # ["1s", "2s", "3s", ...]
    data[:substation_bus_numbers] = substation_bus_numbers  # [1, 2, 3, ...]
    data[:poi_lines] = poi_lines  # ["poi1", "poi2", ...]
    data[:poi_connections] = poi_connections  # [("1s", 1), ("2s", 33), ...]
    
    # For backward compatibility, store first substation as "substationBus"
    # (This will need to be updated when we handle multi-substation optimization)
    data[:substationBus] = 1  # Placeholder - will be updated for multi-POI
end


"""
    parse_network_topology!(data::Dict)

Parse network topology: buses, branches, parent-child relationships.
For multi-POI systems, this includes substation buses and POI connection lines.

Network structure for MPOPF:
- Sset: Substation bus names ["1s", "2s", "3s", "4s", "5s"]
- Nset: ALL buses including substations ["1s", "2s", ..., "5s", 1, 2, ..., 128]
- Nm1set: Only distribution buses [1, 2, 3, ..., 128]
- Lset: ALL lines including POI lines [(1s, 1), (2s, 33), ..., (1, 2), (2, 3), ...]
- L1set: POI connection lines [(1s, 1), (2s, 33), (3s, 120), (4s, 116), (5s, 85)]
- Lm1set: Only distribution lines [(1, 2), (2, 3), ..., (127, 128)]
"""
function parse_network_topology!(data::Dict)
    println("\n--- Parsing Network Topology ---")
    
    # Get all bus names from OpenDSS
    bus_names = OpenDSSDirect.Circuit.AllBusNames()
    
    # Separate substation and distribution buses
    substation_buses = data[:substation_buses]  # Already parsed: ["1s", "2s", ...]
    dist_bus_numbers = Int[]
    
    for name in bus_names
        base_name = split(name, ".")[1]
        if !endswith(base_name, "s")  # Not a substation bus
            bus_num = parse(Int, base_name)
            push!(dist_bus_numbers, bus_num)
        end
    end
    
    dist_bus_numbers = sort(unique(dist_bus_numbers))
    N_dist = maximum(dist_bus_numbers)  # Max distribution bus number
    
    # Create Sset: substation bus names
    Sset = substation_buses  # ["1s", "2s", "3s", "4s", "5s"]
    
    # Create Nset: ALL buses (substations + distribution)
    # Store as Vector{Any} to hold both String and Int
    Nset = Vector{Any}(vcat(substation_buses, dist_bus_numbers))
    
    # Create Nm1set: only distribution buses (excluding substations)
    Nm1set = dist_bus_numbers  # [1, 2, 3, ..., 128]
    
    data[:Sset] = Sset
    data[:Nset] = Nset
    data[:Nm1set] = Nm1set
    data[:N] = length(Nset)  # Total number of buses (substations + distribution)
    data[:N_sub] = length(Sset)  # Number of substation buses
    data[:N_dist] = length(Nm1set)  # Number of distribution buses only
    
    println("  Substation buses (Sset): $(data[:N_sub])")
    println("  Distribution buses (Nm1set): $(data[:N_dist]) ($(minimum(Nm1set)) to $(maximum(Nm1set)))")
    println("  Total buses N = |Nset|: $(data[:N]) (substations + distribution)")
    
    # Get all line names
    line_names = OpenDSSDirect.Lines.AllNames()
    
    # Parse ALL lines (POI + distribution)
    Lset = []  # Will contain Tuple{Any, Any} to handle mixed String/Int bus names
    L1set = []  # POI lines only
    Lm1set = []  # Distribution lines only
    
    parent = Dict{Any, Union{Nothing, Any}}()
    children = Dict{Any, Vector{Any}}()
    
    # Initialize children dict for all buses
    for bus in Nset
        children[bus] = []
    end
    
    # Parse each line
    for line_name in line_names
        OpenDSSDirect.Lines.Name(line_name)
        
        # Get bus names for this line
        buses = OpenDSSDirect.CktElement.BusNames()
        bus1_base = split(buses[1], ".")[1]
        bus2_base = split(buses[2], ".")[1]
        
        # Determine bus types (substation or distribution)
        bus1 = endswith(bus1_base, "s") ? bus1_base : parse(Int, bus1_base)
        bus2 = endswith(bus2_base, "s") ? bus2_base : parse(Int, bus2_base)
        
        # Add to appropriate line set
        line_tuple = (bus1, bus2)
        push!(Lset, line_tuple)
        
        # Classify as POI or distribution line
        is_poi_line = startswith(lowercase(line_name), "poi")
        
        if is_poi_line
            push!(L1set, line_tuple)
            # POI lines: substation -> distribution bus
            # Substation has no parent, distribution bus's parent is substation
            if endswith(string(bus1), "s")
                parent[bus1] = nothing  # Substation has no parent
                parent[bus2] = bus1  # Distribution bus's parent is substation
                push!(children[bus1], bus2)
            else
                parent[bus2] = nothing  # Substation has no parent
                parent[bus1] = bus2  # Distribution bus's parent is substation
                push!(children[bus2], bus1)
            end
        else
            push!(Lm1set, line_tuple)
            # Distribution line: build parent-child relationships
            # Assume bus1 -> bus2 (pointing away from substations)
            if !haskey(parent, bus2) || isnothing(parent[bus2])
                parent[bus2] = bus1
                push!(children[bus1], bus2)
            elseif !haskey(parent, bus1) || isnothing(parent[bus1])
                parent[bus1] = bus2
                push!(children[bus2], bus1)
            else
                # Both have parents, use heuristic (lower bus number is closer to root)
                if isa(bus1, Int) && isa(bus2, Int)
                    if bus1 < bus2
                        parent[bus2] = bus1
                        push!(children[bus1], bus2)
                    else
                        parent[bus1] = bus2
                        push!(children[bus2], bus1)
                    end
                end
            end
        end
    end
    
    data[:Lset] = Lset
    data[:L1set] = L1set  # POI lines
    data[:Lm1set] = Lm1set  # Distribution lines
    data[:parent] = parent
    data[:children] = children
    data[:m] = length(Lset)  # Total number of lines
    data[:m_poi] = length(L1set)  # Number of POI lines
    data[:m_dist] = length(Lm1set)  # Number of distribution lines
    
    println("  POI lines (L1set): $(length(L1set))")
    println("  Distribution lines (Lm1set): $(length(Lm1set))")
    println("  Total lines (Lset): $(length(Lset))")
    
    # Print POI line details for verification
    println("\n  POI Line Details:")
    for (bus1, bus2) in L1set
        println("    $bus1 → $bus2")
    end
end


"""
    parse_line_impedances!(data::Dict, kVA_B::Float64, kV_B::Float64)

Parse line resistance and reactance for all lines (POI + distribution).
Creates unified impedance dictionaries that can handle both String (substation) and Int (distribution) bus names.
"""
function parse_line_impedances!(data::Dict, kVA_B::Float64, kV_B::Float64)
    println("\n--- Parsing Line Impedances ---")
    
    Z_B = kV_B^2 / (kVA_B / 1000)  # Base impedance in Ohms
    
    # Unified dictionaries for ALL lines (Dict{Any, Float64} to handle mixed bus types)
    rdict_pu = Dict{Any, Float64}()
    xdict_pu = Dict{Any, Float64}()
    
    # Separate POI dicts for backward compatibility (keyed by substation bus name)
    poi_rdict = Dict{String, Float64}()
    poi_xdict = Dict{String, Float64}()
    
    line_names = OpenDSSDirect.Lines.AllNames()
    
    for line_name in line_names
        OpenDSSDirect.Lines.Name(line_name)
        
        # Get bus names (e.g., "7.1" -> 7, or "1s.1" -> "1s")
        buses = OpenDSSDirect.CktElement.BusNames()
        bus1_base = split(buses[1], ".")[1]
        bus2_base = split(buses[2], ".")[1]
        
        # Determine bus types (String for substation, Int for distribution)
        bus1 = endswith(bus1_base, "s") ? bus1_base : parse(Int, bus1_base)
        bus2 = endswith(bus2_base, "s") ? bus2_base : parse(Int, bus2_base)
        
        # Get R and X (Ohms)
        R_ohm = OpenDSSDirect.Lines.RMatrix()[1]  # First element of R matrix
        X_ohm = OpenDSSDirect.Lines.XMatrix()[1]  # First element of X matrix
        
        # Convert to per-unit
        R_pu = R_ohm / Z_B
        X_pu = X_ohm / Z_B
        
        # Store in unified dictionaries
        line_tuple = (bus1, bus2)
        rdict_pu[line_tuple] = R_pu
        xdict_pu[line_tuple] = X_pu
        
        # If it's a POI line, also store in POI-specific dict
        if startswith(lowercase(line_name), "poi")
            sub_bus = endswith(string(bus1), "s") ? bus1 : bus2
            poi_rdict[sub_bus] = R_ohm  # Store in Ohms for POI dict
            poi_xdict[sub_bus] = X_ohm
        end
    end
    
    data[:Z_B] = Z_B
    data[:rdict_pu] = rdict_pu
    data[:xdict_pu] = xdict_pu
    data[:poi_rdict] = poi_rdict
    data[:poi_xdict] = poi_xdict
    
    println("  Base impedance Z_B: $(round(Z_B, digits=4)) Ω")
    println("  Total lines with impedances: $(length(rdict_pu))")
    println("  Distribution lines (Lm1set): $(data[:m_dist])")
    println("  POI lines (L1set): $(data[:m_poi])")
    
    # Print POI impedances for verification
    println("\n  POI Line Impedances:")
    for sub_bus in data[:substation_buses]
        R_ohm = get(poi_rdict, sub_bus, NaN)
        X_ohm = get(poi_xdict, sub_bus, NaN)
        R_pu = R_ohm / Z_B
        X_pu = X_ohm / Z_B
        println("    $sub_bus: R=$(round(R_ohm, digits=8)) Ω ($(round(R_pu, digits=10)) pu), X=$(round(X_ohm, digits=8)) Ω ($(round(X_pu, digits=10)) pu)")
    end
end


"""
    parse_loads!(data::Dict, T::Int, LoadShapeLoad::Vector, kVA_B::Float64)

Parse load data from OpenDSS.
"""
function parse_loads!(data::Dict, T::Int, LoadShapeLoad::Vector, kVA_B::Float64)
    println("\n--- Parsing Loads ---")
    
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
    # Use Nm1set (distribution buses only) to determine array size
    Nm1set = data[:Nm1set]
    max_bus = maximum(Nm1set)
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
    
    # Calculate total load
    total_P_kW = sum(values(p_L_R_pu)) * kVA_B
    total_Q_kVAr = sum(values(q_L_R_pu)) * kVA_B
    
    println("  Buses with loads: $(data[:N_L])")
    println("  Total base load: $(round(total_P_kW, digits=1)) kW, $(round(total_Q_kVAr, digits=1)) kVAr")
end


"""
    parse_pv_generators!(data::Dict, T::Int, LoadShapePV::Vector, kVA_B::Float64)

Parse PV/generator data from OpenDSS.
"""
function parse_pv_generators!(data::Dict, T::Int, LoadShapePV::Vector, kVA_B::Float64)
    println("\n--- Parsing PV Systems ---")
    
    gen_names = OpenDSSDirect.PVsystems.AllNames()
    
    Dset = Int[]  # Set of buses with PV
    p_D_R_pu = Dict{Int, Float64}()
    S_D_R = Dict{Int, Float64}()  # Apparent power rating (pu)
    
    # Get max bus number for array sizing (use Nm1set which contains only distribution buses)
    Nm1set = data[:Nm1set]
    max_bus = maximum(Nm1set)
    
    # Initialize 2D array for time-varying PV (indexed by [bus, time])
    p_D_pu = zeros(max_bus, T)
    
    # Only parse if generators exist
    if !isempty(gen_names) && gen_names[1] != "NONE"
        for gen_name in gen_names
            OpenDSSDirect.PVsystems.Name(gen_name)
            
            # Get bus number
            bus_full = OpenDSSDirect.CktElement.BusNames()[1]
            bus_num = parse(Int, split(bus_full, ".")[1])
            
            # Get rated power (kW)
            kW_rated = OpenDSSDirect.PVsystems.kW()
            kVA_rated = OpenDSSDirect.PVsystems.kVARated()
            
            # Convert to per-unit
            p_rated_pu = kW_rated / kVA_B
            S_rated_pu = kVA_rated / kVA_B
            
            # Store (accumulate if multiple PV on same bus)
            if haskey(p_D_R_pu, bus_num)
                p_D_R_pu[bus_num] += p_rated_pu
                S_D_R[bus_num] += S_rated_pu
            else
                p_D_R_pu[bus_num] = p_rated_pu
                S_D_R[bus_num] = S_rated_pu
            end
            
            if !(bus_num in Dset)
                push!(Dset, bus_num)
            end
        end
        
        # Apply time-varying PV shape
        for bus in Dset
            for t in 1:T
                p_D_pu[bus, t] = p_D_R_pu[bus] * LoadShapePV[t]
            end
        end
        
        total_PV_kW = sum(values(p_D_R_pu)) * kVA_B
        println("  Buses with PV: $(length(Dset))")
        println("  Total PV capacity: $(round(total_PV_kW, digits=1)) kW")
    else
        println("  No PV systems found")
    end
    
    data[:Dset] = sort(unique(Dset))
    data[:n_D] = length(data[:Dset])
    data[:p_D_pu] = p_D_pu
    data[:p_D_R_pu] = p_D_R_pu
    data[:S_D_R] = S_D_R
end


"""
    parse_batteries!(data::Dict, T::Int, kVA_B::Float64)

Parse battery/storage data from OpenDSS.
"""
function parse_batteries!(data::Dict, T::Int, kVA_B::Float64)
    println("\n--- Parsing Battery Storage ---")
    
    storage_names = OpenDSSDirect.Storages.AllNames()
    
    Bset = Int[]
    B0_pu = Dict{Int, Float64}()  # Initial SOC (pu)
    B_R = Dict{Int, Float64}()  # Energy capacity (kWh)
    B_R_pu = Dict{Int, Float64}()  # Energy capacity (pu)
    P_B_R = Dict{Int, Float64}()  # Power rating (kW)
    P_B_R_pu = Dict{Int, Float64}()  # Power rating (pu)
    S_B_R_pu = Dict{Int, Float64}()  # Apparent power rating (pu)
    soc_min = Dict{Int, Float64}()
    soc_max = Dict{Int, Float64}()
    
    E_BASE = kVA_B * 1.0  # Energy base: kVA_B * 1 hour = kWh
    
    # Only parse if storage exists
    if !isempty(storage_names) && storage_names[1] != "NONE"
        for storage_name in storage_names
            OpenDSSDirect.Storages.Name(storage_name)
            
            # Get bus number
            bus_full = OpenDSSDirect.CktElement.BusNames()[1]
            bus_num = parse(Int, split(bus_full, ".")[1])
            
            # Get battery parameters
            kWh_stored = OpenDSSDirect.Storages.kWhStored()  # Initial energy
            kWh_rated = OpenDSSDirect.Storages.kWhRated()    # Capacity
            kW_rated = OpenDSSDirect.Storages.kW()           # Power rating
            kVA_rated_storage = OpenDSSDirect.Storages.kVA() # Apparent power
            
            pct_stored = OpenDSSDirect.Storages.pctStored()  # Initial SOC %
            pct_reserve = OpenDSSDirect.Storages.pctReserve() # Min SOC %
            
            # Convert to per-unit
            B0_pu[bus_num] = kWh_stored / E_BASE
            B_R[bus_num] = kWh_rated
            B_R_pu[bus_num] = kWh_rated / E_BASE
            P_B_R[bus_num] = kW_rated
            P_B_R_pu[bus_num] = kW_rated / kVA_B
            S_B_R_pu[bus_num] = kVA_rated_storage / kVA_B
            
            # SOC limits
            soc_min[bus_num] = pct_reserve / 100.0
            soc_max[bus_num] = 1.0
            
            push!(Bset, bus_num)
        end
        
        total_battery_kWh = sum(values(B_R))
        total_battery_kW = sum(values(P_B_R))
        println("  Buses with batteries: $(length(Bset))")
        println("  Total battery capacity: $(round(total_battery_kWh, digits=1)) kWh")
        println("  Total battery power: $(round(total_battery_kW, digits=1)) kW")
    else
        println("  No battery storage found")
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

Parse voltage limits for all buses from OpenDSS.
Uses Dict{Any, Float64} to handle both substation (String) and distribution (Int) buses.
"""
function parse_voltage_limits!(data::Dict)
    println("\n--- Parsing Voltage Limits ---")
    
    # Use Dict{Any, Float64} to handle both String and Int bus names
    Vminpu = Dict{Any, Float64}()
    Vmaxpu = Dict{Any, Float64}()
    
    # Default voltage limits for all buses (substations + distribution)
    for bus in data[:Nset]
        Vminpu[bus] = 0.95
        Vmaxpu[bus] = 1.05
    end
    
    # Also store limits specifically for substation buses (for MPOPF slack bus constraints)
    Vminpu_sub = Dict{String, Float64}()
    Vmaxpu_sub = Dict{String, Float64}()
    for sub_bus in data[:Sset]
        Vminpu_sub[sub_bus] = 0.95
        Vmaxpu_sub[sub_bus] = 1.05
    end
    
    data[:Vminpu] = Vminpu
    data[:Vmaxpu] = Vmaxpu
    data[:Vminpu_sub] = Vminpu_sub
    data[:Vmaxpu_sub] = Vmaxpu_sub
    
    println("  Voltage limits for all buses: [0.95, 1.05] pu")
    println("  Substation buses ($(length(data[:Sset]))): $(data[:Sset])")
    println("  Distribution buses ($(length(data[:Nm1set]))): $(minimum(data[:Nm1set])) to $(maximum(data[:Nm1set]))")
end
