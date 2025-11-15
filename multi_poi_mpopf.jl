# multi_poi_mpopf.jl - Multi-Period Optimal Power Flow for Multi-POI systems
# 
# Self-contained script that includes parser and loads system data from ieee123_5poi_1ph
#
# Onetime script - uncomment if needed, comment once finished
import Pkg
Pkg.activate(joinpath(@__DIR__, "envs", "multi_poi"))
# Pkg.add("Crayons")
# Pkg.add("JuMP")
# Pkg.add("Ipopt")
# Pkg.add("Gurobi")
# Pkg.add("OpenDSSDirect")
# Pkg.instantiate()
# Pkg.precompile()
using JuMP
using Ipopt
using Gurobi
using LinearAlgebra
using Crayons
using Printf
using Statistics
using Plots  # For input curve plotting
using LaTeXStrings  # For LaTeX formatting in plots
using OpenDSSDirect

# Create shared Gurobi environment (suppresses repeated license messages)
const GUROBI_ENV = Gurobi.Env()

# ============================================================================
# PARSER FUNCTION (inline)
# ============================================================================

"""
    parse_system_from_dss_multi_poi(systemName, T; kwargs...)

Parse multi-POI distribution system from OpenDSS files.

Loads IEEE 123-bus 5-POI system with substations at network vertices.

Returns data dictionary with all network sets properly configured for MPOPF:
- Sset: Substation buses ["1s", "2s", "3s", "4s", "5s"]
- Nset: ALL buses (substations + distribution) [133 total]
- Nm1set: Distribution buses only [1, 2, ..., 128]
- Lset: ALL lines (POI + distribution) [132 total]
- L1set: POI connection lines [5 lines]
- Lm1set: Distribution lines only [127 lines]
"""
function parse_system_from_dss_multi_poi(systemName::String, T::Int; kwargs...)
    # Load the system in OpenDSS
    data = Dict{Symbol, Any}()
    
    # Determine path to OpenDSS files
    rawDataFolderPath = get(kwargs, :rawDataFolderPath, nothing)
    if isnothing(rawDataFolderPath)
        rawDataFolderPath = joinpath(@__DIR__, "rawData")
    end
    
    dssFilePath = joinpath(rawDataFolderPath, systemName, "Master.dss")
    
    if !isfile(dssFilePath)
        error("OpenDSS file not found: $dssFilePath")
    end
    
    # Clear any existing circuit and load new one
    OpenDSSDirect.Text.Command("clear")
    OpenDSSDirect.Text.Command("Redirect $(dssFilePath)")
    OpenDSSDirect.Solution.Solve()
    
    println("✓ Loaded system: $systemName from $dssFilePath")
    
    # Store system name and time horizon
    data[:systemName] = systemName
    data[:T] = T
    data[:Tset] = 1:T
    
    # Parse multi-POI elements (substations and POI connections)
    parse_multi_poi_elements!(data)
    
    # Extract parameters from kwargs (need for parsing)
    kVA_B = get(kwargs, :kVA_B, 1000.0)
    kV_B = get(kwargs, :kV_B, 11.547)
    
    # Parse network topology (all buses and lines, stores impedances)
    parse_network_topology!(data, kV_B, kVA_B)
    LoadShapeLoad = get(kwargs, :LoadShapeLoad, ones(T))
    LoadShapePV = get(kwargs, :LoadShapePV, zeros(T))
    LoadShapeCost_dict = get(kwargs, :LoadShapeCost_dict, Dict(i => 0.08 * ones(T) for i in 1:5))
    C_B = get(kwargs, :C_B, 1e-6 * 0.08)
    delta_t_h = get(kwargs, :delta_t_h, 1.0)
    
    data[:kVA_B] = kVA_B
    data[:kV_B] = kV_B
    data[:delta_t_h] = delta_t_h
    
    # Parse line impedances
    parse_line_impedances!(data, kVA_B, kV_B)
    
    # Parse loads
    parse_loads!(data, T, LoadShapeLoad, kVA_B)
    
    # Parse PV generators
    parse_pv_generators!(data, T, LoadShapePV, kVA_B)
    
    # Parse batteries
    parse_batteries!(data, T, kVA_B)
    
    # Parse voltage limits
    parse_voltage_limits!(data)
    
    # Store load shapes and cost data
    data[:LoadShapeLoad] = LoadShapeLoad
    data[:LoadShapePV] = LoadShapePV
    
    # Store cost profiles for each substation
    for (sub_idx, cost_profile) in LoadShapeCost_dict
        data[Symbol("LoadShapeCost_$sub_idx")] = cost_profile
    end
    data[:LoadShapeCost_dict] = LoadShapeCost_dict
    data[:C_B] = C_B
    
    println("\n" * "="^80)
    println("✓ PARSING COMPLETE")
    println("="^80)
    
    return data
end

"""Helper function to identify substation buses and POI lines"""
function parse_multi_poi_elements!(data::Dict)
    println("\n--- Parsing Multi-POI Elements ---")
    
    # Get all bus names
    all_buses = OpenDSSDirect.Circuit.AllBusNames()
    
    # Identify substation buses (buses ending with 's')
    substation_buses = String[]
    substation_bus_numbers = Int[]
    
    for bus_name in all_buses
        base_name = split(bus_name, ".")[1]
        if endswith(base_name, "s")
            push!(substation_buses, base_name)
            # Extract number (e.g., "1s" -> 1)
            bus_num = parse(Int, base_name[1:end-1])
            push!(substation_bus_numbers, bus_num)
        end
    end
    
    # Sort by bus number
    perm = sortperm(substation_bus_numbers)
    substation_buses = substation_buses[perm]
    substation_bus_numbers = substation_bus_numbers[perm]
    
    num_subs = length(substation_buses)
    println("  Number of substations: $num_subs")
    for (i, (bus, num)) in enumerate(zip(substation_buses, substation_bus_numbers))
        println("    Substation $i: Bus $bus (number $num)")
    end
    
    # Identify POI connection lines
    all_lines = OpenDSSDirect.Lines.AllNames()
    poi_lines = String[]
    poi_connections = Tuple{String, Int}[]
    
    for line_name in all_lines
        if startswith(lowercase(line_name), "poi")
            push!(poi_lines, line_name)
            
            # Get connected buses
            OpenDSSDirect.Lines.Name(line_name)
            buses = OpenDSSDirect.CktElement.BusNames()
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
    data[:substation_buses] = substation_buses
    data[:substation_bus_numbers] = substation_bus_numbers
    data[:poi_lines] = poi_lines
    data[:poi_connections] = poi_connections
end

"""Helper function to parse network topology"""
function parse_network_topology!(data::Dict, basekv::Float64, Sbase::Float64)
    println("\n--- Parsing Network Topology ---")
    
    # Get all bus names from OpenDSS
    bus_names = OpenDSSDirect.Circuit.AllBusNames()
    
    # Separate substation and distribution buses
    substation_buses = data[:substation_buses]
    dist_bus_numbers = Int[]
    
    for name in bus_names
        base_name = split(name, ".")[1]
        if !endswith(base_name, "s")
            bus_num = parse(Int, base_name)
            push!(dist_bus_numbers, bus_num)
        end
    end
    
    dist_bus_numbers = sort(unique(dist_bus_numbers))
    
    # Create network sets
    Sset = substation_buses
    Nset = Vector{Any}(vcat(substation_buses, dist_bus_numbers))
    Nm1set = dist_bus_numbers
    
    data[:Sset] = Sset
    data[:Nset] = Nset
    data[:Nm1set] = Nm1set
    data[:N] = length(Nset)
    data[:N_sub] = length(Sset)
    data[:N_dist] = length(Nm1set)
    
    println("  Substation buses (Sset): $(data[:N_sub])")
    println("  Distribution buses (Nm1set): $(data[:N_dist]) ($(minimum(Nm1set)) to $(maximum(Nm1set)))")
    println("  Total buses N = |Nset|: $(data[:N]) (substations + distribution)")
    
    # Parse lines
    line_names = OpenDSSDirect.Lines.AllNames()
    Lset = []
    L1set = []
    Lm1set = []
    
    # Store line impedances for angle computation
    line_impedances = Dict{Tuple, Tuple{Float64, Float64}}()  # line_tuple => (r_pu, x_pu)
    
    parent = Dict{Any, Union{Nothing, Any}}()
    children = Dict{Any, Vector{Any}}()
    
    for bus in Nset
        children[bus] = []
    end
    
    # Build line topology - store as directed tuples and impedances
    for line_name in line_names
        OpenDSSDirect.Lines.Name(line_name)
        buses = OpenDSSDirect.CktElement.BusNames()
        bus1_base = split(buses[1], ".")[1]
        bus2_base = split(buses[2], ".")[1]
        
        bus1 = endswith(bus1_base, "s") ? bus1_base : parse(Int, bus1_base)
        bus2 = endswith(bus2_base, "s") ? bus2_base : parse(Int, bus2_base)
        
        line_tuple = (bus1, bus2)
        push!(Lset, line_tuple)
        
        # Get and store line impedances in per-unit
        r_ohm = OpenDSSDirect.Lines.R1()
        x_ohm = OpenDSSDirect.Lines.X1()
        length_km = OpenDSSDirect.Lines.Length()
        base_impedance = basekv^2 / Sbase
        r_pu = (r_ohm * length_km) / base_impedance
        x_pu = (x_ohm * length_km) / base_impedance
        line_impedances[line_tuple] = (r_pu, x_pu)
        
        is_poi_line = startswith(lowercase(line_name), "poi")
        
        if is_poi_line
            push!(L1set, line_tuple)
        else
            push!(Lm1set, line_tuple)
        end
    end
    
    # Build parent-child relationships via BFS from substations
    # Substations are roots (no parent)
    for s in Sset
        parent[s] = nothing
    end
    
    # BFS to establish topology - use Any type for mixed String/Int buses
    queue = Vector{Any}(collect(Sset))
    visited = Set{Any}(Sset)
    
    # Also track line orientations
    line_orientations = Dict{Tuple{Any,Any}, Tuple{Any,Any}}()  # original → oriented
    
    while !isempty(queue)
        current_bus = popfirst!(queue)
        
        # Find all lines connected to current bus
        for (bus1, bus2) in Lset
            if bus1 == current_bus && !(bus2 in visited)
                # bus1 → bus2, so bus1 is parent of bus2
                parent[bus2] = bus1
                push!(children[bus1], bus2)
                push!(visited, bus2)
                push!(queue, bus2)
                line_orientations[(bus1, bus2)] = (bus1, bus2)  # Already correct orientation
            elseif bus2 == current_bus && !(bus1 in visited)
                # bus2 → bus1, so bus2 is parent of bus1
                parent[bus1] = bus2
                push!(children[bus2], bus1)
                push!(visited, bus1)
                push!(queue, bus1)
                line_orientations[(bus1, bus2)] = (bus2, bus1)  # Needs reversal
            end
        end
    end
    
    # Reorient Lset to match parent→child direction
    Lset_oriented = [haskey(line_orientations, line) ? line_orientations[line] : line for line in Lset]
    L1set_oriented = [haskey(line_orientations, line) ? line_orientations[line] : line for line in L1set]
    Lm1set_oriented = [haskey(line_orientations, line) ? line_orientations[line] : line for line in Lm1set]
    
    data[:Lset] = Lset_oriented
    data[:L1set] = L1set_oriented
    data[:Lm1set] = Lm1set_oriented
    
    # Store undirected versions for topology computation
    data[:Lset_undirected] = Lset  # Original undirected lines
    data[:L1set_undirected] = L1set  # Original undirected POI lines
    data[:Lm1set_undirected] = Lm1set  # Original undirected distribution lines
    
    data[:line_orientations] = line_orientations  # Store for impedance parsing
    data[:line_impedances] = line_impedances  # Store impedances for angle computation
    data[:parent] = parent
    data[:children] = children
    data[:m] = length(Lset_oriented)
    data[:m_poi] = length(L1set_oriented)
    data[:m_dist] = length(Lm1set_oriented)
    
    println("  POI lines (L1set): $(length(L1set_oriented))")
    println("  Distribution lines (Lm1set): $(length(Lm1set_oriented))")
    println("  Total lines (Lset): $(length(Lset_oriented))")
    
    println("\n  POI Line Details:")
    for (bus1, bus2) in L1set_oriented
        println("    $bus1 → $bus2")
    end
end

"""Helper function to parse line impedances"""
function parse_line_impedances!(data::Dict, kVA_B::Float64, kV_B::Float64)
    println("\n--- Parsing Line Impedances ---")
    
    Z_B = kV_B^2 / (kVA_B / 1000)
    
    rdict_pu = Dict{Any, Float64}()
    xdict_pu = Dict{Any, Float64}()
    poi_rdict = Dict{String, Float64}()
    poi_xdict = Dict{String, Float64}()
    
    line_names = OpenDSSDirect.Lines.AllNames()
    
    for line_name in line_names
        OpenDSSDirect.Lines.Name(line_name)
        buses = OpenDSSDirect.CktElement.BusNames()
        bus1_base = split(buses[1], ".")[1]
        bus2_base = split(buses[2], ".")[1]
        
        bus1 = endswith(bus1_base, "s") ? bus1_base : parse(Int, bus1_base)
        bus2 = endswith(bus2_base, "s") ? bus2_base : parse(Int, bus2_base)
        
        R_ohm = OpenDSSDirect.Lines.RMatrix()[1]
        X_ohm = OpenDSSDirect.Lines.XMatrix()[1]
        
        R_pu = R_ohm / Z_B
        X_pu = X_ohm / Z_B
        
        line_tuple_original = (bus1, bus2)
        
        # Use oriented line tuple if available
        line_orientations = get(data, :line_orientations, Dict())
        line_tuple_oriented = get(line_orientations, line_tuple_original, line_tuple_original)
        
        rdict_pu[line_tuple_oriented] = R_pu
        xdict_pu[line_tuple_oriented] = X_pu
        
        if startswith(lowercase(line_name), "poi")
            sub_bus = endswith(string(bus1), "s") ? bus1 : bus2
            poi_rdict[sub_bus] = R_ohm
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
    
    println("\n  POI Line Impedances:")
    for sub_bus in data[:substation_buses]
        R_ohm = get(poi_rdict, sub_bus, NaN)
        X_ohm = get(poi_xdict, sub_bus, NaN)
        R_pu = R_ohm / Z_B
        X_pu = X_ohm / Z_B
        println("    $sub_bus: R=$(round(R_ohm, digits=8)) Ω ($(round(R_pu, digits=10)) pu), X=$(round(X_ohm, digits=8)) Ω ($(round(X_pu, digits=10)) pu)")
    end
end

"""Helper function to find which substation zone a bus belongs to by tracing parent path"""
function find_zone_substation(bus::Int, parent_dict::Dict, Sset)
    current = bus
    # Trace back through parents until we reach a substation
    visited = Set{Any}()
    while true
        if current in visited
            # Cycle detected, default to first substation
            return first(Sset)
        end
        push!(visited, current)
        
        # Check if current node is a substation
        if current in Sset || string(current) in string.(Sset)
            return string(current) in string.(Sset) ? string(current) : string(current) * "s"
        end
        
        # Move to parent
        if haskey(parent_dict, current)
            current = parent_dict[current]
        else
            # No parent found, default to first substation
            return first(Sset)
        end
    end
end

"""Helper function to parse loads"""
function parse_loads!(data::Dict, T::Int, LoadShapeLoad::Vector, kVA_B::Float64)
    println("\n--- Parsing Loads ---")
    
    load_names = OpenDSSDirect.Loads.AllNames()
    NLset = Int[]
    p_L_R_pu = Dict{Int, Float64}()
    q_L_R_pu = Dict{Int, Float64}()
    
    for load_name in load_names
        OpenDSSDirect.Loads.Name(load_name)
        bus_full = OpenDSSDirect.CktElement.BusNames()[1]
        bus_num = parse(Int, split(bus_full, ".")[1])
        
        kW = OpenDSSDirect.Loads.kW()
        kVAR = OpenDSSDirect.Loads.kvar()
        
        p_base_pu = kW / kVA_B
        q_base_pu = kVAR / kVA_B
        
        if haskey(p_L_R_pu, bus_num)
            p_L_R_pu[bus_num] += p_base_pu
            q_L_R_pu[bus_num] += q_base_pu
        else
            p_L_R_pu[bus_num] = p_base_pu
            q_L_R_pu[bus_num] = q_base_pu
            push!(NLset, bus_num)
        end
    end
    
    Nm1set = data[:Nm1set]
    max_bus = maximum(Nm1set)
    p_L_pu = zeros(max_bus, T)
    q_L_pu = zeros(max_bus, T)
    
    # Create spatial load variation: assign each bus to a zone based on parent substation
    # This creates reproducible, systematic spatial diversity
    parent_dict = data[:parent]  # Maps each bus to its parent substation
    
    # Create zone-specific load profiles with phase shifts
    # Zone determined by which substation feeds the bus
    zone_profiles = Dict{String, Vector{Float64}}()
    if T == 1
        for sub in data[:Sset]
            zone_profiles[sub] = [1.0]
        end
    else
        # Diverse daily patterns for each zone: 2 sinusoids, 2 multi-level, 1 ramping
        
        # Zone 1: Sinusoid (smooth residential pattern - peak midday)
        # Range: 0.75 + 0.25*[0,1] = [0.75, 1.00]
        zone_profiles["1s"] = 0.75 .+ 0.25 .* (sin.(range(0, 2π, length=T) .+ 0.0) .+ 1) ./ 2
        
        # Zone 2: Sinusoid (smooth commercial pattern - peak afternoon)
        # Range: 0.75 + 0.25*[0,1] = [0.75, 1.00]
        zone_profiles["2s"] = 0.75 .+ 0.25 .* (sin.(range(0, 2π, length=T) .+ π/3) .+ 1) ./ 2
        
        # Zone 3: Multi-level square (industrial with shifts)
        load_3 = fill(0.70, T)  # Base: low night
        load_3[7:12] .= 1.00    # Morning shift: high
        load_3[13:18] .= 0.90   # Afternoon shift: medium-high
        load_3[19:22] .= 0.80   # Evening: medium
        zone_profiles["3s"] = load_3
        
        # Zone 4: Two-level square (office building - on/off)
        load_4 = fill(0.65, T)  # Base: minimal (off-hours)
        load_4[8:18] .= 1.00    # Business hours: high
        zone_profiles["4s"] = load_4
        
        # Zone 5: Ramping (charging station - gradual increase then drop)
        load_5 = fill(0.75, T)
        load_5[1:8] .= 0.70 .+ (0:7) .* 0.03        # Morning ramp up
        load_5[9:16] .= 0.94 .+ (0:7) .* 0.008      # Peak and continue up slightly
        load_5[17:24] .= 1.00 .- (0:7) .* 0.045     # Evening ramp down
        zone_profiles["5s"] = load_5
    end
    
    # Apply spatially-varying temporal profiles
    for bus in NLset
        # Find which substation zone this bus belongs to
        zone_sub = find_zone_substation(bus, parent_dict, data[:Sset])
        load_profile = zone_profiles[zone_sub]
        
        for t in 1:T
            p_L_pu[bus, t] = p_L_R_pu[bus] * load_profile[t]
            q_L_pu[bus, t] = q_L_R_pu[bus] * load_profile[t]
        end
    end
    
    data[:zone_profiles] = zone_profiles  # Store for reference
    
    data[:NLset] = sort(NLset)
    data[:N_L] = length(NLset)
    data[:p_L_pu] = p_L_pu
    data[:q_L_pu] = q_L_pu
    data[:p_L_R_pu] = p_L_R_pu
    data[:q_L_R_pu] = q_L_R_pu
    
    total_P_kW = sum(values(p_L_R_pu)) * kVA_B
    total_Q_kVAr = sum(values(q_L_R_pu)) * kVA_B
    
    println("  Buses with loads: $(data[:N_L])")
    println("  Total base load: $(round(total_P_kW, digits=1)) kW, $(round(total_Q_kVAr, digits=1)) kVAr")
end

"""Helper function to parse PV generators"""
function parse_pv_generators!(data::Dict, T::Int, LoadShapePV::Vector, kVA_B::Float64)
    println("\n--- Parsing PV Systems ---")
    
    gen_names = OpenDSSDirect.PVsystems.AllNames()
    Dset = Int[]
    p_D_R_pu = Dict{Int, Float64}()
    S_D_R = Dict{Int, Float64}()
    
    Nm1set = data[:Nm1set]
    max_bus = maximum(Nm1set)
    p_D_pu = zeros(max_bus, T)
    
    if !isempty(gen_names) && gen_names[1] != "NONE"
        for gen_name in gen_names
            OpenDSSDirect.PVsystems.Name(gen_name)
            bus_full = OpenDSSDirect.CktElement.BusNames()[1]
            bus_num = parse(Int, split(bus_full, ".")[1])
            
            kW_rated = OpenDSSDirect.PVsystems.kW()
            kVA_rated = OpenDSSDirect.PVsystems.kVARated()
            
            p_rated_pu = kW_rated / kVA_B
            S_rated_pu = kVA_rated / kVA_B
            
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
        
        for bus in Dset
            for t in 1:T
                p_D_pu[bus, t] = p_D_R_pu[bus] * LoadShapePV[t]
            end
        end
        
        println("  PV systems found: $(length(Dset))")
    else
        println("  No PV systems found")
    end
    
    data[:Dset] = sort(Dset)
    data[:n_D] = length(Dset)
    data[:p_D_pu] = p_D_pu
    data[:p_D_R_pu] = p_D_R_pu
    data[:S_D_R] = S_D_R
end

"""Helper function to parse batteries"""
function parse_batteries!(data::Dict, T::Int, kVA_B::Float64)
    println("\n--- Parsing Battery Storage ---")
    
    storage_names = OpenDSSDirect.Storages.AllNames()
    Bset = Int[]
    B0_pu = Dict{Int, Float64}()
    B_R = Dict{Int, Float64}()
    B_R_pu = Dict{Int, Float64}()
    P_B_R = Dict{Int, Float64}()
    P_B_R_pu = Dict{Int, Float64}()
    S_B_R_pu = Dict{Int, Float64}()
    soc_min = Dict{Int, Float64}()
    soc_max = Dict{Int, Float64}()
    
    E_BASE = kVA_B * 1.0
    
    if !isempty(storage_names) && storage_names[1] != "NONE"
        for storage_name in storage_names
            OpenDSSDirect.Storages.Name(storage_name)
            bus_full = OpenDSSDirect.CktElement.BusNames()[1]
            bus_num = parse(Int, split(bus_full, ".")[1])
            
            kWh_stored = OpenDSSDirect.Storages.kWhStored()
            kWh_rated = OpenDSSDirect.Storages.kWhRated()
            kW_rated = OpenDSSDirect.Storages.kW()
            kVA_rated_storage = OpenDSSDirect.Storages.kVA()
            
            B0_pu[bus_num] = kWh_stored / E_BASE
            B_R[bus_num] = kWh_rated
            B_R_pu[bus_num] = kWh_rated / E_BASE
            P_B_R[bus_num] = kW_rated
            P_B_R_pu[bus_num] = kW_rated / kVA_B
            S_B_R_pu[bus_num] = kVA_rated_storage / kVA_B
            
            soc_min[bus_num] = 0.2
            soc_max[bus_num] = 1.0
            
            push!(Bset, bus_num)
        end
        
        println("  Battery storage found: $(length(Bset))")
    else
        println("  No battery storage found")
    end
    
    data[:Bset] = sort(Bset)
    data[:n_B] = length(Bset)
    data[:B0_pu] = B0_pu
    data[:B_R] = B_R
    data[:B_R_pu] = B_R_pu
    data[:P_B_R] = P_B_R
    data[:P_B_R_pu] = P_B_R_pu
    data[:S_B_R_pu] = S_B_R_pu
    data[:soc_min] = soc_min
    data[:soc_max] = soc_max
end

"""Helper function to parse voltage limits"""
function parse_voltage_limits!(data::Dict)
    println("\n--- Parsing Voltage Limits ---")
    
    Vminpu = Dict{Any, Float64}()
    Vmaxpu = Dict{Any, Float64}()
    
    for bus in data[:Nset]
        Vminpu[bus] = 0.95
        # Vminpu[bus] = 0.90  # RELAXED VERSION - commented out
        # Vmaxpu[bus] = 1.05
        Vmaxpu[bus] = 1.10  # RELAXED VERSION - commented out
    end
    
    Vminpu_sub = Dict{String, Float64}()
    Vmaxpu_sub = Dict{String, Float64}()
    for sub_bus in data[:Sset]
        Vminpu_sub[sub_bus] = 0.95
        # Vminpu_sub[sub_bus] = 0.90  # RELAXED VERSION - commented out
        # Vmaxpu_sub[sub_bus] = 1.05
        Vmaxpu_sub[sub_bus] = 1.10  # RELAXED VERSION - commented out
    end
    
    data[:Vminpu] = Vminpu
    data[:Vmaxpu] = Vmaxpu
    data[:Vminpu_sub] = Vminpu_sub
    data[:Vmaxpu_sub] = Vmaxpu_sub
    
    println("  Voltage limits for all buses: [0.95, 1.05] pu")
    println("  Substation buses ($(length(data[:Sset]))): $(data[:Sset])")
    println("  Distribution buses ($(length(data[:Nm1set]))): $(minimum(data[:Nm1set])) to $(maximum(data[:Nm1set]))")
end

# ============================================================================
# SYSTEM PARAMETERS
# ============================================================================

# System identification
systemName = "ieee123_5poi_1ph"

# Simulation parameters
T = 24  # Number of time periods (testing with single period)
delta_t_h = 1.0  # Time step duration in hours
kVA_B = 1000.0
kV_B = 11.547  # 20kV line-to-line -> 11.547kV line-to-neutral

# User-defined parameter: Should substations have same voltage levels?
# - true: Non-slack substations have fixed voltage 
# - false (default): Non-slack substation voltages become decision variables
same_voltage_levels = false

# Plotting options
showPlots = false  # Set to true to display plots in GUI
savePlots = true   # Set to true to save plots to file

# ============================================================================
# TIME HORIZON PARAMETERS
# ============================================================================

# Time-varying load profile (sinusoidal)
if T == 1
    LoadShapeLoad = [1.0]  # Unity load for single period
else
    LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2
end

# Time-varying energy cost ($/kWh) - different for each substation
# Diverse patterns: 2 sinusoids, 2 multi-level squares, 1 flat
if T == 1
    LoadShapeCost_dict = Dict(
        1 => [0.15],
        2 => [0.15],
        3 => [0.15],
        4 => [0.15],
        5 => [0.15]
    )
else
    # Subs 1: Sinusoid (peak early)
    LoadShapeCost_dict = Dict(
        1 => 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .+ π/4) .+ 1) ./ 2,
    )
    
    # Subs 2: Sinusoid (peak late)
    LoadShapeCost_dict[2] = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .- π/3) .+ 1) ./ 2
    
    # Subs 3: Multi-level square wave (3 levels: low-med-high-med)
    cost_3 = fill(0.10, T)  # Base: medium
    cost_3[1:6] .= 0.18      # Hours 1-6: high (night peak)
    cost_3[7:12] .= 0.08     # Hours 7-12: low (morning off-peak)
    cost_3[13:18] .= 0.18    # Hours 13-18: high (afternoon peak)
    cost_3[19:24] .= 0.12    # Hours 19-24: medium (evening)
    LoadShapeCost_dict[3] = cost_3
    
    # Subs 4: Multi-level square wave (2 levels: low-high alternating)
    cost_4 = fill(0.08, T)  # Base: low
    cost_4[9:20] .= 0.20     # Hours 9-20: high (peak period)
    LoadShapeCost_dict[4] = cost_4
    
    # Subs 5: Flat (constant cost)
    LoadShapeCost_dict[5] = fill(0.15, T)  # Constant at mid-level
end

LoadShapePV = zeros(T)  # No PV for now
C_B = 1e-6 * 0.08

# ============================================================================
# PARSE SYSTEM DATA FROM OPENDSS
# ============================================================================

println("\n" * "="^80)
println("PARSING MULTI-POI SYSTEM FROM OPENDSS")
println("="^80)

data = parse_system_from_dss_multi_poi(
    systemName,
    T;
    kVA_B=kVA_B,
    kV_B=kV_B,
    LoadShapeLoad=LoadShapeLoad,
    LoadShapePV=LoadShapePV,
    LoadShapeCost_dict=LoadShapeCost_dict,
    C_B=C_B,
    delta_t_h=delta_t_h
)

# Add user-defined voltage configuration flag
data[:same_voltage_levels] = same_voltage_levels

"""
Export detailed model formulation to text file
"""
function export_model_to_file(model, filename, output_dir=".")
    filepath = joinpath(output_dir, filename)
    open(filepath, "w") do io
        println(io, "="^80)
        println(io, "MPOPF MODEL FORMULATION")
        println(io, "="^80)
        println(io)
        
        # Objective
        println(io, "OBJECTIVE:")
        println(io, objective_function(model))
        println(io)
        println(io, "="^80)
        
        # Variables
        println(io, "\nVARIABLES:")
        for v in all_variables(model)
            println(io, "  $(name(v)): lower=$(has_lower_bound(v) ? lower_bound(v) : "none"), upper=$(has_upper_bound(v) ? upper_bound(v) : "none")")
        end
        println(io)
        println(io, "="^80)
        
        # Constraints
        println(io, "\nCONSTRAINTS:")
        constraint_count = 0
        for (F, S) in list_of_constraint_types(model)
            cons_list = all_constraints(model, F, S)
            for con in cons_list
                constraint_count += 1
                println(io, "\n[$(constraint_count)] $(name(con)):")
                println(io, "  $(constraint_object(con))")
            end
        end
        println(io)
        println(io, "="^80)
        println(io, "\nTotal constraints: $constraint_count")
        println(io, "="^80)
    end
    println("Model exported to: $filename")
end

# Print summary
println("\n" * "="^80)
println("SYSTEM DATA LOADED SUCCESSFULLY")
println("="^80)
println("System: $(data[:systemName])")
println("Substations: $(data[:num_substations])")
println("Substation buses: $(data[:substation_buses])")
println("POI connections:")
for (sub_bus, dist_bus) in data[:poi_connections]
    println("  $sub_bus → Bus $dist_bus")
end
println("\nNetwork:")
println("  Distribution buses: $(data[:N])")
println("  Distribution branches: $(data[:m])")
println("  Total load: $(round(sum(values(data[:p_L_R_pu])) * kVA_B, digits=1)) kW")
println("\nTime horizon: T = $T, Δt = $(delta_t_h) hours")
println("="^80)

# ============================================================================
# MULTI-POI MPOPF SOLVER
# ============================================================================

"""
    solve_multi_poi_mpopf(data; slack_substation="1s", solver=:gurobi)

Solve Multi-Period Optimal Power Flow for multi-POI distribution system.

Arguments:
- data: System data dictionary from parser
- slack_substation: Which substation is slack (default: "1s")
  Valid values: "1s", "2s", "3s", "4s", "5s"
- solver: Optimization solver (:gurobi or :ipopt, default: :gurobi)

Returns:
- result: Dictionary with solution and statistics
"""
function solve_multi_poi_mpopf(data; slack_substation::String="1s", solver::Symbol=:gurobi)
    
    println("\n" * "="^80)
    println("SOLVING MULTI-POI MPOPF")
    println("="^80)
    println("Slack substation: $slack_substation")
    println("Solver: $solver")
    println("="^80)
    
    # ========== 1. UNPACK DATA ==========
    Sset = data[:Sset]  # Substation buses ["1s", "2s", ...]
    Nset = data[:Nset]  # All buses (substations + distribution)
    Nm1set = data[:Nm1set]  # Distribution buses only [1, 2, ..., 128]
    Lset = data[:Lset]  # All lines (POI + distribution)
    L1set = data[:L1set]  # POI connection lines
    Lm1set = data[:Lm1set]  # Distribution lines only
    Tset = data[:Tset]
    T = data[:T]
    
    NLset = data[:NLset]  # Buses with loads
    Dset = data[:Dset]  # Buses with PV
    Bset = data[:Bset]  # Buses with batteries
    
    p_L_pu = data[:p_L_pu]
    q_L_pu = data[:q_L_pu]
    p_L_R_pu = data[:p_L_R_pu]  # Rated load power per bus
    p_D_pu = data[:p_D_pu]
    
    P_B_R_pu = data[:P_B_R_pu]
    B_R_pu = data[:B_R_pu]
    B0_pu = data[:B0_pu]
    soc_min = data[:soc_min]
    soc_max = data[:soc_max]
    S_D_R = data[:S_D_R]
    
    rdict_pu = data[:rdict_pu]
    xdict_pu = data[:xdict_pu]
    Vminpu = data[:Vminpu]
    Vmaxpu = data[:Vmaxpu]
    
    LoadShapeCost_dict = data[:LoadShapeCost_dict]
    C_B = data[:C_B]
    delta_t_h = data[:delta_t_h]
    kVA_B = data[:kVA_B]
    
    parent = data[:parent]
    children = data[:children]
    poi_connections = data[:poi_connections]
    
    Δt = delta_t_h
    P_BASE = kVA_B
    
    # Verify slack substation is valid
    if !(slack_substation in Sset)
        error("Invalid slack substation: $slack_substation. Must be one of: $Sset")
    end
    
    # Create mapping from substation bus name to index
    sub_to_idx = Dict(s => i for (i, s) in enumerate(Sset))
    num_subs = length(Sset)
    
    println("\nSystem overview:")
    println("  Substations: $num_subs")
    println("  Distribution buses: $(length(Nm1set))")
    println("  Total buses: $(length(Nset))")
    println("  POI lines: $(length(L1set))")
    println("  Distribution lines: $(length(Lm1set))")
    println("  Time periods: $T")
    
    # ========== 2. CREATE MODEL ==========
    model = Model()
    if solver == :gurobi
        set_optimizer(model, () -> Gurobi.Optimizer(GUROBI_ENV))
        set_optimizer_attribute(model, "NonConvex", 2)
        set_optimizer_attribute(model, "OutputFlag", 0)  # Silent mode
    else
        set_optimizer(model, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0)  # Silent mode
        set_optimizer_attribute(model, "max_iter", 5000)
        set_optimizer_attribute(model, "tol", 1e-6)
    end
    
    # ========== 3. DEFINE VARIABLES ==========
    
    # Compute total system rated load for substation power limits
    P_L_R_System = sum(values(p_L_R_pu))  # Total rated load in pu
    P_Subs_max = 0.4 * P_L_R_System  # 40% of total system load
    
    println("\nSubstation power limits:")
    println("  Total system rated load: $(round(P_L_R_System, digits=4)) pu = $(round(P_L_R_System * kVA_B, digits=1)) kW")
    println("  Max per substation: $(round(P_Subs_max, digits=4)) pu = $(round(P_Subs_max * kVA_B, digits=1)) kW")
    
    # Substation power injections (indexed by substation name and time)
    # P_Subs_max = Inf  # RELAXED VERSION - commented out to enforce 70% limit
    @variable(model, 0 <= P_Subs[s in Sset, t in Tset] <= P_Subs_max)
    @variable(model, Q_Subs[s in Sset, t in Tset])
    
    # Branch power flows (works with mixed bus types)
    @variable(model, P[line in Lset, t in Tset])
    @variable(model, Q[line in Lset, t in Tset])
    @variable(model, ℓ[line in Lset, t in Tset] >= 0)
    
    # Voltage magnitude squared (for all buses)
    @variable(model, v[bus in Nset, t in Tset])
    
    # Battery variables (indexed by distribution bus number)
    @variable(model, P_B[j in Bset, t in Tset])
    @variable(model, B[j in Bset, t in Tset])
    
    # PV reactive power
    @variable(model, q_D[j in Dset, t in Tset])
    
    println("\nVariables created:")
    println("  P_Subs, Q_Subs: $(num_subs) × $T")
    println("  P, Q, ℓ (branch flows): $(length(Lset)) × $T each")
    println("  v (voltage squared): $(length(Nset)) × $T")
    println("  P_B, B (battery): $(length(Bset)) × $T each")
    println("  q_D (PV reactive): $(length(Dset)) × $T")
    
    # ========== 4. OBJECTIVE FUNCTION ==========
    
    # Energy cost from each substation (different cost profiles)
    @expression(model, energy_cost,
        sum(LoadShapeCost_dict[sub_to_idx[s]][t] * P_Subs[s, t] * P_BASE * Δt 
            for s in Sset, t in Tset))
    
    # Battery degradation cost
    @expression(model, battery_cost,
        sum(C_B * (P_B[j, t] * P_BASE)^2 * Δt for j in Bset, t in Tset))
    
    @objective(model, Min, energy_cost + battery_cost)
    
    println("\nObjective: Minimize energy cost + battery degradation")
    
    # ========== 5. CONSTRAINTS ==========
    
    println("\nAdding constraints...")
    
    # Thermal/ampacity limit constants (defined once outside time loop)
    I_max_pu = 5.0
    ℓ_max = I_max_pu^2  # = 25 pu
    
    for t in Tset
        # ----- 5.1 POWER BALANCE AT SUBSTATIONS -----
        # Each substation injects power through all its connected POI lines
        for s in Sset
            # Find all POI lines connected to this substation
            poi_lines_s = [(i, j) for (i, j) in L1set if i == s]
            
            if isempty(poi_lines_s)
                error("No POI lines found for substation $s in L1set")
            end
            
            # Power balance: Substation power = Sum of power flowing through all POI lines
            @constraint(model, P_Subs[s, t] == sum(P[line, t] for line in poi_lines_s),
                base_name = "RealPowerBalance_Substation_$(s)_t$(t)")
            @constraint(model, Q_Subs[s, t] == sum(Q[line, t] for line in poi_lines_s),
                base_name = "ReactivePowerBalance_Substation_$(s)_t$(t)")
        end
        
        # ----- 5.2 NODAL POWER BALANCE (DISTRIBUTION BUSES) -----
        for j in Nm1set
            # Find all parent lines feeding bus j
            # In general: all lines (i,j_bus) where j_bus == j and i feeds j
            # In radial: parent[j] gives the unique parent, so line is (parent[j], j)
            parent_lines_j = [(i_bus, j_bus) for (i_bus, j_bus) in Lset if j_bus == j && (i_bus, j_bus) in keys(rdict_pu) && i_bus == parent[j]]
            
            # Sum over all child lines emanating from bus j
            sum_Pjk = isempty(children[j]) ? 0.0 : sum(P[(j, k), t] for k in children[j])
            sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])
            
            # Sum over all parent lines (with losses)
            sum_Pij_minus_loss = isempty(parent_lines_j) ? 0.0 : sum(P[line, t] - rdict_pu[line] * ℓ[line, t] for line in parent_lines_j)
            sum_Qij_minus_loss = isempty(parent_lines_j) ? 0.0 : sum(Q[line, t] - xdict_pu[line] * ℓ[line, t] for line in parent_lines_j)
            
            # Nodal injections/withdrawals
            p_L_j_t = (j in NLset) ? p_L_pu[j, t] : 0.0
            q_L_j_t = (j in NLset) ? q_L_pu[j, t] : 0.0
            p_D_j_t = (j in Dset) ? p_D_pu[j, t] : 0.0
            q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0
            P_B_j_t = (j in Bset) ? P_B[j, t] : 0.0
            
            # Real power balance: Σ(P_ij - r_ij*ℓ_ij) - Σ P_jk = p_j
            @constraint(model, sum_Pij_minus_loss - sum_Pjk == p_L_j_t - p_D_j_t - P_B_j_t,
                base_name = "RealPowerBalance_Node$(j)_t$(t)")
            
            # Reactive power balance: Σ(Q_ij - x_ij*ℓ_ij) - Σ Q_jk = q_j
            @constraint(model, sum_Qij_minus_loss - sum_Qjk == q_L_j_t - q_D_j_t,
                base_name = "ReactivePowerBalance_Node$(j)_t$(t)")
        end
        
        # ----- 5.3 VOLTAGE DROP CONSTRAINTS (BFM) -----
        for line in Lset
            (i, j) = line
            r_ij = rdict_pu[line]
            x_ij = xdict_pu[line]
            
            @constraint(model,
                v[j, t] == v[i, t] - 2 * (r_ij * P[line, t] + x_ij * Q[line, t]) + 
                           (r_ij^2 + x_ij^2) * ℓ[line, t],
                base_name = "VoltageDrop_Line_$(i)_$(j)_t$(t)")
        end
        
        # ----- 5.4 CONIC RELAXATION (BCPF) -----
        for line in Lset
            (i, j) = line
            @constraint(model, P[line, t]^2 + Q[line, t]^2 <= v[i, t] * ℓ[line, t],
                base_name = "BCPF_Line_$(i)_$(j)_t$(t)")
        end
        
        # ----- 5.4b THERMAL/AMPACITY LIMITS -----
        # Cap current squared at (I_max)^2 where I_max = 5 pu
        for line in Lset
            (i, j) = line
            @constraint(model, ℓ[line, t] <= ℓ_max,
                base_name = "ThermalLimit_Line_$(i)_$(j)_t$(t)")
        end
        
        # ----- 5.5 VOLTAGE CONSTRAINTS -----
        # Slack substation: fixed voltage
        @constraint(model, v[slack_substation, t] == 1.03^2,
            base_name = "SlackVoltage_$(slack_substation)_t$(t)")
        
        # Non-slack substations: voltage within limits (decision variables)
        for s in Sset
            if s != slack_substation
                @constraint(model, Vminpu[s]^2 <= v[s, t] <= Vmaxpu[s]^2,
                    base_name = "VoltageLimits_Substation_$(s)_t$(t)")
            end
        end
        
        # Distribution buses: voltage within limits
        for j in Nm1set
            @constraint(model, Vminpu[j]^2 <= v[j, t] <= Vmaxpu[j]^2,
                base_name = "VoltageLimits_Node$(j)_t$(t)")
        end
        
        # ----- 5.6 PV REACTIVE POWER LIMITS -----
        for j in Dset
            p_D_val = p_D_pu[j, t]
            S_D_R_val = S_D_R[j]
            q_max_t = sqrt(max(0.0, S_D_R_val^2 - p_D_val^2))
            @constraint(model, -q_max_t <= q_D[j, t] <= q_max_t,
                base_name = "PV_ReactivePowerLimits_Node$(j)_t$(t)")
        end
        
        # ----- 5.7 BATTERY CONSTRAINTS -----
        for j in Bset
            # SOC trajectory
            if t == 1
                @constraint(model, B[j, t] == B0_pu[j] - P_B[j, t] * Δt,
                    base_name = "Battery_SOC_Initial_Node$(j)_t$(t)")
            else
                @constraint(model, B[j, t] == B[j, t-1] - P_B[j, t] * Δt,
                    base_name = "Battery_SOC_Dynamics_Node$(j)_t$(t)")
            end
            
            # SOC limits
            @constraint(model, soc_min[j] * B_R_pu[j] <= B[j, t] <= soc_max[j] * B_R_pu[j],
                base_name = "Battery_SOC_Limits_Node$(j)_t$(t)")
            
            # Power limits
            @constraint(model, -P_B_R_pu[j] <= P_B[j, t] <= P_B_R_pu[j],
                base_name = "Battery_PowerLimits_Node$(j)_t$(t)")
        end
    end
    
    println("Constraints added successfully!")
    
    # ========== 5.8 SAVE MODEL TO FILE ==========
    # Note: results_base_dir will be defined by the calling context
    # For now, we'll export to current directory and move it later
    
    # ========== 6. SOLVE ==========
    println("\n" * "="^80)
    println("SOLVING...")
    println("="^80)
    
    solve_start = time()
    optimize!(model)
    solve_time = time() - solve_start
    
    # ========== 7. EXTRACT RESULTS ==========
    status = termination_status(model)
    
    println("\n" * "="^80)
    if status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL
        println(Crayon(foreground = :green, bold = true)("✓ OPTIMIZATION SUCCESSFUL"))
    else
        println(Crayon(foreground = :red, bold = true)("✗ OPTIMIZATION FAILED: $status"))
    end
    println("="^80)
    
    result = Dict(
        :status => status,
        :solve_time => solve_time,
        :slack_substation => slack_substation
    )
    
    if has_values(model)
        result[:objective] = objective_value(model)
        result[:energy_cost] = value(energy_cost)
        result[:battery_cost] = value(battery_cost)
        
        # Extract all decision variables
        # Substation powers (separate dict for easy access)
        result[:P_Subs] = Dict(s => [value(P_Subs[s, t]) for t in Tset] for s in Sset)
        result[:Q_Subs] = Dict(s => [value(Q_Subs[s, t]) for t in Tset] for s in Sset)
        
        # Branch power flows (P, Q, ℓ)
        result[:P] = Dict(line => [value(P[line, t]) for t in Tset] for line in Lset)
        result[:Q] = Dict(line => [value(Q[line, t]) for t in Tset] for line in Lset)
        result[:ℓ] = Dict(line => [value(ℓ[line, t]) for t in Tset] for line in Lset)
        
        # Voltage magnitude squared (all buses)
        result[:v] = Dict(bus => [value(v[bus, t]) for t in Tset] for bus in Nset)
        result[:v_subs] = Dict(s => [value(v[s, t]) for t in Tset] for s in Sset)  # Substation subset for convenience
        
        # Battery variables (if any)
        result[:P_B] = Dict(j => [value(P_B[j, t]) for t in Tset] for j in Bset)
        result[:B] = Dict(j => [value(B[j, t]) for t in Tset] for j in Bset)
        
        # PV reactive power (if any)
        result[:q_D] = Dict(j => [value(q_D[j, t]) for t in Tset] for j in Dset)
        
        # Print summary
        println("\nSolution Summary:")
        println("  Objective value: \$$(round(result[:objective], digits=2))")
        println("  Energy cost: \$$(round(result[:energy_cost], digits=2))")
        println("  Battery cost: \$$(round(result[:battery_cost], digits=2))")
        println("  Solve time: $(round(solve_time, digits=2)) seconds")
        
        println("\n  Substation Power Dispatch (time-averaged):")
        for s in Sset
            P_avg = mean(result[:P_Subs][s]) * P_BASE
            Q_avg = mean(result[:Q_Subs][s]) * P_BASE
            println("    $s: P = $(round(P_avg, digits=1)) kW, Q = $(round(Q_avg, digits=1)) kVAr")
        end
        
        println("\n  Substation Voltages (time-averaged):")
        for s in Sset
            V_avg = sqrt(mean(result[:v_subs][s]))
            if s == slack_substation
                println("    $s: V = $(round(V_avg, digits=4)) pu (slack - fixed)")
            else
                V_min = sqrt(minimum(result[:v_subs][s]))
                V_max = sqrt(maximum(result[:v_subs][s]))
                println("    $s: V = $(round(V_avg, digits=4)) pu (range: [$(round(V_min, digits=4)), $(round(V_max, digits=4))])")
            end
        end
        
        # PV generation summary
        if !isempty(Dset)
            println("\n  PV Generation (time-averaged):")
            for j in sort(collect(Dset))
                p_D_avg = mean(p_D_pu[j, :]) * P_BASE
                q_D_avg = mean([value(q_D[j, t]) for t in Tset]) * P_BASE
                println("    Bus $j: P_pv = $(round(p_D_avg, digits=1)) kW, Q_pv = $(round(q_D_avg, digits=1)) kVAr")
            end
            total_pv_gen = sum(mean(p_D_pu[j, :]) for j in Dset) * P_BASE
            println("    Total PV: $(round(total_pv_gen, digits=1)) kW (time-averaged)")
        else
            println("\n  No PV systems installed")
        end
        
        # # DIAGNOSTIC OUTPUT - Commented out for production runs
        # # Power balance verification
        # println("\n  Power Balance Verification (t=1):")
        # t = 1
        # total_P_Subs = sum(result[:P_Subs][s][t] for s in Sset) * P_BASE
        # total_load = sum(p_L_pu[j, t] for j in Nm1set) * P_BASE
        # 
        # # Compute losses correctly: P_loss = r * ℓ (both in pu), then scale to kW
        # total_losses_pu = sum(rdict_pu[line] * result[:ℓ][line][t] for line in Lset)
        # total_losses = total_losses_pu * P_BASE
        # 
        # # Debug: Check a few line losses
        # sample_line_losses = []
        # for line in collect(Lset)[1:min(3, length(Lset))]
        #     r_pu = rdict_pu[line]
        #     ℓ_val = result[:ℓ][line][t]
        #     loss_pu = r_pu * ℓ_val
        #     loss_kW = loss_pu * P_BASE
        #     push!(sample_line_losses, (line, r_pu, ℓ_val, loss_kW))
        # end
        # 
        # println("    Total substation power: $(round(total_P_Subs, digits=1)) kW")
        # println("    Total load demand: $(round(total_load, digits=1)) kW")
        # println("    Total network losses: $(round(total_losses, digits=3)) kW ($(round(total_losses_pu, digits=6)) pu)")
        # println("    Sample line losses:")
        # for (line, r, ℓ, loss) in sample_line_losses
        #     println("      Line $line: r=$(round(r, digits=8)) pu, ℓ=$(round(ℓ, digits=6)) pu, loss=$(round(loss, digits=6)) kW")
        # end
        # println("    Balance (P_Subs - Load - Losses): $(round(total_P_Subs - total_load - total_losses, digits=6)) kW")
        # 
        # # Topology connectivity check
        # println("\n  Topology Connectivity Check:")
        # println("    Total buses: $(length(Nset)) ($(length(Sset)) substations + $(length(Nm1set)) distribution)")
        # println("    Total lines: $(length(Lset)) ($(length(L1set)) POI + $(length(Lm1set)) distribution)")
        # 
        # # Check if each distribution bus is reachable from at least one substation (BFS)
        # reachable_buses = Set{Any}(Sset)
        # changed = true
        # while changed
        #     changed = false
        #     for line in Lset
        #         (bus1, bus2) = line
        #         if bus1 in reachable_buses && !(bus2 in reachable_buses)
        #             push!(reachable_buses, bus2)
        #             changed = true
        #         elseif bus2 in reachable_buses && !(bus1 in reachable_buses)
        #             push!(reachable_buses, bus1)
        #             changed = true
        #         end
        #     end
        # end
        # unreachable = setdiff(Set(Nset), reachable_buses)
        # if isempty(unreachable)
        #     println("    ✓ All buses are connected")
        # else
        #     println("    ✗ WARNING: $(length(unreachable)) unreachable buses: $unreachable")
        # end
        # 
        # # Branch flow verification for POI lines
        # println("\n  POI Line Power Flows (t=1):")
        # for (i, j) in L1set
        #     P_val = result[:P][(i, j)][t] * P_BASE
        #     losses_val = rdict_pu[(i, j)] * result[:ℓ][(i, j)][t] * P_BASE
        #     println("    Line ($i → $j): P = $(round(P_val, digits=1)) kW, Losses = $(round(losses_val, digits=3)) kW")
        # end
        # 
        # # Check power balance at first few distribution nodes
        # println("\n  Sample Distribution Node Power Balance (t=1):")
        # sample_nodes = Nm1set[1:min(5, length(Nm1set))]
        # for j in sample_nodes
        #     # Find parent lines (incoming power)
        #     parent_lines_j = [(i_bus, j_bus) for (i_bus, j_bus) in Lset 
        #                       if j_bus == j && (i_bus, j_bus) in keys(rdict_pu) && i_bus == parent[j]]
        #     
        #     # Find child lines (outgoing power)
        #     child_lines_j = [(j, k) for k in children[j] if (j, k) in keys(rdict_pu)]
        #     
        #     # Calculate incoming power (with empty check)
        #     P_in = isempty(parent_lines_j) ? 0.0 : sum(result[:P][line][t] - rdict_pu[line] * result[:ℓ][line][t] for line in parent_lines_j)
        #     
        #     # Calculate outgoing power (with empty check)
        #     P_out = isempty(child_lines_j) ? 0.0 : sum(result[:P][line][t] for line in child_lines_j)
        #     
        #     # Load at this node
        #     P_load = p_L_pu[j, t]
        #     
        #     # Balance
        #     balance = P_in - P_out - P_load
        #     
        #     println("    Node $j: P_in=$(round(P_in*P_BASE, digits=2)) kW, P_out=$(round(P_out*P_BASE, digits=2)) kW, Load=$(round(P_load*P_BASE, digits=2)) kW, Balance=$(round(balance*P_BASE, digits=6)) kW")
        # end
    end
    
    println("="^80)
    
    return result
end

# ============================================================================
# ANGLE COMPUTATION FOR MULTI-POI SYSTEM
# ============================================================================

"""
    compute_angles_from_power_flow(data, result, slack_sub)

Compute bus voltage angles from power flow solution using incidence matrix approach.
Uses linearized power flow: θ = B^T * (B*B^T)^-1 * (x.*P - r.*Q) ./ V

Arguments:
- data: system data dictionary
- result: OPF solution result (with :P_flow, :Q_flow, :v_sq organized by time)
- slack_sub: slack substation (e.g., "1s", "2s", etc.)

Returns Dict with angle statistics.
"""
function compute_angles_from_power_flow(data, result, slack_sub)
    T = data[:T]
    Tset = 1:T
    Nset = data[:Nset]
    Lset = data[:Lset]  # Oriented lines (132 lines)
    
    # Create bus index mapping (133 buses total)
    bus_list = sort(collect(Nset), by = x -> x isa String ? (0, x) : (1, x))
    bus_to_idx = Dict(bus => i for (i, bus) in enumerate(bus_list))
    idx_to_bus = Dict(i => bus for (i, bus) in enumerate(bus_list))
    
    n_buses = length(Nset)  # 133
    n_lines = length(Lset)  # 132
    
    # Build full incidence matrix C (n_lines × n_buses)
    # For each line (i,j): C[line, i] = 1, C[line, j] = -1
    C = zeros(n_lines, n_buses)
    
    for (line_idx, (bus_i, bus_j)) in enumerate(Lset)
        i = bus_to_idx[bus_i]
        j = bus_to_idx[bus_j]
        C[line_idx, i] = 1.0
        C[line_idx, j] = -1.0
    end
    
    # Remove slack bus column to get reduced incidence matrix B
    slack_idx = bus_to_idx[slack_sub]
    non_slack_indices = setdiff(1:n_buses, slack_idx)
    B = C[:, non_slack_indices]  # (n_lines × (n_buses-1))
    
    # Create mapping for non-slack buses
    non_slack_buses = [idx_to_bus[i] for i in non_slack_indices]
    non_slack_to_reduced_idx = Dict(bus => i for (i, bus) in enumerate(non_slack_buses))
    
    # Get line impedances and create vectors
    line_impedances = data[:line_impedances]
    r_vec = zeros(n_lines)
    x_vec = zeros(n_lines)
    
    for (line_idx, line_tuple) in enumerate(Lset)
        # Try both orientations since Lset is oriented but impedances stored with original orientation
        if haskey(line_impedances, line_tuple)
            r_vec[line_idx], x_vec[line_idx] = line_impedances[line_tuple]
        else
            # Try reversed
            line_tuple_rev = (line_tuple[2], line_tuple[1])
            if haskey(line_impedances, line_tuple_rev)
                r_vec[line_idx], x_vec[line_idx] = line_impedances[line_tuple_rev]
            else
                error("Could not find impedance for line $line_tuple or $line_tuple_rev")
            end
        end
    end
    
    # Precompute (B * B^T)^-1 * B for angle computation
    # θ = B^T * (B*B^T)^-1 * (x.*P - r.*Q)
    BBT = B * B'
    BBT_inv = inv(BBT)
    angle_matrix = B' * BBT_inv  # (n_buses-1) × n_lines
    
    # Storage for angles across time
    theta_matrix = Dict{Any, Vector{Float64}}()  # bus => [angle at each t]
    
    # Initialize all buses
    for bus in Nset
        theta_matrix[bus] = zeros(T)
    end
    
    # For each time period, compute angles using incidence matrix
    for t in Tset
        # Extract power flows and voltages at time t
        P_flow = result[:P_flow][t]  # Dict: line_tuple => power flow (pu)
        Q_flow = result[:Q_flow][t]  # Dict: line_tuple => power flow (pu)
        v_sq = result[:v_sq][t]      # Dict: bus => v^2 (pu)
        
        # Build power flow vectors (in line order)
        P_vec = zeros(n_lines)
        Q_vec = zeros(n_lines)
        
        for (line_idx, line_tuple) in enumerate(Lset)
            P_vec[line_idx] = P_flow[line_tuple]
            Q_vec[line_idx] = Q_flow[line_tuple]
        end
        
        # Compute angle differences using linearized power flow
        # δθ = (B^T * (B*B^T)^-1) * (x.*P - r.*Q)
        angle_source = x_vec .* P_vec - r_vec .* Q_vec
        theta_reduced = angle_matrix * angle_source  # Angles for non-slack buses
        
        # Store angles (slack is zero, others from computation)
        for (i, bus) in enumerate(non_slack_buses)
            theta_matrix[bus][t] = theta_reduced[i]
        end
        # Slack bus angle remains zero (already initialized)
    end
    
    # Convert to degrees
    theta_deg = Dict{Any, Vector{Float64}}()
    for (bus, angles_rad) in theta_matrix
        theta_deg[bus] = rad2deg.(angles_rad)
    end
    
    # Compute statistics for all non-slack buses
    Sset = data[:Sset]
    median_angles_deg = Dict{Any, Float64}()
    angle_ranges_deg = Dict{Any, Vector{Float64}}()
    rmse_values_deg = Dict{Any, Float64}()
    
    # Compute stats for all substations except slack
    for sub in Sset
        if sub == slack_sub
            continue  # Skip slack bus (angle is zero by definition)
        end
        
        angles = theta_deg[sub]
        med = median(angles)
        median_angles_deg[sub] = med
        angle_ranges_deg[sub] = [minimum(angles), maximum(angles)]
        rmse_values_deg[sub] = sqrt(mean((angles .- med).^2))
    end
    
    # Total deviation (RMS of all substation RMSEs)
    if length(rmse_values_deg) > 0
        total_deviation_deg = sqrt(sum(v^2 for v in values(rmse_values_deg)) / length(rmse_values_deg))
    else
        total_deviation_deg = 0.0
    end
    
    return Dict(
        :theta_deg => theta_deg,
        :theta_rad => theta_matrix,
        :median_angles_deg => median_angles_deg,
        :angle_ranges_deg => angle_ranges_deg,
        :rmse_values_deg => rmse_values_deg,
        :total_deviation_deg => total_deviation_deg
    )
end

# ============================================================================
# SOLVE MULTI-POI MPOPF FOR ALL SLACK CONFIGURATIONS
# ============================================================================

println("\n\n")
println("╔"*"═"^78*"╗")
println("║"*" "^15*"MULTI-SLACK MPOPF ANALYSIS FOR IEEE 123-BUS 5-POI"*" "^13*"║")
println("╚"*"═"^78*"╝")

# Create results directory
system_folder_name = "$(data[:systemName])_T$(data[:T])"
results_base_dir = joinpath(@__DIR__, "envs", "multi_poi", "processedData", system_folder_name)
plots_dir = joinpath(results_base_dir, "plots")
mkpath(plots_dir)
println("\nResults will be saved to: $results_base_dir")

# Solve for all 5 slack substation configurations
slack_substations = ["1s", "2s", "3s", "4s", "5s"]
results_by_slack = Dict{String, Any}()

println("\n" * "="^80)
println(Crayon(foreground = :cyan, bold = true)("SOLVING MPOPF FOR ALL FIVE SLACK CONFIGURATIONS"))
println("="^80)

for slack_sub in slack_substations
    println("\n" * "-"^80)
    println(Crayon(foreground = :yellow, bold = true)("Solving MPOPF with Substation $slack_sub as slack bus..."))
    println("-"^80)
    
    result = solve_multi_poi_mpopf(data; slack_substation=slack_sub, solver=:gurobi)
    
    if result[:status] == MOI.OPTIMAL || result[:status] == MOI.LOCALLY_SOLVED
        # Reorganize data for angle computation
        # Convert Dict{line => [t values]} to [t => Dict{line => value}]
        T = data[:T]
        P_flow = [Dict(line => result[:P][line][t] for line in data[:Lset]) for t in 1:T]
        Q_flow = [Dict(line => result[:Q][line][t] for line in data[:Lset]) for t in 1:T]
        v_sq = [Dict(bus => result[:v][bus][t] for bus in data[:Nset]) for t in 1:T]
        
        # Package for angle computation
        angle_input = Dict(
            :P_flow => P_flow,
            :Q_flow => Q_flow,
            :v_sq => v_sq
        )
        
        # Compute angles
        println("  Computing voltage angles...")
        angle_result = compute_angles_from_power_flow(data, angle_input, slack_sub)
        
        # Add angle results to main result
        result[:angles] = angle_result
        
        # Compute total losses
        P_loss_total_pu = 0.0
        Q_loss_total_pu = 0.0
        for line in data[:Lset]
            r_pu = data[:rdict_pu][line]
            x_pu = data[:xdict_pu][line]
            for t in data[:Tset]
                P_loss_total_pu += r_pu * result[:ℓ][line][t]
                Q_loss_total_pu += x_pu * result[:ℓ][line][t]
            end
        end
        result[:P_loss_pu] = P_loss_total_pu
        result[:Q_loss_pu] = Q_loss_total_pu
        result[:P_loss_kW] = P_loss_total_pu * data[:kVA_B]
        result[:Q_loss_kVAr] = Q_loss_total_pu * data[:kVA_B]
        
        results_by_slack[slack_sub] = result
        
        # Print summary
        println(Crayon(foreground = :green)("✓ Converged"))
        println("  Objective: \$$(round(result[:objective], digits=2))")
        println("  Solve time: $(round(result[:solve_time], digits=2)) seconds")
        
        # Calculate total power dispatch (time-averaged)
        total_P_kW = sum(mean(result[:P_Subs][s]) for s in data[:Sset]) * data[:kVA_B]
        println("  Total power: $(round(total_P_kW, digits=1)) kW (time-averaged)")
        
        # Print angle statistics
        println("  Total angle deviation: $(round(angle_result[:total_deviation_deg], digits=3))°")
    else
        @warn "MPOPF failed for slack=$slack_sub with status $(result[:status])"
        results_by_slack[slack_sub] = result
    end
end

# ==================================================================================
# PRINT COMPARISON TABLE
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :light_blue, bold = true)("COMPARISON OF SLACK BUS CONFIGURATIONS"))
println("="^80)

# Print header
header_crayon = Crayon(foreground = :white, bold = true)
header = @sprintf("%-15s | %-15s | %-15s | %-12s | %-15s | %-15s", 
                 "Slack Bus", "Operational", "Avg Power", "Solve Time", "Angle Dev", "Total Losses")
header2 = @sprintf("%-15s | %-15s | %-15s | %-12s | %-15s | %-15s",
                  "(Subs)", "Cost (\$)", "(kW)", "(sec)", "(deg)", "(kW)")
separator = "-"^length(header)

println(separator)
println(header_crayon(header))
println(header_crayon(header2))
println(separator)

for slack_sub in slack_substations
    result = results_by_slack[slack_sub]
    
    if result[:status] == MOI.OPTIMAL || result[:status] == MOI.LOCALLY_SOLVED
        total_cost = result[:objective]
        solve_time = result[:solve_time]
        
        # Calculate average total power dispatch
        total_P_kW = sum(mean(result[:P_Subs][s]) for s in data[:Sset]) * data[:kVA_B]
        
        # Get angle deviation
        angle_dev = result[:angles][:total_deviation_deg]
        
        # Get losses
        P_loss_kW = result[:P_loss_kW]
        
        # Format with colors
        slack_idx = parse(Int, slack_sub[1:end-1])
        slack_colors = [:green, :cyan, :light_blue, :light_magenta, :light_yellow]
        color = slack_colors[slack_idx]
        sub_str = Crayon(foreground = color, bold = true)(@sprintf("%-15s", slack_sub))
        
        cost_str = @sprintf("%-15.2f", total_cost)
        power_str = @sprintf("%-15.1f", total_P_kW)
        time_str = @sprintf("%-12.2f", solve_time)
        angle_str = @sprintf("%-15.3f", angle_dev)
        loss_str = @sprintf("%-15.1f", P_loss_kW)
        
        @printf("%s | %s | %s | %s | %s | %s\n", sub_str, cost_str, power_str, time_str, angle_str, loss_str)
    else
        sub_str = Crayon(foreground = :red)(@sprintf("%-15s", slack_sub))
        @printf("%s | %-15s | %-15s | %-12s | %-15s | %-15s\n", sub_str, "FAILED", "N/A", "N/A", "N/A", "N/A")
    end
end

println(separator)

# Print loss statistics summary
println("\n" * "="^80)
println(Crayon(foreground = :light_green, bold = true)("LOSS STATISTICS ACROSS ALL CONFIGURATIONS"))
println("="^80)

loss_values = [results_by_slack[s][:P_loss_kW] for s in slack_substations if haskey(results_by_slack[s], :P_loss_kW)]
if !isempty(loss_values)
    loss_avg = mean(loss_values)
    loss_min = minimum(loss_values)
    loss_max = maximum(loss_values)
    loss_std = std(loss_values)
    
    println("  Average losses: $(round(loss_avg, digits=1)) kW")
    println("  Range: [$(round(loss_min, digits=1)), $(round(loss_max, digits=1))] kW")
    println("  Std deviation: $(round(loss_std, digits=1)) kW")
    println("  Loss percentage (avg): $(round(100 * loss_avg / (sum(sum(results_by_slack["1s"][:P_Subs][s]) for s in data[:Sset]) * data[:kVA_B]), digits=2))%")
    println("\n  Note: Thermal limit constraint added: I_max = 5.0 pu → ℓ_max = 25.0 pu")
end

println(separator)
println("="^80)

# ==================================================================================
# DETAILED RESULTS FOR EACH CONFIGURATION
# ==================================================================================

for slack_sub in slack_substations
    result = results_by_slack[slack_sub]
    
    if !(result[:status] == MOI.OPTIMAL || result[:status] == MOI.LOCALLY_SOLVED)
        continue
    end
    
    println("\n" * "="^80)
    title_crayon = Crayon(foreground = :light_magenta, bold = true)
    println(title_crayon("DETAILED RESULTS: SLACK = $slack_sub"))
    println("="^80)
    
    # Power dispatch summary
    header_crayon = Crayon(foreground = :white, bold = true)
    println(header_crayon("\nSubstation Power Dispatch (total over horizon):"))
    for s in sort(collect(data[:Sset]))
        P_total_kW = sum(result[:P_Subs][s]) * data[:kVA_B]
        Q_total_kVAr = sum(result[:Q_Subs][s]) * data[:kVA_B]
        P_min_kW = minimum(result[:P_Subs][s]) * data[:kVA_B]
        P_max_kW = maximum(result[:P_Subs][s]) * data[:kVA_B]
        Q_min_kVAr = minimum(result[:Q_Subs][s]) * data[:kVA_B]
        Q_max_kVAr = maximum(result[:Q_Subs][s]) * data[:kVA_B]
        
        if s == slack_sub
            slack_marker = Crayon(foreground = :green, bold = true)(" (SLACK)")
            print("  $s: P = $(round(P_total_kW, digits=1)) kW [$(round(P_min_kW, digits=1)), $(round(P_max_kW, digits=1))], Q = $(round(Q_total_kVAr, digits=1)) kVAr [$(round(Q_min_kVAr, digits=1)), $(round(Q_max_kVAr, digits=1))]")
            println(slack_marker)
        else
            println("  $s: P = $(round(P_total_kW, digits=1)) kW [$(round(P_min_kW, digits=1)), $(round(P_max_kW, digits=1))], Q = $(round(Q_total_kVAr, digits=1)) kVAr [$(round(Q_min_kVAr, digits=1)), $(round(Q_max_kVAr, digits=1))]")
        end
    end
    
    # Distribution bus voltage summary
    println(header_crayon("\nVoltage Statistics (distribution buses, time-averaged):"))
    v_vals = [mean(result[:v][bus]) for bus in data[:Nm1set]]
    V_avg = mean(sqrt.(v_vals))
    V_min = minimum(sqrt.(v_vals))
    V_max = maximum(sqrt.(v_vals))
    println("  Average: $(round(V_avg, digits=4)) pu")
    println("  Range: [$(round(V_min, digits=4)), $(round(V_max, digits=4))] pu")
    
    # Line losses - computed from ℓ values (total over horizon)
    # Real power loss on each line: P_loss = r * ℓ (in pu), then convert to kW
    P_loss_total_pu = 0.0
    Q_loss_total_pu = 0.0
    
    # Create 24x1 vector for losses per time period (for debugging)
    P_loss_per_period_pu = zeros(length(data[:Tset]))
    
    for line in data[:Lset]
        r_pu = data[:rdict_pu][line]
        x_pu = data[:xdict_pu][line]
        # Sum over all time periods
        for t in data[:Tset]
            loss_t = r_pu * result[:ℓ][line][t]
            P_loss_total_pu += loss_t
            P_loss_per_period_pu[t] += loss_t
            Q_loss_total_pu += x_pu * result[:ℓ][line][t]
        end
    end
    P_loss_kW = P_loss_total_pu * data[:kVA_B]
    Q_loss_kVAr = Q_loss_total_pu * data[:kVA_B]
    
    println(header_crayon("\nLine Losses (total over horizon):"))
    println("  Σ(r*ℓ) = $(round(P_loss_total_pu, digits=6)) pu")
    println("  Real: $(round(P_loss_kW, digits=1)) kW")
    println("  Reactive: $(round(Q_loss_kVAr, digits=1)) kVAr")
    
    # Print 24x1 vector for slack=1s
    if slack_sub == "1s"
        println(header_crayon("\nLosses per time period (pu) - Σ_lines(r_ij*ℓ_ij) for each t:"))
        for t in data[:Tset]
            println("  t=$t: $(round(P_loss_per_period_pu[t], digits=6)) pu")
        end
        println("  Sum = $(round(sum(P_loss_per_period_pu), digits=6)) pu")
    end
    
    # Power balance verification: Σ P_Subs = Σ p_L + Σ (r*ℓ) - Σ p_D - Σ P_B
    total_P_Subs_pu = sum(sum(result[:P_Subs][s]) for s in data[:Sset])
    total_P_Subs_kW = total_P_Subs_pu * data[:kVA_B]
    
    # Compute total loads (currently allocated to buses)
    total_p_L_pu = sum(sum(data[:p_L_pu][j, t] for t in data[:Tset]) for j in data[:NLset])
    total_p_L_kW = total_p_L_pu * data[:kVA_B]
    
    # Compute total PV generation (zero for now, but ready for future)
    total_p_D_pu = isempty(data[:Dset]) ? 0.0 : sum(sum(data[:p_D_pu][j, t] for t in data[:Tset]) for j in data[:Dset])
    total_p_D_kW = total_p_D_pu * data[:kVA_B]
    
    # Compute total battery power (zero for now, but ready for future)
    total_P_B_pu = 0.0  # isempty(data[:Bset]) ? 0.0 : sum(sum(result[:P_B][j][t] for t in data[:Tset]) for j in data[:Bset]) - Not stored in result yet
    total_P_B_kW = total_P_B_pu * data[:kVA_B]
    
    println(header_crayon("\nPower Balance Check (total over horizon):"))
    println("  Σ P_Subs = $(round(total_P_Subs_pu, digits=6)) pu = $(round(total_P_Subs_kW, digits=1)) kW")
    println("  Σ p_L (loads) = $(round(total_p_L_pu, digits=6)) pu = $(round(total_p_L_kW, digits=1)) kW")
    println("  Σ (r*ℓ) (losses) = $(round(P_loss_total_pu, digits=6)) pu = $(round(P_loss_kW, digits=1)) kW")
    println("  Σ p_D (PV gen) = $(round(total_p_D_kW, digits=1)) kW")
    println("  Σ P_B (battery) = $(round(total_P_B_kW, digits=1)) kW")
    println("  Balance: Σ P_Subs - Σ p_L - Σ (r*ℓ) + Σ p_D + Σ P_B = $(round(total_P_Subs_kW - total_p_L_kW - P_loss_kW + total_p_D_kW + total_P_B_kW, digits=3)) kW (should ≈ 0)")
    
    if slack_sub == "1s"
        println("\n  DEBUG:")
        println("    total_P_Subs_pu calculation: sum over substations of sum over time")
        for s in data[:Sset]
            sum_s = sum(result[:P_Subs][s])
            println("      Substation $s: sum = $(round(sum_s, digits=6)) pu")
        end
        println("    Total = $(round(total_P_Subs_pu, digits=6)) pu")
        println("    Should equal: total_p_L_pu + P_loss_total_pu = $(round(total_p_L_pu, digits=6)) + $(round(P_loss_total_pu, digits=6)) = $(round(total_p_L_pu + P_loss_total_pu, digits=6)) pu")
        
        # Check a few line losses in detail
        println("\n  Checking individual line losses for t=9 (peak loss period):")
        t = 9
        line_count = 0
        total_loss_t9 = 0.0
        max_ℓ_line = nothing
        max_ℓ_val = 0.0
        
        for line in data[:Lset]
            r_pu = data[:rdict_pu][line]
            ℓ_val = result[:ℓ][line][t]
            loss = r_pu * ℓ_val
            total_loss_t9 += loss
            
            if ℓ_val > max_ℓ_val
                max_ℓ_val = ℓ_val
                max_ℓ_line = line
            end
            
            if line_count < 10  # Print first 10 lines
                P_val = result[:P][line][t]
                Q_val = result[:Q][line][t]
                i, j = line
                v_i = result[:v][i][t]
                println("      Line $line: r=$(round(r_pu, sigdigits=4)) pu, ℓ=$(round(ℓ_val, sigdigits=4)) pu, loss=$(round(loss, sigdigits=4)) pu")
                println("        P=$(round(P_val, sigdigits=4)), Q=$(round(Q_val, sigdigits=4)), v_i=$(round(v_i, sigdigits=4))")
                println("        Check: P²+Q² = $(round(P_val^2 + Q_val^2, sigdigits=4)), v*ℓ = $(round(v_i * ℓ_val, sigdigits=4))")
                line_count += 1
            end
        end
        println("      Total loss at t=9: $(round(total_loss_t9, digits=6)) pu")
        println("      (Should match P_loss_per_period_pu[9] = $(round(P_loss_per_period_pu[9], digits=6)) pu)")
        
        if max_ℓ_line !== nothing
            i, j = max_ℓ_line
            P_val = result[:P][max_ℓ_line][t]
            Q_val = result[:Q][max_ℓ_line][t]
            v_i = result[:v][i][t]
            r_pu = data[:rdict_pu][max_ℓ_line]
            println("\n      Line with MAX ℓ at t=9: $max_ℓ_line")
            println("        ℓ = $(round(max_ℓ_val, sigdigits=6)) pu (current squared)")
            println("        I = $(round(sqrt(max_ℓ_val), sigdigits=6)) pu (current magnitude)")
            println("        P=$(round(P_val, sigdigits=4)), Q=$(round(Q_val, sigdigits=4)), v_i=$(round(v_i, sigdigits=4))")
            println("        r=$(round(r_pu, sigdigits=4)) pu, loss=$(round(r_pu * max_ℓ_val, sigdigits=6)) pu")
        end
        
        # Now show ℓ values across all time for the max line
        println("\n  ℓ values for line $max_ℓ_line across all 24 periods:")
        for t in data[:Tset]
            ℓ_val = result[:ℓ][max_ℓ_line][t]
            println("      t=$t: ℓ=$(round(ℓ_val, sigdigits=6)) pu")
        end
    end
    
    # # DETAILED NODAL BALANCE CHECK (only for slack=1s for debugging)
    # # Commented out - power balance verified as correct
    # if slack_sub == "1s"
    #     println(header_crayon("\n" * "="^80))
    #     println(header_crayon("DETAILED NODAL POWER BALANCE CHECK (pu)"))
    #     println(header_crayon("="^80))
    #     println("NRPB equation: Σ P_jk - Σ(P_ij - r_ij*ℓ_ij) = p_j where p_j = P_B_j + p_D_j - p_L_j")
    #     println("For each distribution bus j:")
    #     println("  Σ P_jk       = Total outflow from bus j")
    #     println("  Σ P_ij       = Total inflow to bus j")
    #     println("  Σ (r_ij*ℓ_ij) = Total real power losses on incoming lines")
    #     println("  p_j          = Net generation (P_B + p_D - p_L)")
    #     println("-"^80)
    #     
    #     # Define color scheme
    #     outflow_color = Crayon(foreground = :light_blue, bold = true)
    #     inflow_color = Crayon(foreground = :light_cyan, bold = true)
    #     loss_color = Crayon(foreground = :light_yellow, bold = true)
    #     gen_color = Crayon(foreground = :light_green, bold = true)
    #     balance_color_good = Crayon(foreground = :green, bold = true)
    #     balance_color_bad = Crayon(foreground = :red, bold = true)
    #     
    #     # Get distribution buses (exclude substation buses)
    #     dist_buses = sort(collect(data[:Nm1set]))
    #     
    #     max_imbalance = 0.0
    #     total_imbalance = 0.0
    #     
    #     for j in dist_buses
    #         # Sum over all time periods
    #         sum_Pjk = 0.0      # Outflow
    #         sum_Pij = 0.0      # Inflow
    #         sum_rij_lij = 0.0  # Losses
    #         p_j = 0.0          # Net generation
    #         
    #         for t in data[:Tset]
    #             # Outgoing power (to children): Σ P_jk
    #             if !isempty(data[:children][j])
    #                 for k in data[:children][j]
    #                     child_line = (j, k)
    #                     if child_line in keys(data[:rdict_pu])
    #                         sum_Pjk += result[:P][child_line][t]
    #                     end
    #                 end
    #             end
    #             
    #             # Incoming power and losses (from parent lines)
    #             parent_j = data[:parent][j]
    #             parent_line = (parent_j, j)
    #             if parent_line in keys(data[:rdict_pu])
    #                 P_ij = result[:P][parent_line][t]
    #                 r_ij = data[:rdict_pu][parent_line]
    #                 ℓ_ij = result[:ℓ][parent_line][t]
    #                 
    #                 sum_Pij += P_ij
    #                 sum_rij_lij += r_ij * ℓ_ij
    #             end
    #             
    #             # Net generation: P_B_j + p_D_j - p_L_j
    #             p_L_j_t = (j in data[:NLset]) ? data[:p_L_pu][j, t] : 0.0
    #             p_D_j_t = (j in data[:Dset]) ? data[:p_D_pu][j, t] : 0.0
    #             P_B_j_t = 0.0  # (j in data[:Bset]) ? result[:P_B][j][t] : 0.0
    #             
    #             p_j += P_B_j_t + p_D_j_t - p_L_j_t
    #         end
    #         
    #         # Verify: Σ(r_ij*ℓ_ij) should equal Σ P_ij - Σ P_jk + p_j
    #         # From NRPB: Σ P_jk - Σ P_ij + Σ(r_ij*ℓ_ij) = p_j
    #         # Therefore: Σ(r_ij*ℓ_ij) = Σ P_ij - Σ P_jk + p_j
    #         sum_rij_lij_computed = sum_Pij - sum_Pjk + p_j
    #         
    #         # Balance check: Σ P_jk - Σ(P_ij - r_ij*ℓ_ij) - p_j should ≈ 0
    #         # Which is: Σ P_jk - Σ P_ij + Σ(r_ij*ℓ_ij) - p_j
    #         balance = sum_Pjk - sum_Pij + sum_rij_lij - p_j
    #         abs_balance = abs(balance)
    #         
    #         max_imbalance = max(max_imbalance, abs_balance)
    #         total_imbalance += abs_balance
    #         
    #         # Color code the balance
    #         balance_crayon = abs_balance < 1e-6 ? balance_color_good : balance_color_bad
    #         
    #         # Print: Σ P_jk - p_j, Σ P_ij, Σ(r_ij*ℓ_ij) [where third = first - second]
    #         print("Bus $(rpad(j, 4)): ")
    #         print(outflow_color("Σ P_jk-p_j = $(rpad(round(sum_Pjk - p_j, digits=6), 11))"))
    #         print(" | ")
    #         print(inflow_color("Σ P_ij = $(rpad(round(sum_Pij, digits=6), 11))"))
    #         print(" | ")
    #         print(loss_color("Σ r*ℓ = $(rpad(round(sum_rij_lij, digits=6), 11))"))
    #         print(" (comp: $(rpad(round(sum_rij_lij_computed, digits=6), 11)))")
    #         print(" | ")
    #         println(balance_crayon("Bal = $(round(balance, sigdigits=3))"))
    #     end
    #     
    #     println("-"^80)
    #     println("Summary:")
    #     println("  Max absolute imbalance: $(round(max_imbalance, sigdigits=4)) pu")
    #     println("  Total absolute imbalance: $(round(total_imbalance, sigdigits=4)) pu")
    #     println(header_crayon("="^80))
    # end
    
    # Angle statistics
    println(header_crayon("\nVoltage Angle Statistics (substations):"))
    angle_results = result[:angles]
    for s in sort(collect(data[:Sset]))
        if s == slack_sub
            slack_marker = Crayon(foreground = :green, bold = true)(" (SLACK)")
            print("  $s: θ = 0.0° (slack reference)")
            println(slack_marker)
        else
            if haskey(angle_results[:median_angles_deg], s)
                median_angle = angle_results[:median_angles_deg][s]
                angle_range = angle_results[:angle_ranges_deg][s]
                rmse_angle = angle_results[:rmse_values_deg][s]
                println("  $s: θ_med = $(round(median_angle, digits=3))°, range = [$(round(angle_range[1], digits=3))°, $(round(angle_range[2], digits=3))°], RMSE = $(round(rmse_angle, digits=3))°")
            end
        end
    end
    println("  Total angle deviation (RMS): $(round(angle_results[:total_deviation_deg], digits=3))°")
    
    # Voltage magnitude statistics for substations
    println(header_crayon("\nVoltage Magnitude Statistics (substations):"))
    for s in sort(collect(data[:Sset]))
        V_vals = sqrt.(result[:v][s])
        V_avg = mean(V_vals)
        V_min = minimum(V_vals)
        V_max = maximum(V_vals)
        
        if s == slack_sub
            slack_marker = Crayon(foreground = :green, bold = true)(" (SLACK)")
            print("  $s: V_avg = $(round(V_avg, digits=4)) pu, range = [$(round(V_min, digits=4)), $(round(V_max, digits=4))] pu")
            println(slack_marker)
        else
            println("  $s: V_avg = $(round(V_avg, digits=4)) pu, range = [$(round(V_min, digits=4)), $(round(V_max, digits=4))] pu")
        end
    end
    
    println("\nObjective Value: \$$(round(result[:objective], digits=2))")
    println("Solve Time: $(round(result[:solve_time], digits=2)) seconds")
end

println("\n" * "="^80)
success_crayon = Crayon(foreground = :light_green, bold = true)
info_crayon = Crayon(foreground = :light_blue)
println(success_crayon("✓ MULTI-SLACK MPOPF ANALYSIS COMPLETE!"))
println(info_crayon("ℹ Analyzed all $(length(slack_substations)) slack bus configurations."))
println("ℹ Results saved to: $results_base_dir")
println("="^80)

# ==================================================================================
# GENERATE COMPARISON PLOTS
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :cyan, bold = true)("GENERATING COMPARISON PLOTS"))
println("="^80)

# Define colors and markers for all plots (once, outside the loop)
slack_colors = [RGB(0x00/255, 0x72/255, 0xB2/255),    # #0072B2 deep blue
                RGB(0xE6/255, 0x9F/255, 0x00/255),    # #E69F00 amber
                RGB(0x00/255, 0x9E/255, 0x73/255),    # #009E73 teal
                RGB(0xCC/255, 0x79/255, 0xA7/255),    # #CC79A7 magenta
                RGB(0xD5/255, 0x5E/255, 0x00/255)]    # #D55E00 rust red
marker_shapes = [:circle, :square, :diamond, :utriangle, :star5]

# Extract data for all successful configurations
Tset = 1:data[:T]
substations = sort(collect(data[:Sset]))
kVA_B = data[:kVA_B]

# For each slack configuration, extract time series data
for slack_sub in slack_substations
    result = results_by_slack[slack_sub]
    
    if !(result[:status] == MOI.OPTIMAL || result[:status] == MOI.LOCALLY_SOLVED)
        continue
    end
    
    # Prepare data for stacked bar plots
    # P_Subs: Real power dispatch (kW) for each substation at each time step
    P_matrix = zeros(length(substations), data[:T])
    Q_matrix = zeros(length(substations), data[:T])
    Cost_matrix = zeros(length(substations), data[:T])
    
    for (i, s) in enumerate(substations)
        P_kW = result[:P_Subs][s] .* kVA_B
        Q_kVAr = result[:Q_Subs][s] .* kVA_B
        
        P_matrix[i, :] = P_kW
        Q_matrix[i, :] = Q_kVAr
        
        # Cost calculation: use substation-specific cost profile
        # Extract substation index (e.g., "1s" -> 1)
        s_idx = parse(Int, string(s)[1:end-1])
        energy_rate_per_t = data[:LoadShapeCost_dict][s_idx]  # $/kWh for this substation
        Δt = data[:delta_t_h]  # hours
        Cost_matrix[i, :] = P_kW .* energy_rate_per_t .* Δt
    end
    
    # Smart x-tick spacing for clarity
    xtick_sparse = if T <= 24
        collect(1:2:T)
    else
        step = max(2, div(T, 8))
        collect(1:step:T)
    end
    
    # Plot 1: Line plot of P_Subs (Real Power) for all substations
    p1 = plot(xlabel = L"Time Period $t$ (hour)", 
              ylabel = L"Real Power $P_{\mathrm{Subs}}$ (kW)",
              title = "Active Power Dispatch (Slack: Subs $slack_sub)",
              legend = :topright,
              legend_columns = 2,
              legendfontsize = 8,
              legend_background_color = RGBA(1,1,1,0.7),
              size = (710, 460),
              theme = :dao,
              dpi = 400,
              tickfont = font(9, "Computer Modern"),
              guidefont = font(10, "Computer Modern"),
              titlefont = font(11, "Computer Modern"),
              xticks = (xtick_sparse, string.(xtick_sparse)),
              grid = true,
              gridlinewidth = 0.5,
              gridalpha = 0.25,
              gridstyle = :solid,
              framestyle = :box,
              left_margin = 3Plots.mm,
              right_margin = 3Plots.mm,
              top_margin = 2Plots.mm,
              bottom_margin = 4Plots.mm)
    for (i, s) in enumerate(substations)
        plot!(p1, Tset, P_matrix[i, :], 
              label = L"Subs %$i",
              color = slack_colors[i],
              markershape = marker_shapes[i],
              markersize = 5,
              markerstrokewidth = 1.5,
              markerstrokecolor = :black,
              markeralpha = 0.9,
              linewidth = 2.2)
    end
    
    # Plot 2: Line plot of Q_Subs (Reactive Power) for all substations
    p2 = plot(xlabel = L"Time Period $t$ (hour)", 
              ylabel = L"Reactive Power $Q_{\mathrm{Subs}}$ (kVAr)",
              title = "Reactive Power Dispatch (Slack: Subs $slack_sub)",
              legend = :topright,
              legend_columns = 2,
              legendfontsize = 8,
              legend_background_color = RGBA(1,1,1,0.7),
              size = (710, 460),
              theme = :dao,
              dpi = 400,
              tickfont = font(9, "Computer Modern"),
              guidefont = font(10, "Computer Modern"),
              titlefont = font(11, "Computer Modern"),
              xticks = (xtick_sparse, string.(xtick_sparse)),
              grid = true,
              gridlinewidth = 0.5,
              gridalpha = 0.25,
              gridstyle = :solid,
              framestyle = :box,
              left_margin = 3Plots.mm,
              right_margin = 3Plots.mm,
              top_margin = 2Plots.mm,
              bottom_margin = 4Plots.mm)
    for (i, s) in enumerate(substations)
        plot!(p2, Tset, Q_matrix[i, :], 
              label = L"Subs %$i",
              color = slack_colors[i],
              markershape = marker_shapes[i],
              markersize = 5,
              markerstrokewidth = 1.5,
              markerstrokecolor = :black,
              markeralpha = 0.9,
              linewidth = 2.2)
    end
    
    # Plot 3: Line plot of Cost for all substations
    p3 = plot(xlabel = L"Time Period $t$ (hour)", 
              ylabel = L"Energy Cost $C$ (\$)",
              title = "Energy Cost per Period (Slack: Subs $slack_sub)",
              legend = :topright,
              legend_columns = 2,
              legendfontsize = 8,
              legend_background_color = RGBA(1,1,1,0.7),
              size = (710, 460),
              theme = :dao,
              dpi = 400,
              tickfont = font(9, "Computer Modern"),
              guidefont = font(10, "Computer Modern"),
              titlefont = font(11, "Computer Modern"),
              xticks = (xtick_sparse, string.(xtick_sparse)),
              grid = true,
              gridlinewidth = 0.5,
              gridalpha = 0.25,
              gridstyle = :solid,
              framestyle = :box,
              left_margin = 3Plots.mm,
              right_margin = 3Plots.mm,
              top_margin = 2Plots.mm,
              bottom_margin = 4Plots.mm)
    for (i, s) in enumerate(substations)
        plot!(p3, Tset, Cost_matrix[i, :], 
              label = L"Subs %$i",
              color = slack_colors[i],
              markershape = marker_shapes[i],
              markersize = 5,
              markerstrokewidth = 1.5,
              markerstrokecolor = :black,
              markeralpha = 0.9,
              linewidth = 2.2)
    end
    
    # Combined Plot: Voltage angles (top) and magnitudes (bottom) - matching small3poi style
    angle_matrix = zeros(length(substations), T)
    voltage_matrix = zeros(length(substations), T)
    median_angles = Dict{String, Float64}()
    
    for (i, s) in enumerate(substations)
        angle_matrix[i, :] = rad2deg.(result[:angles][:theta_rad][s])
        voltage_matrix[i, :] = sqrt.(result[:v][s])
        
        # Store median angles for non-slack substations
        if s != slack_sub && haskey(result[:angles][:median_angles_deg], s)
            median_angles[s] = result[:angles][:median_angles_deg][s]
        end
    end
    
    # Calculate y-axis limits for angles
    all_angles = vec(angle_matrix)
    θ_min = minimum(all_angles)
    θ_max = maximum(all_angles)
    θ_range = θ_max - θ_min
    θ_margin = max(0.5, 0.15 * θ_range)
    angle_ylims = (θ_min - θ_margin, θ_max + θ_margin)
    
    # Calculate y-axis limits for voltages
    all_voltages = vec(voltage_matrix)
    v_min = minimum(all_voltages)
    v_max = maximum(all_voltages)
    v_range = v_max - v_min
    v_margin = max(0.01, 0.1 * v_range)
    voltage_ylims = (v_min - v_margin, v_max + v_margin)
    
    # Lighter colors for voltage panel (bottom)
    voltage_colors = [RGB(0x66/255, 0xB2/255, 0xFF/255),  # sky blue
                      RGB(0xFF/255, 0xD5/255, 0x80/255),  # light amber
                      RGB(0x7F/255, 0xD1/255, 0xAE/255),  # mint green
                      RGB(0xE6/255, 0x9F/255, 0xCC/255),  # light pink
                      RGB(0xFF/255, 0xB3/255, 0x8C/255)]  # peach
    
    # TOP PANEL: Angles with median reference lines
    p_angle = plot(xlabel = "",
                   ylabel = L"Angle $\theta$ [°]",
                   title = "Substation Voltage Profile vs Time (Slack: Subs $slack_sub)",
                   legend = :outertop,
                   legend_columns = 5,
                   legendfontsize = 7,
                   legend_background_color = RGBA(1,1,1,0.9),
                   size = (710, 460),
                   theme = :dao,
                   dpi = 400,
                   tickfont = font(9, "Computer Modern"),
                   guidefont = font(10, "Computer Modern"),
                   titlefont = font(11, "Computer Modern"),
                   xticks = (xtick_sparse, string.(xtick_sparse)),
                   xformatter = _->"",
                   ylims = angle_ylims,
                   grid = true,
                   gridlinewidth = 0.5,
                   gridalpha = 0.25,
                   gridstyle = :solid,
                   framestyle = :box,
                   left_margin = 3Plots.mm,
                   right_margin = 3Plots.mm,
                   top_margin = 2Plots.mm,
                   bottom_margin = 0Plots.mm)
    
    # Plot median reference lines first (dashed, behind trajectories)
    for (i, s) in enumerate(substations)
        if haskey(median_angles, s)
            δ_med = median_angles[s]
            hline!(p_angle, [δ_med],
                   linestyle = :dash,
                   linewidth = 2.0,
                   linecolor = slack_colors[i],
                   alpha = 0.8,
                   label = L"\delta_{%$i} = %$(round(δ_med, digits=2))°")
        end
    end
    
    # Plot angle trajectories
    for (i, s) in enumerate(substations)
        plot!(p_angle, Tset, angle_matrix[i, :], 
              label = L"\theta^t_{%$i}",
              color = slack_colors[i],
              markershape = marker_shapes[i],
              markersize = 5,
              markerstrokewidth = 1.5,
              markerstrokecolor = :black,
              markeralpha = 0.9,
              linewidth = 2.2)
    end
    
    # BOTTOM PANEL: Voltages (lighter colors, dashed lines)
    p_voltage = plot(xlabel = L"Time Period $t$ (hour)",
                     ylabel = L"Voltage $V$ [pu]",
                     legend = false,
                     size = (710, 460),
                     theme = :dao,
                     dpi = 400,
                     tickfont = font(9, "Computer Modern"),
                     guidefont = font(10, "Computer Modern"),
                     titlefont = font(11, "Computer Modern"),
                     xticks = (xtick_sparse, string.(xtick_sparse)),
                     ylims = voltage_ylims,
                     grid = true,
                     gridlinewidth = 0.5,
                     gridalpha = 0.25,
                     gridstyle = :solid,
                     framestyle = :box,
                     left_margin = 3Plots.mm,
                     right_margin = 3Plots.mm,
                     top_margin = 0Plots.mm,
                     bottom_margin = 4Plots.mm)
    
    # Plot voltage trajectories (lighter colors, dashed, square markers)
    for (i, s) in enumerate(substations)
        plot!(p_voltage, Tset, voltage_matrix[i, :], 
              label = L"V_{\mathrm{Subs},%$i}",
              color = voltage_colors[i],
              markershape = :square,
              markersize = 4,
              markerstrokewidth = 1.5,
              markerstrokecolor = :black,
              markeralpha = 0.9,
              linewidth = 2.0,
              linestyle = :dash)
    end
    
    # Stack the two panels vertically
    p4 = plot(p_angle, p_voltage,
              layout = (2,1),
              size = (710, 540),
              plot_title = "",
              link = :x)
    
    # Save plots
    plots_dir = joinpath(results_base_dir, "plots")
    mkpath(plots_dir)
    
    savefig(p1, joinpath(plots_dir, "P_Subs_timeseries_slack_$(slack_sub).png"))
    savefig(p2, joinpath(plots_dir, "Q_Subs_timeseries_slack_$(slack_sub).png"))
    savefig(p3, joinpath(plots_dir, "Cost_timeseries_slack_$(slack_sub).png"))
    savefig(p4, joinpath(plots_dir, "angle_voltage_slack$(slack_sub).png"))
    
    println("  ✓ Saved plots for slack=$slack_sub")
end

println("\n✓ All plots saved to: $(joinpath(results_base_dir, "plots"))")
println("="^80)

# ==================================================================================
# PLOT INPUT CURVES (matching small3poi style)
# ==================================================================================
println("\n" * "="^80)
println(Crayon(foreground = :magenta, bold = true)("GENERATING INPUT CURVE PLOTS"))
println("="^80)

# Smart x-tick spacing for clarity
xtick_sparse = if T <= 24
    collect(1:2:T)
else
    step = max(2, div(T, 8))
    collect(1:step:T)
end

# Smart time step label
Δt = data[:delta_t_h]
dt_label = if Δt >= 1.0
    "Δt = $(Δt) h"
else
    dt_minutes = round(Int, Δt * 60)
    "Δt = $(dt_minutes) min"
end

# Color palette (matching substation plots)
cost_colors = slack_colors  # Reuse the same colors for consistency

# ==== TOP PANEL: LOAD PROFILES FOR ALL 5 ZONES ====
p_load = plot(
    dpi=400,
    ylabel="Normalized Load",
    legend=:topright,
    legend_columns=2,
    legendfontsize=8,
    legend_background_color=RGBA(1,1,1,0.7),
    framestyle=:box,
    ylims=(0.60, 1.05),
    xticks=(xtick_sparse, string.(xtick_sparse)),
    xformatter=_->"",  # Hide x-labels on top panel
    grid=true,
    gridlinewidth=0.5,
    gridalpha=0.25,
    gridstyle=:solid,
    titlefont=font(11, "Computer Modern"),
    guidefont=font(10, "Computer Modern"),
    tickfont=font(9, "Computer Modern"),
    left_margin=3Plots.mm,
    right_margin=3Plots.mm,
    top_margin=2Plots.mm,
    bottom_margin=0Plots.mm,
    title="Input Curves: Zone Load Profiles and Substation Costs ($(dt_label))"
)

# Plot each zone's load profile
zone_profiles = data[:zone_profiles]
for (i, sub) in enumerate(data[:Sset])
    plot!(
        p_load, Tset, zone_profiles[sub],
        color=slack_colors[i],
        lw=2.2,
        markershape=marker_shapes[i],
        markersize=5,
        markerstrokewidth=1.5,
        markerstrokecolor=:black,
        markeralpha=0.9,
        label=L"\lambda^t_{\mathrm{Zone}\,%$i}"
    )
end

# ==== BOTTOM PANEL: COST CURVES ====
# Extract cost data
LoadShapeCost_dict = data[:LoadShapeCost_dict]
num_subs = length(data[:Sset])
all_costs = vcat([LoadShapeCost_dict[i] .* 100 for i in 1:num_subs]...)  # Convert to cents/kWh
cost_min = floor(0.95 * minimum(all_costs))
cost_max = ceil(1.05 * maximum(all_costs))

p_cost = plot(
    dpi=400,
    xlabel="Time Period (t)",
    ylabel="Cost [cents/kWh]",
    legend=:topright,
    legend_columns=2,
    legendfontsize=8,
    legend_background_color=RGBA(1,1,1,0.7),
    framestyle=:box,
    ylims=(cost_min, cost_max),
    xticks=(xtick_sparse, string.(xtick_sparse)),
    grid=true,
    gridlinewidth=0.5,
    gridalpha=0.25,
    gridstyle=:solid,
    guidefont=font(10, "Computer Modern"),
    tickfont=font(9, "Computer Modern"),
    left_margin=3Plots.mm,
    right_margin=3Plots.mm,
    top_margin=0Plots.mm,
    bottom_margin=4Plots.mm
)

# Plot each substation's cost curve
for i in 1:num_subs
    cost_cents = LoadShapeCost_dict[i] .* 100
    plot!(
        p_cost, Tset, cost_cents,
        color=cost_colors[i],
        lw=2.2,
        markershape=marker_shapes[i],
        markersize=5,
        markerstrokewidth=1.5,
        markerstrokecolor=:black,
        markeralpha=0.9,
        label=L"C^t_%$i"
    )
end

# Stack the two panels vertically with IEEE 2-column dimensions
p_input = plot(p_load, p_cost, 
         layout=(2,1), 
         size=(710, 460),
         plot_title="",
         link=:x)

savefig(p_input, joinpath(plots_dir, "input_curves.png"))
println("  ✓ Saved input curves plot")

println("="^80)

#=
# ============================================================================
# OLD 3-POI FUNCTIONS AND CODE - COMMENTED OUT
# ============================================================================
# These functions are specific to the simple 3-POI 4-bus star topology

"""
    compute_network_matrices(data, slack_node)

Compute the incidence matrices B, B_T (tree), B_⊥ (chord) for the 4-bus system
with the specified slack bus removed.

Topology:
- Node 0 (substation 1): can be slack
- Node 1 (substation 2): can be slack
- Node 2 (substation 3): can be slack
- Node 3 (load): never slack
- Link 1: connects node 0 → node 3
- Link 2: connects node 1 → node 3
- Link 3: connects node 2 → node 3

Arguments:
- data: system data dictionary
- slack_node: which node is slack (0, 1, or 2)

Returns a Dict with:
- :B - reduced incidence matrix (m × n) after removing slack node
- :B_T - tree part of B (n × n)
- :B_⊥ - chord part of B ((m-n) × n) - empty for radial network
- :spanning_tree_links - indices of tree links
- :chord_links - indices of chord links (empty for radial)
"""
function compute_network_matrices(data, slack_node::Int)
    # Network topology for 4-bus system:
    # Nodes: 0 (substation 1), 1 (substation 2), 2 (substation 3), 3 (load)
    # Links: e=1 (0→3), e=2 (1→3), e=3 (2→3)
    
    # Full incidence matrix C (before removing slack)
    # Rows: nodes 0, 1, 2, 3
    # Cols: links 1, 2, 3
    # C[i,e] = +1 if link e leaves node i
    #        = -1 if link e enters node i
    #        =  0 otherwise
    
    # Link 1: 0 → 3 (leaves node 0, enters node 3)
    # Link 2: 1 → 3 (leaves node 1, enters node 3)
    # Link 3: 2 → 3 (leaves node 2, enters node 3)
    
    C_full = [
        1   0   0;   # node 0: link 1 leaves, links 2,3 not incident
        0   1   0;   # node 1: link 2 leaves, links 1,3 not incident  
        0   0   1;   # node 2: link 3 leaves, links 1,2 not incident
       -1  -1  -1    # node 3: all links enter
    ]
    
    # Remove slack node row and transpose to get B
    if slack_node == 0
        # Remove node 0 (row 1), keep nodes 1, 2, 3
        B = C_full[2:end, :]'
        n_nodes_reduced = 3
    elseif slack_node == 1
        # Remove node 1 (row 2), keep nodes 0, 2, 3
        B = C_full[[1,3,4], :]'
        n_nodes_reduced = 3
    elseif slack_node == 2
        # Remove node 2 (row 3), keep nodes 0, 1, 3
        B = C_full[[1,2,4], :]'
        n_nodes_reduced = 3
    else
        error("Invalid slack node: $slack_node. Must be 0, 1, or 2.")
    end
    
    m_links = 3
    
    # For a radial network, all links are in the spanning tree
    spanning_tree_links = [1, 2, 3]  # All three links are in the tree
    chord_links = Int[]               # No chord links
    
    # B_T: tree part (should be n × n and invertible)
    B_T = B[spanning_tree_links, :]
    
    # Check if B_T is invertible
    det_B_T = det(B_T)
    
    # B_⊥: chord part ((m-n) × n matrix) - empty for radial network
    B_chord = zeros(0, n_nodes_reduced)
    
    println("  Network with slack=$slack_node: det(B_T) = $(round(det_B_T, digits=3))")
    
    return Dict(
        :B => B,
        :B_T => B_T,
        :B_chord => B_chord,
        :n_nodes => n_nodes_reduced,
        :m_links => m_links,
        :spanning_tree_links => spanning_tree_links,
        :chord_links => chord_links,
        :C_full => C_full,
        :slack_node => slack_node
    )
end


"""
    compute_beta_angles(data, opf_result, slack_node, v_nonslack=nothing)

Compute the angle difference β_ij = ∠(v_i - z*_ij S_ij) for each link (i,j) ∈ E
after solving the OPF, as per equation (27) in the paper.

Arguments:
- data: system data dictionary
- opf_result: Dict with OPF solution (:P_1j, :Q_1j, :P_2j, :Q_2j, :P_3j, :Q_3j, :v_j, etc.)
- slack_node: which node is slack (0, 1, or 2)
- v_nonslack: Dict with voltage magnitudes squared for non-slack substations (if variable)
              Keys are substation indices (1, 2, or 3)

Returns Dict with:
- :beta - vector of angle differences for all links (in radians)
- :beta_deg - vector of angle differences for all links (in degrees)
"""
function compute_beta_angles(data, opf_result, slack_node::Int, v_nonslack=nothing)
    # Compute β_ij = ∠(v_i - z*_ij S_ij) for each link (i,j) ∈ E
    # as per equation (27) in the paper
    
    # Extract voltage magnitudes
    same_voltage_levels = data[:same_voltage_levels]
    
    if same_voltage_levels
        # Fixed voltages
        V_0_pu = data[:V_1_pu]  # Node 0 (substation 1)
        V_1_pu = data[:V_2_pu]  # Node 1 (substation 2)
        V_2_pu = data[:V_3_pu]  # Node 2 (substation 3)
    else
        # Variable voltages for non-slack substations
        # v_nonslack is a Dict with keys being substation indices
        if slack_node == 0
            V_0_pu = data[:V_1_pu]  # Slack (fixed)
            V_1_pu = sqrt(v_nonslack[2])  # Non-slack (variable)
            V_2_pu = sqrt(v_nonslack[3])  # Non-slack (variable)
        elseif slack_node == 1
            V_0_pu = sqrt(v_nonslack[1])  # Non-slack (variable)
            V_1_pu = data[:V_2_pu]  # Slack (fixed)
            V_2_pu = sqrt(v_nonslack[3])  # Non-slack (variable)
        else  # slack_node == 2
            V_0_pu = sqrt(v_nonslack[1])  # Non-slack (variable)
            V_1_pu = sqrt(v_nonslack[2])  # Non-slack (variable)
            V_2_pu = data[:V_3_pu]  # Slack (fixed)
        end
    end
    
    V_3_pu = sqrt(opf_result[:v_j])  # Node 3 (load) - from OPF
    
    # Set angles for each node (slack = 0, others use input angles for β computation)
    if slack_node == 0
        δ_0 = 0.0
        δ_1 = data[:delta_2_rad]
        δ_2 = data[:delta_3_rad]
    elseif slack_node == 1
        δ_0 = data[:delta_1_rad]
        δ_1 = 0.0
        δ_2 = data[:delta_3_rad]
    else  # slack_node == 2
        δ_0 = data[:delta_1_rad]
        δ_1 = data[:delta_2_rad]
        δ_2 = 0.0
    end
    
    # Extract power flows
    P_03 = opf_result[:P_1j]  # Power on link 1: node 0 → node 3
    Q_03 = opf_result[:Q_1j]
    P_13 = opf_result[:P_2j]  # Power on link 2: node 1 → node 3
    Q_13 = opf_result[:Q_2j]
    P_23 = opf_result[:P_3j]  # Power on link 3: node 2 → node 3
    Q_23 = opf_result[:Q_3j]
    
    # Extract line impedances
    r_1 = data[:r1_pu]
    x_1 = data[:x1_pu]
    z_1 = r_1 + im*x_1
    
    r_2 = data[:r2_pu]
    x_2 = data[:x2_pu]
    z_2 = r_2 + im*x_2
    
    r_3 = data[:r3_pu]
    x_3 = data[:x3_pu]
    z_3 = r_3 + im*x_3
    
    # Compute β for each link: β_ij = ∠(v_i - z*_ij S_ij)
    # Link 1: from node 0 to node 3
    S_03 = P_03 + im*Q_03
    V_0_complex = V_0_pu * exp(im * δ_0)
    V_0_minus_zS = V_0_complex - conj(z_1) * S_03
    β_1 = angle(V_0_minus_zS)
    
    # Link 2: from node 1 to node 3
    S_13 = P_13 + im*Q_13
    V_1_complex = V_1_pu * exp(im * δ_1)
    V_1_minus_zS = V_1_complex - conj(z_2) * S_13
    β_2 = angle(V_1_minus_zS)
    
    # Link 3: from node 2 to node 3
    S_23 = P_23 + im*Q_23
    V_2_complex = V_2_pu * exp(im * δ_2)
    V_2_minus_zS = V_2_complex - conj(z_3) * S_23
    β_3 = angle(V_2_minus_zS)
    
    β_vec = [β_1, β_2, β_3]
    
    return Dict(
        :beta => β_vec,
        :beta_deg => rad2deg.(β_vec)
    )
end


"""
    solve_three_poi_mpopf(data, slack_node)

Solve multi-period optimal power flow for 3-POI system with specified slack bus.

Arguments:
- data: system data dictionary
- slack_node: which node is slack (0, 1, or 2, corresponding to substations 1, 2, and 3)

Returns Dict with time-indexed arrays for all decision variables.
"""
function solve_three_poi_mpopf(data, slack_node::Int)
    # Extract time set and parameters
    Tset = data[:Tset]
    T = data[:T]
    Δt = data[:delta_t_h]
    P_BASE = data[:kVA_B]
    same_voltage_levels = data[:same_voltage_levels]
    
    # Compute network topology matrices for this slack configuration
    println("  Computing network matrices for slack bus = node $slack_node...")
    network_matrices = compute_network_matrices(data, slack_node)
    
    # Create model
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "max_iter", 3000)
    set_optimizer_attribute(model, "tol", 1e-6)
    
    # ========== VARIABLES (time-indexed) ==========
    @variable(model, 0 <= P_1j[t in Tset] <= data[:P_1_max_pu])
    @variable(model, Q_1j[t in Tset])
    @variable(model, l_1j[t in Tset] >= 0)
    @variable(model, 0 <= P_2j[t in Tset] <= data[:P_2_max_pu])
    @variable(model, Q_2j[t in Tset])
    @variable(model, l_2j[t in Tset] >= 0)
    @variable(model, 0 <= P_3j[t in Tset] <= data[:P_3_max_pu])
    @variable(model, Q_3j[t in Tset])
    @variable(model, l_3j[t in Tset] >= 0)
    @variable(model, v_j[t in Tset])
    
    # Voltage variables for non-slack substations (if same_voltage_levels = false)
    if !same_voltage_levels
        if slack_node == 0
            # Node 0 is slack, nodes 1 and 2 voltages become decision variables
            @variable(model, data[:Vminpu]^2 <= v_1[t in Tset] <= data[:Vmaxpu]^2)
            @variable(model, data[:Vminpu]^2 <= v_2[t in Tset] <= data[:Vmaxpu]^2)
        elseif slack_node == 1
            # Node 1 is slack, nodes 0 and 2 voltages become decision variables
            @variable(model, data[:Vminpu]^2 <= v_0[t in Tset] <= data[:Vmaxpu]^2)
            @variable(model, data[:Vminpu]^2 <= v_2[t in Tset] <= data[:Vmaxpu]^2)
        else  # slack_node == 2
            # Node 2 is slack, nodes 0 and 1 voltages become decision variables
            @variable(model, data[:Vminpu]^2 <= v_0[t in Tset] <= data[:Vmaxpu]^2)
            @variable(model, data[:Vminpu]^2 <= v_1[t in Tset] <= data[:Vmaxpu]^2)
        end
    end
    
    # Extract parameters (time-independent)
    r1j = data[:r1_pu]
    x1j = data[:x1_pu]
    r2j = data[:r2_pu]
    x2j = data[:x2_pu]
    r3j = data[:r3_pu]
    x3j = data[:x3_pu]
    Vminpu = data[:Vminpu]
    Vmaxpu = data[:Vmaxpu]
    
    # Set substation voltages based on same_voltage_levels flag
    if same_voltage_levels
        # Fixed voltages
        v_1 = data[:V_1_pu]^2
        v_2 = data[:V_2_pu]^2
        v_3 = data[:V_3_pu]^2
    else
        # Voltage depends on slack configuration - will be used in constraints below
        # (handled differently in constraint section)
    end

    # ========== CONSTRAINTS (for each time period) ==========
    for t in Tset
        # Load at time t
        P_L_t = data[:P_L_pu][t]
        Q_L_t = data[:Q_L_pu][t]
        
        # 1. Real power balance at load
        @constraint(model, (P_1j[t] - r1j * l_1j[t]) + (P_2j[t] - r2j * l_2j[t]) + (P_3j[t] - r3j * l_3j[t]) == P_L_t)
        
        # 2. Reactive power balance at load
        @constraint(model, (Q_1j[t] - x1j * l_1j[t]) + (Q_2j[t] - x2j * l_2j[t]) + (Q_3j[t] - x3j * l_3j[t]) == Q_L_t)
        
        # 3-5. KVL constraints for lines 1, 2, 3 (node 0/1/2 → node 3)
        # 6-8. SOC relaxation for lines 1, 2, 3
        if same_voltage_levels
            # Fixed voltages
            v_0_val = data[:V_1_pu]^2
            v_1_val = data[:V_2_pu]^2
            v_2_val = data[:V_3_pu]^2
            
            @constraint(model, v_j[t] == v_0_val - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
            @constraint(model, v_j[t] == v_1_val - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
            @constraint(model, v_j[t] == v_2_val - 2 * (r3j * P_3j[t] + x3j * Q_3j[t]) + (r3j^2 + x3j^2) * l_3j[t])
            @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_0_val * l_1j[t])
            @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_1_val * l_2j[t])
            @constraint(model, P_3j[t]^2 + Q_3j[t]^2 >= v_2_val * l_3j[t])
        else
            # Variable voltages for non-slack substations
            if slack_node == 0
                # Node 0 (Sub 1) is slack with fixed voltage
                # Nodes 1 (Sub 2) and 2 (Sub 3) have variable voltages v_1[t], v_2[t]
                v_0_val = data[:V_1_pu]^2  # Slack voltage (fixed)
                
                @constraint(model, v_j[t] == v_0_val - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
                @constraint(model, v_j[t] == v_1[t] - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
                @constraint(model, v_j[t] == v_2[t] - 2 * (r3j * P_3j[t] + x3j * Q_3j[t]) + (r3j^2 + x3j^2) * l_3j[t])
                @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_0_val * l_1j[t])
                @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_1[t] * l_2j[t])
                @constraint(model, P_3j[t]^2 + Q_3j[t]^2 >= v_2[t] * l_3j[t])
            elseif slack_node == 1
                # Node 1 (Sub 2) is slack with fixed voltage
                # Nodes 0 (Sub 1) and 2 (Sub 3) have variable voltages v_0[t], v_2[t]
                v_1_val = data[:V_2_pu]^2  # Slack voltage (fixed)
                
                @constraint(model, v_j[t] == v_0[t] - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
                @constraint(model, v_j[t] == v_1_val - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
                @constraint(model, v_j[t] == v_2[t] - 2 * (r3j * P_3j[t] + x3j * Q_3j[t]) + (r3j^2 + x3j^2) * l_3j[t])
                @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_0[t] * l_1j[t])
                @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_1_val * l_2j[t])
                @constraint(model, P_3j[t]^2 + Q_3j[t]^2 >= v_2[t] * l_3j[t])
            else  # slack_node == 2
                # Node 2 (Sub 3) is slack with fixed voltage
                # Nodes 0 (Sub 1) and 1 (Sub 2) have variable voltages v_0[t], v_1[t]
                v_2_val = data[:V_3_pu]^2  # Slack voltage (fixed)
                
                @constraint(model, v_j[t] == v_0[t] - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
                @constraint(model, v_j[t] == v_1[t] - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
                @constraint(model, v_j[t] == v_2_val - 2 * (r3j * P_3j[t] + x3j * Q_3j[t]) + (r3j^2 + x3j^2) * l_3j[t])
                @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_0[t] * l_1j[t])
                @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_1[t] * l_2j[t])
                @constraint(model, P_3j[t]^2 + Q_3j[t]^2 >= v_2_val * l_3j[t])
            end
        end
        
        # 9. Voltage limits (load bus always constrained)
        @constraint(model, Vminpu^2 <= v_j[t] <= Vmaxpu^2)
    end

    # ========== OBJECTIVE ==========
    # Total energy cost over all time periods
    # Cost is time-varying and different for each substation
    LoadShapeCost_1 = data[:LoadShapeCost_1]
    LoadShapeCost_2 = data[:LoadShapeCost_2]
    LoadShapeCost_3 = data[:LoadShapeCost_3]
    
    @expression(model, energy_cost_total,
        sum(LoadShapeCost_1[t] * P_1j[t] * P_BASE * Δt + 
            LoadShapeCost_2[t] * P_2j[t] * P_BASE * Δt +
            LoadShapeCost_3[t] * P_3j[t] * P_BASE * Δt for t in Tset))
    
    @objective(model, Min, energy_cost_total)

    # ========== SOLVE ==========
    println("\nSolving MPOPF with T=$(T) time periods...")
    solve_start = time()
    optimize!(model)
    solve_time_val = time() - solve_start
    
    # Check termination status
    status = termination_status(model)
    
    if status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL
        status_crayon = Crayon(foreground = :green, bold = true)
        println(status_crayon("✓ MPOPF Solved: $status"))
        println("Solve time: $(round(solve_time_val, digits=2)) seconds")
    else
        status_crayon = Crayon(foreground = :red, bold = true)
        println(status_crayon("✗ MPOPF Failed: $status"))
    end

    # ========== COLLECT RESULTS ==========
    result = Dict(
        :status => status,
        :objective => has_values(model) ? objective_value(model) : NaN,
        :solve_time => solve_time_val,
        :slack_node => slack_node,
        :same_voltage_levels => same_voltage_levels,
        # Store time-indexed arrays
        :P_1j => [value(P_1j[t]) for t in Tset],
        :Q_1j => [value(Q_1j[t]) for t in Tset],
        :l_1j => [value(l_1j[t]) for t in Tset],
        :P_2j => [value(P_2j[t]) for t in Tset],
        :Q_2j => [value(Q_2j[t]) for t in Tset],
        :l_2j => [value(l_2j[t]) for t in Tset],
        :P_3j => [value(P_3j[t]) for t in Tset],
        :Q_3j => [value(Q_3j[t]) for t in Tset],
        :l_3j => [value(l_3j[t]) for t in Tset],
        :v_j => [value(v_j[t]) for t in Tset],
        # Store network matrices for angle computation
        :network_matrices => network_matrices
    )
    
    # Store non-slack substation voltages if variable
    if !same_voltage_levels
        result[:v_nonslack] = Dict()
        if slack_node == 0
            # Nodes 1 and 2 voltages were variable
            result[:v_nonslack][2] = [value(v_1[t]) for t in Tset]
            result[:v_nonslack][3] = [value(v_2[t]) for t in Tset]
        elseif slack_node == 1
            # Nodes 0 and 2 voltages were variable
            result[:v_nonslack][1] = [value(v_0[t]) for t in Tset]
            result[:v_nonslack][3] = [value(v_2[t]) for t in Tset]
        else  # slack_node == 2
            # Nodes 0 and 1 voltages were variable
            result[:v_nonslack][1] = [value(v_0[t]) for t in Tset]
            result[:v_nonslack][2] = [value(v_1[t]) for t in Tset]
        end
    end
    
    return result
end

"""
    analyze_slack_configuration_mpopf(slack_node, data)

Solve MPOPF with specified slack bus and compute angle statistics across time.

Returns Dict with:
- :objective - total cost over all time periods
- :angles_deg - matrix of angles [substation × time]
- :median_angles_deg - median angle for each non-slack substation
- :angle_ranges_deg - [min, max] angle for each non-slack substation
- :rmse_values_deg - RMSE (temporal variation) for each non-slack substation
- :total_deviation_deg - grand RMS of all RMSE values
- Other power flow results
"""
function analyze_slack_configuration_mpopf(slack_node::Int, data)
    # Solve MPOPF
    result = solve_three_poi_mpopf(data, slack_node)
    
    if !(result[:status] == MOI.OPTIMAL || result[:status] == MOI.LOCALLY_SOLVED)
        @warn "MPOPF failed for slack_node=$slack_node"
        return result
    end
    
    # Extract time parameters
    T = data[:T]
    Tset = data[:Tset]
    
    # Compute angles for each time period
    # angles_matrix[substation_idx, time_idx]
    angles_matrix = Dict()
    
    for t in Tset
        # Package results for this time period
        opf_result_t = Dict(
            :P_1j => result[:P_1j][t],
            :Q_1j => result[:Q_1j][t],
            :P_2j => result[:P_2j][t],
            :Q_2j => result[:Q_2j][t],
            :P_3j => result[:P_3j][t],
            :Q_3j => result[:Q_3j][t],
            :v_j => result[:v_j][t]
        )
        
        # Get non-slack voltages if variable (Dict with substation indices as keys)
        v_nonslack_t = nothing
        if haskey(result, :v_nonslack)
            v_nonslack_t = Dict()
            for (sub_idx, v_vec) in result[:v_nonslack]
                v_nonslack_t[sub_idx] = v_vec[t]
            end
        end
        
        # Compute β angles
        beta_result = compute_beta_angles(data, opf_result_t, slack_node, v_nonslack_t)
        β = beta_result[:beta]
        
        # Get network matrices
        B_T = result[:network_matrices][:B_T]
        
        # Solve for angles: θ_★ = B_T^(-1) * β
        θ_star = B_T \ β
        
        # Project to (-π, π]
        θ_star = mod.(θ_star .+ π, 2π) .- π
        
        # Map θ_star to substation indices based on which nodes remain after slack removal
        # θ_star has 3 elements (one for each non-slack node after removing slack)
        if slack_node == 0
            # Removed node 0 (Sub 1), kept nodes 1, 2, 3
            # θ_star[1] → node 1 (Sub 2)
            # θ_star[2] → node 2 (Sub 3)
            # θ_star[3] → node 3 (Load) - not stored
            node_to_sub = [2, 3]
        elseif slack_node == 1
            # Removed node 1 (Sub 2), kept nodes 0, 2, 3
            # θ_star[1] → node 0 (Sub 1)
            # θ_star[2] → node 2 (Sub 3)
            # θ_star[3] → node 3 (Load) - not stored
            node_to_sub = [1, 3]
        else  # slack_node == 2
            # Removed node 2 (Sub 3), kept nodes 0, 1, 3
            # θ_star[1] → node 0 (Sub 1)
            # θ_star[2] → node 1 (Sub 2)
            # θ_star[3] → node 3 (Load) - not stored
            node_to_sub = [1, 2]
        end
        
        # Store angles for non-slack substations (first 2 elements of θ_star)
        for (i, sub_idx) in enumerate(node_to_sub)
            angle_deg = rad2deg(θ_star[i])
            if !haskey(angles_matrix, sub_idx)
                angles_matrix[sub_idx] = Float64[]
            end
            push!(angles_matrix[sub_idx], angle_deg)
        end
    end
    
    # Compute statistics across time for each non-slack substation
    median_angles_deg = Dict()
    angle_ranges_deg = Dict()
    rmse_values_deg = Dict()
    
    for (sub, angles_vec) in angles_matrix
        # Median
        median_angles_deg[sub] = median(angles_vec)
        
        # Range [min, max]
        angle_ranges_deg[sub] = [minimum(angles_vec), maximum(angles_vec)]
        
        # RMSE (temporal variation around median)
        rmse_values_deg[sub] = sqrt(mean((angles_vec .- median_angles_deg[sub]).^2))
    end
    
    # Total angle deviation (grand RMS)
    if length(rmse_values_deg) > 0
        total_deviation = sqrt(sum(v^2 for v in values(rmse_values_deg)) / length(rmse_values_deg))
    else
        total_deviation = 0.0
    end
    
    # Add statistics to result
    result[:angles_matrix] = angles_matrix
    result[:median_angles_deg] = median_angles_deg
    result[:angle_ranges_deg] = angle_ranges_deg
    result[:rmse_values_deg] = rmse_values_deg
    result[:total_deviation_deg] = total_deviation
    
    return result
end


println("\n" * "="^80)
println("NETWORK TOPOLOGY")
println("="^80)
println("4-bus system: 3 substations (nodes 0, 1, 2) + 1 load (node 3)")
println("Links: 0→3 (link 1), 1→3 (link 2), 2→3 (link 3)")
println("Radial tree structure (no loops)")
println("="^80)

# ==================================================================================
# SOLVE MPOPF FOR BOTH SLACK BUS CONFIGURATIONS
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :cyan, bold = true)("SOLVING MULTI-PERIOD OPF FOR ALL THREE SLACK CONFIGURATIONS"))
println("="^80)

# Create results directory structure: processedData/<systemName>_T<horizon>/<voltage_config>/
system_folder_name = "$(data[:systemName])_T$(data[:T])"
voltage_config_folder = data[:same_voltage_levels] ? "symmetrical_voltages" : "variable_voltages"
results_base_dir = joinpath(@__DIR__, "envs", "multi_poi", "processedData", system_folder_name, voltage_config_folder)
plots_dir = joinpath(results_base_dir, "plots")
mkpath(plots_dir)

println("\nResults will be saved to: $results_base_dir")
println("="^80)

# Analyze each slack bus configuration
results_by_slack = Dict()

for slack_sub in [1, 2, 3]
    slack_node = slack_sub - 1  # Convert substation number to node index
    println("\n" * "-"^80)
    println(Crayon(foreground = :yellow, bold = true)("Solving MPOPF with Substation $slack_sub (node $slack_node) as slack bus..."))
    println("-"^80)
    results_by_slack[slack_sub] = analyze_slack_configuration_mpopf(slack_node, data)
end

# ==================================================================================
# PRINT COMPARISON TABLES (TABLE II FIRST, THEN TABLE I)
# ==================================================================================

# ==================================================================================
# DETAILED ANGLE COORDINATION (TABLE II) FOR EACH CONFIGURATION
# ==================================================================================

for slack_sub in sort(collect(keys(results_by_slack)))
    result = results_by_slack[slack_sub]
    
    println("\n" * "="^80)
    title2_crayon = Crayon(foreground = :light_magenta, bold = true)
    println(title2_crayon("TABLE II: ANGLE COORDINATION DETAILS (SLACK: SUBSTATION $slack_sub)"))
    println("="^80)
    
    header2_crayon = Crayon(foreground = :white, bold = true)
    header_detail = @sprintf("%-15s | %-20s | %-25s | %-15s", 
                            "Substation", "Median", "Angle Range", "Temporal RMSE")
    header_detail2 = @sprintf("%-15s | %-20s | %-25s | %-15s",
                             "(Non-Slack)", "Angle (°)", "(°)", "(°)")
    sep_detail = "-"^length(header_detail)
    
    println(sep_detail)
    println(header2_crayon(header_detail))
    println(header2_crayon(header_detail2))
    println(sep_detail)
    
    # Get non-slack substations
    non_slack_subs = sort(collect(keys(result[:angles_matrix])))
    
    for sub in non_slack_subs
        median_angle = result[:median_angles_deg][sub]
        angle_range = result[:angle_ranges_deg][sub]
        rmse = result[:rmse_values_deg][sub]
        
        sub_str = Crayon(foreground = :yellow)(@sprintf("%-15d", sub))
        median_str = @sprintf("θ%-2d,%-2d = %7.3f", sub, slack_sub, median_angle)
        range_str = @sprintf("[%7.3f, %7.3f]", angle_range[1], angle_range[2])
        rmse_str = @sprintf("%-15.3f", rmse)
        
        @printf("%s | %-20s | %-25s | %s\n", sub_str, median_str, range_str, rmse_str)
    end
    
    println(sep_detail)
    
    # Print computed substation angles (showing median and range)
    println()
    println("Computed Substation Angles (across time):")
    
    # Create angle arrays with slack = 0 and non-slack from results
    θ_medians = Dict{Int, Float64}()
    θ_ranges = Dict{Int, Vector{Float64}}()
    
    for sub in 1:3
        if sub == slack_sub
            θ_medians[sub] = 0.0
            θ_ranges[sub] = [0.0, 0.0]
        else
            θ_medians[sub] = result[:median_angles_deg][sub]
            θ_ranges[sub] = result[:angle_ranges_deg][sub]
        end
    end
    
    # Print each substation's angle
    theta_colors = [:green, :yellow, :cyan]
    for sub in 1:3
        color = theta_colors[sub]
        if sub == slack_sub
            theta_str = Crayon(foreground = color, bold = true)("θ$sub (Sub $sub, slack)")
            @printf("  %s = %7.3f° (fixed)\n", theta_str, θ_medians[sub])
        else
            theta_str = Crayon(foreground = color, bold = true)("θ$sub (Sub $sub)")
            @printf("  %s = %7.3f° (median), range [%7.3f°, %7.3f°]\n",
                    theta_str, θ_medians[sub], θ_ranges[sub][1], θ_ranges[sub][2])
        end
    end
    
    # Print angle differences
    println()
    delta_str = Crayon(foreground = :magenta, bold = true)("Angle differences (median):")
    println(delta_str)
    @printf("  Δθ₂₁ (θ₂ - θ₁) = %7.3f°\n", θ_medians[2] - θ_medians[1])
    @printf("  Δθ₃₁ (θ₃ - θ₁) = %7.3f°\n", θ_medians[3] - θ_medians[1])
    @printf("  Δθ₃₂ (θ₃ - θ₂) = %7.3f°\n", θ_medians[3] - θ_medians[2])
    
    # Print bus voltages (time-averaged)
    println()
    v_crayon = Crayon(foreground = :light_yellow)
    v_var_crayon = Crayon(foreground = :light_magenta, bold = true)
    println("Bus Voltages (time-averaged):")
    
    same_v_levels = data[:same_voltage_levels]  # Local copy to avoid soft scope warning
    slack_node_val = result[:slack_node]  # Get slack node from result
    
    v_3_pu = mean(sqrt.(result[:v_j]))  # Load bus (node 3) - averaged over time
    
    if same_v_levels
        # Fixed voltages
        v_0_pu = data[:V_1_pu]  # Substation 1 (node 0) - fixed
        v_1_pu = data[:V_2_pu]  # Substation 2 (node 1) - fixed
        v_2_pu = data[:V_3_pu]  # Substation 3 (node 2) - fixed
        
        println("  ", v_crayon("Node 0 (Sub 1):"), " V = $(round(v_0_pu, digits=4)) pu (fixed)")
        println("  ", v_crayon("Node 1 (Sub 2):"), " V = $(round(v_1_pu, digits=4)) pu (fixed)")
        println("  ", v_crayon("Node 2 (Sub 3):"), " V = $(round(v_2_pu, digits=4)) pu (fixed)")
        println("  ", v_crayon("Node 3 (Load): "), " V = $(round(v_3_pu, digits=4)) pu (avg)")
    else
        # Variable voltages for non-slack substations
        if slack_node_val == 0
            # Node 0 is slack (fixed), nodes 1 and 2 are variable
            v_0_pu = data[:V_1_pu]
            v_1_pu_avg = mean(sqrt.(result[:v_nonslack][2]))
            v_1_pu_min = minimum(sqrt.(result[:v_nonslack][2]))
            v_1_pu_max = maximum(sqrt.(result[:v_nonslack][2]))
            v_2_pu_avg = mean(sqrt.(result[:v_nonslack][3]))
            v_2_pu_min = minimum(sqrt.(result[:v_nonslack][3]))
            v_2_pu_max = maximum(sqrt.(result[:v_nonslack][3]))
            
            println("  ", v_crayon("Node 0 (Sub 1):"), " V = $(round(v_0_pu, digits=4)) pu (slack, fixed)")
            println("  ", v_var_crayon("Node 1 (Sub 2):"), " V = $(round(v_1_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_1_pu_min, digits=4)), $(round(v_1_pu_max, digits=4))] pu")
            println("  ", v_var_crayon("Node 2 (Sub 3):"), " V = $(round(v_2_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_2_pu_min, digits=4)), $(round(v_2_pu_max, digits=4))] pu")
        elseif slack_node_val == 1
            # Node 1 is slack (fixed), nodes 0 and 2 are variable
            v_0_pu_avg = mean(sqrt.(result[:v_nonslack][1]))
            v_0_pu_min = minimum(sqrt.(result[:v_nonslack][1]))
            v_0_pu_max = maximum(sqrt.(result[:v_nonslack][1]))
            v_1_pu = data[:V_2_pu]
            v_2_pu_avg = mean(sqrt.(result[:v_nonslack][3]))
            v_2_pu_min = minimum(sqrt.(result[:v_nonslack][3]))
            v_2_pu_max = maximum(sqrt.(result[:v_nonslack][3]))
            
            println("  ", v_var_crayon("Node 0 (Sub 1):"), " V = $(round(v_0_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_0_pu_min, digits=4)), $(round(v_0_pu_max, digits=4))] pu")
            println("  ", v_crayon("Node 1 (Sub 2):"), " V = $(round(v_1_pu, digits=4)) pu (slack, fixed)")
            println("  ", v_var_crayon("Node 2 (Sub 3):"), " V = $(round(v_2_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_2_pu_min, digits=4)), $(round(v_2_pu_max, digits=4))] pu")
        else  # slack_node_val == 2
            # Node 2 is slack (fixed), nodes 0 and 1 are variable
            v_0_pu_avg = mean(sqrt.(result[:v_nonslack][1]))
            v_0_pu_min = minimum(sqrt.(result[:v_nonslack][1]))
            v_0_pu_max = maximum(sqrt.(result[:v_nonslack][1]))
            v_1_pu_avg = mean(sqrt.(result[:v_nonslack][2]))
            v_1_pu_min = minimum(sqrt.(result[:v_nonslack][2]))
            v_1_pu_max = maximum(sqrt.(result[:v_nonslack][2]))
            v_2_pu = data[:V_3_pu]
            
            println("  ", v_var_crayon("Node 0 (Sub 1):"), " V = $(round(v_0_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_0_pu_min, digits=4)), $(round(v_0_pu_max, digits=4))] pu")
            println("  ", v_var_crayon("Node 1 (Sub 2):"), " V = $(round(v_1_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_1_pu_min, digits=4)), $(round(v_1_pu_max, digits=4))] pu")
            println("  ", v_crayon("Node 2 (Sub 3):"), " V = $(round(v_2_pu, digits=4)) pu (slack, fixed)")
        end
        
        println("  ", v_crayon("Node 3 (Load): "), " V = $(round(v_3_pu, digits=4)) pu (avg)")
    end
    
    # Print power dispatch summary (time-averaged)
    println()
    p1_crayon = Crayon(foreground = :light_green)
    p2_crayon = Crayon(foreground = :light_cyan)
    p3_crayon = Crayon(foreground = :light_blue)
    
    P_BASE = data[:kVA_B]
    println("Power Dispatches (time-averaged):")
    p1_avg = mean(result[:P_1j]) * P_BASE
    q1_avg = mean(result[:Q_1j]) * P_BASE
    p2_avg = mean(result[:P_2j]) * P_BASE
    q2_avg = mean(result[:Q_2j]) * P_BASE
    p3_avg = mean(result[:P_3j]) * P_BASE
    q3_avg = mean(result[:Q_3j]) * P_BASE
    
    println("  ", p1_crayon("Substation 1:"), " P = $(round(p1_avg, digits=1)) kW, Q = $(round(q1_avg, digits=1)) kVAr")
    println("  ", p2_crayon("Substation 2:"), " P = $(round(p2_avg, digits=1)) kW, Q = $(round(q2_avg, digits=1)) kVAr")
    println("  ", p3_crayon("Substation 3:"), " P = $(round(p3_avg, digits=1)) kW, Q = $(round(q3_avg, digits=1)) kVAr")
    
    cost_crayon = Crayon(foreground = :light_red, bold = true)
    println()
    println("  ", cost_crayon("Total Energy Cost:"), " \$$(round(result[:objective], digits=2))")
    println("  Solve time: $(round(result[:solve_time], digits=2)) seconds")
end

# ==================================================================================
# SUMMARY COMPARISON TABLE (TABLE I)
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :light_blue, bold = true)("TABLE I: COMPARISON OF SLACK BUS CONFIGURATIONS"))
println("="^80)

# Print header
header_crayon = Crayon(foreground = :white, bold = true)
header = @sprintf("%-15s | %-20s | %-20s", "Slack Bus", "Total Angle", "Operational")
header2 = @sprintf("%-15s | %-20s | %-20s", "(Substation)", "Deviation* (°)", "Cost (\$)")
separator = "-"^length(header)

println(separator)
println(header_crayon(header))
println(header_crayon(header2))
println(separator)

for slack_sub in sort(collect(keys(results_by_slack)))
    result = results_by_slack[slack_sub]
    
    total_cost = result[:objective]
    total_dev = result[:total_deviation_deg]
    
    # Format with colors
    slack_colors = [:green, :cyan, :light_blue]
    color = slack_colors[slack_sub]
    sub_str = Crayon(foreground = color, bold = true)(@sprintf("%-15d", slack_sub))
    
    dev_str = @sprintf("%-20.3f", total_dev)
    cost_str = @sprintf("%-20.2f", total_cost)
    
    @printf("%s | %s | %s\n", sub_str, dev_str, cost_str)
end

println(separator)
println("*RMS of temporal angle variations across all non-slack substations")
println("="^80)

println("\n" * "="^80)
success_crayon = Crayon(foreground = :light_green, bold = true)
info_crayon = Crayon(foreground = :light_blue)
println(success_crayon("✓ MPOPF SCRIPT COMPLETE!"), " Analyzed all three slack bus configurations.")
println(info_crayon("ℹ Tables I & II show angle coordination and cost comparison."))
println("="^80)

# ==================================================================================
# PLOT INPUT CURVES
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :magenta, bold = true)("GENERATING INPUT CURVE PLOTS"))
println("="^80)

# Load plotter
include(joinpath(@__DIR__, "envs", "multi_poi", "Plotter.jl"))

# Generate input curves plot
plot_filename = joinpath(plots_dir, "mpopf_input_curves.png")
plot_input_curves(data, showPlots=showPlots, savePlots=savePlots, filename=plot_filename)
print_curve_statistics(data)

# Generate angle and voltage trajectory plots for each slack configuration
println("\n" * "="^80)
println(Crayon(foreground = :magenta, bold = true)("GENERATING ANGLE AND VOLTAGE TRAJECTORY PLOTS"))
println("="^80)

for slack_sub in sort(collect(keys(results_by_slack)))
    result = results_by_slack[slack_sub]
    
    println("\nGenerating plot for Slack: Substation $(slack_sub)...")
    traj_filename = joinpath(plots_dir, "angle_voltage_slack$(slack_sub).png")
    plot_angle_voltage_trajectories(result, data, slack_sub, 
                                     showPlots=showPlots, 
                                     savePlots=savePlots, 
                                     filename=traj_filename)
    if savePlots
        println("  ✓ Plot saved: $(traj_filename)")
    end
end

println("="^80)

=#
# End of commented-out 3-POI code