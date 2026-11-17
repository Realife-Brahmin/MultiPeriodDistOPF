"""
Shared OpenDSS network parser for MPOPF algorithms.
Numerical parsing is unchanged; callers provide their own Julia environment.
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
        rawDataFolderPath = @__DIR__
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
    impedance_scale=1.0,
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
    parse_line_impedances!(data, kVA_B, kV_B; impedance_scale=impedance_scale)
    
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
    # Decks with named buses (e.g. ieee9500xf, SMART-DS) are numbered here; decks
    # with numbered buses keep their numbers and the code below, unchanged.
    all(n -> tryparse(Int, split(n, ".")[1]) !== nothing, bus_names) || return parse_named_topology!(data)
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
    parse_named_topology!(data::Dict)

Topology of a deck whose buses have names. The deck must have exactly one
enabled source; the enabled lines from it must form a tree (a closed tie that
makes a loop is an error). Bus 1 is the source bus, the others are numbered
breadth first from it; every branch is (parent, child). data[:busName][i] is
the name of bus i, data[:bus_id] maps a (lower-case) name to its number and
data[:lineName] a branch to its OpenDSS line. Buses that no enabled line
reaches (the middle buses of open ties) are left out.
"""
function parse_named_topology!(data::Dict)
    bname(s) = lowercase(String(split(s, ".")[1]))
    src = nothing
    for n in OpenDSSDirect.Vsources.AllNames()
        OpenDSSDirect.Circuit.SetActiveElement("Vsource.$n")
        OpenDSSDirect.CktElement.Enabled() || continue
        isnothing(src) || error("deck $(data[:systemName]): more than one enabled source; the model needs one substation")
        src = bname(OpenDSSDirect.CktElement.BusNames()[1])
    end
    isnothing(src) && error("deck $(data[:systemName]): no enabled source")
    adj = Dict{String, Vector{Tuple{String, String}}}()
    nlines = 0
    for ln in OpenDSSDirect.Lines.AllNames()
        OpenDSSDirect.Circuit.SetActiveElement("Line.$ln")
        OpenDSSDirect.CktElement.Enabled() || continue
        bs = OpenDSSDirect.CktElement.BusNames()
        a, b = bname(bs[1]), bname(bs[2])
        a == b && error("line $ln connects bus $a to itself")
        push!(get!(adj, a, Tuple{String, String}[]), (b, ln)); push!(get!(adj, b, Tuple{String, String}[]), (a, ln))
        nlines += 1
    end
    id = Dict(src => 1); order = [src]; via = Dict{String, String}(); par = Dict{String, String}()
    q = 1
    while q <= length(order)
        x = order[q]; q += 1
        for (y, ln) in get(adj, x, Tuple{String, String}[])
            ln == get(via, x, "") && continue
            haskey(id, y) && error("line $ln closes a loop at bus $y: the model needs a radial network")
            push!(order, y); id[y] = length(order); via[y] = ln; par[y] = x
        end
    end
    length(order) - 1 == nlines || error("deck $(data[:systemName]): $(nlines - length(order) + 1) enabled lines are not connected to the source")
    N = length(order)
    data[:N] = N
    data[:Nset] = collect(1:N)
    data[:m] = nlines
    data[:busName] = order
    data[:bus_id] = id
    parent = Dict{Int, Union{Nothing, Int}}(1 => nothing)
    children = Dict{Int, Vector{Int}}(i => Int[] for i in 1:N)
    Lset = Tuple{Int, Int}[]; lineName = Dict{Tuple{Int, Int}, String}()
    for y in order[2:end]
        i, j = id[par[y]], id[y]
        parent[j] = i; push!(children[i], j); push!(Lset, (i, j)); lineName[(i, j)] = via[y]
    end
    data[:Lset] = Lset
    data[:parent] = parent
    data[:children] = children
    data[:lineName] = lineName
    data[:L1set] = [(i, j) for (i, j) in Lset if i == 1]
    data[:Lm1set] = [(i, j) for (i, j) in Lset if i != 1]
    data[:Nm1set] = collect(2:N)
    println("Named buses: $(N) buses from source bus $(src), $(nlines) lines")
end

# Number of the bus of an element's bus string ("7.1" -> 7; a named bus via data[:bus_id]).
_bus_number(data::Dict, s::AbstractString) =
    haskey(data, :bus_id) ? data[:bus_id][lowercase(String(split(s, ".")[1]))] : parse(Int, split(s, ".")[1])


"""
    parse_line_impedances!(data::Dict, kVA_B::Float64, kV_B::Float64)

Parse line resistance and reactance.
"""
function parse_line_impedances!(data::Dict, kVA_B::Float64, kV_B::Float64;
    impedance_scale::Float64=1.0)
    Z_B = kV_B^2 / (kVA_B / 1000)  # Base impedance in Ohms
    
    rdict = Dict{Tuple{Int,Int}, Float64}()
    xdict = Dict{Tuple{Int,Int}, Float64}()
    rdict_pu = Dict{Tuple{Int,Int}, Float64}()
    xdict_pu = Dict{Tuple{Int,Int}, Float64}()
    
    line_names = OpenDSSDirect.Lines.AllNames()
    named = haskey(data, :lineName)
    named && (line_names = [data[:lineName][e] for e in data[:Lset]])   # enabled lines, as (parent, child)

    for (k, line_name) in enumerate(line_names)
        OpenDSSDirect.Lines.Name(line_name)

        # Get bus numbers from bus names (e.g., "7.1" -> 7)
        buses = OpenDSSDirect.CktElement.BusNames()
        i, j = named ? data[:Lset][k] : (parse(Int, split(buses[1], ".")[1]), parse(Int, split(buses[2], ".")[1]))
        
        # Get R and X (Ohms)
        R_ohm = impedance_scale * OpenDSSDirect.Lines.RMatrix()[1]
        X_ohm = impedance_scale * OpenDSSDirect.Lines.XMatrix()[1]
        
        rdict[(i,j)] = R_ohm
        xdict[(i,j)] = X_ohm
        rdict_pu[(i,j)] = R_ohm / Z_B
        xdict_pu[(i,j)] = X_ohm / Z_B
    end
    
    data[:Z_B] = Z_B
    data[:impedance_scale] = impedance_scale
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
        bus_num = _bus_number(data, bus_full)
        
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
            bus_num = _bus_number(data, bus_full)
            haskey(data, :bus_id) && haskey(p_D_R_pu, bus_num) && error("two PV units on bus $(bus_full): the model has one per bus")
            
            # Get generator rating: Pmpp, the rated maximum power. Not
            # PVsystems.kW(), which is the unit's PRESENT output: for a unit
            # with a daily shape it is zero at the hour the circuit was
            # compiled. That read every large10k unit (Daily=
            # LoadShapePVDefault, 350 kW each) as 0 kW until 2026-10-05; units
            # without a daily shape (ieee123, ieee2522) were unaffected.
            kW_rated = parse(Float64, OpenDSSDirect.Text.Command("? PVSystem.$gen_name.Pmpp"))
            kW_rated > 0 || error("PV unit $gen_name has no rated power (Pmpp = $kW_rated)")
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
            bus_num = _bus_number(data, bus_full)
            haskey(data, :bus_id) && haskey(B_R, bus_num) && error("two storage units on bus $(bus_full): the model has one per bus")
            
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
