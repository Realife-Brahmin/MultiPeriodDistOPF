module parseBranchData

export parse_branch_data

include("../helperFunctions.jl")
using .helperFunctions: myprintln

using DataStructures
using Parameters

#region parse_branch_data
"""
    parse_branch_data(systemName::String; kVA_B=1000, kV_B=2.4018,  verbose::Bool=false)

Parse branch data from the BranchData.dss file for a given system.

This function reads and parses the BranchData.dss file for the specified system, extracting relevant branch data and initializing various parameters. 
It handles the extraction of branch properties such as resistance, reactance, and connectivity information.

# Arguments
- `systemName::String`: The name of the system for which to parse branch data.
- `kVA_B::Float64`: The base power in kVA for per-unit calculations (default: 1000).
- `kV_B::Float64`: The base voltage in kV for per-unit calculations (default: 2.4018).
- `Z_B::Union{Nothing, Float64}`: The base impedance for per-unit calculations. If not provided, it is calculated using `kVA_B` and `kV_B` (default: 5.768643240000001).
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `branchData::Dict`: A dictionary containing the parsed branch data, including:
    - `Nset`: Set of all bus numbers.
    - `Lset`: Set of branches (edges).
    - `rdict`: Resistance of each branch.
    - `xdict`: Reactance of each branch.
    - `rdict_pu`: Per-unit resistance of each branch.
    - `xdict_pu`: Per-unit reactance of each branch.
    - `parent`: Parent node of each node.
    - `children`: Children nodes of each node.
    - `N`: Total number of buses.
    - `m`: Total number of branches.
    - `N1set`: Substation node (bus 1).
    - `Nm1set`: Buses not including substation bus (1).
    - `Nc1set`: Buses connected to substation bus (1).
    - `Nnc1set`: Buses not connected to substation bus (1).
    - `L1set`: Branches where one node is the substation (1).
    - `Lm1set`: Branches where no node is the substation (1).
    - `N1`: Number of substation nodes.
    - `Nm1`: Number of buses not including substation bus.
    - `Nc1`: Number of buses connected to substation bus.
    - `Nnc1`: Number of buses not connected to substation bus.
    - `m1`: Number of branches where one node is the substation.
    - `mm1`: Number of branches where no node is the substation.
"""
function parse_branch_data(systemName::String;
    kVA_B=1000,
    kV_B=2.4018,
    verbose::Bool=false)

    ## Now let's take a small detour to ensure that we know exactly what base impedance to compute pu values of impedances:

    MVA_B = kVA_B / 1000
    Z_B = (kV_B^2) / MVA_B

    verbose = true
    myprintln(verbose, "Final Z_B = $Z_B")

    # Todo: Ensure that substation bus being equal to 1 is not taken for granted, have some kwarg or something to ensure that even bus 153 can be the substation bus

    wd = @__DIR__
    # Construct the file path using wd
    filename = joinpath(wd, "..", "..", "rawData", systemName, "BranchData.dss")

    # Initialize data structures
    Nset = Set{Int}()                             # Set of all bus numbers
    Lset = Set{Tuple{Int,Int}}()                  # Set of branches (edges)
    rdict = Dict{Tuple{Int,Int},Float64}()     # Resistance of each branch
    xdict = Dict{Tuple{Int,Int},Float64}()     # Reactance of each branch

    rdict_pu = Dict{Tuple{Int,Int},Float64}()  # Per-unit resistance of each branch
    xdict_pu = Dict{Tuple{Int,Int},Float64}()  # Per-unit reactance of each branch

    # Define the parent dictionary to hold Int or nothing
    parent = Dict{Int,Union{Int,Nothing}}()    # Parent node of each node
    children = Dict{Int,Vector{Int}}()         # Children nodes of each node

    # Initialize additional sets and parameters
    N1set = Set{Int}()                            # Substation node (bus 1)
    Nm1set = Set{Int}()                           # Buses not including substation bus (1)
    Nc1set = Set{Int}()                           # Buses connected to substation bus (1)
    Nnc1set = Set{Int}()                          # Buses not connected to substation bus (1)
    L1set = Set{Tuple{Int,Int}}()                 # Branches where one node is the substation (1)
    Lm1set = Set{Tuple{Int,Int}}()                # Branches where no node is the substation (1)

    # Regular expression to match key=value pairs with optional spaces
    kv_pattern = r"(\w+)\s*=\s*([\S]+)"

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Remove comments and strip whitespace
            line = split(line, "!")[1]
            line = strip(line)

            # Skip empty lines
            if isempty(line)
                continue
            end

            # Parse lines starting with "New Line."
            if startswith(line, "New Line.")
                # Extract parameters from the line using regular expressions
                branch_info = Dict{String,String}()
                for m in eachmatch(kv_pattern, line)
                    key = strip(m.captures[1])
                    value = strip(m.captures[2])
                    branch_info[key] = value
                end

                # Extract from_bus and to_bus
                if haskey(branch_info, "Bus1") && haskey(branch_info, "Bus2")
                    # Extract bus numbers
                    bus1_str = branch_info["Bus1"]
                    bus2_str = branch_info["Bus2"]

                    # Extract the integer part before any decimal
                    bus1_parts = split(bus1_str, ".")
                    bus2_parts = split(bus2_str, ".")

                    from_bus = parse(Int, bus1_parts[1])
                    to_bus = parse(Int, bus2_parts[1])

                    # Add buses to Nset
                    Nset = union(Nset, [from_bus, to_bus])

                    # Add to branch set
                    Lset = union(Lset, [(from_bus, to_bus)])

                    # Update parent and children dictionaries
                    parent[to_bus] = from_bus
                    if haskey(children, from_bus)
                        push!(children[from_bus], to_bus)
                    else
                        children[from_bus] = [to_bus]
                    end

                    # Ensure to_bus is a key in the children dictionary, even if it might not have any child nodes
                    get!(children, to_bus, Int[])

                    # Ensure from_bus is a key in the parent dictionary, even if it might not have any parent node (only substation bus)
                    get!(parent, from_bus, nothing)

                    # Update N1set, Nm1set, L1set, Lm1set
                    if from_bus == 1
                        push!(N1set, from_bus)
                        push!(Nc1set, to_bus)
                        push!(L1set, (from_bus, to_bus))
                    elseif to_bus == 1
                        push!(N1set, from_bus)
                        push!(Nc1set, from_bus)
                        push!(L1set, (from_bus, to_bus))
                    else
                        Nm1set = union(Nm1set, [from_bus, to_bus])
                        push!(Lm1set, (from_bus, to_bus))
                    end

                    # # Nm1set  will be all nodes except substation node
                    Nm1set = setdiff(Nset, N1set)

                    # # Nnc1set will be all nodes that are not connected to the substation
                    Nnc1set = setdiff(Nm1set, Nc1set)

                else
                    error("Bus1 or Bus2 not specified for a line in BranchData.dss")
                end

                # Extract resistance (r1) and reactance (x1), and calculate per-unit values
                rdict[(from_bus, to_bus)] = haskey(branch_info, "r1") ? parse(Float64, branch_info["r1"]) : 0.0
                xdict[(from_bus, to_bus)] = haskey(branch_info, "x1") ? parse(Float64, branch_info["x1"]) : 0.0
            end
        end
    end

    baseValuesDict = get_base_units(systemName, Nset, Lset, kVA_B=kVA_B, kV_B=kV_B, verbose=verbose)

    Z_B_dict = Dict{Tuple{Int,Int},Float64}() # Base impedance for each branch
    for (from_bus, to_bus) ∈ Lset
        Z_B_dict[(from_bus, to_bus)] = baseValuesDict[:kV_B_dict][(from_bus, to_bus)]^2 / baseValuesDict[:MVA_B_dict][(from_bus, to_bus)]
        # Calculate per-unit values for resistance and reactance
        rdict_pu[(from_bus, to_bus)] = rdict[(from_bus, to_bus)] / Z_B_dict[(from_bus, to_bus)]
        xdict_pu[(from_bus, to_bus)] = xdict[(from_bus, to_bus)] / Z_B_dict[(from_bus, to_bus)]
    end

    # Calculate total number of buses and branches
    N = length(Nset)
    m = length(Lset)

    # Calculate additional parameters
    N1 = length(N1set)
    Nm1 = length(Nm1set)
    Nc1 = length(Nc1set)
    Nnc1 = length(Nnc1set)
    m1 = length(L1set)
    mm1 = length(Lm1set)

    # Sort each data structure in place or by reassigning to the same variable
    Nset = sort(collect(Nset))
    Lset = sort(collect(Lset))
    N1set = sort(collect(N1set))
    Nm1set = sort(collect(Nm1set))
    Nc1set = sort(collect(Nc1set))
    Nnc1set = sort(collect(Nnc1set))
    L1set = sort(collect(L1set))
    Lm1set = sort(collect(Lm1set))

    @unpack kVA_B_dict, kV_B_dict, MVA_B_dict = baseValuesDict
    # Create a dictionary with all outputs
    branchData = Dict(
        :Nset => Nset,
        :Lset => Lset,
        :rdict => rdict,
        :xdict => xdict,
        :rdict_pu => rdict_pu,
        :xdict_pu => xdict_pu,
        :parent => parent,
        :children => children,
        :kV_B_dict => kV_B_dict,
        :kVA_B_dict => kVA_B_dict,
        :MVA_B_dict => MVA_B_dict,
        :N => N,
        :m => m,
        :N1set => N1set,
        :Nm1set => Nm1set,
        :Nc1set => Nc1set,
        :Nnc1set => Nnc1set,
        :L1set => L1set,
        :Lm1set => Lm1set,
        :N1 => N1,
        :Nm1 => Nm1,
        :Nc1 => Nc1,
        :Nnc1 => Nnc1,
        :m1 => m1,
        :mm1 => mm1,
        :Z_B_dict => Z_B_dict,
    )

    return branchData

end
#endregion

function parse_transformers(file_path::String)
    transformers = Dict{Int,Float64}()
    open(file_path, "r") do file
        current_line = ""
        for line in eachline(file)
            line = strip(line)
            if startswith(line, "New Transformer")
                current_line = line
            elseif startswith(line, "~")
                current_line *= " " * line
            else
                if !isempty(current_line)
                    # Use regular expressions to extract the relevant information
                    bus1_match = match(r"wdg=1 bus=(\d+)", current_line)
                    bus2_match = match(r"wdg=2 bus=(\d+)", current_line)
                    kv2_match = match(r"wdg=2 bus=\d+ conn=\w+ kv=(\d+\.\d+)", current_line)

                    if bus1_match !== nothing && bus2_match !== nothing && kv2_match !== nothing
                        bus1 = parse(Int, bus1_match.captures[1])
                        bus2 = parse(Int, bus2_match.captures[1])
                        kv2 = parse(Float64, kv2_match.captures[1])
                        transformers[bus2] = kv2
                    end
                    current_line = ""
                end
            end
        end
        # Handle the last transformer if the file does not end with a new transformer line
        if !isempty(current_line)
            bus1_match = match(r"wdg=1 bus=(\d+)", current_line)
            bus2_match = match(r"wdg=2 bus=(\d+)", current_line)
            kv2_match = match(r"wdg=2 bus=\d+ conn=\w+ kv=(\d+\.\d+)", current_line)

            if bus1_match !== nothing && bus2_match !== nothing && kv2_match !== nothing
                bus1 = parse(Int, bus1_match.captures[1])
                bus2 = parse(Int, bus2_match.captures[1])
                kv2 = parse(Float64, kv2_match.captures[1])
                transformers[bus2] = kv2
            end
        end
    end
    @show transformers
    return transformers
end

function map_voltage_levels(Nset, Lset, transformers; kV_B=12.66)
    voltage_levels = Dict()
    visited = Set{Int}()
    stack = Stack{Int}()

    # Initialize the stack with the substation bus (assuming it's bus 1)
    push!(stack, 1)
    voltage_levels[1] = kV_B  # Assuming the substation bus is at 12.66 kV

    while !isempty(stack)
        bus = pop!(stack)
        visited = union(visited, Set([bus]))

        for (i, j) in Lset
            if i == bus && !(j in visited)
                if j in keys(transformers)
                    voltage_levels[j] = transformers[j]
                else
                    voltage_levels[j] = voltage_levels[i]
                end
                push!(stack, j)
            elseif j == bus && !(i in visited)
                if i in transformers
                    voltage_levels[i] = transformers[i]
                else
                    voltage_levels[i] = voltage_levels[j]
                end
                push!(stack, i)
            end
        end
    end

    return voltage_levels
end

function get_base_units(systemName::String, Nset, Lset; kVA_B=1000, kV_B=2.4018, verbose::Bool=false)
    kVA_B_dict = Dict()
    MVA_B_dict = Dict()
    kV_B_dict = Dict()

    if systemName == "ieee730_1ph" || systemName == "ieee729_1ph"
        file_path = joinpath(@__DIR__, "..", "..", "rawData", systemName, "LoadXfmrs.dss")
        transformers = parse_transformers(file_path)
        @show transformers
        voltage_levels = map_voltage_levels(Nset, Lset, transformers; kV_B=kV_B)

        for bus in Nset
            kV_B_dict[bus] = voltage_levels[bus]
            kVA_B_dict[bus] = kVA_B
            MVA_B_dict[bus] = kVA_B / 1000
        end

        for (i, j) in Lset
            kV_B_dict[(i, j)] = voltage_levels[j]  # Use secondary voltage for branches
            kVA_B_dict[(i, j)] = kVA_B
            MVA_B_dict[(i, j)] = kVA_B / 1000
        end
    else
        for bus in Nset
            kV_B_dict[bus] = kV_B
            kVA_B_dict[bus] = kVA_B
            MVA_B_dict[bus] = kVA_B / 1000
        end

        for (i, j) in Lset
            kV_B_dict[(i, j)] = kV_B
            kVA_B_dict[(i, j)] = kVA_B
            MVA_B_dict[(i, j)] = kVA_B / 1000
        end
    end

    baseValuesDict = Dict(
        :kVA_B_dict => kVA_B_dict,
        :MVA_B_dict => MVA_B_dict,
        :kV_B_dict => kV_B_dict
    )
    return baseValuesDict
end

end # module parseBranchData