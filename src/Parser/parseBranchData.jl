# parseBranchData.jl

module parseBranchData

export parse_branch_data

include("../helperFunctions.jl")
using .helperFunctions: myprintln

function parse_branch_data(systemName::String;
    kVA_B = 1000,
    kV_B = 2.4018,
    Z_B = 5.768643240000001,
    verbose::Bool = false)

    ## Now let's take a small detour to ensure that we know exactly what base impedance to compute pu values of impedances:
    
    MVA_B = kVA_B/1000,

    # Rule: If Z_B is not provided and kVA_B is specified but not kV_B, throw an error
    if isnothing(Z_B) && !isnothing(kVA_B) && isnothing(kV_B)
        error("Error: You must specify both kV_B and kVA_B to calculate Z_B, or provide Z_B directly.")
    end

    # Rule: If Z_B is not specified, and kV_B is provided, compute Z_B using the default or given kVA_B
    if isnothing(Z_B) && !isnothing(kV_B)
        Z_B = (kV_B^2) / MVA_B
        myprintln(verbose, "Computed Z_B = (kV_B^2) / MVA_B = $Z_B using kV_B = $kV_B and kVA_B = $kVA_B")
    end

    # Rule: If Z_B is provided and kV_B/kVA_B are not specified, just use the provided Z_B
    if !isnothing(Z_B)
        myprintln(verbose, "Using user-specified Z_B = $Z_B")
    end

    # Rule: If neither Z_B nor kV_B/kVA_B are provided, use the default value of Z_B
    if isnothing(Z_B) && isnothing(kV_B) && isnothing(kVA_B)
        Z_B = 5.768643240000001  # Default value
        myprintln(verbose, "Using default Z_B = $Z_B")
    end

    # Now Z_B is guaranteed to be set to a valid value at this point
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

                # Calculate per-unit values for resistance and reactance
                rdict_pu[(from_bus, to_bus)] = rdict[(from_bus, to_bus)] / Z_B
                xdict_pu[(from_bus, to_bus)] = xdict[(from_bus, to_bus)] / Z_B

            end
        end
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
        :mm1 => mm1
    )

    return branchData

end


end # module
