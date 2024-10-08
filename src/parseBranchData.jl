# parseBranchData.jl

module parseBranchData

export parse_branch_data

function parse_branch_data(systemName::String)

    # Todo: Ensure that substation bus being equal to 1 is not taken for granted, have some kwarg or something to ensure that even bus 153 can be the substation bus
    
    wd = @__DIR__
    # Construct the file path using wd
    filename = joinpath(wd, "..", "rawData", systemName, "BranchData.dss")

    # Initialize data structures
    Nset = Set{Int}()                         # Set of all bus numbers
    Lset = Set{Tuple{Int,Int}}()             # Set of branches (edges)
    rdict = Dict{Tuple{Int,Int},Float64}()  # Resistance of each branch
    xdict = Dict{Tuple{Int,Int},Float64}()  # Reactance of each branch
    parent = Dict{Int,Int}()                 # parent node of each node
    children = Dict{Int,Vector{Int}}()       # children nodes of each node

    # Initialize additional sets and parameters
    N1set = Set{Int}()                       # Buses connected to substation bus (1)
    Nm1set = Set{Int}()                      # Buses not connected to substation bus (1)
    L1set = Set{Tuple{Int,Int}}()           # Branches where one node is 1
    Lm1set = Set{Tuple{Int,Int}}()          # Branches where no node is 1

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

                    # Update N1set, Nm1set, L1set, Lm1set
                    if from_bus == 1
                        push!(N1set, to_bus)
                        push!(L1set, (from_bus, to_bus))
                    elseif to_bus == 1
                        push!(N1set, from_bus)
                        push!(L1set, (from_bus, to_bus))
                    else
                        Nm1set = union(Nm1set, [from_bus, to_bus])
                        push!(Lm1set, (from_bus, to_bus))
                    end
                else
                    error("Bus1 or Bus2 not specified for a line in BranchData.dss")
                end

                # Extract resistance (r1) and reactance (x1)
                if haskey(branch_info, "r1")
                    rdict[(from_bus, to_bus)] = parse(Float64, branch_info["r1"])
                else
                    rdict[(from_bus, to_bus)] = 0.0  # Default to zero if not specified
                end

                if haskey(branch_info, "x1")
                    xdict[(from_bus, to_bus)] = parse(Float64, branch_info["x1"])
                else
                    xdict[(from_bus, to_bus)] = 0.0  # Default to zero if not specified
                end
            end
        end
    end

    # Calculate total number of buses and branches
    N = length(Nset)
    m = length(Lset)

    # Calculate additional parameters
    N1 = length(N1set)
    Nm1 = length(Nm1set)
    m1 = length(L1set)
    mm1 = length(Lm1set)

    # Create a dictionary with all outputs
    branchData = Dict(
        :Nset => Nset,
        :Lset => Lset,
        :rdict => rdict,
        :xdict => xdict,
        :parent => parent,
        :children => children,
        :N => N,
        :m => m,
        :N1set => N1set,
        :Nm1set => Nm1set,
        :L1set => L1set,
        :Lm1set => Lm1set,
        :N1 => N1,
        :Nm1 => Nm1,
        :m1 => m1,
        :mm1 => mm1
    )

    return branchData

end


end # module
