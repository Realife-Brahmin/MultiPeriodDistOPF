# parseBranchData.jl

module parseBranchData

export parse_branch_data

function parse_branch_data(filename::String)
    Lset = Set{Tuple{Int, Int}}()
    r = Dict{Tuple{Int, Int}, Float64}()
    x = Dict{Tuple{Int, Int}, Float64}()
    Parent = Dict{Int, Int}()
    Children = Dict{Int, Vector{Int}}()

    Nodes = Set{Int}()

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Skip comments and empty lines
            if isempty(strip(line)) || startswith(strip(line), "!") || startswith(strip(line), "//")
                continue
            end
            # Parse branch data
            # Example line: "Line.1   bus1=1   bus2=2   r=0.01   x=0.02"
            tokens = split(strip(line))
            if startswith(tokens[1], "Line.")
                branch_info = Dict{String, Any}()
                for token in tokens[2:end]
                    if occursin("=", token)
                        key, value = split(token, "=")
                        branch_info[key] = value
                    end
                end
                # Get bus indices
                i = parse(Int, branch_info["bus1"])
                j = parse(Int, branch_info["bus2"])
                Nodes = union(Nodes, [i, j])
                Lset = union(Lset, [(i, j)])
                # Get resistance and reactance
                r[(i, j)] = parse(Float64, branch_info["r"])
                x[(i, j)] = parse(Float64, branch_info["x"])
                # Build Parent and Children mappings
                Parent[j] = i
                if haskey(Children, i)
                    push!(Children[i], j)
                else
                    Children[i] = [j]
                end
            end
        end
    end

    # Ensure all nodes have an entry in Children
    for node in Nodes
        if !haskey(Children, node)
            Children[node] = Int[]
        end
    end

    # Define L1set and Lm1set
    SubstationNode = 1  # Assuming node 1 is the substation
    L1set = Set{Tuple{Int, Int}}()
    for (i, j) in Lset
        if i == SubstationNode
            L1set = union(L1set, [(i, j)])
        end
    end
    Lm1set = setdiff(Lset, L1set)

    return Lset, L1set, Lm1set, r, x, Parent, Children
end

end # module
