# parseLoadData.jl

module parseLoadData

export parse_load_data

function parse_load_data(filename::String, T::Int)
    Nset = Set{Int}()
    p_L = Dict{Int,Vector{Float64}}()
    q_L = Dict{Int,Vector{Float64}}()

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Skip comments and empty lines
            if isempty(strip(line)) || startswith(strip(line), "!") || startswith(strip(line), "//")
                continue
            end
            # Parse load data
            # Example line: "New Load.Load1   Bus1=2   kW=100   kvar=50"
            tokens = split(strip(line))
            if occursin("Load.", tokens[2]) || occursin("Load.", tokens[1])
                load_info = Dict{String,Any}()
                for token in tokens
                    if occursin("=", token)
                        key, value = split(token, "=")
                        load_info[key] = value
                    end
                end
                # Get bus number
                if haskey(load_info, "Bus1")
                    bus = load_info["Bus1"]
                    # Extract bus number (assuming format like "2" or "Bus2")
                    if occursin(".", bus)
                        # If format is "Bus.2", extract the number
                        bus_id = parse(Int, split(bus, ".")[2])
                    else
                        bus_id = parse(Int, bus)
                    end
                    Nset = union(Nset, [bus_id])

                    # Extract kW and kvar
                    p_kw = parse(Float64, load_info["kW"])
                    q_kvar = parse(Float64, load_info["kvar"])

                    # Assume load is constant over time (or modify as needed)
                    p_L[bus_id] = [p_kw for t in 1:T]
                    q_L[bus_id] = [q_kvar for t in 1:T]
                end
            end
        end
    end

    # Identify the substation node (assuming node 1)
    SubstationNode = 1
    Nset = union(Nset, [SubstationNode])

    N = maximum(Nset)
    Nm1set = setdiff(Nset, [SubstationNode])

    # Ensure all nodes have p_L and q_L (zero if not specified)
    for node in Nset
        if !haskey(p_L, node)
            p_L[node] = [0.0 for t in 1:T]
            q_L[node] = [0.0 for t in 1:T]
        end
    end

    return N, Nset, Nm1set, p_L, q_L
end

end # module
