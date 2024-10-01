# parseBusData.jl

module parseBusData

export parse_bus_data

function parse_bus_data(filename::String)
    # Initialize variables
    Nset = Set{Int}()
    bus_data = Dict{Int, Dict{String, Any}}()

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Skip comments and empty lines
            if isempty(strip(line)) || startswith(strip(line), "!") || startswith(strip(line), "//")
                continue
            end
            # Parse bus data
            # Example line: "Bus.1   kv=12.47   pu=1.0"
            tokens = split(strip(line))
            if startswith(tokens[1], "Bus.")
                bus_id = parse(Int, split(tokens[1], ".")[2])
                Nset = union(Nset, [bus_id])
                # Extract other bus parameters if needed
                bus_info = Dict{String, Any}()
                for token in tokens[2:end]
                    if occursin("=", token)
                        key, value = split(token, "=")
                        bus_info[key] = parse(Float64, value)
                    end
                end
                bus_data[bus_id] = bus_info
            end
        end
    end

    N = maximum(Nset)
    Nm1set = setdiff(Nset, [1])  # Assuming bus 1 is the substation
    return N, Nset, Nm1set, bus_data
end

end # module
