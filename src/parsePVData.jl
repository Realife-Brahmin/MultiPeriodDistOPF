# parsePVData.jl

module parsePVData

export parse_pv_data

function parse_pv_data(filename::String, T::Int)
    Dset = Set{Int}()
    p_D = Dict{Int, Vector{Float64}}()

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Skip comments and empty lines
            if isempty(strip(line)) || startswith(strip(line), "!") || startswith(strip(line), "//")
                continue
            end
            # Parse PV system data
            # Example line: "PVSystem.1   bus1=4   pmpp=50"
            tokens = split(strip(line))
            if startswith(tokens[1], "PVSystem.")
                pv_info = Dict{String, Any}()
                for token in tokens[2:end]
                    if occursin("=", token)
                        key, value = split(token, "=")
                        pv_info[key] = value
                    end
                end
                j = parse(Int, pv_info["bus1"])
                Dset = union(Dset, [j])
                # Extract pmpp (maximum power point)
                pmpp = parse(Float64, pv_info["pmpp"])
                # Assume PV generation profile (you can replace this with actual data)
                # For now, we assume a simple sinusoidal generation over the day
                p_D[j] = [pmpp * sin(π * (t / T)) for t in 1:T]
            end
        end
    end

    return Dset, p_D
end

end # module
