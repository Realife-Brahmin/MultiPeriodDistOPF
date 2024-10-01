# parseBatteryData.jl

module parseBatteryData

export parse_battery_data

function parse_battery_data(filename::String)
    Bset = Set{Int}()
    battery_params = Dict{Int, Dict{String, Any}}()

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Skip comments and empty lines
            if isempty(strip(line)) || startswith(strip(line), "!") || startswith(strip(line), "//")
                continue
            end
            # Parse battery data
            # Example line: "Storage.1   bus1=5   kW=100   kWh=200"
            tokens = split(strip(line))
            if startswith(tokens[1], "Storage.")
                storage_info = Dict{String, Any}()
                for token in tokens[2:end]
                    if occursin("=", token)
                        key, value = split(token, "=")
                        storage_info[key] = value
                    end
                end
                j = parse(Int, storage_info["bus1"])
                Bset = union(Bset, [j])
                battery_params[j] = storage_info
                # You can extract other parameters like capacity, max charge/discharge rates, etc.
            end
        end
    end

    return Bset, battery_params
end

end # module
