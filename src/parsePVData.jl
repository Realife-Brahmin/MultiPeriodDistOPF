# parsePVData.jl

module parsePVData

export parse_pv_data

function parse_pv_data(systemName::String, T::Int; LoadShape=nothing)

    # get wd: the path of <this> file
    wd = @__DIR__
    # Construct the file path using wd
    filename = joinpath(wd, "..", "rawData", systemName, "PVSystem.dss")

    # Initialize data structures
    Dset = Set{Int}()                    # Set of nodes with PVs
    p_D_R = Dict{Int,Float64}()        # Rated PV active power (Pmpp, kW)
    S_D_R = Dict{Int,Float64}()        # Rated PV apparent power (kVA)
    irrad = Dict{Int,Float64}()        # Irradiance values for each PV
    p_D = Dict{Int,Vector{Float64}}()  # PV active power profile over time

    LoadShapePV = []  # This will hold the T-values after scaling

    # LoadShape default file if LoadShape is not provided
    if LoadShape === nothing
        loadshape_filename = joinpath(wd, "..", "rawData", systemName, "LoadShapePVDefault.dss")
        defaultLoadShape = []

        # Read LoadShapePVDefault.dss
        open(loadshape_filename, "r") do file
            for line in eachline(file)
                # Remove comments and strip whitespace
                line = split(line, "!")[1]
                line = strip(line)
                if isempty(line)
                    continue
                end
                # Parse the load shape values
                if occursin("mult=[", line)
                    shape_str = strip(split(line, "mult=")[2], ['[', ']'])
                    defaultLoadShape = [parse(Float64, val) for val in split(shape_str)]
                end
            end
        end


        # Scale or subsample the load shape based on T
        if T < 24
            # Take the middle T values
            start_idx = Int(floor((24 - T) / 2)) + 1
            LoadShapePV = defaultLoadShape[start_idx:(start_idx+T-1)]
        elseif T > 24
            # Supersample to generate T values
            for i in 1:T
                # Use a basic supersampling method by linearly interpolating between adjacent values
                scaled_idx = (i - 1) * (24 - 1) / (T - 1) + 1
                lower_idx = Int(floor(scaled_idx))
                upper_idx = min(lower_idx + 1, 24)
                frac = scaled_idx - lower_idx
                interpolated_value = (1 - frac) * defaultLoadShape[lower_idx] + frac * defaultLoadShape[upper_idx]
                push!(LoadShapePV, interpolated_value)
            end
        else
            # T == 24, use the default load shape
            LoadShapePV = defaultLoadShape
        end

        # Now set LoadShape to this scaled result for the rest of the function
        LoadShape = LoadShapePV
    end

    # Open and read the PVSystem.dss file
    open(filename, "r") do file
        for line in eachline(file)
            # Remove comments and strip whitespace
            line = split(line, "!")[1]
            line = strip(line)

            # Skip empty lines
            if isempty(line)
                continue
            end

            # Parse lines starting with "New PVSystem."
            if startswith(line, "New PVSystem.")
                # Extract parameters from the line
                tokens = split(line)
                pv_info = Dict{String,String}()

                for token in tokens[2:end]
                    if occursin("=", token)
                        key, value = split(token, "=")
                        key = strip(key)
                        value = strip(value)
                        pv_info[key] = value
                    end
                end

                # Extract bus number (e.g., Bus1=6.1)
                if haskey(pv_info, "Bus1")
                    bus_str = pv_info["Bus1"]
                    # Extract the integer part before the decimal point
                    bus_parts = split(bus_str, ".")
                    j = parse(Int, bus_parts[1])  # Node number
                    Dset = union(Dset, [j])
                else
                    error("Bus1 not specified for a PV in PVSystem.dss")
                end

                # Extract rated active power (Pmpp)
                if haskey(pv_info, "Pmpp")
                    p_D_R[j] = parse(Float64, pv_info["Pmpp"])
                else
                    p_D_R[j] = 0.0  # Default to zero if not specified
                end

                # Extract rated apparent power (kVA)
                if haskey(pv_info, "kVA")
                    S_D_R[j] = parse(Float64, pv_info["kVA"])
                else
                    S_D_R[j] = 0.0  # Default to zero if not specified
                end

                # Extract irradiance (irrad)
                if haskey(pv_info, "irrad")
                    irrad[j] = parse(Float64, pv_info["irrad"])
                else
                    irrad[j] = 1.0  # Default irradiance if not specified
                end

                # Calculate active power profile over time using LoadShape
                p_D[j] = [p_D_R[j] * LoadShape[t] for t in 1:T]
            end
        end
    end

    n_D = length(Dset)

    pvData = Dict(
        :n_D => n_D,
        :Dset => Dset,
        :p_D_R => p_D_R,
        :S_D_R => S_D_R,
        :irrad => irrad,
        :p_D => p_D,
        :LoadShapePV => LoadShapePV  # Add the scaled load shape to the output
    )

    # Return the extracted data as a dictionary
    return pvData

end

end # module
