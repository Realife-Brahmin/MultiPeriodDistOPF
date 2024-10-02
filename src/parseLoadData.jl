# parseLoadData.jl

module parseLoadData

export parse_load_data

function parse_load_data(systemName::String, T::Int)
    # Construct the file path
    filename = joinpath("..", "rawData", systemName, "Loads.dss")

    # Initialize data structures
    Nset = Set{Int}()                 # Set of nodes with loads
    p_L_R = Dict{Int,Float64}()      # Rated active power for each load (kW)
    q_L_R = Dict{Int,Float64}()      # Rated reactive power for each load (kVAR)
    V_minpu = Dict{Int,Float64}()    # Minimum per-unit voltage for each node
    V_maxpu = Dict{Int,Float64}()    # Maximum per-unit voltage for each node

    # Placeholder for load shapes (to be updated later)
    p_L = Dict{Int,Vector{Float64}}()  # Active power profile over time
    q_L = Dict{Int,Vector{Float64}}()  # Reactive power profile over time

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

            # Parse lines starting with "New load."
            if startswith(line, "New load.")
                # Extract parameters from the line
                tokens = split(line)
                load_info = Dict{String,String}()

                for token in tokens[2:end]
                    if occursin("=", token)
                        key, value = split(token, "=")
                        key = strip(key)
                        value = strip(value)
                        load_info[key] = value
                    end
                end

                # Extract bus number (e.g., Bus=2.1)
                if haskey(load_info, "Bus")
                    bus_str = load_info["Bus"]
                    # Extract the integer part before the decimal point
                    bus_parts = split(bus_str, ".")
                    j = parse(Int, bus_parts[1])  # Node number
                    Nset = union(Nset, [j])
                else
                    error("Bus not specified for a load in Loads.dss")
                end

                # Extract rated active power (kw)
                if haskey(load_info, "kw")
                    p_L_R[j] = parse(Float64, load_info["kw"])
                else
                    p_L_R[j] = 0.0  # Default to zero if not specified
                end

                # Extract rated reactive power (kvar)
                if haskey(load_info, "kvar")
                    q_L_R[j] = parse(Float64, load_info["kvar"])
                else
                    q_L_R[j] = 0.0  # Default to zero if not specified
                end

                # Extract minimum per-unit voltage (Vminpu)
                if haskey(load_info, "Vminpu")
                    V_minpu[j] = parse(Float64, load_info["Vminpu"])
                else
                    V_minpu[j] = 0.95  # Default value if not specified
                end

                # Extract maximum per-unit voltage (Vmaxpu)
                if haskey(load_info, "Vmaxpu")
                    V_maxpu[j] = parse(Float64, load_info["Vmaxpu"])
                else
                    V_maxpu[j] = 1.05  # Default value if not specified
                end

                # Initialize load profiles with rated values (to be updated with actual load shapes)
                p_L[j] = [p_L_R[j] for t in 1:T]
                q_L[j] = [q_L_R[j] for t in 1:T]
            end
        end
    end

    # Return the extracted data
    return Nset, p_L_R, q_L_R, V_minpu, V_maxpu, p_L, q_L
end

end # module
