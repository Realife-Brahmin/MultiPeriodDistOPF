module parseLoadData

export parse_load_data

include("helperFunctions.jl")
using .helperFunctions: generateLoadShape

function parse_load_data(systemName::String, T::Int;
    LoadShape=nothing, 
    filenameLoadShape::String="LoadShapeDefault.dss")

    # get wd: the path of <this> file
    wd = @__DIR__
    # Construct the file path for Loads.dss using wd
    filename_load = joinpath(wd, "..", "rawData", systemName, "Loads.dss")

    # Initialize data structures
    NLset = Set{Int}()                 # Set of nodes with loads
    p_L_R = Dict{Int,Float64}()      # Rated active power for each load (kW)
    q_L_R = Dict{Int,Float64}()      # Rated reactive power for each load (kVAR)
    Vminpu_L = Dict{Int,Float64}()    # Minimum per-unit voltage for each node
    Vmaxpu_L = Dict{Int,Float64}()    # Maximum per-unit voltage for each node

    # Placeholder for load shapes (to be updated later)
    p_L = Dict{Int,Vector{Float64}}()  # Active power profile over time
    q_L = Dict{Int,Vector{Float64}}()  # Reactive power profile over time

    # If user does not provide a LoadShape, generate one using the helper function
    if LoadShape === nothing
        LoadShape = generateLoadShape(T, filenameLoadShape=filenameLoadShape)
    end

    # Open and read the Loads.dss file
    open(filename_load, "r") do file
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
                    NLset = union(NLset, [j])
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
                    Vminpu_L[j] = parse(Float64, load_info["Vminpu"])
                else
                    Vminpu_L[j] = 0.95  # Default value if not specified
                end

                # Extract maximum per-unit voltage (Vmaxpu)
                if haskey(load_info, "Vmaxpu")
                    Vmaxpu_L[j] = parse(Float64, load_info["Vmaxpu"])
                else
                    Vmaxpu_L[j] = 1.05  # Default value if not specified
                end

                # Initialize load profiles using the provided or generated LoadShape
                p_L[j] = [p_L_R[j] * LoadShape[t] for t in 1:T]
                q_L[j] = [q_L_R[j] * LoadShape[t] for t in 1:T]
            end
        end
    end

    N_L = length(NLset)

    loadData = Dict(
        :N_L => N_L,
        :NLset => NLset,
        :p_L_R => p_L_R,
        :q_L_R => q_L_R,
        :Vminpu_L => Vminpu_L,
        :Vmaxpu_L => Vmaxpu_L,
        :p_L => p_L,
        :q_L => q_L,
        :LoadShape => LoadShape  # Store the LoadShape used
    )

    # Return the extracted data as a dictionary
    return loadData
end

end
