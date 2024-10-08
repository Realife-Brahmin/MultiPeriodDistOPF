# parsePVData.jl

module parsePVData

export parse_pv_data

include("helperFunctions.jl")
using .helperFunctions: generateLoadShape

function parse_pv_data(systemName::String, T::Int; LoadShape=nothing, filenameLoadShape=nothing)
    # get wd: the path of <this> file
    wd = @__DIR__
    # Construct the file path using wd
    filename = joinpath(wd, "..", "rawData", systemName, "PVSystem.dss")

    # Initialize data structures
    Dset = Set{Int}()                     # Set of nodes with PVs
    p_D_R = Dict{Int,Float64}()         # Rated PV active power (Pmpp, kW)
    S_D_R = Dict{Int,Float64}()         # Rated PV apparent power (kVA)
    irrad = Dict{Int,Float64}()         # Irradiance values for each PV
    p_D = Dict{Int,Vector{Float64}}()   # PV active power profile over time
    Vminpu_D = Dict{Int,Float64}()   # PV bus voltage lower limit (artifically created data entry)
    Vmaxpu_D = Dict{Int,Float64}()   # PV bus voltage upper limit (artifically created data entry)

    LoadShapePV = LoadShape 

    # If the user doesn't provide a LoadShape, generate it using the helper function
    if LoadShape === nothing
        LoadShapePV = generateLoadShape(T, filenameLoadShape=filenameLoadShape)
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

                # Force bus voltage limits (Vminpu)
                if haskey(pv_info, "Vminpu") # will never be read, as it doesn't exist for the OpenDSS class PVSystem
                    Vminpu_D[j] = parse(Float64, pv_info["Vminpu"])
                else
                    Vminpu_D[j] = 0.95  # Default to 0.95 if not specified
                end

                # Force bus voltage limits (Vmaxpu)
                if haskey(pv_info, "Vmaxpu") # will never be read, as it doesn't exist for the OpenDSS class PVSystem
                    Vmaxpu_D[j] = parse(Float64, pv_info["Vmaxpu"])
                else
                    Vmaxpu_D[j] = 1.05  # Default to 1.05 if not specified
                end

                # Calculate active power profile over time using the user-provided or generated LoadShapePV
                p_D[j] = [p_D_R[j] * LoadShapePV[t] for t in 1:T]
            end
        end
    end

    n_D = length(Dset)

    # Return the extracted data as a dictionary
    pvData = Dict(
        :n_D => n_D,
        :Dset => Dset,
        :p_D_R => p_D_R,
        :S_D_R => S_D_R,
        :irrad => irrad,
        :p_D => p_D,
        :Vminpu_D => Vminpu_D,
        :Vmaxpu_D => Vmaxpu_D,
        :LoadShapePV => LoadShapePV  # Store the LoadShapePV used (either user-supplied or generated)
    )

    return pvData
end

end # module
