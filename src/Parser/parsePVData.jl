module parsePVData

export parse_pv_data

include("../helperFunctions.jl")
using .helperFunctions: generateLoadShape

# Note: PVsystem and NOT PVSystem will be read (note the lowercase third letter)

#region parse_pv_data
"""
    parse_pv_data(systemName::String, T::Int; N_L=nothing, kVA_B=1000, LoadShape=nothing, filenameLoadShape=nothing)

Parse PV data from the PVsystem.dss file for a given system.

This function reads and parses the PVsystem.dss file for the specified system, extracting relevant PV data and initializing various parameters. 
It handles the extraction of PV properties such as rated active and apparent power, voltage limits, and PV profiles over time.

# Arguments
- `systemName::String`: The name of the system for which to parse PV data.
- `T::Int`: The number of time steps for the PV profile.
- `N_L::Union{Nothing, Int}`: The number of loads, used to calculate the percentage of nodes with PVs (default: nothing).
- `kVA_B::Float64`: The base power in kVA for per-unit calculations (default: 1000).
- `LoadShape::Union{Nothing, Vector{Float64}}`: An optional vector specifying the PV load shape over time. If not provided, it is generated using the `filenameLoadShape` (default: nothing).
- `filenameLoadShape::Union{Nothing, String}`: The filename for the default load shape (default: nothing).

# Returns
- `pvData::Dict`: A dictionary containing the parsed PV data, including:
    - `n_D`: Number of nodes with PVs.
    - `DER_percent`: Percentage of nodes with PVs.
    - `Dset`: Set of nodes with PVs.
    - `p_D_R`: Rated PV active power for each node (kW).
    - `p_D_R_pu`: Rated PV active power for each node in per-unit.
    - `S_D_R`: Rated PV apparent power for each node (kVA).
    - `S_D_R_pu`: Rated PV apparent power for each node in per-unit.
    - `irrad`: Irradiance values for each PV.
    - `p_D`: PV active power profile over time.
    - `p_D_pu`: PV active power profile over time in per-unit.
    - `Vminpu_D`: Minimum per-unit voltage for each node.
    - `Vmaxpu_D`: Maximum per-unit voltage for each node.
    - `LoadShapePV`: The load shape used for the PV profiles.

# Steps
1. **File Path Construction**: Constructs the file path for PVsystem.dss using the system name.
2. **Data Initialization**: Initializes various parameters and dictionaries to store the extracted data.
3. **Load Shape Generation**: Generates a load shape if not provided by the user.
4. **File Reading**: Reads the PVsystem.dss file for the specified system.
5. **Data Extraction**: Extracts PV properties such as rated active and apparent power, voltage limits, and PV profiles over time.
6. **Return Data**: Returns a dictionary containing the parsed PV data.
"""
function parse_pv_data(systemName::String, T::Int;
    N_L = nothing,
    kVA_B = 1000,
    LoadShape=nothing,
    filenameLoadShape=nothing)
    # get wd: the path of <this> file
    wd = @__DIR__
    # Construct the file path for PVsystem.dss using wd
    filename_pv = joinpath(wd, "..", "..", "rawData", systemName, "PVsystem.dss")

    # Initialize data structures for PV systems
    Dset = Set{Int}()                     # Set of nodes with PVs
    p_D_R = Dict{Int,Float64}()            # Rated PV active power (Pmpp, kW)
    p_D_R_pu = Dict{Int,Float64}()         # Rated PV active power [pu]
    S_D_R = Dict{Int,Float64}()            # Rated PV apparent power (kVA)
    S_D_R_pu = Dict{Int,Float64}()         # Rated PV apparent power [pu]
    irrad = Dict{Int,Float64}()            # Irradiance values for each PV
    p_D = Dict{Tuple{Int,Int},Float64}()   # PV active power profile over time
    p_D_pu = Dict{Tuple{Int,Int},Float64}()# PV active power profile over time [pu]
    Vminpu_D = Dict{Int,Float64}()         # PV bus voltage lower limit (artificially created data entry)
    Vmaxpu_D = Dict{Int,Float64}()         # PV bus voltage upper limit (artificially created data entry)

    LoadShapePV = LoadShape 

    # If the user doesn't provide a LoadShape, generate it using the helper function
    if LoadShape === nothing
        LoadShapePV = generateLoadShape(T, filenameLoadShape=filenameLoadShape)
    end

    # Open and read the PVsystem.dss file
    open(filename_pv, "r") do file
        for line in eachline(file)
            # Remove comments and strip whitespace
            line = split(line, "!")[1]
            line = strip(line)

            # Skip empty lines
            if isempty(line)
                continue
            end

            # Parse lines starting with "New PVsystem."
            if startswith(line, "New PVsystem.")
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
                    error("Bus not specified for a PV in PVsystem.dss")
                end

                # Extract rated active power (Pmpp, kW)
                if haskey(pv_info, "Pmpp")
                    p_D_R[j] = parse(Float64, pv_info["Pmpp"])
                    p_D_R_pu[j] = p_D_R[j]/kVA_B
                else
                    p_D_R[j] = 0.0  # Default to zero if not specified
                    p_D_R_pu[j] = p_D_R[j] / kVA_B
                end

                # Extract rated apparent power (kVA)
                if haskey(pv_info, "kVA")
                    S_D_R[j] = parse(Float64, pv_info["kVA"])
                    S_D_R_pu[j] = S_D_R[j]/kVA_B
                else
                    S_D_R[j] = 0.0  # Default to zero if not specified
                    S_D_R_pu[j] = S_D_R[j] / kVA_B
                end

                # Extract irradiance values
                if haskey(pv_info, "irradiance")
                    irrad[j] = parse(Float64, pv_info["irradiance"])
                else
                    irrad[j] = 1.0  # Default to 1.0 if not specified
                end

                # Extract minimum per-unit voltage (Vminpu)
                if haskey(pv_info, "Vminpu")
                    Vminpu_D[j] = parse(Float64, pv_info["Vminpu"])
                else
                    Vminpu_D[j] = 0.95  # Default value if not specified
                end

                # Extract maximum per-unit voltage (Vmaxpu)
                if haskey(pv_info, "Vmaxpu")
                    Vmaxpu_D[j] = parse(Float64, pv_info["Vmaxpu"])
                else
                    Vmaxpu_D[j] = 1.05  # Default value if not specified
                end

                # Initialize PV profiles using the provided or generated LoadShapePV
                for t in 1:T
                    p_D[(j, t)] = p_D_R[j] * LoadShapePV[t]
                    p_D_pu[(j, t)] = p_D_R_pu[j] * LoadShapePV[t]
                end
            end
        end
    end

    n_D = length(Dset)
    DER_percent = 0

    # If the user doesn't provide a N_L, set DER_percent to 100, else compute an actual percentage
    if N_L === nothing
        DER_percent = 100
    else # 1% if it is less than that (but nonzero)
        DER_percent = Int(ceil(n_D/N_L * 100))
    end

    Dset = sort(collect(Dset))
    
    pvData = Dict(
        :n_D => n_D,
        :DER_percent => DER_percent,
        :Dset => Dset,
        :p_D_R => p_D_R,
        :p_D_R_pu => p_D_R_pu,
        :S_D_R => S_D_R,
        :S_D_R_pu => S_D_R_pu,
        :irrad => irrad,
        :p_D => p_D,
        :p_D_pu => p_D_pu,
        :Vminpu_D => Vminpu_D,
        :Vmaxpu_D => Vmaxpu_D,
        :LoadShapePV => LoadShapePV  # Store the LoadShapePV used
    )

    # Return the extracted data as a dictionary
    return pvData
end
#endregion

end # module parsePVData
