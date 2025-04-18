module parseBatteryData

export parse_battery_data

using Parameters

include("../helperFunctions.jl")
import .helperFunctions as HF

#region parse_battery_data
"""
    parse_battery_data(systemName::String; kVA_B=1000, N_L=nothing, verbose::Bool=false)

Parse battery data from the Storage.dss file for a given system.

This function reads and parses the Storage.dss file for the specified system, extracting relevant battery data and initializing various parameters. 
It handles the extraction of battery properties such as initial state of charge (SOC), rated power, efficiencies, and voltage limits.

# Arguments
- `systemName::String`: The name of the system for which to parse battery data.
- `kVA_B::Float64`: The base power in kVA for per-unit calculations (default: 1000).
- `N_L::Union{Nothing, Int}`: The number of loads, used to calculate the percentage of buses with batteries (default: nothing).
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `storageData::Dict`: A dictionary containing the parsed battery data, including:
    - `Bset`: Set of bus numbers with batteries.
    - `n_B`: Number of buses with batteries.
    - `Batt_percent`: Percentage of buses with batteries.
    - `B0`: Initial SOC (kWhr) for each battery.
    - `B0_pu`: Initial SOC (kWhr) for each battery in per-unit.
    - `Bref`: Reference SOC (kWhr) for each battery.
    - `Bref_percent`: Reference SOC as a percentage of rated capacity.
    - `Bref_pu`: Reference SOC (kWhr) for each battery in per-unit.
    - `B_R`: Rated storage capacity for each battery (kWhr).
    - `B_R_pu`: Rated storage capacity for each battery in per-unit.
    - `eta_C`: Charging efficiency for each battery.
    - `eta_D`: Discharging efficiency for each battery.
    - `P_B_R`: Rated charging/discharging power (kW).
    - `P_B_R_pu`: Rated charging/discharging power in per-unit.
    - `S_B_R`: Rated inverter apparent power (kVA).
    - `S_B_R_pu`: Rated inverter apparent power in per-unit.
    - `Vminpu_B`: Minimum voltage for each battery (pu).
    - `Vmaxpu_B`: Maximum voltage for each battery (pu).
    - `soc_min`: Minimum SOC (% reserve).
    - `soc_max`: Maximum SOC (% maximum allowed).
    - `soc_0`: Initial SOC (% stored).

# Steps
1. **File Reading**: Reads the Storage.dss file for the specified system.
2. **Data Extraction**: Extracts battery properties such as initial SOC, rated power, efficiencies, and voltage limits.
3. **Data Initialization**: Initializes various parameters and dictionaries to store the extracted data.
4. **Verbose Output**: Optionally prints detailed information about the parsing process if `verbose` is true.
5. **Return Data**: Returns a dictionary containing the parsed battery data.
"""
function parse_battery_data(systemName::String;
    N_L = nothing,
    verbose::Bool=false,
    gedDict_ud=nothing,
    baseValuesDict=nothing)

    if isnothing(baseValuesDict)
        error("baseValuesDict must be provided")
    else
        @unpack kVA_B_dict, kV_B_dict, MVA_B_dict = baseValuesDict;
    end
    # Get the working directory of this script
    wd = @__DIR__

    if isnothing(gedDict_ud)
        gedDict_ud = Dict(:DER_Percent_ud => 20, :DER_Rating_factor_ud => 1, :Batt_Percent_ud => 30, :Batt_Rating_factor_ud => 1)
    end
    @unpack DER_Percent_ud, DER_Rating_factor_ud, Batt_Percent_ud, Batt_Rating_factor_ud = gedDict_ud

    # If the user doesn't provide a N_L, set Batt_percent to 100, else compute an actual percentage
    if N_L === nothing
        Batt_percent = 100
    else # 1% if it is less than that (but nonzero)
        n_B = Int(ceil(Batt_Percent_ud / 100 * N_L))
        Batt_percent = Int(ceil(n_B / N_L * 100))
    end

    Batt_Rating_factor_str = HF.trim_number_for_printing(Batt_Rating_factor_ud * 100, sigdigits=4)

    # Construct the full file path for the Storage.dss file
    filename_batt = joinpath(wd, "..", "..", "rawData", systemName, "Storage_$(Batt_percent)_$(Batt_Rating_factor_str).dss")

    HF.myprintln(verbose, "Reading Storage file from: $filename_batt")

    # Initialize output data structures
    Bset = Set{Int}()  # Set of bus numbers with batteries
    B0 = Dict{Int,Float64}()  # Initial SOC (kWhr) for each battery
    B0_pu = Dict{Int,Float64}()  # Initial SOC (kWhr) for each battery in [pu h]
    B_R = Dict{Int,Float64}()  # Rated storage capacity for each battery (kWhr)
    B_R_pu = Dict{Int,Float64}()  # Rated storage capacity for each battery [pu h]
    eta_C = Dict{Int,Float64}()  # Charging efficiency for each battery
    eta_D = Dict{Int,Float64}()  # Discharging efficiency for each battery
    P_B_R = Dict{Int,Float64}()  # Rated charging/discharging power (kW)
    P_B_R_pu = Dict{Int,Float64}()  # Rated charging/discharging power [pu]
    S_B_R = Dict{Int,Float64}()  # Rated inverter apparent power (kVA)
    S_B_R_pu = Dict{Int,Float64}()  # Rated inverter apparent power (pu)
    Vminpu_B = Dict{Int,Float64}()  # Minimum voltage for each battery (pu)
    Vmaxpu_B = Dict{Int,Float64}()  # Maximum voltage for each battery (pu)
    soc_min = Dict{Int,Float64}()  # Minimum SOC (%reserve)
    soc_max = Dict{Int,Float64}()  # Maximum SOC (%maximum allowed, (not a real thing in OpenDSS Storage modelling))
    soc_0 = Dict{Int,Float64}() # Initial SOC (%stored)

    # Regular expression to match key-value pairs
    kv_pattern = r"(\w+)\s*=\s*([\S]+)"

    # Read and parse the Storage.dss file
    open(filename_batt, "r") do file
        for line in eachline(file)
            # Remove comments and strip whitespace
            line = split(line, "!")[1]
            line = strip(line)

            # Print the line being processed
            HF.myprintln(verbose, "Processing line: $line")

            # Skip empty lines
            if isempty(line)
                HF.myprintln(verbose, "Skipping empty line.")
                continue
            end

            # Parse lines starting with "New Storage."
            if startswith(line, "New Storage.")
                HF.myprintln(verbose, "Found a New Storage entry.")

                storage_info = Dict{String,String}()
                for m in eachmatch(kv_pattern, line)
                    key = strip(m.captures[1])
                    value = strip(m.captures[2])
                    storage_info[key] = value
                end

                # Print the parsed storage information
                HF.myprintln(verbose, "Parsed storage info: $storage_info")

                # Extract bus number
                if haskey(storage_info, "Bus1")
                    bus_str = storage_info["Bus1"]
                    bus_parts = split(bus_str, ".")
                    bus = parse(Int, bus_parts[1])  # Extract the integer part
                    Bset = union(Bset, [bus])
                    HF.myprintln(verbose, "Battery located at bus $bus.")
                else
                    error("Bus1 not specified for a storage in Storage.dss")
                end

                # Extract rated power and energy values
                P_B_R[bus] = haskey(storage_info, "kWrated") ? parse(Float64, storage_info["kWrated"]) : 0.0
                P_B_R_pu[bus] = P_B_R[bus]/kVA_B_dict[bus]
                S_B_R[bus] = haskey(storage_info, "kVA") ? parse(Float64, storage_info["kVA"]) : 0.0
                S_B_R_pu[bus] = S_B_R[bus]/kVA_B_dict[bus]
                B_R[bus] = haskey(storage_info, "kWhrated") ? parse(Float64, storage_info["kWhrated"]) : 0.0
                B_R_pu[bus] = B_R[bus]/kVA_B_dict[bus]
                HF.myprintln(verbose, "P_B_R: $(P_B_R[bus]), S_B_R: $(S_B_R[bus]), B_R: $(B_R[bus])")

                # Extract efficiencies
                eta_C[bus] = haskey(storage_info, "EffCharge") ? parse(Float64, storage_info["EffCharge"]) / 100 : 0.95
                eta_D[bus] = haskey(storage_info, "EffDischarge") ? parse(Float64, storage_info["EffDischarge"]) / 100 : 0.95
                HF.myprintln(verbose, "eta_C: $(eta_C[bus]), eta_D: $(eta_D[bus])")

                # Extract minimum and maximum SOC
                soc_max[bus] = 0.95
                soc_min[bus] = haskey(storage_info, "reserve") ? parse(Float64, storage_info["reserve"]) / 100 : 0.3
                HF.myprintln(verbose, "soc_max: $(soc_max[bus]), soc_min: $(soc_min[bus])")

                # Extract initial soc percentage (soc_0)
                soc_0[bus] = haskey(storage_info, "stored") ? parse(Float64, storage_info["stored"]) / 100 : 0.625

                # Compute initial SOC (B0)
                B0[bus] = soc_0[bus] * B_R[bus]
                B0_pu[bus] = B0[bus]/kVA_B_dict[bus]
                HF.myprintln(verbose, "Initial SOC B0: $(B0[bus])")

                # Extract voltage limits
                Vminpu_B[bus] = haskey(storage_info, "Vminpu") ? parse(Float64, storage_info["Vminpu"]) : 0.90
                Vmaxpu_B[bus] = haskey(storage_info, "Vmaxpu") ? parse(Float64, storage_info["Vmaxpu"]) : 1.10
                HF.myprintln(verbose, "Vminpu: $(Vminpu_B[bus]), Vmaxpu: $(Vmaxpu_B[bus])")
            end
        end
    end

    # Compute the cardinality of Bset
    n_B = length(Bset)

    # By default, ensuring that batteries always get back to their original SOCs at the end of the optimization horizon
    Bref = B0 
    Bref_pu = B0_pu
    Bref_percent = Dict(j => Bref_pu[j] / B_R_pu[j] for j ∈ Bset)

    Bset = sort(collect(Bset))

    # Create a dictionary with all outputs
    storageData = Dict(
        :Bset => Bset,
        :n_B => n_B,
        :Batt_percent => Batt_percent,
        :B0 => B0,
        :B0_pu => B0_pu,
        :Bref => Bref,
        :Bref_percent => Bref_percent,
        :Bref_pu => Bref_pu,
        :B_R => B_R,
        :B_R_pu => B_R_pu,
        :eta_C => eta_C,
        :eta_D => eta_D,
        :P_B_R => P_B_R,
        :P_B_R_pu => P_B_R_pu,
        :S_B_R => S_B_R,
        :S_B_R_pu => S_B_R_pu,
        :Vminpu_B => Vminpu_B,
        :Vmaxpu_B => Vmaxpu_B,
        :soc_min => soc_min,
        :soc_max => soc_max,
        :soc_0 => soc_0
    )

    HF.myprintln(verbose, "Final parsed storage data: $storageData")

    return storageData
end
#endregion

end # module parseBatteryData
