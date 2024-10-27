# parseBatteryData.jl

module parseBatteryData

export parse_battery_data

include("../helperFunctions.jl")
import .helperFunctions: myprintln

function parse_battery_data(systemName::String;
    kVA_B = 1000,
    N_L = nothing,
    verbose::Bool=false)
    # Get the working directory of this script
    wd = @__DIR__

    # Construct the full file path for the Storage.dss file
    filename = joinpath(wd, "..", "..", "rawData", systemName, "Storage.dss")

    myprintln(verbose, "Reading Storage file from: $filename")

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
    open(filename, "r") do file
        for line in eachline(file)
            # Remove comments and strip whitespace
            line = split(line, "!")[1]
            line = strip(line)

            # Print the line being processed
            myprintln(verbose, "Processing line: $line")

            # Skip empty lines
            if isempty(line)
                myprintln(verbose, "Skipping empty line.")
                continue
            end

            # Parse lines starting with "New Storage."
            if startswith(line, "New Storage.")
                myprintln(verbose, "Found a New Storage entry.")

                storage_info = Dict{String,String}()
                for m in eachmatch(kv_pattern, line)
                    key = strip(m.captures[1])
                    value = strip(m.captures[2])
                    storage_info[key] = value
                end

                # Print the parsed storage information
                myprintln(verbose, "Parsed storage info: $storage_info")

                # Extract bus number
                if haskey(storage_info, "Bus1")
                    bus_str = storage_info["Bus1"]
                    bus_parts = split(bus_str, ".")
                    bus = parse(Int, bus_parts[1])  # Extract the integer part
                    Bset = union(Bset, [bus])
                    myprintln(verbose, "Battery located at bus $bus.")
                else
                    error("Bus1 not specified for a storage in Storage.dss")
                end

                # Extract rated power and energy values
                P_B_R[bus] = haskey(storage_info, "kWrated") ? parse(Float64, storage_info["kWrated"]) : 0.0
                P_B_R_pu[bus] = P_B_R[bus]/kVA_B
                S_B_R[bus] = haskey(storage_info, "kVA") ? parse(Float64, storage_info["kVA"]) : 0.0
                S_B_R_pu[bus] = S_B_R[bus]/kVA_B
                B_R[bus] = haskey(storage_info, "kWhrated") ? parse(Float64, storage_info["kWhrated"]) : 0.0
                B_R_pu[bus] = B_R[bus]/kVA_B
                myprintln(verbose, "P_B_R: $(P_B_R[bus]), S_B_R: $(S_B_R[bus]), B_R: $(B_R[bus])")

                # Extract efficiencies
                eta_C[bus] = haskey(storage_info, "EffCharge") ? parse(Float64, storage_info["EffCharge"]) / 100 : 0.95
                eta_D[bus] = haskey(storage_info, "EffDischarge") ? parse(Float64, storage_info["EffDischarge"]) / 100 : 0.95
                myprintln(verbose, "eta_C: $(eta_C[bus]), eta_D: $(eta_D[bus])")

                # Extract minimum and maximum SOC
                soc_max[bus] = 0.95
                soc_min[bus] = haskey(storage_info, "reserve") ? parse(Float64, storage_info["reserve"]) / 100 : 0.3
                myprintln(verbose, "soc_max: $(soc_max[bus]), soc_min: $(soc_min[bus])")

                # Extract initial soc percentage (soc_0)
                soc_0[bus] = haskey(storage_info, "stored") ? parse(Float64, storage_info["stored"]) / 100 : 0.625

                # Compute initial SOC (B0)
                B0[bus] = soc_0[bus] * B_R[bus]
                B0_pu[bus] = B0[bus]/kVA_B
                myprintln(verbose, "Initial SOC B0: $(B0[bus])")

                # Extract voltage limits
                Vminpu_B[bus] = haskey(storage_info, "Vminpu") ? parse(Float64, storage_info["Vminpu"]) : 0.95
                Vmaxpu_B[bus] = haskey(storage_info, "Vmaxpu") ? parse(Float64, storage_info["Vmaxpu"]) : 1.05
                myprintln(verbose, "Vminpu: $(Vminpu_B[bus]), Vmaxpu: $(Vmaxpu_B[bus])")
            end
        end
    end

    # Compute the cardinality of Bset
    n_B = length(Bset)
    # If the user doesn't provide a N_L, Set Batt_percent to 100, else compute an actual percentage
    if N_L === nothing
        Batt_percent = 100
    else # 1% if it is less than that (but nonzero)
        Batt_percent = Int(ceil(n_B / N_L * 100))
    end

    # By default, ensuring that batteries always get back to their original SOCs at the end of the optimization horizon
    Bref = B0 
    Bref_pu = B0_pu

    Bset = sort(collect(Bset))

    # Create a dictionary with all outputs
    storageData = Dict(
        :Bset => Bset,
        :n_B => n_B,
        :Batt_percent => Batt_percent,
        :B0 => B0,
        :B0_pu => B0_pu,
        :Bref => Bref,
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

    myprintln(verbose, "Final parsed storage data: $storageData")

    return storageData
end

end # module
