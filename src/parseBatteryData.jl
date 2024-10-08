# parseBatteryData.jl

module parseBatteryData

export parse_battery_data

function parse_battery_data(systemName::String)
    # Get the working directory of this script
    wd = @__DIR__

    # Construct the full file path for the Storage.dss file
    filename = joinpath(wd, "..", "rawData", systemName, "Storage.dss")

    # Initialize output data structures
    Bset = Set{Int}()  # Set of bus numbers with batteries
    B0 = Dict{Int,Float64}()  # Initial SOC (kWhr) for each battery
    B_R = Dict{Int,Float64}()  # Rated storage capacity for each battery (kWhr)
    eta_C = Dict{Int,Float64}()  # Charging efficiency for each battery
    eta_D = Dict{Int,Float64}()  # Discharging efficiency for each battery
    P_B_R = Dict{Int,Float64}()  # Rated charging/discharging power (kW)
    S_B_R = Dict{Int,Float64}()  # Rated inverter apparent power (kVA)
    Vminpu_B = Dict{Int,Float64}()  # Minimum voltage for each battery (pu)
    Vmaxpu_B = Dict{Int,Float64}()  # Maximum voltage for each battery (pu)
    soc_min = Dict{Int,Float64}()  # Minimum SOC (%reserve)
    soc_max = Dict{Int,Float64}()  # Maximum SOC (%stored)

    # Regular expression to match key-value pairs
    kv_pattern = r"(\w+)\s*=\s*([\S]+)"

    # Read and parse the Storage.dss file
    open(filename, "r") do file
        for line in eachline(file)
            # Remove comments and strip whitespace
            line = split(line, "!")[1]
            line = strip(line)

            # Skip empty lines
            if isempty(line)
                continue
            end

            # Parse lines starting with "New Storage."
            if startswith(line, "New Storage.")
                storage_info = Dict{String,String}()
                for m in eachmatch(kv_pattern, line)
                    key = strip(m.captures[1])
                    value = strip(m.captures[2])
                    storage_info[key] = value
                end

                # Extract bus number
                if haskey(storage_info, "Bus1")
                    bus_str = storage_info["Bus1"]
                    bus_parts = split(bus_str, ".")
                    bus = parse(Int, bus_parts[1])  # Extract the integer part
                    Bset = union(Bset, [bus])
                else
                    error("Bus1 not specified for a storage in Storage.dss")
                end

                # Extract rated power and energy values
                P_B_R[bus] = haskey(storage_info, "kWrated") ? parse(Float64, storage_info["kWrated"]) : 0.0
                S_B_R[bus] = haskey(storage_info, "kVA") ? parse(Float64, storage_info["kVA"]) : 0.0
                B_R[bus] = haskey(storage_info, "kWhrated") ? parse(Float64, storage_info["kWhrated"]) : 0.0

                # Extract efficiencies
                eta_C[bus] = haskey(storage_info, "%EffCharge") ? parse(Float64, storage_info["%EffCharge"]) / 100 : 1.0
                eta_D[bus] = haskey(storage_info, "%EffDischarge") ? parse(Float64, storage_info["%EffDischarge"]) / 100 : 1.0

                # Extract minimum and maximum SOC
                soc_max[bus] = haskey(storage_info, "%stored") ? parse(Float64, storage_info["%stored"]) / 100 : 1.0
                soc_min[bus] = haskey(storage_info, "%reserve") ? parse(Float64, storage_info["%reserve"]) / 100 : 0.0

                # Compute initial SOC (B0)
                B0[bus] = soc_min[bus] * B_R[bus]

                # Extract voltage limits
                Vminpu_B[bus] = haskey(storage_info, "Vminpu") ? parse(Float64, storage_info["Vminpu"]) : 0.95
                Vmaxpu_B[bus] = haskey(storage_info, "Vmaxpu") ? parse(Float64, storage_info["Vmaxpu"]) : 1.05
            end
        end
    end

    # Compute the cardinality of Bset
    n_B = length(Bset)

    # Create a dictionary with all outputs
    storageData = Dict(
        :Bset => Bset,
        :n_B => n_B,
        :B0 => B0,
        :B_R => B_R,
        :eta_C => eta_C,
        :eta_D => eta_D,
        :P_B_R => P_B_R,
        :S_B_R => S_B_R,
        :Vminpu_B => Vminpu_B,
        :Vmaxpu_B => Vmaxpu_B,
        :soc_min => soc_min,
        :soc_max => soc_max
    )

    return storageData
end


end # module
