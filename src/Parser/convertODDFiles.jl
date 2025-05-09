module openDSSFilesGenerator

export
    convert_files_loads_dss,
    generate_battery_dss_from_loads_dss,
    generate_pvsystem_dss_from_loads_dss

include("../helperFunctions.jl")
import .helperFunctions as HF

using CSV
using DataFrames

function convert_files_loads_dss(input_file_path::String, output_file_path::String; removeEmptyLoads::Bool=true)
    # Open the input file for reading
    input_file = open(input_file_path, "r")
    # Open the output file for writing
    output_file = open(output_file_path, "w")

    # Process each line in the input file
    for line in eachline(input_file)
        # Skip empty lines
        if isempty(line)
            continue
        end

        # Parse the line to extract relevant information
        pattern = r"New Load\.L(\d+) bus1=(\d+) model=(\d+) phases=(\d+) conn=wye kv=(\d+\.\d+) kW=(\d+\.\d+) kvar=(\d+\.\d+) Vminpu=(\d+\.\d+) Vmaxpu=(\d+\.\d+)"
        match_data = match(pattern, line)
        if match_data !== nothing
            bus_number = match_data.captures[2]
            kV = match_data.captures[5]
            kW = parse(Float64, match_data.captures[6])
            kvar = parse(Float64, match_data.captures[7])

            # Check if the load is non-empty or if removeEmptyLoads is disabled
            if !removeEmptyLoads || (kW != 0.0 || kvar != 0.0)
                # Create the new line in the desired format
                new_line = "New load.S$(bus_number)P Bus=$(bus_number).1 Phases=1 Model=1 kV=$(kV) kW=$(kW) kVar=$(kvar) Vminpu=0.95 Vmaxpu=1.05 Daily=LoadShapeLoadDefault"
                println(output_file, new_line)
            end
        end
    end

    # Close the files
    close(input_file)
    close(output_file)
end

username = ENV["USERNAME"]
# # Example usage
# input_file_path = "C:/Users/" * username * "/Documents/documents_general/MultiPeriodDistOPF/rawData/OpenDSS730node/Loads.dss"

# output_file_path = "C:/Users/" * username * "/Documents/documents_general/MultiPeriodDistOPF/rawData/ieee730_1ph/Loads.dss"
# convert_files_loads_dss(input_file_path, output_file_path)

function extract_parameter(parts::Vector{String}, param::String)
    for part in parts
        if startswith(part, param)
            return split(part, "=")[2]
        end
    end
    return nothing
end

function generate_battery_dss_from_loads_dss(loads_file_path::String; Batt_percent=30.0, Batt_rating_factor=1.0)
    # Read the Loads.dss file
    lines = readlines(loads_file_path)

    # Extract non-zero loads
    load_buses = []
    for line in lines
        if startswith(line, "New load.")
            parts = split(line)
            parts = map(String, parts)  # Convert SubString elements to String
            bus = extract_parameter(parts, "Bus")
            bus = split(bus, ".")[1]  # Take the integer part of the bus number
            kV = parse(Float64, extract_parameter(parts, "kV"))
            kW = parse(Float64, extract_parameter(parts, "kW"))
            if kW > 0
                push!(load_buses, (bus, kV, kW))
            end
        end
    end

    # Calculate the number of batteries
    N_L = length(load_buses)
    num_batteries = ceil(Int, Batt_percent / 100 * N_L)
    actual_Batt_percent = ceil(Int, num_batteries / N_L * 100)

    # Select evenly spaced buses for battery placement
    selected_buses = [load_buses[i] for i in 1:div(N_L, num_batteries):N_L][1:num_batteries]

    # Generate the Storage<BattPercent>.dss file
    batt_rating_factor_str = HF.trim_number_for_printing(Batt_rating_factor * 100, sigdigits=4)
    storage_file_path = joinpath(dirname(loads_file_path), "Storage_$(actual_Batt_percent)_$(batt_rating_factor_str).dss")
    open(storage_file_path, "w") do file
        for (i, (bus, kV, kW)) in enumerate(selected_buses)
            kWrated = Batt_rating_factor * kW
            kVA = 1.2 * kWrated  # 120% of the load
            kWhrated = 4 * kVA
            println(file, "New Storage.Battery$(bus) Phases=1 Bus1=$(bus) kV=$(kV) kVA=$(round(kVA, digits=4)) kWrated=$(round(kWrated, digits=4)) kWhrated=$(round(kWhrated, digits=4)) %stored=62.5 %reserve=30 %EffCharge=95 %EffDischarge=95 %idlingkW=0 Vminpu=0.95 Vmaxpu=1.05 DispMode=External")
        end
    end

    println("Generated $storage_file_path with $num_batteries batteries.")
end

function generate_pvsystem_dss_from_loads_dss(loads_file_path::String; DER_Percent=20.0, DER_rating_factor=1.0)
    # Read the Loads.dss file
    lines = readlines(loads_file_path)

    # Extract non-zero loads
    load_buses = []
    for line in lines
        if startswith(line, "New load.")
            parts = split(line)
            parts = map(String, parts)  # Convert SubString elements to String
            bus = extract_parameter(parts, "Bus")
            bus = split(bus, ".")[1]  # Take the integer part of the bus number
            kV = parse(Float64, extract_parameter(parts, "kV"))
            kW = parse(Float64, extract_parameter(parts, "kW"))
            if kW > 0
                push!(load_buses, (bus, kV, kW))
            end
        end
    end

    # Calculate the number of PV systems
    N_L = length(load_buses)
    num_pvsystems = ceil(Int, DER_Percent / 100 * N_L)
    actual_DER_percent = ceil(Int, num_pvsystems / N_L * 100)

    # Select evenly spaced buses for PV system placement
    selected_buses = [load_buses[i] for i in 1:div(N_L, num_pvsystems):N_L][1:num_pvsystems]

    # Generate the PVSystem<DER_Percent>.dss file
    der_rating_factor_str = HF.trim_number_for_printing(DER_rating_factor * 100, sigdigits=2)
    pvsystem_file_path = joinpath(dirname(loads_file_path), "PVSystem_$(actual_DER_percent)_$(der_rating_factor_str).dss")
    open(pvsystem_file_path, "w") do file
        for (i, (bus, kV, kW)) in enumerate(selected_buses)
            Pmpp = DER_rating_factor * kW
            kVA = 1.2 * Pmpp
            println(file, "New PVsystem.PV$(bus) irrad=1.0 Phases=1 Bus1=$(bus).1 kV=$(kV) kVA=$(round(kVA, digits=4)) Pmpp=$(round(Pmpp, digits=4)) %cutin=0.001 %cutout=0.001")
        end
    end

    println("Generated $pvsystem_file_path with $num_pvsystems PV systems.")
end

# Example usage:
generate_battery_dss_from_loads_dss("C:/Users/" * username * "/Documents/documents_general/MultiPeriodDistOPF/rawData/ieee123_1ph/Loads.dss", Batt_percent=50, Batt_rating_factor=1.0)
generate_pvsystem_dss_from_loads_dss("C:/Users/" * username * "/Documents/documents_general/MultiPeriodDistOPF/rawData/ieee123_1ph/Loads.dss", DER_Percent=40, DER_rating_factor=1.0)

end # module openDSSFilesGenerator