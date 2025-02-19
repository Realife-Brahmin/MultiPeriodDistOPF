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

# Example usage
input_file_path = "C:/Users/aryan/Documents/documents_general/MultiPeriodDistOPF/rawData/OpenDSS730node/Loads.dss"
output_file_path = "C:/Users/aryan/Documents/documents_general/MultiPeriodDistOPF/rawData/ieee730_1ph/Loads.dss"
convert_files_loads_dss(input_file_path, output_file_path)