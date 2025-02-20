#!/usr/bin/env julia

"""
Reads `linedata.txt` and `linenew.dss` to compare R/X values in ohms vs. per-unit,
then outputs a DataFrame with columns: i, j, r1_factor, x1_factor.

r1_factor = r_pu / r_ohm
x1_factor = x_pu / x_ohm
"""

using Printf
using DataFrames

# Helper function to ensure (bus1, bus2) are in ascending numerical order
function sort_bus_pair(bus1::AbstractString, bus2::AbstractString)
    b1 = parse(Int, String(bus1))
    b2 = parse(Int, String(bus2))
    return b1 < b2 ? (bus1, bus2) : (bus2, bus1)
end

function parse_linedata(linedata_path::String)
    d = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    for line in eachline(linedata_path)
        l = strip(line)
        if isempty(l) || startswith(l, "//")
            continue
        end
        parts = split(l)
        if length(parts) < 4
            continue
        end
        i_bus, j_bus = parts[1], parts[2]
        r_pu = parse(Float64, parts[3])
        x_pu = parse(Float64, parts[4])
        bus_pair = sort_bus_pair(i_bus, j_bus)
        d[bus_pair] = (r_pu, x_pu)
    end
    return d
end

function parse_linenew(linenew_path::String)
    pattern = r"bus1\s*=\s*(\S+).*?bus2\s*=\s*(\S+).*?r1\s*=\s*([0-9.]+).*?x1\s*=\s*([0-9.]+)"
    d = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    for line in eachline(linenew_path)
        l = strip(line)
        if !startswith(l, "New Line.")
            continue
        end
        m = match(pattern, l)
        if m !== nothing
            i_bus = m.captures[1]
            j_bus = m.captures[2]
            r_ohm = parse(Float64, m.captures[3])
            x_ohm = parse(Float64, m.captures[4])
            bus_pair = sort_bus_pair(i_bus, j_bus)
            d[bus_pair] = (r_ohm, x_ohm)
        end
    end
    return d
end

function main()
    linedata_file = "./rawData/powerflowcomparision/Node730DG10/BFM/linedata.txt"
    linenew_file = "./rawData/powerflowcomparision/Node730DG10/Opendss/linenew.dss"

    d_data = parse_linedata(linedata_file)
    d_new = parse_linenew(linenew_file)

    all_keys = union(keys(d_data), keys(d_new))

    # Sort by numerical bus ID
    function bus_key(x)
        (parse(Int, x[1]), parse(Int, x[2]))
    end
    sorted_keys = sort(collect(all_keys), by=bus_key)

    # Build a DataFrame
    results = DataFrame(i=String[], j=String[], r1_factor=Float64[], x1_factor=Float64[])
    for (i, j) in sorted_keys
        if haskey(d_data, (i, j)) && haskey(d_new, (i, j))
            r_pu, x_pu = d_data[(i, j)]
            r_ohm, x_ohm = d_new[(i, j)]
            if abs(r_ohm) > 1e-12 && abs(x_ohm) > 1e-12
                push!(results, (i, j, r_pu / r_ohm, x_pu / x_ohm))
            end
        end
    end

    # Display the DataFrame
    show(results, allrows=true, allcols=true)
    return results
end

comparison_table = main()

vscodedisplay(comparison_table)