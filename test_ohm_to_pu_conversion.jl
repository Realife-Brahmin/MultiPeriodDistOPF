#!/usr/bin/env julia

"""
Reads `linedata.txt` and `linenew.dss` to compare R/X values in ohms vs. per-unit,
then outputs a DataFrame with columns:
1) i
2) j
3) r_pu
4) x_pu
5) r_ohm
6) x_ohm
7) r1_factor = r_pu / r_ohm
8) x1_factor = x_pu / x_ohm

Now includes lines from `linenew.dss` even if they're commented out with "//" 
(but still contain "New Line." somewhere).
Also adds a "Transformer" column with Yes/No based on whether the line has '!' in linenew.dss.
"""

using Printf
using DataFrames

# Helper function to ensure (bus1, bus2) are in ascending numerical order
function sort_bus_pair(bus1::AbstractString, bus2::AbstractString)
    b1 = parse(Int, String(bus1))
    b2 = parse(Int, String(bus2))
    return b1 < b2 ? (bus1, bus2) : (bus2, bus1)
end

# Classify each line based on r1_factor (ratio) and x_ohm
function get_remark(r_ratio::Float64, x_ohm::Float64)::String
    if r_ratio <= 1 && x_ohm >= 0.005
        return "12.66kV to 12.66kV"
    elseif r_ratio >= 1 && x_ohm >= 0.005
        return "12.66kV to 0.24kV"
    elseif r_ratio > 1 && x_ohm <= 0.005
        return "0.24kV to 0.24kV"
    else
        @error "Unknown classification for r_ratio=$r_ratio and x_ohm=$x_ohm"
        return "Unknown"
    end
end

function parse_linedata(linedata_path::String)
    # Returns a dict of ((bus1, bus2) => (r_pu, x_pu))
    d = Dict{Tuple{String,String}, Tuple{Float64,Float64}}()
    for line in eachline(linedata_path)
        l = strip(line)
        # Skip empty lines or any line starting with "//"
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
    # Matches lines containing (bus1=..., bus2=..., r1=..., x1=...)
    # even if there's "//" in front or not, as long as "New Line." is present.
    # Also checks for '!' in the original line for the Transformer column.
    pattern = r"bus1\s*=\s*(\S+).*?bus2\s*=\s*(\S+).*?r1\s*=\s*([0-9.]+).*?x1\s*=\s*([0-9.]+)"
    d = Dict{Tuple{String,String}, Tuple{Float64,Float64,String}}()
    for raw_line in eachline(linenew_path)
        l = strip(raw_line)
        if occursin(r"(?i)new\s+line\.", l)
            m = match(pattern, l)
            if m !== nothing
                i_bus = m.captures[1]
                j_bus = m.captures[2]
                r_ohm = parse(Float64, m.captures[3])
                x_ohm = parse(Float64, m.captures[4])
                bus_pair = sort_bus_pair(i_bus, j_bus)
                transformer_flag = occursin('!', raw_line) ? "Yes" : "No"
                d[bus_pair] = (r_ohm, x_ohm, transformer_flag)
            end
        end
    end
    return d
end

function main()
    linedata_file = "./rawData/powerflowcomparision/Node730DG10/BFM/linedata.txt"
    linenew_file = "./rawData/powerflowcomparision/Node730DG10/Opendss/linenew.dss"

    d_data = parse_linedata(linedata_file)
    d_new  = parse_linenew(linenew_file)

    all_keys = union(keys(d_data), keys(d_new))

    # Sort by numerical bus IDs
    function bus_key(x)
        (parse(Int, x[1]), parse(Int, x[2]))
    end
    sorted_keys = sort(collect(all_keys), by=bus_key)

    # Build a DataFrame with columns in the specified order:
    # i, j, x_ohm, x_pu, x_factor, Transformer, r_ohm, r_pu, r_factor
    results = DataFrame(
        i=String[],
        j=String[],
        x_ohm=Float64[],
        x_pu=Float64[],
        x_factor=Float64[],
        Transformer=String[],
        r_ohm=Float64[],
        r_pu=Float64[],
        r_factor=Float64[]
    )

    for (i, j) in sorted_keys
        if haskey(d_data, (i, j)) && haskey(d_new, (i, j))
            (r_pu, x_pu) = d_data[(i, j)]
            (r_ohm, x_ohm, transformer_flag) = d_new[(i, j)]
            if abs(r_ohm) > 1e-12 && abs(x_ohm) > 1e-12
                x_factor = x_pu / x_ohm
                r_factor = r_pu / r_ohm
                push!(
                    results,
                    (
                        i,
                        j,
                        x_ohm,
                        x_pu,
                        x_factor,
                        transformer_flag,
                        r_ohm,
                        r_pu,
                        r_factor
                    )
                )
            end
        end
    end

    # Display the DataFrame
    show(results, allrows=true, allcols=true)
    return results
end

comparison_table = main()

vscodedisplay(comparison_table)