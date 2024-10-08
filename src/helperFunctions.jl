# helperFunctions.jl
module helperFunctions

export generateLoadShape, myprintln

function generateLoadShape(T::Int; filenameLoadShape=nothing)
    wd = @__DIR__

    defaultLoadShape = []
    if filenameLoadShape === nothing
        # Default to LoadShapePVDefault.dss if no filenameLoadShape is provided
        filenameLoadShape = "LoadShapePVDefault.dss"
    end

    # Construct the full file path using the provided filenameLoadShape
    loadshape_filepath = joinpath(wd, "..", "rawData", filenameLoadShape)

    # Read the LoadShape file
    open(loadshape_filepath, "r") do file
        for line in eachline(file)
            line = split(line, "!")[1]  # Remove comments
            line = strip(line)
            if isempty(line)
                continue
            end
            # Parse the load shape values
            if occursin("mult=[", line)
                shape_str = strip(split(line, "mult=")[2], ['[', ']'])
                defaultLoadShape = [parse(Float64, val) for val in split(shape_str)]
            end
        end
    end

    # Generate load shape values based on T
    LoadShape = []
    if T < 24
        # Take the middle T values
        start_idx = Int(floor((24 - T) / 2)) + 1
        LoadShape = defaultLoadShape[start_idx:(start_idx+T-1)]
    elseif T > 24
        # Supersample the values
        for i in 1:T
            scaled_idx = (i - 1) * (24 - 1) / (T - 1) + 1
            lower_idx = Int(floor(scaled_idx))
            upper_idx = min(lower_idx + 1, 24)
            frac = scaled_idx - lower_idx
            interpolated_value = (1 - frac) * defaultLoadShape[lower_idx] + frac * defaultLoadShape[upper_idx]
            push!(LoadShape, interpolated_value)
        end
    else
        # T == 24, use the default load shape
        LoadShape = defaultLoadShape
    end

    return LoadShape
end


function myprintln(msg::String, verbose::Bool=true)
    if verbose
        println(msg)
    end
end

end
