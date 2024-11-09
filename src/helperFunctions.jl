# helperFunctions.jl
module helperFunctions

export generateLoadShape, generateBinaryLoadShape, myprintln

function myprintln(verbose::Bool, msg::String)
    if verbose
        println(msg)
    end
end

function generateBinaryLoadShape(T::Int; filenameLoadShape::String="LoadShapePSubsCostDefault.dss",
    hi::Union{Float64,Nothing}=nothing,
    lo::Union{Float64,Nothing}=nothing,
    peakHoursFraction::Float64=0.3,
    verbose::Bool=false)
    # Get the working directory of this script
    wd = @__DIR__

    # Construct the full file path using the filename
    filepath = joinpath(wd, "..", "rawData", filenameLoadShape)

    # Print the file path to confirm it's the correct file
    myprintln(verbose, "Attempting to read LoadShape from: $filepath")

    # Initialize a vector for prices
    prices = Float64[]

    # Read the LoadShapePSubsCostDefault.dss file
    myprintln(verbose, "Reading file content:")
    open(filepath, "r") do file
        for line in eachline(file)
            line = strip(line)
            # Print each line read from the file
            myprintln(verbose, "Read line: $line")
            # Skip empty lines and comments
            if isempty(line) || startswith(line, "!")
                myprintln(verbose, "Skipping comment or empty line")
                continue
            end
            # Check if the line contains the 'mult=' string
            if occursin("mult=", line)
                # Extract the values inside the brackets
                price_str = replace(line, r".*mult=\[(.*)\]" => s"\1")
                # Parse the prices from the string
                prices = [parse(Float64, p) for p in split(price_str)]
                myprintln(verbose, "Parsed prices: $prices")
            else
                myprintln(verbose, "Skipping line, does not contain mult=")
            end
        end
    end

    # Check if prices is empty
    if isempty(prices)
        error("No data found in $filenameLoadShape. The prices array is empty.")
    end

    # Compute hi and lo based on the user input or default to the max/min of lambdaVals
    hi = isnothing(hi) ? maximum(prices) : hi
    lo = isnothing(lo) ? minimum(prices) : lo

    myprintln(verbose, "Computed hi: $hi, lo: $lo")

    # Subsampling or supersampling the 24-step price array to get lambdaVals0
    if T < 24
        # Subsample the middle T points
        start_idx = Int(floor((24 - T) / 2)) + 1
        lambdaVals0 = prices[start_idx:(start_idx+T-1)]
        myprintln(verbose, "Subsampling: lambdaVals0 = $lambdaVals0")
    elseif T > 24
        # Supersample using interpolation to create T values from the 24-step array
        lambdaVals0 = Float64[]
        myprintln(verbose, "Supersampling to $T values:")
        for i in 1:T
            scaled_idx = (i - 1) * (24 - 1) / (T - 1) + 1
            lower_idx = Int(floor(scaled_idx))
            upper_idx = min(lower_idx + 1, 24)
            frac = scaled_idx - lower_idx
            interpolated_value = (1 - frac) * prices[lower_idx] + frac * prices[upper_idx]
            push!(lambdaVals0, interpolated_value)
        end
        myprintln(verbose, "Supersampled: lambdaVals0 = $lambdaVals0")
    else
        # If T == 24, use the original prices
        lambdaVals0 = prices
        myprintln(verbose, "T == 24, using original prices: lambdaVals0 = $lambdaVals0")
    end

    # Now that we have lambdaVals0, compute the number of peak hours
    num_peak_hours = max(1, floor(Int, peakHoursFraction * T))

    myprintln(verbose, "Number of peak hours: $num_peak_hours")

    # Initialize the cost array to the low value (lo)
    costArray = fill(lo, T)
    myprintln(verbose, "Initialized costArray with low values: $costArray")

    # Sort lambdaVals0 in descending order and get the sorted indices
    sortedIndices = sortperm(lambdaVals0, rev=true)
    myprintln(verbose, "Sorted indices for peak hours: $sortedIndices")

    # Assign the high value (hi) to the peak hours
    costArray[sortedIndices[1:num_peak_hours]] .= hi
    myprintln(verbose, "Updated costArray with peak values: $costArray")

    # Create the output dictionary
    costData = Dict(
        :LoadShapeCost => costArray,
        :peakCost => hi,
        :offPeakCost => lo,
        :peakHoursFraction => peakHoursFraction
    )

    return costData
end

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
    LoadShape = Float64[]
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

    LoadShapeLoad = [Float64(x) for x in LoadShape] # a patch for T>24 where the resultant returned vector is Vector{Any} instead of Vector{Float64}
    return LoadShapeLoad
end

end # module
