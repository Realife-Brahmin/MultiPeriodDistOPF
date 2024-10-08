# parseSystemSimulationData.jl

module parseSystemSimulationData

export parse_system_simulation_data

function parse_system_simulation_data(systemName::String)
    # Initialize parameters with default values
    substationBus = 1       # Default substation bus number
    V_Subs = 1.0            # Default per-unit voltage at substation
    V_base = 1.0            # Default base voltage in kV (line-to-ground)
    Δt = 1.0                # Default time step in hours

    wd = @__DIR__
    # Construct the file path using wd
    filename = joinpath(wd, "..", "rawData", systemName, "SysSim.dss")

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Remove comments and strip whitespace
            line = split(line, "!")[1]
            line = strip(line)

            # Skip empty lines
            if isempty(line)
                continue
            end

            # Parse the 'Edit "Vsource.source"' line
            if startswith(line, "Edit \"Vsource.source\"")
                # Extract parameters from the line
                tokens = split(line)
                for token in tokens[2:end]
                    if occursin("=", token)
                        key, value = split(token, "=")
                        key = strip(key)
                        value = strip(value)
                        if key == "bus1"
                            substationBus = parse(Int, value)
                        elseif key == "pu"
                            V_Subs = parse(Float64, value)
                        elseif key == "basekv"
                            V_base = parse(Float64, value)
                        end
                    end
                end
                # Parse the 'Set stepsize' line
            elseif startswith(line, "Set stepsize")
                # Extract the stepsize value
                pattern = r"Set stepsize\s*=\s*(.*)"
                m = match(pattern, line)
                if m !== nothing
                    stepsize_str = strip(m.captures[1])
                    # Handle units (assumed to be hours 'h')
                    if occursin("h", stepsize_str)
                        Δt = parse(Float64, replace(stepsize_str, "h" => ""))
                    elseif occursin("m", stepsize_str)
                        Δt = parse(Float64, replace(stepsize_str, "m" => "")) / 60.0
                    elseif occursin("s", stepsize_str)
                        Δt = parse(Float64, replace(stepsize_str, "s" => "")) / 3600.0
                    else
                        # Default to hours if no unit is specified
                        Δt = parse(Float64, stepsize_str)
                    end
                end
            end
        end
    end

    sysSimData = Dict(
        :substationBus => substationBus,
        :V_Subs => V_Subs,
        :V_base => V_base,
        :delta_t => Δt
    )

    # Return the extracted parameters as a dictionary
    return sysSimData

end

end # module
