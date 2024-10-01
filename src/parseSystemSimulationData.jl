# parseSystemSimulationData.jl

module parseSystemSimulationData

export parse_system_simulation_data

function parse_system_simulation_data(filename::String)
    # Initialize variables
    T = 24  # Default to 24 time periods if not specified
    C = Dict{Int, Float64}()  # Cost data
    η_C = 0.9  # Charging efficiency
    η_D = 0.9  # Discharging efficiency
    V_base = 1.0  # Base voltage (per unit)
    V_min = 0.95
    V_max = 1.05

    # Open and read the file
    open(filename, "r") do file
        for line in eachline(file)
            # Skip comments and empty lines
            if isempty(strip(line)) || startswith(strip(line), "!") || startswith(strip(line), "//")
                continue
            end
            # Parse system simulation data
            # Example lines:
            # "Set T=24"
            # "Set Cost=50"
            # "Set η_C=0.9"
            # "Set η_D=0.9"
            # "Set V_base=1.0"
            tokens = split(strip(line))
            if startswith(tokens[1], "Set")
                key_value = split(tokens[2], "=")
                key = key_value[1]
                value = key_value[2]
                if key == "T"
                    T = parse(Int, value)
                elseif key == "Cost"
                    for t in 1:T
                        C[t] = parse(Float64, value)
                    end
                elseif key == "η_C"
                    η_C = parse(Float64, value)
                elseif key == "η_D"
                    η_D = parse(Float64, value)
                elseif key == "V_base"
                    V_base = parse(Float64, value)
                elseif key == "V_min"
                    V_min = parse(Float64, value)
                elseif key == "V_max"
                    V_max = parse(Float64, value)
                end
            end
        end
    end

    V_Subs = 1.03
    Tset = 1:T
    v_min = (V_min * V_base)^2
    v_max = (V_max * V_base)^2
    return T, Tset, C, η_C, η_D, V_base, V_Subs, v_min, v_max
end

end # module
