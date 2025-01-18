module parseSystemSimulationData

export 
    parse_system_simulation_data, 
    post_process_data

using Parameters: @unpack, @pack!

#region parse_system_simulation_data
"""
    parse_system_simulation_data(systemName::String, T::Int; numAreas=1, alpha=1e-3, gamma=1e-3, kVA_B=1000, objfun0="genCostMin", objfun2="scd", temporal_decmp=false, PSubsMax_kW=Inf, inputForecastDescription="nonspecific", solver="Ipopt", tSOC_hard=true)

Parse system simulation data from the SysSim.dss file for a given system.

This function reads and parses the SysSim.dss file for the specified system, extracting relevant simulation parameters and initializing various settings. 
It handles the extraction of parameters such as substation bus, voltage, base voltage, and time step size.

# Arguments
- `systemName::String`: The name of the system for which to parse data.
- `T::Int`: The number of time steps for the simulation.
- `numAreas::Int`: The number of areas in the system (default: 1).
- `alpha::Float64`: The alpha parameter for the objective function (default: 1e-3).
- `gamma::Float64`: The gamma parameter for the objective function (default: 1e-3).
- `kVA_B::Float64`: The base power in kVA for per-unit calculations (default: 1000).
- `objfun0::String`: The primary objective function type (default: "genCostMin").
- `objfun2::String`: The secondary objective function type (default: "scd").
- `temporal_decmp::Bool`: A flag to enable temporal decomposition (default: false).
- `PSubsMax_kW::Float64`: The maximum substation power in kW (default: Inf).
- `inputForecastDescription::String`: A description of the input forecast (default: "nonspecific").
- `solver::String`: The solver to use for optimization (default: "Ipopt").
- `tSOC_hard::Bool`: A flag to enable hard terminal SOC constraints (default: true).

# Returns
- `sysSimData::Dict`: A dictionary containing the parsed system simulation data.

# Steps
1. **Parameter Initialization**: Initializes parameters with default values.
2. **File Path Construction**: Constructs the file path for SysSim.dss using the system name.
3. **File Reading**: Reads the SysSim.dss file for the specified system.
4. **Data Extraction**: Extracts simulation parameters such as substation bus, voltage, base voltage, and time step size.
5. **Data Initialization**: Initializes various parameters and dictionaries to store the extracted data.
6. **Return Data**: Returns a dictionary containing the parsed system simulation data.
"""
function parse_system_simulation_data(systemName::String, T::Int;
    numAreas=1,
    kVA_B = 1000,
    objfun0 = "subsPowerCostMin",
    objfun2 = "scd",
    temporal_decmp = false,
    PSubsMax_kW = Inf,
    inputForecastDescription = "nonspecific",
    solver = "Ipopt",
    tSOC_hard = false,
    relax_terminal_soc_constraint = false
    )

    # Initialize parameters with default values
    substationBus = 1       # Default substation bus number
    V_Subs = 1.0            # Default per-unit voltage at substation
    kV_B = 1.0            # Default base voltage in kV (line-to-ground)
    Δt = min(1.0, 24.0/T) # Currently has no effect as what's read from SysSim.dss takes precedence (which is always 1h as of now)

    wd = @__DIR__
    # Construct the file path using wd
    filename = joinpath(wd, "..", "..",  "rawData", systemName, "SysSim.dss")

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
                            kV_B = parse(Float64, value)
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

    MVA_B = kVA_B/1000
    Z_B = (kV_B)^2/MVA_B
    
    Tset = collect(1:T)
    Tset = sort(collect(Tset))

    if objfun0 == "subsPowerCostMin"
        objfunString = "Cost of Substation Power"
        objfunSense = "Min"
        objfunPrefix = "subsPowerCost_min"
        objfunUnit = "\$"
    elseif objfun0 == "lineLossMin"
        objfunString = "Line Losses"
        objfunSense = "Min"
        objfunPrefix = "lineLoss_min"
        objfunUnit = "kW"
    elseif objfun0 == "subsPowerMin"
        objfunString = "Substation Power"
        objfunSense = "Min"
        objfunPrefix = "subsPower_min"
        objfunUnit = "kW"
    else
        objfunString = "unknown objective";
        objfunSense = "Min"
        objfunPrefix = "unknown_obj" 
    end

    if objfun2 == "scd"
        objfunAppendix = "with_scd"
    else
        objfunAppendix = ""
    end

    objfunConciseDescription = objfunPrefix*"_"*objfunAppendix 

    if temporal_decmp == true
        temporalDecmpString = "Temporally Decomposed via DDP"
        temporalDecmpAppendix = "tmprl_dcmpsd"
    elseif temporal_decmp == false
        temporalDecmpString = "Temporally Brute-forced"
        temporalDecmpAppendix = "tmprl_bruteforced"
    else
        @error "floc"
    end

    if numAreas > 1
        spatialDecString = "Spatially Decomposed into $(numAreas) areas"
        spatialDecAppendix = "spat_dcmpsd_$(numAreas)_areas"
    elseif numAreas == 1
        spatialDecString = "Spatially Centralized"
        spatialDecAppendix = "spat_centr_system"
    else
        @error "floc"
    end
    machine_ID = gethostname()

    simNatureString = temporalDecmpString * ", " * spatialDecString
    simNatureAppendix = temporalDecmpAppendix * "_" * spatialDecAppendix

    macroItrsCompleted = 0
    solution_time = -1 

    sysSimData = Dict(
        :inputForecastDescription => inputForecastDescription,
        :machine_ID => machine_ID,
        :macroItrsCompleted => macroItrsCompleted,
        :systemName => systemName, # user input
        :numAreas => numAreas, # user input
        :solution_time => solution_time,
        :substationBus => substationBus,
        :V_Subs => V_Subs,
        :kV_B => kV_B,
        :kVA_B => kVA_B,
        :Z_B => Z_B,
        :delta_t => Δt,
        :objfun0 => objfun0, # user input
        :objfun2 => objfun2, # user input
        :objfunString => objfunString,
        :objfunSense => objfunSense,
        :objfunPrefix => objfunPrefix,
        :objfunAppendix => objfunAppendix,
        :objfunConciseDescription => objfunConciseDescription,
        :objfunUnit => objfunUnit,
        :PSubsMax_kW => PSubsMax_kW,
        :relax_terminal_soc_constraint => relax_terminal_soc_constraint,
        :simNatureAppendix => simNatureAppendix,
        :simNatureString => simNatureString,
        :solver => solver,
        :spatialDecAppendix => spatialDecAppendix,
        :spatialDecString => spatialDecString,
        :tSOC_hard => tSOC_hard,
        :T => T, # user input
        :temporal_decmp => temporal_decmp,
        :temporalDecmpString => temporalDecmpString,
        :temporalDecmpAppendix => temporalDecmpAppendix,
        :Tset => Tset
    )

    # Return the extracted parameters as a dictionary
    return sysSimData

end
#endregion

#region post_process_data
"""
    post_process_data(data)

Post-process the parsed data to add additional information.

This function performs post-processing on the parsed data to add additional information such as DER and battery percentages.

# Arguments
- `data::Dict`: A dictionary containing the parsed data.

# Returns
- `data::Dict`: The updated dictionary with additional post-processed information.
"""
function post_process_data(data)
    
    @unpack DER_percent, Batt_percent = data;
    gedString = "$(DER_percent)% PVs and $(Batt_percent)% Batteries"
    gedAppendix = "pv_$(DER_percent)_batt_$(Batt_percent)"
    @pack! data = gedAppendix, gedString

    return data
end
#endregion

end # module parseSystemSimulationData
