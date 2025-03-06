module openDSSValidator

export 
    compute_highest_allTime_voltage_discrepancy,
    export_validation_decision_variables, 
    get_battery_powers_opendss_powerflow_for_timestep_t,
    get_load_powers_opendss_powerflow_for_timestep_t,
    get_pv_powers_opendss_powerflow_for_timestep_t,
    get_source_bus, 
    get_substation_lines, 
    get_terminal_soc_values_opendss_powerflow,
    get_voltages_opendss_powerflow_for_timestep_t,
    set_battery_controls_opendss_powerflow_for_timestep_t,
    set_custom_load_shape!, 
    set_pv_controls_opendss_powerflow_for_timestep_t,
    validate_opf_against_opendss

using CSV
using DataFrames
using JuMP: value
using OpenDSSDirect
using Parameters: @unpack, @pack!

include("helperFunctions.jl")
import .helperFunctions as HF

#region compute_highest_allTime_voltage_discrepancy
"""
    compute_highest_allTime_voltage_discrepancy(modelDict, valdVals)

Compute the highest voltage discrepancy over all time steps.

This function calculates the maximum voltage discrepancy between the optimization results in `modelDict` and the powerflow results in `valdVals` over all time steps.
"""
function compute_highest_allTime_voltage_discrepancy(modelDict, valdVals)
    @unpack modelVals, data = modelDict
    # Initialize the maximum voltage discrepancy
    disc_voltage_all_time_pu = 0.0

    # Retrieve the model voltage variable `v`
    v = modelVals[:v]
    @unpack Tset = data;
    # Iterate over each timestep
    for t in Tset
        # Retrieve the dictionary of bus voltages for the current timestep
        vald_voltages_dict = valdVals[:vald_voltages_vs_t_1toT_pu][t]

        # Iterate over each bus in the dictionary
        for (bus_index, vald_voltage) in vald_voltages_dict
            # Compute the model voltage for the bus at the current timestep
            model_voltage = sqrt(v[bus_index, t])

            # Calculate the discrepancy
            discrepancy = abs(vald_voltage - model_voltage)

            # Update the maximum discrepancy if needed
            if discrepancy > disc_voltage_all_time_pu
                disc_voltage_all_time_pu = discrepancy
            end
        end
    end

    return disc_voltage_all_time_pu
end
#endregion

#region get_battery_powers_opendss_powerflow_for_timestep_t
"""
    get_battery_powers_opendss_powerflow_for_timestep_t(; verbose::Bool=false)

Retrieve battery real and reactive powers for a given timestep from OpenDSS powerflow simulation.

This function calculates the total real and reactive power for all battery systems at a given timestep from the OpenDSS powerflow simulation.
"""
function get_battery_powers_opendss_powerflow_for_timestep_t(;
    verbose::Bool=false)
    # Initialize total battery power values
    vald_battery_real_power_t_kW = 0.0
    vald_battery_reactive_power_t_kVAr = 0.0
    vald_battery_real_power_transaction_magnitude_t_kW = 0.0
    vald_battery_reactive_power_transaction_magnitude_t_kVAr = 0.0

    # Retrieve all battery system names
    battery_names = OpenDSSDirect.Storages.AllNames()

    # Iterate through each battery system to calculate total real and reactive power
    for battery_name in battery_names
        OpenDSSDirect.Circuit.SetActiveElement("Storage.$battery_name")
        battery_powers = OpenDSSDirect.CktElement.Powers()

        # Retrieve and sum real and reactive power components
        real_power = -real(battery_powers[1])  # Negate for power injection convention
        reactive_power = -imag(battery_powers[1])

        # Accumulate total battery power for this timestep
        vald_battery_real_power_t_kW += real_power
        vald_battery_reactive_power_t_kVAr += reactive_power

        # Track transaction magnitude for each component
        vald_battery_real_power_transaction_magnitude_t_kW += abs(real_power)
        vald_battery_reactive_power_transaction_magnitude_t_kVAr += abs(reactive_power)

        HF.myprintln(verbose, "Battery $battery_name | Real Power: $real_power kW, Reactive Power: $reactive_power kVAr")
    end

    # Store results in a dictionary and return
    batteryPowersDict_t = Dict(
        :vald_battery_real_power_t_kW => vald_battery_real_power_t_kW,
        :vald_battery_reactive_power_t_kVAr => vald_battery_reactive_power_t_kVAr,
        :vald_battery_real_power_transaction_magnitude_t_kW => vald_battery_real_power_transaction_magnitude_t_kW,
        :vald_battery_reactive_power_transaction_magnitude_t_kVAr => vald_battery_reactive_power_transaction_magnitude_t_kVAr
    )

    return batteryPowersDict_t
end
#endregion

#region get_load_powers_opendss_powerflow_for_timestep_t
"""
    get_load_powers_opendss_powerflow_for_timestep_t()

Retrieve load real and reactive powers for a given timestep from OpenDSS powerflow simulation.

This function calculates the total real and reactive power for all loads at a given timestep from the OpenDSS powerflow simulation.
"""
function get_load_powers_opendss_powerflow_for_timestep_t(;
    t=nothing,
    LoadShapeLoad=nothing,
    verbose::Bool=false)
    # Initialize total load power values
    total_load_t_kW = 0.0
    total_load_t_kVAr = 0.0

    if isnothing(t)
        # If no specific timestep is provided, use the current time step
        t = OpenDSSDirect.Solution.Step()
        myprintln(verbose, "t not provided, using current time step from OpenDSS: $t")
    end

    if isnothing(LoadShapeLoad)
        # If no load shape is provided, use a default value
        myprintln(verbose, "LoadShapeLoad not provided, so cannot check for correctness of load dispatch")
    end
    # Retrieve all load names
    load_names = OpenDSSDirect.Loads.AllNames()
    OpenDSSDirect.Loads.First()
    # Iterate through each load to calculate total real and reactive power
    for (count, load_name) in enumerate(load_names)
        OpenDSSDirect.Circuit.SetActiveElement("Load.$load_name")
        load_powers = OpenDSSDirect.CktElement.Powers()
        total_load_t_kW += real(load_powers[1])
        total_load_t_kVAr += imag(load_powers[1])

        # Print details only for the first 5 loads
        if verbose && count <= 50 && !isnothing(LoadShapeLoad)
            # Retrieve the rated kW from the load
            rated_kW = OpenDSSDirect.Loads.kW()
            # Retrieve the load shape name from OpenDSS (if any)
            shape_name = OpenDSSDirect.Loads.Daily()
            OpenDSSDirect.Loads.Next()
            # HF.myprintln(verbose,
            #     "Load $load_name => Intended  Power: $(rated_kW*LoadShapeLoad[t]) kW,
            #     Actual Power: $(real(load_powers[1])) kW"
            # )
        end
    end

    # Store results in a dictionary and return
    loadPowersDict_t = Dict(
        :total_load_t_kW => total_load_t_kW,
        :total_load_t_kVAr => total_load_t_kVAr
    )

    return loadPowersDict_t
end
#endregion

#region get_pv_powers_opendss_powerflow_for_timestep_t
"""
    get_pv_powers_opendss_powerflow_for_timestep_t()

Retrieve PV real and reactive powers for a given timestep from OpenDSS powerflow simulation.

This function calculates the total real and reactive power for all PV systems at a given timestep from the OpenDSS powerflow simulation.
"""
function get_pv_powers_opendss_powerflow_for_timestep_t()
    # Initialize total PV power values
    total_pv_t_kW = 0.0
    total_pv_t_kVAr = 0.0

    # Retrieve all PV system names
    pv_names = OpenDSSDirect.PVsystems.AllNames()

    # Iterate through each PV system to calculate total real and reactive power
    for pv_name in pv_names
        OpenDSSDirect.Circuit.SetActiveElement("PVSystem.$pv_name")
        pv_powers = OpenDSSDirect.CktElement.Powers()
        total_pv_t_kW -= real(pv_powers[1])  # Negate since power is injected
        total_pv_t_kVAr -= imag(pv_powers[1])
    end

    # Store results in a dictionary and return
    pvPowersDict_t = Dict(
        :total_pv_t_kW => total_pv_t_kW,
        :total_pv_t_kVAr => total_pv_t_kVAr
    )

    return pvPowersDict_t
end
#endregion

#region get_source_bus
"""
    get_source_bus()

Retrieve the source bus name from OpenDSS.

This function retrieves the name of the source bus from the OpenDSS simulation.
"""
function get_source_bus()
    vsource_element = OpenDSSDirect.Vsources.First()
    vsource_name = OpenDSSDirect.Vsources.Name()
    OpenDSSDirect.Circuit.SetActiveElement("Vsource.$vsource_name")
    source_bus = OpenDSSDirect.CktElement.BusNames()[1]
    return source_bus
end
#endregion

#region get_substation_lines
"""
    get_substation_lines(substation_bus::String)

Retrieve the lines connected to the substation bus from OpenDSS.

This function retrieves the names of all lines connected to the specified substation bus from the OpenDSS simulation.
"""
function get_substation_lines(substation_bus::String)
    substation_lines = []

    # Iterate over all lines in the circuit
    line_id = Lines.First()
    while line_id > 0
        line_name = Lines.Name()
        OpenDSSDirect.Circuit.SetActiveElement("Line.$line_name")
        bus1 = OpenDSSDirect.CktElement.BusNames()[1]
        bus2 = OpenDSSDirect.CktElement.BusNames()[2]

        if bus1 == substation_bus || bus2 == substation_bus
            push!(substation_lines, line_name)
        end

        line_id = Lines.Next()
    end

    return substation_lines
end
#endregion

#region get_substation_powers_opendss_powerflow_for_timestep_t
"""
    get_substation_powers_opendss_powerflow_for_timestep_t(data; useVSourcePower::Bool=true)

Retrieve substation real and reactive powers for a given timestep from OpenDSS powerflow simulation.

This function calculates the total real and reactive power at the substation for a given timestep, either using VSource power or individual line powers.
"""
function get_substation_powers_opendss_powerflow_for_timestep_t(data; useVSourcePower::Bool=true)
    # Unpack the substation bus from data
    @unpack substationBus = data

    # Initialize variables for substation power totals
    P_substation_total_t_kW = 0.0
    Q_substation_total_t_kVAr = 0.0

    if useVSourcePower
        # Directly use VSource power if the kwarg is true
        OpenDSSDirect.Circuit.SetActiveElement("Vsource.source")
        vsource_powers = -OpenDSSDirect.CktElement.Powers()
        P_substation_total_t_kW = real(vsource_powers[1])
        Q_substation_total_t_kVAr = imag(vsource_powers[1])
    else
        # Use individual line powers if the kwarg is false
        substation_lines = get_substation_lines(substationBus)

        for line in substation_lines
            OpenDSSDirect.Circuit.SetActiveElement("Line.$line")
            line_powers = OpenDSSDirect.CktElement.Powers()
            P_line = sum(real(line_powers[1]))
            Q_line = sum(imag(line_powers[1]))

            P_substation_total_t_kW += P_line
            Q_substation_total_t_kVAr += Q_line
        end
    end

    # Return the results as a dictionary
    substationPowersDict_t = Dict(
        :P_substation_total_t_kW => P_substation_total_t_kW,
        :Q_substation_total_t_kVAr => Q_substation_total_t_kVAr
    )

    return substationPowersDict_t
end
#endregion

#region get_terminal_soc_values_opendss_powerflow
"""
    get_terminal_soc_values_opendss_powerflow(data)

Retrieve terminal state of charge (SOC) values from OpenDSS powerflow simulation.

This function calculates the SOC values and SOC violations for all storage elements at the terminal timestep from the OpenDSS powerflow simulation.
"""
function get_terminal_soc_values_opendss_powerflow(data)
    # Unpack required fields from data
    @unpack Bref_percent, B_R = data

    # Initialize output dictionary
    terminalSOCDict = Dict(
        :vald_Bj_T_pu => Dict{Int, Float64}(),          # Maps storage_number to SOC in pu
        :soc_violation_j_kWh => Dict{Int, Float64}(),   # Maps storage_number to SOC violation in kWh
        :vald_terminal_soc_violation_kWh => 0.0         # Cumulative SOC violation in kWh
    )

    # Initialize storage_id to start iterating over Storages
    storage_id = OpenDSSDirect.Storages.First()

    # Iterate over all storage elements
    while storage_id > 0
        # Get storage name and determine storage number assuming 'batteryX' naming convention
        storage_name = OpenDSSDirect.Storages.Name()
        storage_number = parse(Int, split(storage_name, "battery")[2])

        # Retrieve the SOC in per-unit for this battery
        vald_Bj_T_pu = OpenDSSDirect.Storages.puSOC()
        terminalSOCDict[:vald_Bj_T_pu][storage_number] = vald_Bj_T_pu  # Store SOC in dict with storage_number as key

        # Calculate SOC violation in kWh
        soc_violation = abs(vald_Bj_T_pu - Bref_percent[storage_number]) * B_R[storage_number]
        terminalSOCDict[:soc_violation_j_kWh][storage_number] = soc_violation  # Store SOC violation in dict

        # Accumulate total SOC violation in kWh
        terminalSOCDict[:vald_terminal_soc_violation_kWh] += soc_violation

        # Move to the next storage element
        storage_id = OpenDSSDirect.Storages.Next()
    end

    return terminalSOCDict
end
#endregion

#region get_voltages_opendss_powerflow_for_timestep_t
"""
    get_voltages_opendss_powerflow_for_timestep_t()

Retrieve bus voltage magnitudes for a given timestep from OpenDSS powerflow simulation.

This function retrieves the voltage magnitudes for all buses at a given timestep from the OpenDSS powerflow simulation.
"""
function get_voltages_opendss_powerflow_for_timestep_t()
    # Initialize a dictionary to store voltages with integer bus numbers as keys
    vald_voltage_dict_t_pu = Dict{Int,Float64}()

    # Get the bus names and corresponding voltage magnitudes
    bus_names = OpenDSSDirect.Circuit.AllBusNames()
    bus_voltages = OpenDSSDirect.Circuit.AllBusMagPu()

    # Populate the dictionary with integer keys
    for (i, bus_name) in enumerate(bus_names)
        # Assuming bus names are integers in string format, like "1", "2", etc.
        bus_number = parse(Int, bus_name)
        vald_voltage_dict_t_pu[bus_number] = bus_voltages[i]
    end

    return vald_voltage_dict_t_pu
end
#endregion

#region set_custom_load_shape!
"""
    set_custom_load_shape!(LoadShapeArray::Vector{Float64}; verbose::Bool=false)

Set a custom load shape in OpenDSS.

This function defines a custom load shape in OpenDSS using the provided load shape array and applies it to all loads in the system.
"""
function set_custom_load_shape!(LoadShapeArray::Vector{Float64};
    verbose::Bool=false)
    # Define LoadShapeLoad in OpenDSS with the provided LoadShapeArray array
    loadshape_command = "New Loadshape.LoadShapeLoad npts = $(length(LoadShapeArray)) interval = 1 mult = [" *
                        join(LoadShapeArray, " ") * "]"
    OpenDSSDirect.Text.Command(loadshape_command)
    HF.myprintln(verbose, "Defined LoadShapeLoad with provided LoadShapeArray")

    # Apply this load shape to all loads in the system
    load_id = OpenDSSDirect.Loads.First()
    while load_id > 0
        OpenDSSDirect.Loads.Daily("LoadShapeLoad")
        load_id = OpenDSSDirect.Loads.Next()
    end
    HF.myprintln(verbose, "Applied LoadShapeLoad to all loads")
end
#endregion

#region set_battery_controls_opendss_powerflow_for_timestep_t
"""
    set_battery_controls_opendss_powerflow_for_timestep_t(modelDict, t; verbose=false)

Set battery control actions for a given timestep in OpenDSS powerflow simulation.

This function sets the charging, discharging, and reactive power levels for all batteries at a given timestep in the OpenDSS powerflow simulation.
"""
function set_battery_controls_opendss_powerflow_for_timestep_t(modelDict, t; verbose=false)
    @unpack modelVals, data = modelDict
    # Unpack necessary data
    P_c = modelVals[:P_c]
    P_d = modelVals[:P_d]
    q_B = modelVals[:q_B]
    @unpack kVA_B_dict = data

    # Set battery power levels
    storage_id = OpenDSSDirect.Storages.First()
    while storage_id > 0
        storage_name = OpenDSSDirect.Storages.Name()
        storage_number = parse(Int, split(storage_name, "battery")[2])

        # Calculate power levels based on optimization model variables
        charge_power_kW = P_c[storage_number, t] * kVA_B_dict[storage_number]
        discharge_power_kW = P_d[storage_number, t] * kVA_B_dict[storage_number]
        net_power_kW = discharge_power_kW - charge_power_kW
        reactive_power_kVAr = q_B[storage_number, t] * kVA_B_dict[storage_number]

        # Command to set battery power levels in OpenDSS
        command_str = "Edit Storage.Battery$(storage_number) kW=$(net_power_kW) kvar=$(reactive_power_kVAr)"
        OpenDSSDirect.Text.Command(command_str)

        # Optionally print command for verification
        verbose = true
        if verbose
            println("Time Step $t: Setting battery $storage_number with command: $command_str")
        end

        # Move to the next storage element
        storage_id = OpenDSSDirect.Storages.Next()
    end
end
#endregion

#region set_pv_controls_opendss_powerflow_for_timestep_t
"""
    set_pv_controls_opendss_powerflow_for_timestep_t(modelDict, t; verbose::Bool=false)

Set PV control actions for a given timestep in OpenDSS powerflow simulation.

This function sets the real and reactive power levels for all PV systems at a given timestep in the OpenDSS powerflow simulation.
"""
function set_pv_controls_opendss_powerflow_for_timestep_t(modelDict, t; verbose::Bool=false)
    @unpack modelVals, data = modelDict;
    # Unpack necessary data fields from `data`
    @unpack kVA_B_dict, p_D_pu = data
    q_D = modelVals[:q_D]  # Access q_D from the model

    # Set power levels for PV systems at each time step
    pv_id = OpenDSSDirect.PVsystems.First()
    while pv_id > 0
        pv_name = OpenDSSDirect.PVsystems.Name()
        pv_number = parse(Int, split(pv_name, "pv")[2])

        # Retrieve real and reactive power setpoints for this PV system and timestep
        p_D_t_kW = p_D_pu[(pv_number, t)] * kVA_B_dict[pv_number]
        q_D_t_kVAr = q_D[pv_number, t] * kVA_B_dict[pv_number]

        # Set real and reactive power for the PV system
        OpenDSSDirect.PVsystems.Pmpp(p_D_t_kW)
        OpenDSSDirect.PVsystems.kvar(q_D_t_kVAr)

        # verbose = true
        if verbose
            println("Setting PV for bus $(pv_number) at t = $(t): p_D_t_kW = $(p_D_t_kW), q_D_t_kVAr = $(q_D_t_kVAr)")
        end

        # Move to the next PV system
        pv_id = OpenDSSDirect.PVsystems.Next()
    end
end
#endregion

#region validate_opf_against_opendss
"""
    validate_opf_against_opendss(modelDict; verbose::Bool=false)

Validate the OPF results against OpenDSS powerflow simulation.

This function compares the optimization results stored in `modelDict` with the powerflow results from OpenDSS stored in `valdVals`.
It performs the following steps:
1. Initializes the OpenDSS simulation environment.
2. Sets custom load shapes and control actions for PV and battery systems.
3. Runs the powerflow simulation for each timestep.
4. Retrieves and stores various powerflow results, including:
    - Circuit losses
    - Substation real and reactive powers
    - Load real and reactive powers
    - PV real and reactive powers
    - Battery real and reactive powers
    - Voltage magnitudes for all buses
5. Aggregates generation data and calculates discrepancies between the optimization and powerflow results.
6. Computes terminal state of charge (SOC) values and SOC violations for storage elements.
7. Returns the updated `modelDict` with validation results.

The function ensures that the state variables from the powerflow simulation match those solved by the optimizer.
"""
function validate_opf_against_opendss(modelDict; verbose::Bool=false)
    @unpack modelVals, data = modelDict
    # Initialize valdVals dictionary for timestep-specific and cumulative results
    @unpack T, LoadShapePV, Dset, Bset = data
    LoadShapeLoad = data[:LoadShapeLoad]

    valdVals = Dict(
        # Timestep-specific fields
        :vald_PLoss_vs_t_1toT_kW => zeros(T),
        :vald_PSubs_vs_t_1toT_kW => zeros(T),
        :vald_QLoss_vs_t_1toT_kVAr => zeros(T),
        :vald_QSubs_vs_t_1toT_kVAr => zeros(T),
        :vald_load_real_power_vs_t_1toT_kW => zeros(T),
        :vald_load_reactive_power_vs_t_1toT_kVAr => zeros(T),
        :vald_pv_real_power_vs_t_1toT_kW => zeros(T),
        :vald_pv_reactive_power_vs_t_1toT_kVAr => zeros(T),
        :vald_battery_real_power_vs_t_1toT_kW => zeros(T),
        :vald_battery_reactive_power_vs_t_1toT_kVAr => zeros(T),
        :vald_static_cap_reactive_power_vs_t_1toT_kVAr => zeros(T),
        :vald_total_gen_reactive_power_vs_t_1toT_kVAr => zeros(T),
        :vald_total_gen_real_power_vs_t_1toT_kW => zeros(T),
        :vald_voltages_vs_t_1toT_pu => Vector{Dict{Int,Float64}}(undef, T),
        :vald_PSubsCost_vs_t_1toT_dollar => zeros(T),
        :vald_battery_real_power_transaction_magnitude_vs_t_1toT_kW => zeros(T),
        :vald_battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr => zeros(T),

        # Cumulative fields
        :vald_battery_reactive_power_allT_kVAr => 0.0,
        :vald_battery_reactive_power_transaction_magnitude_allT_kVAr => 0.0,
        :vald_battery_real_power_allT_kW => 0.0,
        :vald_battery_real_power_transaction_magnitude_allT_kW => 0.0,
        :vald_pv_reactive_power_allT_kVAr => 0.0,
        :vald_pv_real_power_allT_kW => 0.0,
        :vald_PLoss_allT_kW => 0.0,
        :vald_PSubs_allT_kW => 0.0,
        :vald_QLoss_allT_kVAr => 0.0,
        :vald_QSubs_allT_kVAr => 0.0,
        :vald_static_cap_reactive_power_allT_kVAr => 0.0,
        :vald_total_gen_real_power_allT_kW => 0.0,
        :vald_total_gen_reactive_power_allT_kVAr => 0.0,
        :vald_load_real_power_allT_kW => 0.0,
        :vald_load_reactive_power_allT_kVAr => 0.0,
        :vald_PSubsCost_allT_dollar => 0.0,
        :vald_solution_time => 0.0,
        :vald_static_cap_reactive_power_allT_kVAr => 0.0,
        :vald_substation_real_power_peak_allT_kW => 0.0,
        :vald_terminal_soc_violation_kWh => 0.0
    )

    useVSourcePower = true
    # Set paths for DSS files
    @unpack systemName, rawDataFolderPath = data;
    dss_dir = joinpath(rawDataFolderPath, systemName)
    # dss_dir = joinpath(@__DIR__, "rawData", system_name)
    dss_file = joinpath(dss_dir, "Master.dss")
    HF.myprintln(verbose, "Master.dss file path: $(dss_file)")

    # Initialize OpenDSS
    OpenDSSDirect.Text.Command("Clear")
    OpenDSSDirect.Text.Command("Redirect \"$dss_file\"")

    # Set custom load shape
    set_custom_load_shape!(LoadShapeLoad)

    @show LoadShapeLoad

    # Loop through each timestep to perform power flow and store valdVals
    for t in 1:T
        set_pv_controls_opendss_powerflow_for_timestep_t(modelDict, t)
        set_battery_controls_opendss_powerflow_for_timestep_t(modelDict, t)
        Solution.Solve()

        # Retrieve circuit losses
        totalLosses_t_kVA = Circuit.Losses() ./ 1000
        totalLosses_t_kW, totalLosses_t_kVAr = real(totalLosses_t_kVA), -imag(totalLosses_t_kVA)
        valdVals[:vald_PLoss_vs_t_1toT_kW][t] = totalLosses_t_kW
        valdVals[:vald_PLoss_allT_kW] += totalLosses_t_kW
        valdVals[:vald_QLoss_vs_t_1toT_kVAr][t] = totalLosses_t_kVAr
        valdVals[:vald_QLoss_allT_kVAr] += totalLosses_t_kVAr

        # Retrieve substation real and reactive powers post powerflow for this timestep
        substationPowersDict_t = get_substation_powers_opendss_powerflow_for_timestep_t(data, useVSourcePower=true)
        @unpack P_substation_total_t_kW, Q_substation_total_t_kVAr = substationPowersDict_t

        valdVals[:vald_PSubs_vs_t_1toT_kW][t] = P_substation_total_t_kW
        valdVals[:vald_PSubs_allT_kW] += P_substation_total_t_kW
        valdVals[:vald_QSubs_vs_t_1toT_kVAr][t] = Q_substation_total_t_kVAr
        valdVals[:vald_QSubs_allT_kVAr] += Q_substation_total_t_kVAr
        valdVals[:vald_substation_real_power_peak_allT_kW] = max(valdVals[:vald_substation_real_power_peak_allT_kW], P_substation_total_t_kW)

        # Cost calculation
        valdVals[:vald_PSubsCost_vs_t_1toT_dollar][t] = data[:LoadShapeCost][t] * P_substation_total_t_kW * data[:delta_t]
        valdVals[:vald_PSubsCost_allT_dollar] += valdVals[:vald_PSubsCost_vs_t_1toT_dollar][t]

        # Retrieve load real and reactive powers
        if t == 1
            verboseLoad = true
        else
            verboseLoad = false
        end
        loadPowersDict_t = get_load_powers_opendss_powerflow_for_timestep_t(verbose=verboseLoad, t=t, LoadShapeLoad=LoadShapeLoad)
        @unpack total_load_t_kW, total_load_t_kVAr = loadPowersDict_t

        valdVals[:vald_load_real_power_vs_t_1toT_kW][t] = total_load_t_kW
        valdVals[:vald_load_reactive_power_vs_t_1toT_kVAr][t] = total_load_t_kVAr
        valdVals[:vald_load_real_power_allT_kW] += total_load_t_kW
        valdVals[:vald_load_reactive_power_allT_kVAr] += total_load_t_kVAr

        # PV and battery power calculations
        pvPowersDict_t = get_pv_powers_opendss_powerflow_for_timestep_t()
        @unpack total_pv_t_kW, total_pv_t_kVAr = pvPowersDict_t

        valdVals[:vald_pv_real_power_vs_t_1toT_kW][t] = total_pv_t_kW
        valdVals[:vald_pv_real_power_allT_kW] += total_pv_t_kW
        valdVals[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] = total_pv_t_kVAr
        valdVals[:vald_pv_reactive_power_allT_kVAr] += total_pv_t_kVAr

        batteryPowersDict_t = get_battery_powers_opendss_powerflow_for_timestep_t(verbose=verbose)
        @unpack vald_battery_real_power_t_kW, vald_battery_reactive_power_t_kVAr, vald_battery_real_power_transaction_magnitude_t_kW, vald_battery_reactive_power_transaction_magnitude_t_kVAr = batteryPowersDict_t

        valdVals[:vald_battery_reactive_power_vs_t_1toT_kVAr][t] = vald_battery_reactive_power_t_kVAr
        valdVals[:vald_battery_reactive_power_allT_kVAr] += vald_battery_reactive_power_t_kVAr
        valdVals[:vald_battery_reactive_power_transaction_magnitude_t_kVAr] = vald_battery_reactive_power_transaction_magnitude_t_kVAr
        valdVals[:vald_battery_reactive_power_transaction_magnitude_allT_kVAr] += vald_battery_reactive_power_transaction_magnitude_t_kVAr
        valdVals[:vald_battery_real_power_vs_t_1toT_kW][t] = vald_battery_real_power_t_kW
        valdVals[:vald_battery_real_power_allT_kW] += vald_battery_real_power_t_kW
        valdVals[:vald_battery_real_power_transaction_magnitude_t_kW] = vald_battery_real_power_transaction_magnitude_t_kW
        valdVals[:vald_battery_real_power_transaction_magnitude_allT_kW] += vald_battery_real_power_transaction_magnitude_t_kW

        # Aggregate generation data
        valdVals[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t] = valdVals[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] + valdVals[:vald_battery_reactive_power_vs_t_1toT_kVAr][t]
        valdVals[:vald_total_gen_reactive_power_allT_kVAr] += valdVals[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t]

        valdVals[:vald_total_gen_real_power_vs_t_1toT_kW][t] = valdVals[:vald_pv_real_power_vs_t_1toT_kW][t] + valdVals[:vald_battery_real_power_vs_t_1toT_kW][t]
        valdVals[:vald_total_gen_real_power_allT_kW] += valdVals[:vald_total_gen_real_power_vs_t_1toT_kW][t]

        # Retrieve voltage magnitudes
        valdVals[:vald_voltages_vs_t_1toT_pu][t] = get_voltages_opendss_powerflow_for_timestep_t()
    end

    # Battery Terminal SOC Checking
    terminalSOCDict = get_terminal_soc_values_opendss_powerflow(data)
    @unpack vald_terminal_soc_violation_kWh = terminalSOCDict
    @pack! valdVals = vald_terminal_soc_violation_kWh

    # Discrepancy Calculations
    disc_voltage_all_time_pu = compute_highest_allTime_voltage_discrepancy(modelDict, valdVals)
    line_loss_discrepancies = abs.(valdVals[:vald_PLoss_vs_t_1toT_kW] .- data[:PLoss_vs_t_1toT_kW])
    disc_line_loss_all_time_kW = maximum(line_loss_discrepancies)

    disc_PSubs_vs_t_1toT_kW = abs.(valdVals[:vald_PSubs_vs_t_1toT_kW] .- data[:PSubs_vs_t_1toT_kW])
    disc_PSubs_all_time_kW = maximum(disc_PSubs_vs_t_1toT_kW)

    disc_QSubs_vs_t_1toT_kVAr = abs.(valdVals[:vald_QSubs_vs_t_1toT_kVAr] .- data[:QSubs_vs_t_1toT_kVAr])
    disc_QSubs_all_time_kVAr = maximum(disc_QSubs_vs_t_1toT_kVAr)

    @pack! valdVals = disc_voltage_all_time_pu, disc_line_loss_all_time_kW, disc_PSubs_all_time_kW, disc_QSubs_all_time_kVAr

    @pack! modelDict = valdVals
    # return valdVals
    return modelDict
end
#endregion

end # module openDSSValidator
