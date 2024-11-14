module openDSSValidator

export export_validation_decision_variables, get_source_bus, get_substation_lines, set_custom_load_shape!, set_battery_controls_for_timestep_t, set_pv_controls_for_timestep_t,
validate_opf_against_opendss

using CSV
using DataFrames
using JuMP: value
using OpenDSSDirect
using Parameters: @unpack

include("helperFunctions.jl")
using .helperFunctions: myprintln

function validate_opf_against_opendss(model, data; filename="validation_results.csv")
    # Set paths for DSS files
    system_name = data[:systemName]
    dss_dir = joinpath(dirname(@__DIR__), "rawData", system_name)
    dss_file = joinpath(dss_dir, "Master.dss")
    println(dss_file)
    # Initialize OpenDSS
    OpenDSSDirect.Text.Command("Clear")
    # OpenDSSDirect.Text.Command("Redirect $dss_file")
    OpenDSSDirect.Text.Command("Redirect '$dss_file'")


    # List all components in the circuit
    component_names = OpenDSSDirect.Circuit.AllElementNames()

    # Display component names and types
    for component in component_names
        println("Component: $component")
    end

    # Unpack data
    T = data[:T]
    kVA_B = data[:kVA_B]
    load_shape_pv = data[:LoadShapePV]
    @unpack Dset, Bset = data  # PV and Battery bus sets

    # Extract battery charge (P_c) and discharge (P_d) from the model
    P_c = model[:P_c]
    P_d = model[:P_d]

    # Initialize results DataFrame
    results = DataFrame(
        t=1:T,
        PLoss_kW=zeros(T),
        PSubs_kW=zeros(T),
        QSubs_kVAr=zeros(T),
        Voltages=Vector{Vector{Float64}}(undef, T)
    )

    for t in 1:T
        # Set power levels for PVs based on current timestep
        for pv_bus in Dset
            OpenDSSDirect.Circuit.SetActiveElement("PVSystem.pv$(pv_bus)")
            OpenDSSDirect.PVSystem.kW(pv_bus, load_shape_pv[t] * kVA_B)
        end

        # # Set battery power for each battery bus based on P_c and P_d values
        # for battery_bus in Bset
        #     charge_power_kW = value(P_d[battery_bus, t]) * kVA_B
        #     discharge_power_kW = value(P_c[battery_bus, t]) * kVA_B

        #     # Calculate net power for the battery (discharge - charge)
        #     net_power_kW = charge_power_kW - discharge_power_kW

        #     # Set battery power at this timestep
        #     OpenDSSDirect.Circuit.SetActiveElement("Storage.$battery_bus")
        #     OpenDSSDirect.Storage.kW(battery_bus, net_power_kW)
        # end

        # Solve the power flow
        OpenDSSDirect.Solution.Solve()

        # Retrieve circuit losses
        total_losses = OpenDSSDirect.Circuit.Losses() ./ 1000  # Convert W to kW
        results.PLoss_kW[t] = total_losses[1]

        # Retrieve substation real and reactive power
        OpenDSSDirect.Circuit.SetActiveElement("Line.L1")
        substation_powers = OpenDSSDirect.CktElement.Powers()
        results.PSubs_kW[t] = sum(substation_powers[1:2:end])  # Summing real power across phases
        results.QSubs_kVAr[t] = sum(substation_powers[2:2:end])  # Summing reactive power across phases

        # Capture voltage magnitudes at all buses
        results.Voltages[t] = OpenDSSDirect.Circuit.AllBusVmagPu()

        # Print the key results for this timestep
        println("Time: $t, PLoss: $(results.PLoss_kW[t]) kW, PSubs: $(results.PSubs_kW[t]) kW, QSubs: $(results.QSubs_kVAr[t]) kVAr")
    end

    # Save the results
    CSV.write(filename, results)
    println("Validation results written to $filename")
end

# Todo: This function should be in exporter instead of this mod
function export_validation_decision_variables(vald, data; verbose::Bool=false)

    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunConciseDescription, processedDataFolderPath, simNatureAppendix, solver = data
    base_dir = joinpath(processedDataFolderPath, systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

    # Create the directory if it doesn't exist
    if !isdir(base_dir)
        if verbose
            println("Creating directory: $base_dir")
        end
        mkpath(base_dir)
    end

    # Define the filename with the appropriate structure
    filename = joinpath(base_dir, "Horizon_$(T)_$(machine_ID)_$(solver)_validationDecisionVariables_$(gedAppendix)_for_$(objfunConciseDescription)_via_$(simNatureAppendix).txt")

    # Write the vald dictionary to a CSV file
    CSV.write(filename, vald)

    # Print confirmation if verbose is enabled
    if verbose
        println("Validation decision variables written to $filename")
    end
end

function get_source_bus()
    vsource_element = OpenDSSDirect.Vsources.First()
    vsource_name = OpenDSSDirect.Vsources.Name()
    OpenDSSDirect.Circuit.SetActiveElement("Vsource.$vsource_name")
    source_bus = OpenDSSDirect.CktElement.BusNames()[1]
    return source_bus
end

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

function set_custom_load_shape!(LoadShapeArray::Vector{Float64};
    verbose::Bool=false)
    # Define LoadShapeLoad in OpenDSS with the provided LoadShapeArray array
    loadshape_command = "New Loadshape.LoadShapeLoad npts = $(length(LoadShapeArray)) interval = 1 mult = [" *
                        join(LoadShapeArray, " ") * "]"
    OpenDSSDirect.Text.Command(loadshape_command)
    myprintln(verbose, "Defined LoadShapeLoad with provided LoadShapeArray")

    # Apply this load shape to all loads in the system
    load_id = OpenDSSDirect.Loads.First()
    while load_id > 0
        OpenDSSDirect.Loads.Daily("LoadShapeLoad")
        load_id = OpenDSSDirect.Loads.Next()
    end
    myprintln(verbose, "Applied LoadShapeLoad to all loads")
end

function set_battery_controls_for_timestep_t(model, data, t; verbose=false)
    # Unpack necessary data
    P_c = model[:P_c]
    P_d = model[:P_d]
    q_B = model[:q_B]
    @unpack kVA_B = data

    # Set battery power levels
    storage_id = OpenDSSDirect.Storages.First()
    while storage_id > 0
        storage_name = OpenDSSDirect.Storages.Name()
        storage_number = parse(Int, split(storage_name, "battery")[2])

        # Calculate power levels based on optimization model variables
        charge_power_kW = value(P_c[storage_number, t]) * kVA_B
        discharge_power_kW = value(P_d[storage_number, t]) * kVA_B
        net_power_kW = discharge_power_kW - charge_power_kW
        reactive_power_kVAr = value(q_B[storage_number, t]) * kVA_B

        # Command to set battery power levels in OpenDSS
        command_str = "Edit Storage.Battery$(storage_number) kW=$(net_power_kW) kvar=$(reactive_power_kVAr)"
        OpenDSSDirect.Text.Command(command_str)

        # Optionally print command for verification
        if verbose
            println("Time Step $t: Setting battery $storage_number with command: $command_str")
        end

        # Move to the next storage element
        storage_id = OpenDSSDirect.Storages.Next()
    end
end

function set_pv_controls_for_timestep_t(model, data, t; verbose::Bool=false)
    # Unpack necessary data fields from `data`
    @unpack kVA_B, p_D_pu = data
    q_D = model[:q_D]  # Access q_D from the model

    # Set power levels for PV systems at each time step
    pv_id = OpenDSSDirect.PVsystems.First()
    while pv_id > 0
        pv_name = OpenDSSDirect.PVsystems.Name()
        pv_number = parse(Int, split(pv_name, "pv")[2])

        # Retrieve real and reactive power setpoints for this PV system and timestep
        p_D_t_kW = p_D_pu[pv_number][t] * kVA_B
        q_D_t_kVAr = value(q_D[pv_number, t]) * kVA_B

        # Set real and reactive power for the PV system
        OpenDSSDirect.PVsystems.Pmpp(p_D_t_kW)
        OpenDSSDirect.PVsystems.kvar(q_D_t_kVAr)

        if verbose
            println("Setting PV for bus $(pv_number) at t = $(t): p_D_t_kW = $(p_D_t_kW), q_D_t_kVAr = $(q_D_t_kVAr)")
        end

        # Move to the next PV system
        pv_id = OpenDSSDirect.PVsystems.Next()
    end
end



end # module openDSSValidator
