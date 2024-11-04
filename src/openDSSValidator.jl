module openDSSValidator

using CSV
using DataFrames
using OpenDSSDirect
using Parameters: @unpack
export get_source_bus, get_substation_lines, validate_opf_against_opendss

# function validate_opf_against_opendss(model, data)
#     # Extract systemName from data dictionary
#     systemName = data[:systemName]
    
#     # Build the path to the Master.dss file
#     filename = joinpath(dirname(@__DIR__), "rawData", systemName, "Master.dss")
    
#     # Set up and run the OpenDSS commands
#     dss("""
#         clear
#         redirect "$filename"
#         // solve
#     """)

#     # Initialize load aggregation variables
#     loadnumber = Loads.First()
#     kWsum = 0.0
#     kvarsum = 0.0
    
#     # Sum up kW and kvar across all loads in OpenDSS
#     while loadnumber > 0
#         kWsum += Loads.kW()
#         kvarsum += Loads.kvar()
#         loadnumber = Loads.Next()
#     end

#     # Return the aggregated results
#     return kWsum, kvarsum
# end

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

        # Set battery power for each battery bus based on P_c and P_d values
        for battery_bus in Bset
            charge_power_kW = value(P_d[battery_bus, t]) * kVA_B
            discharge_power_kW = value(P_c[battery_bus, t]) * kVA_B

            # Calculate net power for the battery (discharge - charge)
            net_power_kW = charge_power_kW - discharge_power_kW

            # Set battery power at this timestep
            OpenDSSDirect.Circuit.SetActiveElement("Storage.$battery_bus")
            OpenDSSDirect.Storage.kW(battery_bus, net_power_kW)
        end

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

end # module openDSSValidator
