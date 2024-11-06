using OpenDSSDirect
using CSV, DataFrames
using Parameters: @unpack
using JuMP: value

# Set paths for DSS files
system_name = data[:systemName]
dss_dir = joinpath(@__DIR__, "rawData", system_name)
dss_file = joinpath(dss_dir, "Master.dss")
println(dss_file)

# Initialize OpenDSS
OpenDSSDirect.Text.Command("Clear")
OpenDSSDirect.Text.Command("Redirect \"$dss_file\"")

# Unpack data
@unpack T, kVA_B, LoadShapePV, Dset, Bset = data
LoadShapeLoad = data[:LoadShapeLoad];

# Extract battery charge (P_c) and discharge (P_d) from the model
P_c = model[:P_c];
P_d = model[:P_d];
q_D = model[:q_D];
q_B = model[:q_B];
@unpack p_D_pu = data;

# Set the custom load shape before each power flow solution
set_custom_load_shape!(LoadShapeLoad)

# Initialize results DataFrame
results = DataFrame(
    t=1:T,
    PLoss_kW=zeros(T),
    PSubs_kW=zeros(T),
    QSubs_kVAr=zeros(T),
    TotalLoad_kW=zeros(T),
    TotalLoad_kVAr=zeros(T),
    TotalPV_kW=zeros(T),
    TotalPV_kVAr=zeros(T),
    TotalBattery_kW=zeros(T),
    TotalBattery_kVAr=zeros(T),
    Voltages=Vector{Vector{Float64}}(undef, T)
)

for t in 1:T

    # Set power levels for PV systems at each time step
    pv_id = PVsystems.First()
    while pv_id > 0
        # Get the PV name and the bus number (assuming naming convention like 'pv1', 'pv2', etc.)
        pv_name = PVsystems.Name()
        pv_number = parse(Int, split(pv_name, "pv")[2])

        # Define the real and reactive power settings for the PV system
        p_D_t_kW = p_D_pu[pv_number][t] * kVA_B
        q_D_t_kVAr = value(q_D[pv_number, t]) * kVA_B

        println("Setting PV for bus $(pv_number) at t = $(t): p_D_t_kW = $(p_D_t_kW), q_D_t_kVAr = $(q_D_t_kVAr)")

        # Attempt to set the real and reactive power for the PV system
        PVsystems.Pmpp(p_D_t_kW)
        PVsystems.kvar(q_D_t_kVAr)


        # Move to the next PV system
        pv_id = PVsystems.Next()
    end

    # Iterate over each battery in OpenDSSDirect
    storage_id = Storages.First()
    while storage_id > 0
        # Get the storage element's name
        storage_name = Storages.Name()

        # Calculate the charge and discharge powers based on optimization values
        # Assuming the storage name has a numeric suffix that corresponds to the index in the optimization variables
        storage_number = parse(Int, split(storage_name, "battery")[2])

        charge_power_kW = value(P_d[storage_number, t]) * kVA_B
        discharge_power_kW = value(P_c[storage_number, t]) * kVA_B

        # Calculate net power for the battery (discharge - charge)
        net_power_kW = discharge_power_kW - charge_power_kW
        reactive_power_kVAr = value(q_B[storage_number, t]) * kVA_B

        # Construct the OpenDSS command to set the battery's active and reactive power
        command_str = "Edit Storage.$storage_name kW=$net_power_kW kvar=$reactive_power_kVAr"
        println(command_str)
        println("Executing command for storage $storage_name at t = $t: $command_str")

        # Execute the command
        OpenDSSDirect.Text.Command(command_str)

        # Move to the next storage element
        storage_id = Storages.Next()
    end

    # Solve the power flow
    Solution.Solve()

    # Retrieve circuit losses
    total_losses = Circuit.Losses() ./ 1000
    results.PLoss_kW[t] = real(total_losses)

    # Get substation bus and lines connected to it
    substation_bus = get_source_bus()
    substation_lines = get_substation_lines(substation_bus)

    # Calculate the total substation power by summing power from each line
    P_substation_total_kW = 0.0
    Q_substation_total_kVAr = 0.0

    for line in substation_lines
        Circuit.SetActiveElement("Line.$line")
        line_powers = CktElement.Powers()
        P_line = sum(real(line_powers[1]))
        Q_line = sum(imag(line_powers[1]))

        P_substation_total_kW += P_line
        Q_substation_total_kVAr += Q_line
    end

    # Also retrieve the VSource substation power for comparison
    Circuit.SetActiveElement("Vsource.source")
    vsource_powers = -CktElement.Powers()
    P_vsource_kW = real(vsource_powers[1])
    Q_vsource_kVAr = imag(vsource_powers[1])

    # Store total substation power based on line summation
    results.PSubs_kW[t] = P_substation_total_kW
    results.QSubs_kVAr[t] = Q_substation_total_kVAr

    # Calculate and store total load, PV power, and battery power

    # Initialize total load values for each timestep
    total_load_kW = 0.0
    total_load_kVAr = 0.0

    # Retrieve and sum up load outputs directly after power flow solution
    load_names = Loads.AllNames()  # Retrieve all load names for iteration

    for load_name in load_names
        OpenDSSDirect.Circuit.SetActiveElement("Load.$load_name")  # Set each load as the active element
        load_powers = CktElement.Powers()  # Retrieve the power output for the load

        # Extract real and reactive power components
        actual_load_kW = real(load_powers[1])
        actual_load_kVAr = imag(load_powers[1])

        # Accumulate the total load power
        total_load_kW += actual_load_kW
        total_load_kVAr += actual_load_kVAr
    end

    println("Total Load Power after power flow solution: kW = $(total_load_kW), kvar = $(total_load_kVAr)")

    total_pv_kW = 0.0
    total_pv_kVAr = 0.0
    # Retrieve and sum up PV system outputs after power flow solution
    pv_names = PVsystems.AllNames()
    for pv_name in pv_names
        OpenDSSDirect.Circuit.SetActiveElement("PVSystem.$pv_name")
        actual_p_D_kW = -real(CktElement.Powers()[1])
        actual_q_D_kVAr = -imag(CktElement.Powers()[1])
        total_pv_kW += actual_p_D_kW
        total_pv_kVAr += actual_q_D_kVAr
    end

    total_battery_kW = 0.0
    total_battery_kVAr = 0.0
    # Sum up the battery storage outputs after power flow
    battery_names = Storages.AllNames()
    for battery in battery_names
        Circuit.SetActiveElement("Storage.$battery")
        battery_powers = -CktElement.Powers()
        total_battery_kW += real(battery_powers[1])
        total_battery_kVAr += imag(battery_powers[1])
    end

    # Store the computed values in results DataFrame
    results.TotalLoad_kW[t] = total_load_kW
    results.TotalLoad_kVAr[t] = total_load_kVAr
    results.TotalPV_kW[t] = total_pv_kW
    results.TotalPV_kVAr[t] = total_pv_kVAr
    results.TotalBattery_kW[t] = total_battery_kW
    results.TotalBattery_kVAr[t] = total_battery_kVAr

    # Capture voltage magnitudes at all buses
    results.Voltages[t] = Circuit.AllBusMagPu()

    # Print the key results for this timestep
    println("\n" * "*"^30)
    println("   Time Step: $t")
    println("*"^30)
    println("   Power Loss              : $(results.PLoss_kW[t]) kW")
    println("   Substation Power (VSource): $P_vsource_kW kW")
    println("   Reactive Power (VSource) : $Q_vsource_kVAr kVAr")
    println("   Total Load Power        : $(results.TotalLoad_kW[t]) kW")
    println("   Total Load Reactive Power: $(results.TotalLoad_kVAr[t]) kVAr")
    println("   Total PV Power          : $(results.TotalPV_kW[t]) kW")
    println("   Total PV Reactive Power : $(results.TotalPV_kVAr[t]) kVAr")
    println("   Total Battery Power     : $(results.TotalBattery_kW[t]) kW")
    println("   Total Battery Reactive Power: $(results.TotalBattery_kVAr[t]) kVAr")
    println("*"^30 * "\n")
end

# Define the path and filename based on the specified structure
@unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunConciseDescription, simNatureAppendix = data
base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

# Create the directory if it doesn't exist
if !isdir(base_dir)
    println("Creating directory: $base_dir")
    mkpath(base_dir)
end

# Define the filename with the appropriate structure
filename = joinpath(base_dir, "Horizon_$(T)_$(machine_ID)_postsimValidation_$(gedAppendix)_for_$(objfunConciseDescription)_via_$(simNatureAppendix).txt")

CSV.write(filename, results)
println("Validation results written to $filename")
