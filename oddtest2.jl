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
LoadShapeSim = data[:LoadShape];

# Extract battery charge (P_c) and discharge (P_d) from the model
P_c = model[:P_c];
P_d = model[:P_d];
q_D = model[:q_D];
q_B = model[:q_B];
@unpack p_D_pu = data;

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

    # Set power levels for PV systems
    pv_id = PVsystems.First()
    while pv_id > 0
        pv_name = PVsystems.Name()
        pv_number = parse(Int, split(pv_name, "pv")[2])

        # Set real and reactive power for the PV system
        p_D_t_kW = p_D_pu[pv_number][t] * kVA_B
        q_D_t_kVAr = value(q_D[pv_number, t]) * kVA_B

        println("Setting PV for bus $(pv_number) at t = $(t): p_D_t_kW = $(p_D_t_kW), q_D_t_kVAr = $(q_D_t_kVAr)")

        # Apply the settings
        PVsystems.Pmpp() = p_D_t_kW
        PVsystems.kvar() = q_D_t_kVAr

        # Re-select the element to force an update
        OpenDSSDirect.Circuit.SetActiveElement("PVSystem.$pv_name")
        println("Re-selected PVSystem.$pv_name to verify settings.")

        # Fetch and print to confirm
        actual_p_D_kW = PVsystems.Pmpp()
        actual_q_D_kVAr = PVsystems.kvar()
        println("Actual PV values after setting for bus $(pv_number) at t = $(t): kW = $(actual_p_D_kW), kvar = $(actual_q_D_kVAr)")

        # Move to the next PV
        pv_id = PVsystems.Next()
    end

    # Set battery power for each battery bus based on P_c and P_d values
    storage_id = Storages.First()
    while storage_id > 0
        storage_name = Storages.Name()
        storage_number = parse(Int, split(storage_name, "battery")[2])

        charge_power_kW = value(P_d[storage_number, t]) * kVA_B
        discharge_power_kW = value(P_c[storage_number, t]) * kVA_B
        Pdc_t_kW = charge_power_kW - discharge_power_kW
        q_B_t_kVAr = value(q_B[storage_number, t]) * kVA_B

        OpenDSSDirect.Text.Command("Edit Storage.$storage_name kW=$(Pdc_t_kW) kvar=$(q_B_t_kVAr)")
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
    total_load_kW = 0.0
    total_load_kVAr = 0.0
    total_pv_kW = 0.0
    total_pv_kVAr = 0.0
    total_battery_kW = 0.0
    total_battery_kVAr = 0.0

    # Sum up the loads
    load_id = Loads.First()
    while load_id > 0
        total_load_kW += Loads.kW() * LoadShapeSim[t]
        total_load_kVAr += Loads.kvar() * LoadShapeSim[t]
        load_id = Loads.Next()
    end

    # After solving the power flow
    Solution.Solve()

    # Retrieve all PV system names
    pv_names = PVsystems.AllNames()

    # Verify post-powerflow PV values manually
    for pv_name in pv_names
        OpenDSSDirect.Circuit.SetActiveElement("PVSystem.$pv_name")
        actual_p_D_kW = -real(CktElement.Powers()[1])   # Retrieve the real power output (kW) directly
        actual_q_D_kVAr = -imag(CktElement.Powers()[1]) # Retrieve the reactive power output (kVAr) directly
        println("Post-powerflow PV values for $pv_name: kW = $actual_p_D_kW, kvar = $actual_q_D_kVAr")
        total_pv_kW += actual_p_D_kW
        total_pv_kVAr += actual_q_D_kVAr
    end

    println("total_pv_kW_ODD = $(total_pv_kW) kW")
    println("total_pv_kVAr_ODD for t = $t = $(total_pv_kVAr) kVAr")

    # Sum up the battery storage based on power flow
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

# Save the results
filename = "validation_results.csv"
CSV.write(filename, results)
println("Validation results written to $filename")
