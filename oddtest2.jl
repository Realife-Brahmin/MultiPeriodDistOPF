using OpenDSSDirect
using CSV, DataFrames
using Parameters: @unpack

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

# Extract battery charge (P_c) and discharge (P_d) from the model
P_c = model[:P_c]
P_d = model[:P_d]
q_D = model[:q_D]
q_B = model[:q_B]
@unpack p_D_pu = data

# Initialize results DataFrame
results = DataFrame(
    t=1:T,
    PLoss_kW=zeros(T),
    PSubs_kW=zeros(T),
    QSubs_kVAr=zeros(T),
    Voltages=Vector{Vector{Float64}}(undef, T)
)

for t in 1:T
    # Set power levels for PV systems
    pv_id = PVsystems.First()
    while pv_id > 0
        pv_name = PVsystems.Name()
        pv_number = parse(Int, split(pv_name, "pv")[2])

        p_D_t_kW = p_D_pu[pv_number] * kVA_B * LoadShapePV[t]
        PVsystems.kW() = p_D_t_kW
        PVsystems.kvar() = q_D[pv_number, t]
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
    OpenDSSDirect.Solution.Solve()

    # Retrieve circuit losses
    total_losses = OpenDSSDirect.Circuit.Losses() ./ 1000
    results.PLoss_kW[t] = real(total_losses)

    # Get substation bus and lines connected to it
    substation_bus = get_source_bus()
    substation_lines = get_substation_lines(substation_bus)

    # Calculate the total substation power by summing power from each line
    P_substation_total_kW = 0.0
    Q_substation_total_kVAr = 0.0

    for line in substation_lines
        OpenDSSDirect.Circuit.SetActiveElement("Line.$line")
        line_powers = OpenDSSDirect.CktElement.Powers()
        P_line = sum(real(line_powers[1]))
        Q_line = sum(imag(line_powers[1]))

        P_substation_total_kW += P_line
        Q_substation_total_kVAr += Q_line

        # println("Line: $line, P_line: $P_line kW, Q_line: $Q_line kVAr")
    end

    # Also retrieve the VSource substation power for comparison
    OpenDSSDirect.Circuit.SetActiveElement("Vsource.source")
    vsource_powers = -OpenDSSDirect.CktElement.Powers()
    P_vsource_kW = real(vsource_powers[1])
    Q_vsource_kVAr = imag(vsource_powers[1])

    # Store total substation power based on line summation
    results.PSubs_kW[t] = P_substation_total_kW
    results.QSubs_kVAr[t] = Q_substation_total_kVAr

    # Capture voltage magnitudes at all buses
    results.Voltages[t] = OpenDSSDirect.Circuit.AllBusMagPu()

    # Print the key results for this timestep
    println("\n" * "*"^30)
    println("   Time Step: $t")
    println("*"^30)
    println("   Power Loss           : $(results.PLoss_kW[t]) kW")
    println("   Substation Power (Lines): $P_substation_total_kW kW")
    println("   Reactive Power (Lines)  : $Q_substation_total_kVAr kVAr")
    println("   Substation Power (VSource): $P_vsource_kW kW")
    println("   Reactive Power (VSource) : $Q_vsource_kVAr kVAr")
    println("*"^30 * "\n")
end

# Save the results
filename = "validation_results.csv"
CSV.write(filename, results)
println("Validation results written to $filename")
