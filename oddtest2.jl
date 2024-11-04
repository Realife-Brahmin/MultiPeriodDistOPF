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
        net_power_kW = charge_power_kW - discharge_power_kW

        # Use `Text.Command` to set battery values directly
        OpenDSSDirect.Text.Command("Edit Storage.$storage_name kW=$net_power_kW kvar=0")
        OpenDSSDirect.Text.Command("Edit Storage.$storage_name State=" * (net_power_kW >= 0 ? "2" : "1"))

        storage_id = Storages.Next()
    end

    # Solve the power flow
    OpenDSSDirect.Solution.Solve()

    # Retrieve circuit losses
    total_losses = OpenDSSDirect.Circuit.Losses() ./ 1000
    # println(total_losses)
    results.PLoss_kW[t] = real(total_losses)

    # Retrieve substation real and reactive power
    OpenDSSDirect.Circuit.SetActiveElement("Line.L1")
    substation_powers = OpenDSSDirect.CktElement.Powers()
    results.PSubs_kW[t] = real(substation_powers[1])
    results.QSubs_kVAr[t] = imag(substation_powers[1])

    # Capture voltage magnitudes at all buses
    results.Voltages[t] = OpenDSSDirect.Circuit.AllBusMagPu()

    # Print the key results for this timestep
    println("\n" * "*"^30)
    println("   Time Step: $t")
    println("*"^30)
    println("   Power Loss      : $(results.PLoss_kW[t]) kW")
    println("   Substation Power: $(results.PSubs_kW[t]) kW")
    println("   Reactive Power  : $(results.QSubs_kVAr[t]) kVAr")
    println("*"^30 * "\n")
end

# Save the results
filename = "validation_results.csv"
CSV.write(filename, results)
println("Validation results written to $filename")
