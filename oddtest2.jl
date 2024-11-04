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
@unpack T, kVA_B, LoadShapePV, Dset, Bset = data;
# T = data[:T]
# kVA_B = data[:kVA_B]
# load_shape_pv = data[:LoadShapePV]
# @unpack Dset, Bset = data  # PV and Battery bus sets

# Extract battery charge (P_c) and discharge (P_d) from the model
P_c = model[:P_c]
P_d = model[:P_d]
q_D = model[:q_D]
# p_D_pu = data[:p_D_pu]
@unpack p_D_pu = data;

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
        pv_number = parse(Int, split(pv_name, "pv")[2])  # Extract numeric part from name

        p_D_t_kW = p_D_pu[pv_number] * kVA_B * LoadShapePV[t] 
        # pv_kW = value(p_D[pv_number, t]) * kVA_B
        PVsystems.kW() = p_D_t_kW  # Set real power output for the PV system
        PVsystems.kvar() = q_D[pv_number, t]  # Set reactive power output for the PV system if needed
        pv_id = PVsystems.Next()
    end

    # Set battery power for each battery bus based on P_c and P_d values
    storage_id = Storages.First()
    while storage_id > 0
        storage_name = Storages.Name()
        storage_number = parse(Int, split(storage_name, "battery")[2])  # Extract numeric part from name

        charge_power_kW = value(P_d[storage_number, t]) * kVA_B
        discharge_power_kW = value(P_c[storage_number, t]) * kVA_B
        net_power_kW = charge_power_kW - discharge_power_kW

        Storages.kW() = net_power_kW  # Set net real power for the storage
        Storages.kvar(0.0)         # Set reactive power for the storage if needed
        storage_id = Storages.Next()
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
