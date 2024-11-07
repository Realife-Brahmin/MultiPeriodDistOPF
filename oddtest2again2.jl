using OpenDSSDirect
using CSV, DataFrames
using Parameters: @unpack
using JuMP: value

include("src/helperFunctions.jl")
using .helperFunctions: myprintln

include("src/openDSSValidator.jl")
using .openDSSValidator: get_source_bus, get_substation_lines

verbose = false

# Set paths for DSS files
system_name = data[:systemName]
dss_dir = joinpath(@__DIR__, "rawData", system_name)
dss_file = joinpath(dss_dir, "Master.dss")
myprintln(verbose, "Master.dss file path: $(dss_file)")

# Initialize OpenDSS
OpenDSSDirect.Text.Command("Clear")
OpenDSSDirect.Text.Command("Redirect \"$dss_file\"")

# Unpack data
@unpack T, kVA_B, LoadShapePV, Dset, Bset = data
LoadShapeLoad = data[:LoadShapeLoad]

# Extract controllables from the solved optimization model
P_c = model[:P_c]
P_d = model[:P_d]
q_D = model[:q_D]
q_B = model[:q_B]
@unpack p_D_pu, delta_t, objfun0, objfun2, alpha, LoadShapeCost = data;

# Set custom load shape
set_custom_load_shape!(LoadShapeLoad)

# Initialize Dictionary for timestep-specific and cumulative results
results = Dict(
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
    :vald_voltages_vs_t_1toT_pu => Vector{Vector{Float64}}(undef, T),
    :vald_PSubsCost_vs_t_1toT_dollar => zeros(T),
    :vald_battery_real_power_transaction_magnitude_vs_t_1toT_kW => zeros(T),
    :vald_battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr => zeros(T),

    # Cumulative fields
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

# Loop through each timestep to perform power flow and store results
for t in 1:T
    # Set power levels for PV systems at each time step
    pv_id = PVsystems.First()
    while pv_id > 0
        pv_name = PVsystems.Name()
        pv_number = parse(Int, split(pv_name, "pv")[2])

        p_D_t_kW = p_D_pu[pv_number][t] * kVA_B
        q_D_t_kVAr = value(q_D[pv_number, t]) * kVA_B
        PVsystems.Pmpp(p_D_t_kW)
        PVsystems.kvar(q_D_t_kVAr)

        pv_id = PVsystems.Next()
    end

    # Set battery power levels
    storage_id = Storages.First()
    while storage_id > 0
        storage_name = Storages.Name()
        storage_number = parse(Int, split(storage_name, "battery")[2])

        charge_power_kW = value(P_c[storage_number, t]) * kVA_B
        discharge_power_kW = value(P_d[storage_number, t]) * kVA_B
        net_power_kW = discharge_power_kW - charge_power_kW
        reactive_power_kVAr = value(q_B[storage_number, t]) * kVA_B

        OpenDSSDirect.Text.Command("Edit Storage.Battery$(storage_number) kW=$(net_power_kW) kvar=$(reactive_power_kVAr)")
        storage_id = Storages.Next()
    end

    Solution.Solve()

    # Retrieve circuit losses
    total_losses = Circuit.Losses() ./ 1000
    results[:vald_PLoss_vs_t_1toT_kW][t] = real(total_losses)
    results[:vald_PLoss_allT_kW] += real(total_losses)
    results[:vald_QLoss_vs_t_1toT_kVAr][t] = -imag(total_losses) # maybe this is cheating but it is not a top priority for me to investigate reactive power losses in the grid
    results[:vald_QLoss_allT_kVAr] += -imag(total_losses)

    # Substation power calculations
    substation_bus = get_source_bus()
    substation_lines = get_substation_lines(substation_bus)

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

    Circuit.SetActiveElement("Vsource.source")
    vsource_powers = -CktElement.Powers()
    P_vsource_kW = real(vsource_powers[1])
    Q_vsource_kVAr = imag(vsource_powers[1])

    results[:vald_PSubs_vs_t_1toT_kW][t] = P_substation_total_kW
    results[:vald_PSubs_allT_kW] += P_substation_total_kW
    results[:vald_QSubs_vs_t_1toT_kVAr][t] = Q_substation_total_kVAr
    results[:vald_QSubs_allT_kVAr] += Q_substation_total_kVAr
    results[:vald_substation_real_power_peak_allT_kW] = max(results[:vald_substation_real_power_peak_allT_kW], P_vsource_kW)

    # Cost calculation
    results[:vald_PSubsCost_vs_t_1toT_dollar][t] = kVA_B * LoadShapeCost[t] * P_substation_total_kW * delta_t
    results[:vald_PSubsCost_allT_dollar] += results[:vald_PSubsCost_vs_t_1toT_dollar][t]

    # Reactive and real power totals
    total_load_kW = 0.0
    total_load_kVAr = 0.0
    load_names = Loads.AllNames()
    for load_name in load_names
        OpenDSSDirect.Circuit.SetActiveElement("Load.$load_name")
        load_powers = CktElement.Powers()
        total_load_kW += real(load_powers[1])
        total_load_kVAr += imag(load_powers[1])
    end

    results[:vald_load_real_power_vs_t_1toT_kW][t] = total_load_kW
    results[:vald_load_reactive_power_vs_t_1toT_kVAr][t] = total_load_kVAr
    results[:vald_load_real_power_allT_kW] += total_load_kW
    results[:vald_load_reactive_power_allT_kVAr] += total_load_kVAr

    # PV power calculations
    total_pv_kW = 0.0
    total_pv_kVAr = 0.0
    pv_names = PVsystems.AllNames()
    for pv_name in pv_names
        OpenDSSDirect.Circuit.SetActiveElement("PVSystem.$pv_name")
        total_pv_kW -= real(CktElement.Powers()[1])
        total_pv_kVAr -= imag(CktElement.Powers()[1])
    end

    results[:vald_pv_real_power_vs_t_1toT_kW][t] = total_pv_kW
    results[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] = total_pv_kVAr

    # Battery power calculations
    total_battery_kW = 0.0
    total_battery_kVAr = 0.0
    battery_names = Storages.AllNames()
    for battery_name in battery_names
        Circuit.SetActiveElement("Storage.$battery_name")
        total_battery_kW += real(-CktElement.Powers()[1])
        total_battery_kVAr += imag(-CktElement.Powers()[1])
    end

    results[:vald_battery_real_power_vs_t_1toT_kW][t] = total_battery_kW
    results[:vald_battery_reactive_power_vs_t_1toT_kVAr][t] = total_battery_kVAr

    results[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t] = results[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] + results[:vald_battery_reactive_power_vs_t_1toT_kVAr][t] + results[:vald_static_cap_reactive_power_vs_t_1toT_kVAr][t]

    results[:vald_total_gen_reactive_power_allT_kVAr] += results[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t]

    results[:vald_total_gen_real_power_vs_t_1toT_kW][t] = results[:vald_pv_real_power_vs_t_1toT_kW][t] + results[:vald_battery_real_power_vs_t_1toT_kW][t]

    # Other cumulative battery metrics
    results[:vald_battery_real_power_transaction_magnitude_vs_t_1toT_kW][t] = abs(total_battery_kW)
    results[:vald_battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr][t] = abs(total_battery_kVAr)

    # Capture voltage magnitudes at all buses
    results[:vald_voltages_vs_t_1toT_pu][t] = Circuit.AllBusMagPu()

    # Print key results for this timestep
    println("\n" * "*"^30)
    println("   Time Step: $t")
    println("*"^30)
    println("   Power Loss              : $(results[:vald_PLoss_vs_t_1toT_kW][t]) kW")
    println("   Substation Power (VSource): $P_vsource_kW kW")
    println("   Reactive Power (VSource) : $Q_vsource_kVAr kVAr")
    println("   Total Load Power        : $(results[:vald_load_real_power_vs_t_1toT_kW][t]) kW")
    println("   Total Load Reactive Power: $(results[:vald_load_reactive_power_vs_t_1toT_kVAr][t]) kVAr")
    println("   Total PV Power          : $(results[:vald_pv_real_power_vs_t_1toT_kW][t]) kW")
    println("   Total PV Reactive Power : $(results[:vald_pv_reactive_power_vs_t_1toT_kVAr][t]) kVAr")
    println("   Total Battery Power     : $(results[:vald_battery_real_power_vs_t_1toT_kW][t]) kW")
    println("   Total Battery Reactive Power: $(results[:vald_battery_reactive_power_vs_t_1toT_kVAr][t]) kVAr")
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
myprintln(verbose, "Validation results written to $filename")
