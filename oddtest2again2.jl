using OpenDSSDirect
using CSV, DataFrames
using Parameters: @unpack, @pack!
using JuMP: value
using Revise

include("src/helperFunctions.jl")
using .helperFunctions: myprintln

# include("src/openDSSValidator.jl")
includet("src/openDSSValidator.jl")
using .openDSSValidator: export_validation_decision_variables, 
    get_source_bus, 
    get_substation_lines,
    get_voltages_opendss_powerflow_for_timestep_t, 
    set_battery_controls_opendss_powerflow_for_timestep_t, 
    set_custom_load_shape!, 
    set_pv_controls_opendss_powerflow_for_timestep_t
# using .openDSSValidator

include("src/exporter.jl")
using .Exporter: export_key_validation_results

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
vald = Dict(
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

# Loop through each timestep to perform power flow and store vald
for t in 1:T

    set_pv_controls_opendss_powerflow_for_timestep_t(model, data, t)

    set_battery_controls_opendss_powerflow_for_timestep_t(model, data, t)

    Solution.Solve()

    # Retrieve circuit losses
    total_losses = Circuit.Losses() ./ 1000
    vald[:vald_PLoss_vs_t_1toT_kW][t] = real(total_losses)
    vald[:vald_PLoss_allT_kW] += real(total_losses)
    vald[:vald_QLoss_vs_t_1toT_kVAr][t] = -imag(total_losses) # maybe this is cheating but it is not a top priority for me to investigate reactive power losses in the grid
    vald[:vald_QLoss_allT_kVAr] += -imag(total_losses)

    # Substation power calculations
    substation_bus = get_source_bus()
    substation_lines = get_substation_lines(substation_bus)

    P_substation_total_t_kW = 0.0
    Q_substation_total_t_kVAr = 0.0

    for line in substation_lines
        Circuit.SetActiveElement("Line.$line")
        line_powers = CktElement.Powers()
        P_line = sum(real(line_powers[1]))
        Q_line = sum(imag(line_powers[1]))

        P_substation_total_t_kW += P_line
        Q_substation_total_t_kVAr += Q_line
    end

    Circuit.SetActiveElement("Vsource.source")
    vsource_powers = -CktElement.Powers()
    P_vsource_kW = real(vsource_powers[1])
    Q_vsource_kVAr = imag(vsource_powers[1])

    vald[:vald_PSubs_vs_t_1toT_kW][t] = P_substation_total_t_kW
    vald[:vald_PSubs_allT_kW] += P_substation_total_t_kW
    vald[:vald_QSubs_vs_t_1toT_kVAr][t] = Q_substation_total_t_kVAr
    vald[:vald_QSubs_allT_kVAr] += Q_substation_total_t_kVAr
    vald[:vald_substation_real_power_peak_allT_kW] = max(vald[:vald_substation_real_power_peak_allT_kW], P_substation_total_t_kW)

    # Cost calculation
    vald[:vald_PSubsCost_vs_t_1toT_dollar][t] = LoadShapeCost[t] * P_substation_total_t_kW * delta_t
    vald[:vald_PSubsCost_allT_dollar] += vald[:vald_PSubsCost_vs_t_1toT_dollar][t]

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

    vald[:vald_load_real_power_vs_t_1toT_kW][t] = total_load_kW
    vald[:vald_load_reactive_power_vs_t_1toT_kVAr][t] = total_load_kVAr
    vald[:vald_load_real_power_allT_kW] += total_load_kW
    vald[:vald_load_reactive_power_allT_kVAr] += total_load_kVAr

    # PV power calculations
    total_pv_kW = 0.0
    total_pv_kVAr = 0.0
    pv_names = PVsystems.AllNames()
    for pv_name in pv_names
        OpenDSSDirect.Circuit.SetActiveElement("PVSystem.$pv_name")
        total_pv_kW -= real(CktElement.Powers()[1])
        total_pv_kVAr -= imag(CktElement.Powers()[1])
    end

    vald[:vald_pv_real_power_vs_t_1toT_kW][t] = total_pv_kW

    vald[:vald_pv_real_power_allT_kW] += total_pv_kW

    vald[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] = total_pv_kVAr

    vald[:vald_pv_reactive_power_allT_kVAr] += total_pv_kVAr

    # Battery power calculations
    vald_battery_real_power_t_kW = 0.0
    vald_battery_reactive_power_t_kVAr = 0.0
    battery_names = Storages.AllNames()
    for battery_name in battery_names
        Circuit.SetActiveElement("Storage.$battery_name")
        vald_battery_real_power_t_kW += real(-CktElement.Powers()[1])
        vald_battery_reactive_power_t_kVAr += imag(-CktElement.Powers()[1])

        vald[:vald_battery_real_power_transaction_magnitude_vs_t_1toT_kW][t] += abs(vald_battery_real_power_t_kW)
        vald[:vald_battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr][t] += abs(vald_battery_reactive_power_t_kVAr)

    end

    vald[:vald_battery_real_power_vs_t_1toT_kW][t] = vald_battery_real_power_t_kW

    vald[:vald_battery_real_power_allT_kW] += vald_battery_real_power_t_kW

    vald[:vald_battery_reactive_power_vs_t_1toT_kVAr][t] = vald_battery_reactive_power_t_kVAr

    vald[:vald_battery_reactive_power_allT_kVAr] += vald_battery_reactive_power_t_kVAr

    vald[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t] = vald[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] + vald[:vald_battery_reactive_power_vs_t_1toT_kVAr][t] + vald[:vald_static_cap_reactive_power_vs_t_1toT_kVAr][t]

    vald[:vald_total_gen_reactive_power_allT_kVAr] += vald[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t]

    vald[:vald_total_gen_real_power_vs_t_1toT_kW][t] = vald[:vald_pv_real_power_vs_t_1toT_kW][t] + vald[:vald_battery_real_power_vs_t_1toT_kW][t]

    vald[:vald_total_gen_real_power_allT_kW] += vald[:vald_total_gen_real_power_vs_t_1toT_kW][t]

    vald_voltage_dict_t_pu = get_voltages_opendss_powerflow_for_timestep_t()
    # Store the dictionary in the results
    vald[:vald_voltages_vs_t_1toT_pu][t] = vald_voltage_dict_t_pu

end

# Battery Terminal SOC Checking
@unpack Bref_percent, B_R = data
vald[:vald_terminal_soc_violation_kWh] = 0.0

# Initialize storage_id with the first storage element
global storage_id = Storages.First() # It is weird that I have to specify it as a global variable here and I don't have a definite explanation on why so. I think this will get resolved once this script is in its own function

# Begin the while loop
while storage_id > 0
    storage_name = Storages.Name()
    storage_number = parse(Int, split(storage_name, "battery")[2])  # assuming 'batteryX' naming
    # Retrieve the SOC in per-unit for this battery
    Bj_T = Storages.puSOC()

    # Calculate SOC violation in kWh
    soc_violation_j_kWh = abs(Bj_T - Bref_percent[storage_number]) * B_R[storage_number]
    vald[:vald_terminal_soc_violation_kWh] += soc_violation_j_kWh

    # Move to the next storage element
    global storage_id = Storages.Next()  # Only use Storages.Next() to advance to the next storage
end

# Initialize the global discrepancy variable for voltage
global disc_voltage_all_time_pu = 0.0
v = model[:v]

# Iterate over each timestep
for t in 1:T
    # Retrieve the dictionary of bus voltages for the current timestep
    vald_voltages_dict = vald[:vald_voltages_vs_t_1toT_pu][t]

    # Iterate over each bus in the dictionary
    for (bus_index, vald_voltage) in vald_voltages_dict

        # Compute the model voltage for the bus at the current timestep
        model_voltage = sqrt(value(v[bus_index, t]))

        # Calculate the discrepancy and update the maximum discrepancy if needed
        discrepancy = abs(vald_voltage - model_voltage)
        global disc_voltage_all_time_pu
        if discrepancy > disc_voltage_all_time_pu
            disc_voltage_all_time_pu = discrepancy
        end
    end
end

println("Maximum All Time Voltage Discrepancy: ", disc_voltage_all_time_pu, " pu")

line_loss_discrepancies = abs.(vald[:vald_PLoss_vs_t_1toT_kW] .- data[:PLoss_vs_t_1toT_kW])
disc_line_loss_all_time_kW = maximum(line_loss_discrepancies)
println("Maximum All Time Line Loss Discrepancy: ", disc_line_loss_all_time_kW, " kW")

disc_PSubs_vs_t_1toT_kW = abs.(vald[:vald_PSubs_vs_t_1toT_kW] .- data[:PSubs_vs_t_1toT_kW])
disc_PSubs_all_time_kW = maximum(disc_PSubs_vs_t_1toT_kW)
println("Maximum All Time Substation Borrowed Real Power Discrepancy: ", disc_PSubs_all_time_kW, " kW")

disc_QSubs_vs_t_1toT_kVAr = abs.(vald[:vald_QSubs_vs_t_1toT_kVAr] .- data[:QSubs_vs_t_1toT_kVAr])
disc_QSubs_all_time_kVAr = maximum(disc_QSubs_vs_t_1toT_kVAr)
println("Maximum All Time Substation Borrowed Reactive Power Discrepancy: ", disc_QSubs_all_time_kVAr, " kVAr")

@pack! vald = disc_voltage_all_time_pu, disc_line_loss_all_time_kW, disc_PSubs_all_time_kW, disc_QSubs_all_time_kVAr;

export_validation_decision_variables(vald, data, verbose=true)

export_key_validation_results(vald, data)