using OpenDSSDirect
using CSV, DataFrames
using Parameters: @unpack, @pack!
using JuMP: value
using Revise

include("src/helperFunctions.jl")
using .helperFunctions: myprintln

# include("src/openDSSValidator.jl")
includet("src/openDSSValidator.jl")
    using .openDSSValidator: compute_highest_allTime_voltage_discrepancy,
    export_validation_decision_variables,
    get_battery_powers_opendss_powerflow_for_timestep_t,
    get_load_powers_opendss_powerflow_for_timestep_t,
    get_pv_powers_opendss_powerflow_for_timestep_t,
    get_source_bus, 
    get_substation_lines,
    get_substation_powers_opendss_powerflow_for_timestep_t,
    get_terminal_soc_values_opendss_powerflow,
    get_voltages_opendss_powerflow_for_timestep_t, 
    set_battery_controls_opendss_powerflow_for_timestep_t, 
    set_custom_load_shape!, 
    set_pv_controls_opendss_powerflow_for_timestep_t
# using .openDSSValidator

include("src/exporter.jl")
using .Exporter: 
    export_validation_key_results

verbose = false
useVSourcePower = true
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
    totalLosses_t_kVA = Circuit.Losses() ./ 1000
    totalLosses_t_kW, totalLosses_t_kVAr = real(totalLosses_t_kVA), -imag(totalLosses_t_kVA) # maybe this is cheating but it is not a top priority for me to investigate reactive power losses in the grid
    vald[:vald_PLoss_vs_t_1toT_kW][t] = totalLosses_t_kW
    vald[:vald_PLoss_allT_kW] += totalLosses_t_kW
    vald[:vald_QLoss_vs_t_1toT_kVAr][t] = totalLosses_t_kVAr 
    vald[:vald_QLoss_allT_kVAr] += totalLosses_t_kVAr

    # Retrieve substation real and reactive powers post powerflow for this timestep
    substationPowersDict_t = get_substation_powers_opendss_powerflow_for_timestep_t(data, useVSourcePower=useVSourcePower)
    @unpack P_substation_total_t_kW, Q_substation_total_t_kVAr = substationPowersDict_t;

    vald[:vald_PSubs_vs_t_1toT_kW][t] = P_substation_total_t_kW
    vald[:vald_PSubs_allT_kW] += P_substation_total_t_kW
    vald[:vald_QSubs_vs_t_1toT_kVAr][t] = Q_substation_total_t_kVAr
    vald[:vald_QSubs_allT_kVAr] += Q_substation_total_t_kVAr
    vald[:vald_substation_real_power_peak_allT_kW] = max(vald[:vald_substation_real_power_peak_allT_kW], P_substation_total_t_kW)

    # Cost calculation
    vald[:vald_PSubsCost_vs_t_1toT_dollar][t] = LoadShapeCost[t] * P_substation_total_t_kW * delta_t
    vald[:vald_PSubsCost_allT_dollar] += vald[:vald_PSubsCost_vs_t_1toT_dollar][t]

    # Retrieve load real and reactive powers post powerflow for this timestep
    loadPowersDict_t = get_load_powers_opendss_powerflow_for_timestep_t()
    @unpack total_load_t_kW, total_load_t_kVAr = loadPowersDict_t;

    vald[:vald_load_real_power_vs_t_1toT_kW][t] = total_load_t_kW
    vald[:vald_load_reactive_power_vs_t_1toT_kVAr][t] = total_load_t_kVAr
    vald[:vald_load_real_power_allT_kW] += total_load_t_kW
    vald[:vald_load_reactive_power_allT_kVAr] += total_load_t_kVAr

    pvPowersDict_t = get_pv_powers_opendss_powerflow_for_timestep_t()
    @unpack total_pv_t_kW, total_pv_t_kVAr = pvPowersDict_t;

    vald[:vald_pv_real_power_vs_t_1toT_kW][t] = total_pv_t_kW
    vald[:vald_pv_real_power_allT_kW] += total_pv_t_kW
    vald[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] = total_pv_t_kVAr
    vald[:vald_pv_reactive_power_allT_kVAr] += total_pv_t_kVAr

    batteryPowersDict_t = get_battery_powers_opendss_powerflow_for_timestep_t(verbose=verbose)

    @unpack vald_battery_real_power_t_kW, vald_battery_reactive_power_t_kVAr, vald_battery_real_power_transaction_magnitude_t_kW, vald_battery_reactive_power_transaction_magnitude_t_kVAr = batteryPowersDict_t;

    vald[:vald_battery_reactive_power_vs_t_1toT_kVAr][t] = vald_battery_reactive_power_t_kVAr
    vald[:vald_battery_reactive_power_allT_kVAr] += vald_battery_reactive_power_t_kVAr
    vald[:vald_battery_reactive_power_transaction_magnitude_t_kVAr] = vald_battery_reactive_power_transaction_magnitude_t_kVAr
    vald[:vald_battery_reactive_power_transaction_magnitude_allT_kVAr] += vald_battery_reactive_power_transaction_magnitude_t_kVAr
    vald[:vald_battery_real_power_vs_t_1toT_kW][t] = vald_battery_real_power_t_kW
    vald[:vald_battery_real_power_allT_kW] += vald_battery_real_power_t_kW
    vald[:vald_battery_real_power_transaction_magnitude_t_kW] = vald_battery_real_power_transaction_magnitude_t_kW
    vald[:vald_battery_real_power_transaction_magnitude_allT_kW] += vald_battery_real_power_transaction_magnitude_t_kW

    vald[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t] = vald[:vald_pv_reactive_power_vs_t_1toT_kVAr][t] + vald[:vald_battery_reactive_power_vs_t_1toT_kVAr][t] + vald[:vald_static_cap_reactive_power_vs_t_1toT_kVAr][t]

    vald[:vald_total_gen_reactive_power_allT_kVAr] += vald[:vald_total_gen_reactive_power_vs_t_1toT_kVAr][t]

    vald[:vald_total_gen_real_power_vs_t_1toT_kW][t] = vald[:vald_pv_real_power_vs_t_1toT_kW][t] + vald[:vald_battery_real_power_vs_t_1toT_kW][t]

    vald[:vald_total_gen_real_power_allT_kW] += vald[:vald_total_gen_real_power_vs_t_1toT_kW][t]

    vald_voltage_dict_t_pu = get_voltages_opendss_powerflow_for_timestep_t()
    # Store the dictionary in the results
    vald[:vald_voltages_vs_t_1toT_pu][t] = vald_voltage_dict_t_pu

end

# Battery Terminal SOC Checking
terminalSOCDict = get_terminal_soc_values_opendss_powerflow(data)
@unpack vald_terminal_soc_violation_kWh = terminalSOCDict
@pack! vald = vald_terminal_soc_violation_kWh

# Checking discrepancies in value of state/output variables between OpenDSS and optsim model
disc_voltage_all_time_pu = compute_highest_allTime_voltage_discrepancy(model, data, vald)
myprintln(verbose, "Maximum voltage discrepancy across all timesteps and buses: $disc_voltage_all_time_pu pu")

line_loss_discrepancies = abs.(vald[:vald_PLoss_vs_t_1toT_kW] .- data[:PLoss_vs_t_1toT_kW])
disc_line_loss_all_time_kW = maximum(line_loss_discrepancies)
myprintln(verbose, "Maximum All Time Line Loss Discrepancy: $(disc_line_loss_all_time_kW) kW")

disc_PSubs_vs_t_1toT_kW = abs.(vald[:vald_PSubs_vs_t_1toT_kW] .- data[:PSubs_vs_t_1toT_kW])
disc_PSubs_all_time_kW = maximum(disc_PSubs_vs_t_1toT_kW)
myprintln(verbose, "Maximum All Time Substation Borrowed Real Power Discrepancy: $(disc_PSubs_all_time_kW) kW")

disc_QSubs_vs_t_1toT_kVAr = abs.(vald[:vald_QSubs_vs_t_1toT_kVAr] .- data[:QSubs_vs_t_1toT_kVAr])
disc_QSubs_all_time_kVAr = maximum(disc_QSubs_vs_t_1toT_kVAr)
myprintln(verbose, "Maximum All Time Substation Borrowed Reactive Power Discrepancy: $(disc_QSubs_all_time_kVAr) kVAr")

@pack! vald = disc_voltage_all_time_pu, disc_line_loss_all_time_kW, disc_PSubs_all_time_kW, disc_QSubs_all_time_kVAr;

export_validation_decision_variables(vald, data, verbose=true)

export_validation_key_results(vald, data)