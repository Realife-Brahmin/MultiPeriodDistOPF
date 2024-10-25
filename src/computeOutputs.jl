module computeOutputs

include("./functionRetriever.jl")
using .functionRetriever

import JuMP: value

export compute_output_values
using Parameters: @unpack, @pack!

function compute_output_values(model, data)

    # Initialize arrays of size T
    @unpack T = data;
    fval_vs_t_1toT = zeros(Float64, T)

    # Scalars to store cumulative values
    fval_allT = 0.0

    PLoss_vs_t_1toT_kW = get_loss_real_power(model, data, horizon="1toT")
    PLoss_allT_kW = get_loss_real_power(model, data, horizon="allT")

    PSubs_vs_t_1toT_kW = get_substation_real_power(model, data, horizon="1toT")
    PSubs_allT_kW = get_substation_real_power(model, data, horizon="allT")

    PSubsCost_vs_t_1toT_dollar = get_substation_power_cost(model, data, horizon="1toT")
    PSubsCost_allT_dollar = get_substation_power_cost(model, data, horizon="allT")

    scd_vs_t_1toT_kW = get_scd(model, data; horizon="1toT")
    scd_allT_kW = get_scd(model, data; horizon="allT")

    terminal_soc_violation_kWh = get_terminal_SOC_violation(model, data)

    QLoss_vs_t_1toT_kVAr = get_loss_reactive_power(model, data, horizon="1toT")
    QLoss_allT_kVAr = get_loss_reactive_power(model, data, horizon="allT")

    battery_real_power_vs_t_1toT_kW = get_battery_real_power(model, data, horizon="1toT")
    battery_real_power_allT_kW = get_battery_real_power(model, data, horizon="allT")

    battery_reactive_power_vs_t_1toT_kVAr = get_battery_reactive_power(model, data, horizon="1toT")
    battery_reactive_power_allT_kVAr = get_battery_reactive_power(model, data, horizon="allT")

    pv_real_power_vs_t_1toT_kW = get_pv_real_power(model, data, horizon="1toT")
    pv_real_power_allT_kW = get_pv_real_power(model, data, horizon="allT")

    pv_reactive_power_vs_t_1toT_kVAr = get_pv_reactive_power(model, data, horizon="1toT")
    pv_reactive_power_allT_kVAr = get_pv_reactive_power(model, data, horizon="allT")

    load_real_power_vs_t_1toT_kW = get_load_real_power(model, data, horizon="1toT")
    load_real_power_allT_kW = get_load_real_power(model, data, horizon="allT")

    load_reactive_power_vs_t_1toT_kVAr = get_load_reactive_power(model, data, horizon="1toT")
    load_reactive_power_allT_kVAr = get_load_reactive_power(model, data, horizon="allT")

    # Loop over time steps to compute all required values
    @unpack Tset = data;
    for t in Tset
        @unpack kVA_B = data;

        PLoss_vs_t_1toT = PLoss_vs_t_1toT_kW./kVA_B
        PSubsCost_vs_t_1toT = PSubsCost_vs_t_1toT_dollar./kVA_B
        scd_vs_t_1toT = scd_vs_t_1toT_kW./kVA_B

        @unpack objfun0, objfun2 = data;
        if objfun0 == "genCostMin"
            fval_vs_t_1toT[t] = PSubsCost_vs_t_1toT[t]
        elseif objfun0 == "lineLossMin"
            fval_vs_t_1toT[t] = PLoss_vs_t_1toT[t]
        else
            @error "floc"
        end

        if objfun2 == "scd"
            fval_vs_t_1toT[t] += scd_vs_t_1toT[t]
        elseif objfun2 == "none"
            # do nothing
        else
            @error "floc"
        end

        fval_allT += fval_vs_t_1toT[t]

    end


    @pack! data =
        battery_real_power_allT_kW,
        battery_real_power_vs_t_1toT_kW,
        battery_reactive_power_allT_kVAr,
        battery_reactive_power_vs_t_1toT_kVAr,
    fval_allT,
        load_reactive_power_vs_t_1toT_kVAr,
        load_reactive_power_allT_kVAr,
        load_real_power_vs_t_1toT_kW,
        load_real_power_allT_kW,
        PSubsCost_allT_dollar,
        PSubsCost_vs_t_1toT_dollar,
    PSubs_allT_kW,
        PSubs_vs_t_1toT_kW,
        pv_real_power_allT_kW,
        pv_real_power_vs_t_1toT_kW,
        pv_reactive_power_allT_kVAr,
        pv_reactive_power_vs_t_1toT_kVAr,
        QLoss_allT_kVAr,
        QLoss_vs_t_1toT_kVAr,
        scd_allT_kW,
    scd_vs_t_1toT_kW,
    scd_allT_kW,
    terminal_soc_violation_kWh,
    QLoss_vs_t_1toT_kVAr,
    QLoss_allT_kVAr,
    battery_reactive_power_vs_t_1toT_kVAr,
    battery_reactive_power_allT_kVAr,
    battery_real_power_vs_t_1toT_kW,
    battery_real_power_allT_kW,
    pv_reactive_power_vs_t_1toT_kVAr,
    pv_reactive_power_allT_kVAr,
    pv_real_power_vs_t_1toT_kW,
    pv_real_power_allT_kW

    return data  # Return the updated data dictionary
end

end  # module computeOutputs
