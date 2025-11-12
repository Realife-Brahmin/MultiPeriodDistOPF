module computeOutputs

include("./functionRetriever.jl")
using .functionRetriever

# import JuMP: value, solve_time
using JuMP

export compute_output_values
using Parameters: @unpack, @pack!

#region compute_output_values
"""
    compute_output_values(modelDict; verbose::Bool=false)

Compute and return various output values from the model.

This function extracts data from `modelDict` and computes several output variables post-optimization. 
The computed output variables are appended to the data in `modelDict`.
"""
function compute_output_values(modelDict;
    verbose::Bool=false,
    forwardPass::Bool=false)

    # @unpack model, data = modelDict;
    @unpack data = modelDict

    # Initialize arrays of size T
    @unpack T = data;
    fval_vs_t_1toT = zeros(Float64, T)

    # Scalars to store cumulative values
    fval_allT = 0.0

    PLoss_vs_t_1toT_kW = get_loss_real_power(modelDict, horizon="1toT")
    PLoss_allT_kW = get_loss_real_power(modelDict, horizon="allT")

    QSubs_vs_t_1toT_kVAr = get_substation_reactive_power(modelDict, horizon="1toT")
    QSubs_allT_kVAr = get_substation_reactive_power(modelDict, horizon="allT")

    PSubs_vs_t_1toT_kW = get_substation_real_power(modelDict, horizon="1toT")
    PSubs_allT_kW = get_substation_real_power(modelDict, horizon="allT")

    PSubsCost_vs_t_1toT_dollar = get_substation_power_cost(modelDict, horizon="1toT")
    PSubsCost_allT_dollar = get_substation_power_cost(modelDict, horizon="allT")

    scd_vs_t_1toT_kW = get_scd(modelDict; horizon="1toT")
    scd_allT_kW = get_scd(modelDict; horizon="allT")

    terminal_soc_violation_kWh = get_terminal_SOC_violation(modelDict)

    QLoss_vs_t_1toT_kVAr = get_loss_reactive_power(modelDict, horizon="1toT")
    QLoss_allT_kVAr = get_loss_reactive_power(modelDict, horizon="allT")

    battery_real_power_vs_t_1toT_kW = get_battery_real_power(modelDict, horizon="1toT")
    battery_real_power_allT_kW = get_battery_real_power(modelDict, horizon="allT")

    battery_reactive_power_vs_t_1toT_kVAr = get_battery_reactive_power(modelDict, horizon="1toT")
    battery_reactive_power_allT_kVAr = get_battery_reactive_power(modelDict, horizon="allT")

    pv_real_power_vs_t_1toT_kW = get_pv_real_power(modelDict, horizon="1toT")
    pv_real_power_allT_kW = get_pv_real_power(modelDict, horizon="allT")

    pv_reactive_power_vs_t_1toT_kVAr = get_pv_reactive_power(modelDict, horizon="1toT")
    pv_reactive_power_allT_kVAr = get_pv_reactive_power(modelDict, horizon="allT")

    load_real_power_vs_t_1toT_kW = get_load_real_power(data, horizon="1toT")
    load_real_power_allT_kW = get_load_real_power(data, horizon="allT")

    load_reactive_power_vs_t_1toT_kVAr = get_load_reactive_power(data, horizon="1toT")
    load_reactive_power_allT_kVAr = get_load_reactive_power(data, horizon="allT")

    static_cap_reactive_power_vs_t_1toT_kVAr = get_static_capacitor_reactive_power(modelDict, horizon="1toT")
    static_cap_reactive_power_allT_kVAr = get_static_capacitor_reactive_power(modelDict, horizon="allT")

    substation_real_power_peak_allT_kW = get_substation_real_power_peak(modelDict)

    total_gen_real_power_vs_t_1toT_kW = get_total_generation_real_power(modelDict, horizon="1toT")
    total_gen_real_power_allT_kW = get_total_generation_real_power(modelDict, horizon="allT")

    total_gen_reactive_power_vs_t_1toT_kVAr = get_total_generation_reactive_power(modelDict, horizon="1toT")
    total_gen_reactive_power_allT_kVAr = get_total_generation_reactive_power(modelDict, horizon="allT")

    battery_real_power_transaction_magnitude_vs_t_1toT_kW = get_battery_real_power_transaction_magnitude(modelDict, horizon="1toT")
    battery_real_power_transaction_magnitude_allT_kW = get_battery_real_power_transaction_magnitude(modelDict, horizon="allT")
    battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr = get_battery_reactive_power_transaction_magnitude(modelDict, horizon="1toT")
    battery_reactive_power_transaction_magnitude_allT_kVAr = get_battery_reactive_power_transaction_magnitude(modelDict, horizon="allT")
    
    modelSize = get_model_size(modelDict)
    @unpack num_decvars, num_lincons, num_nonlincons = modelSize;

    if !forwardPass
        solution_time = get_solution_time(modelDict)
    else 
        solution_time = 0.0
    end
    # Loop over time steps to compute all required values
    @unpack Tset = data;
    for t in Tset
        @unpack kVA_B = data;

        PLoss_vs_t_1toT = PLoss_vs_t_1toT_kW./kVA_B
        PSubsCost_vs_t_1toT = PSubsCost_vs_t_1toT_dollar./kVA_B
        scd_vs_t_1toT = scd_vs_t_1toT_kW./kVA_B

        @unpack objfun0, objfun2 = data;
        if objfun0 == "subsPowerCostMin"
            fval_vs_t_1toT[t] = PSubsCost_vs_t_1toT[t]
        elseif objfun0 == "lineLossMin"
            fval_vs_t_1toT[t] = PLoss_vs_t_1toT[t]
        elseif objfun0 == "subsPowerMin"
            fval_vs_t_1toT[t] = PSubs_vs_t_1toT_kW[t]
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
    
    # println("***********************************************")
    # println("Substation power cost: ", PSubsCost_allT_dollar)
    # println("************************************************")

    @pack! data =
        battery_reactive_power_allT_kVAr,
        battery_reactive_power_transaction_magnitude_allT_kVAr,
        battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr,
        battery_reactive_power_vs_t_1toT_kVAr,
        battery_real_power_allT_kW,
        battery_real_power_transaction_magnitude_allT_kW,
        battery_real_power_transaction_magnitude_vs_t_1toT_kW,
        battery_real_power_vs_t_1toT_kW,
        fval_allT,
        fval_vs_t_1toT,
        load_reactive_power_vs_t_1toT_kVAr,
        load_reactive_power_allT_kVAr,
        load_real_power_vs_t_1toT_kW,
        load_real_power_allT_kW,
        num_decvars,
        num_lincons,
        num_nonlincons,
        PSubsCost_allT_dollar,
        PSubsCost_vs_t_1toT_dollar,
        PSubs_allT_kW,
        PSubs_vs_t_1toT_kW,
        PLoss_allT_kW,
        PLoss_vs_t_1toT_kW,
        pv_real_power_allT_kW,
        pv_real_power_vs_t_1toT_kW,
        pv_reactive_power_allT_kVAr,
        pv_reactive_power_vs_t_1toT_kVAr,
        QLoss_allT_kVAr,
        QLoss_vs_t_1toT_kVAr,
        QSubs_allT_kVAr,
        QSubs_vs_t_1toT_kVAr,
        scd_allT_kW,
        scd_vs_t_1toT_kW,
        solution_time,
        static_cap_reactive_power_allT_kVAr,
        static_cap_reactive_power_vs_t_1toT_kVAr,
        substation_real_power_peak_allT_kW,
        terminal_soc_violation_kWh,
        total_gen_reactive_power_allT_kVAr,
        total_gen_reactive_power_vs_t_1toT_kVAr,
        total_gen_real_power_allT_kW,
        total_gen_real_power_vs_t_1toT_kW


    @pack! modelDict = data;
    return modelDict
end
#endregion

end  # module computeOutputs
