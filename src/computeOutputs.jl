module computeOutputs

include("./functionRetriever.jl")
using .functionRetriever

import JuMP: value

export compute_output_values
using Parameters: @unpack, @pack!

function compute_output_values(model, data)

    # Unpack relevant data
    @unpack Tset, Lset, Bset = data;
    Tset = sort(collect(Tset))
    Lset = sort(collect(Lset))
    Bset = sort(collect(Bset))
    @unpack LoadShapeCost, rdict_pu, eta_C, eta_D, delta_t, kVA_B = data;
    
    C, r, η_C, η_D = LoadShapeCost, rdict_pu, eta_C, eta_D

    P_Subs = model[:P_Subs]
    P_c = model[:P_c]
    P_d = model[:P_d]
    P = model[:P]
    l = model[:l]

    # Initialize arrays of size T
    @unpack T = data;
    # PSubs_vs_t_1toT = zeros(Float64, T)
    PSubs_vs_t_1toT_kW = zeros(Float64, T)
    # PSubsCost_vs_t_1toT = zeros(Float64, T)
    PSubsCost_vs_t_1toT_dollar = zeros(Float64, T)
    # PLoss_vs_t_1toT = zeros(Float64, T)
    PLoss_vs_t_1toT_kW = zeros(Float64, T)
    fval_vs_t_1toT = zeros(Float64, T)
    # scd_vs_t_1toT = zeros(Float64, T)
    scd_vs_t_1toT_kW = zeros(Float64, T)

    # Scalars to store cumulative values
    # PSubs_allT = 0.0
    PSubs_allT_kW = 0.0
    # PSubsCost_allT = 0.0
    PSubsCost_allT_dollar = 0.0
    # PLoss_allT = 0.0
    PLoss_allT_kW = 0.0
    fval_allT = 0.0
    # scd_allT = 0.0
    scd_allT_kW = 0.0
    terminal_soc_violation_kWh = 0.0

    PLoss_vs_t_1toT_kW = get_real_power_loss(model, data, horizon="1toT")
    PLoss_allT_kW = get_real_power_loss(model, data, horizon="allT")

    PSubs_vs_t_1toT_kW = get_substation_power(model, data, horizon="1toT")
    PSubs_allT_kW = get_substation_power(model, data, horizon="allT")

    PSubsCost_vs_t_1toT_dollar = get_substation_power_cost(model, data, horizon="1toT")
    PSubsCost_allT_dollar = get_substation_power_cost(model, data, horizon="allT")

    scd_vs_t_1toT_kW = get_scd(model, data; horizon="1toT")
    scd_allT_kW = get_scd(model, data; horizon="allT")

    terminal_soc_violation_kWh = get_terminal_SOC_violation(model, data)

    QLoss_vs_t_1toT_kVAr = get_reactive_power_loss(model, data, horizon="1toT")
    QLoss_allT_kVAr = get_reactive_power_loss(model, data, horizon="allT")

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

        # fval_vs_t_1toT[t] = PSubsCost_vs_t_1toT[t] + scd_vs_t_1toT[t]
        fval_allT += fval_vs_t_1toT[t]

        # # Compute Terminal SOC Constraint Violation
        # @unpack T, kVA_B, Bset, Bref_pu = data;
        # B = model[:B]
        # terminal_soc_violation_kWh = kVA_B * sum(abs(value(B[j, T]) - Bref_pu[j]) for j in Bset)
    end

    # Store the computed arrays and scalars in `data` for later export
    @pack! data = fval_vs_t_1toT, fval_allT, PLoss_vs_t_1toT_kW, PLoss_allT_kW, PSubs_vs_t_1toT_kW, PSubs_allT_kW, PSubsCost_vs_t_1toT_dollar, PSubsCost_allT_dollar, scd_vs_t_1toT_kW, scd_allT_kW, terminal_soc_violation_kWh,
    QLoss_vs_t_1toT_kVAr, QLoss_allT_kVAr

    # data[:PSubs_vs_t_1toT] = PSubs_vs_t_1toT
    # data[:PSubsCost_vs_t_1toT] = PSubsCost_vs_t_1toT
    # data[:PLoss_vs_t_1toT] = PLoss_vs_t_1toT
    # data[:fval_vs_t_1toT] = fval_vs_t_1toT
    # data[:scd_vs_t_1toT] = scd_vs_t_1toT

    # data[:PSubs_allT] = PSubs_allT
    # data[:PSubsCost_allT] = PSubsCost_allT
    # data[:PLoss_allT] = PLoss_allT
    # data[:fval_allT] = fval_allT
    # data[:scd_allT] = scd_allT

    return data  # Return the updated data dictionary
end

end  # module computeOutputs
