module computeOutputs

import JuMP: value

export compute_output_values
using Parameters: @unpack 

function compute_output_values(model, data)

    # Unpack relevant data
    @unpack Tset, Lset, Bset = data;
    Tset = sort(collect(Tset))
    Lset = sort(collect(Lset))
    Bset = sort(collect(Bset))
    @unpack LoadShapeCost, rdict_pu, eta_C, eta_D, delta_t = data;
    
    C, r, η_C, η_D = LoadShapeCost, rdict_pu, eta_C, eta_D

    P_Subs = model[:P_Subs]
    P_c = model[:P_c]
    P_d = model[:P_d]
    P = model[:P]
    l = model[:l]

    # Initialize arrays of size T
    @unpack T = data;
    PSubs_vs_t_1toT = zeros(Float64, T)
    PSubsCost_vs_t_1toT = zeros(Float64, T)
    PLoss_vs_t_1toT = zeros(Float64, T)
    fval_vs_t_1toT = zeros(Float64, T)
    scd_vs_t_1toT = zeros(Float64, T)

    # Scalars to store cumulative values
    PSubs_allT = 0.0
    PSubsCost_allT = 0.0
    PLoss_allT = 0.0
    fval_allT = 0.0
    scd_allT = 0.0

    # Loop over time steps to compute all required values
    for t in Tset
        # 1. Compute P_Subs over time and total PSubs_allT
        PSubs_vs_t_1toT[t] = value(P_Subs[t])
        PSubs_allT += PSubs_vs_t_1toT[t]

        # 2. Compute PSubsCost over time and total cost
        PSubsCost_vs_t_1toT[t] = C[t] * value(P_Subs[t]) * delta_t
        PSubsCost_allT += PSubsCost_vs_t_1toT[t]

        # 3. Compute Power Loss (PLoss) over time
        PLoss_vs_t_1toT[t] = sum(r[i, j] * value(l[(i, j), t]) for (i, j) in Lset)
        PLoss_allT += PLoss_vs_t_1toT[t]

        # 4. Compute SCD over time using the provided formula
        scd_vs_t_1toT[t] = sum(max(value(P_c[j, t]), value(P_d[j, t])) - abs(value(P_c[j, t]) - value(P_d[j, t])) for j in Bset)
        scd_allT += scd_vs_t_1toT[t]

        # 5. Compute Objective Function (fval) over time
        fval_vs_t_1toT[t] = PSubsCost_vs_t_1toT[t] + scd_vs_t_1toT[t]
        fval_allT += fval_vs_t_1toT[t]
    end

    # Store the computed arrays and scalars in `data` for later export
    data[:PSubs_vs_t_1toT] = PSubs_vs_t_1toT
    data[:PSubsCost_vs_t_1toT] = PSubsCost_vs_t_1toT
    data[:PLoss_vs_t_1toT] = PLoss_vs_t_1toT
    data[:fval_vs_t_1toT] = fval_vs_t_1toT
    data[:scd_vs_t_1toT] = scd_vs_t_1toT

    data[:PSubs_allT] = PSubs_allT
    data[:PSubsCost_allT] = PSubsCost_allT
    data[:PLoss_allT] = PLoss_allT
    data[:fval_allT] = fval_allT
    data[:scd_allT] = scd_allT

    return data  # Return the updated data dictionary
end

end  # module computeOutputs
