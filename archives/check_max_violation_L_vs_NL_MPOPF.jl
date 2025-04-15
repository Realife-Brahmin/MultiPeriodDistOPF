function check_max_violation(modelVals, data)
    @unpack q_B, P_c, P_d = modelVals
    @unpack S_B_R_pu, kVA_B, Bset, Tset = data

    max_violation = 0.0

    for j in Bset
        for t in Tset
            P_B = abs(P_d[j, t] - P_c[j, t])
            q_B_val = q_B[j, t]
            discrepancy = (S_B_R_pu[j]^2 - (P_B^2 + q_B_val^2))
            if discrepancy < 0
                println("j = ", j, ", t = ", t, ", discrepancy = ", discrepancy)
                max_violation = min(max_violation, discrepancy)
            else
                inactiveness = kVA_B * sqrt(discrepancy)
                println("j = $j t = $t inactiveness = $inactiveness kW")
            end


        end
    end

    println("Maximum violation: ", max_violation)
    return max_violation
end

# Example usage:
# modelVals = ... # Your modelVals dictionary with solved values
# data = ... # Your data dictionary with necessary parameters
check_max_violation(modelVals, data)