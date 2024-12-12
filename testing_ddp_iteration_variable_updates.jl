# Give me a quick script which compares modelVals_ddp_vs_t_vs_k for the same t_ddp for all k_ddp(s) for the variable B
# For example, modelVals_ddp_vs_t_vs_k[1, 1][:B] gives the values for variable B for the first t_ddp and the first k_ddp. Use a for loop so that for all k_ddp for 1 to 3, the same index in Bset is used, so something like modelVals_ddp_vs_t_vs_k[1, 1][:B][j], modelVals_ddp_vs_t_vs_k[1, 2][:B][j] and modelVals_ddp_vs_t_vs_k[1, 3][:B][j] are compared for the same j ∈ Bset.

# everything is already loaded, just create the for loop of comparison
ddpModel = modelDict
@unpack data, modelVals, modelVals_ddp_vs_t_vs_k = ddpModel
@unpack Bset, Compset = data

t_ddp = 3
for j in Compset
    values = []
    for k_ddp in 1:3
        push!(values, modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:v][j, t_ddp])
    end
    println("j = $j: ", join(values, ", "))
end
