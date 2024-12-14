# everything is already loaded, just create the for loop of comparison

@unpack data = modelDict
@unpack temporal_decmp = data

if !temporal_decmp
    @error "Cannot test for DDP as it wasn't run"
    return
end

ddpModel = modelDict
@unpack data, modelVals, mu,modelVals_ddp_vs_t_vs_k = ddpModel
@unpack Bset, Compset = data

μ = mu
t_ddp = rand(1:T)
println("t_ddp = $t_ddp")
for j in Bset
    values_mu = []
    for k_ddp in 1:3
        push!(values_mu, trim_number_for_printing(μ[j, t_ddp, k_ddp]))
    end
    println("j = $j: μ = ", join(values_mu, ", "))

end
