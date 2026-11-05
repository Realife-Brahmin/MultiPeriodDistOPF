using Serialization, LinearAlgebra, Printf, TOML

reference, candidate = deserialize.(ARGS[1:2])
result = Dict{String,Any}("reference_status"=>reference[:status],
    "candidate_status"=>candidate[:status],
    "reference_iterations"=>reference[:iterations],
    "candidate_iterations"=>candidate[:iterations])
same_iterations = reference[:iterations] == candidate[:iterations]
strict = reference[:status] == candidate[:status] == 0 && same_iterations
for key in (:x, :u)
    @assert length(reference[key]) == length(candidate[key])
    worst_scaled = 0.0
    worst_absolute = 0.0
    for (a,b) in zip(reference[key],candidate[key])
        @assert size(a) == size(b) && all(isfinite,a) && all(isfinite,b)
        absolute = norm(a-b,Inf)
        scaled = absolute/max(1,norm(a,Inf),norm(b,Inf))
        worst_scaled=max(worst_scaled,scaled)
        worst_absolute=max(worst_absolute,absolute)
    end
    result["$(key)_scaled_inf_difference"]=worst_scaled
    result["$(key)_absolute_inf_difference"]=worst_absolute
    strict &= worst_scaled <= 1e-10
end
result["strict_existing_match_gate"]=strict
open(ARGS[3],"w") do io; TOML.print(io,result); end
println("FULL_STUDY_SOLUTION_COMPARISON ",result)
# Preserve a failed equivalence test as a reported scientific result. Each arm
# must independently pass the original numerical/physical/OpenDSS validators.
