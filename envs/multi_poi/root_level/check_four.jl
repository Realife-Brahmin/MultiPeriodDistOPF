# The four-substation regions (multi_windows.csv, stages four_*) against the predicted tetrahedron, face by
# face. A face is the set of measured boundary points x = r u whose predicted limiting substation is k.
# For each face: a plane a.x = 1 fitted to its points (least squares), the RMS distance of the points from it
# relative to its distance from the origin (0 = planar), the angle between its normal and the predicted
# face's normal, and the ratio of the two planes' distances from the origin (1 = on the prediction, >1 =
# outward). The predicted normal is the daily-mean sensitivity row for cyclic batteries; for idle batteries
# the face is the hour that binds most often among that face's directions.
#   julia --project=envs/msopf envs/multi_poi/root_level/check_four.jl
using Printf, LinearAlgebra, Statistics
const RES = joinpath(@__DIR__, "..", "results")
rows(f) = [split(l, ',') for l in readlines(f)[2:end]]

function sens_for(stage)
    f = first(filter(p -> occursin("multi_sens_$(stage)_T24.csv", p), [joinpath(r, x) for (r, _, fs) in walkdir(joinpath(RES, "kamiak", "four_oct20")) for x in fs]))
    P = Dict{String,Vector{Float64}}(); K = Dict{String,Vector{Vector{Float64}}}()
    for r in rows(f)
        push!(get!(P, r[1], Float64[]), parse(Float64, r[3])); push!(get!(K, r[1], Vector{Float64}[]), parse.(Float64, r[4:end]))
    end
    return P, K
end

const STAGE = Dict("subs1:subs2:subs3:subs4" => "four_ieee123a", "subs1:subs3:subs4:subs5" => "four_ieee123b", "source:subs2:subs3:subs4" => "four_ieee2522")
data = rows(joinpath(RES, "multiperiod", "multi_windows.csv"))
for key in unique((r[1], r[2]) for r in data), case in ("nobatt", "cyclic")
    P, K = sens_for(STAGE[key[2]])
    sel = [r for r in data if (r[1], r[2]) == key && r[4] == case && r[12] == "false"]
    println("\n", key[1], " ", key[2], "  ", case)
    for k in sort(unique(r[11] for r in sel))
        pts = [parse(Float64, r[9]) .* parse.(Float64, r[6:8]) for r in sel if r[11] == k]
        length(pts) >= 4 || (println(@sprintf("  face %-6s %d points: too few to fit", k, length(pts))); continue)
        X = reduce(hcat, pts)'; a = X \ ones(length(pts))
        d = 1 / norm(a); rms = sqrt(mean(((X * a) .- 1) .^ 2)) / norm(a)
        if case == "cyclic"
            n_pred = -sum(K[k]) ; c_pred = sum(P[k])                     # face: c + n.x = 0 with n = -sum_t K
        else
            # the binding hour of this face: the most frequent arg-min hour along its directions
            hrs = [argmin([K[k][t]' * (x / norm(x)) < 0 ? P[k][t] / -(K[k][t]' * (x / norm(x))) : Inf for t in eachindex(P[k])]) for x in pts]
            t = argmax([count(==(h), hrs) for h in 1:length(P[k])])
            n_pred = -K[k][t]; c_pred = P[k][t]
        end
        d_pred = c_pred / norm(n_pred)
        ang = acosd(clamp(dot(a, n_pred) / (norm(a) * norm(n_pred)), -1, 1))
        @printf("  face %-6s %2d points: planar to %.1f%% (RMS / distance), normal %.1f deg from the predicted, distance %.3f vs %.3f predicted (x%.2f)\n",
                k, length(pts), 100rms / d, ang, d, d_pred, d / d_pred)
    end
end
