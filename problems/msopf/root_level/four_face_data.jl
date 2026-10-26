# Figure data for the four-substation faces seen straight on (paper figure mp_faces): for each system, each
# face of the predicted tetrahedron with batteries dispatched (daily-mean sensitivities, n_k . x = c_k with
# n_k = -sum_t K_k(t), c_k = sum_t P_k(t)), its three corners and the measured boundary points predicted on
# it, in the face's own plane coordinates (degrees, origin at the foot of the plane), with what bound each
# point and its measured / predicted reach. Also prints the rows of the paper's cross-system table.
# Same basis as the viewer's Face-on view (four_regions.html).
#   julia --project=problems/msopf problems/msopf/root_level/four_face_data.jl
using Printf, LinearAlgebra, Statistics
const RES = joinpath(@__DIR__, "..", "results")
rows(f) = [split(l, ',') for l in readlines(f)[2:end]]

# (stage, system, sources, output label)
const CASES = [("four_ieee123a", "ieee123_5poi_1ph", "subs1:subs2:subs3:subs4", "ieee123"),
               ("four_ieee9500", "ieee9500_1ph", "subs1:subs2:subs3:subs4", "ieee2705"),
               ("four_ieee2522", "ieee2522C_1ph", "source:subs2:subs3:subs4", "ieee2522"),
               # decks with service transformers, band at customer buses only (2026-10-25 campaign, regions_oct25)
               ("four_ieee9500xf", "ieee9500xf_1ph", "subs1:subs2:subs3:subs4", "ieee2705xf"),
               ("four_smartds_p1u_typeiii", "smartds_p1u_xf_typeIII", "subs1:subs2:subs3:subs4", "sfp1u")]

function sens(stage)
    f = first(p for (r, _, fs) in walkdir(joinpath(RES, "kamiak")) for x in fs for p in (joinpath(r, x),) if endswith(p, "multi_sens_$(stage)_T24.csv"))
    names = String[]; P = Dict{String,Vector{Float64}}(); K = Dict{String,Vector{Vector{Float64}}}()
    for r in rows(f)
        k = String(r[1]); k in names || push!(names, k)
        push!(get!(P, k, Float64[]), parse(Float64, r[3])); push!(get!(K, k, Vector{Float64}[]), parse.(Float64, r[4:end]))
    end
    return names, P, K
end

basis(n) = (u = normalize(n); a = abs(u[3]) < 0.9 ? [0.0, 0, 1] : [1.0, 0, 0]; e1 = normalize(cross(u, a)); (u, e1, cross(u, e1)))

data = rows(joinpath(RES, "multiperiod", "multi_windows_fib78.csv"))
for (stage, sys, srcs, lab) in CASES
    names, P, K = sens(stage)
    N = [-sum(K[k]) for k in names]; C = [sum(P[k]) for k in names]
    meas = [r for r in data if r[1] == sys && r[2] == srcs && r[12] == "false"]
    out = ["kind,face,case,u_deg,v_deg,bound,ratio"]
    for k in 1:4
        n, e1, e2 = basis(N[k]); x0 = N[k] * C[k] / dot(N[k], N[k])
        o = setdiff(1:4, k)
        for (i, j) in ((o[1], o[2]), (o[1], o[3]), (o[2], o[3]))
            v = [N[k] N[i] N[j]]' \ [C[k], C[i], C[j]]
            push!(out, @sprintf("vertex,%d,,%.5f,%.5f,,", k - 1, dot(v - x0, e1), dot(v - x0, e2)))
        end
        for r in meas
            findfirst(==(r[11]), names) == k || continue
            x = parse(Float64, r[9]) .* parse.(Float64, r[6:8])
            b = String(r[19]); bi = b in names ? findfirst(==(b), names) - 1 : (b == "voltage" ? -1 : -2)
            push!(out, @sprintf("point,%d,%s,%.5f,%.5f,%d,%.5f", k - 1, r[4], dot(x - x0, e1), dot(x - x0, e2), bi,
                                parse(Float64, r[9]) / parse(Float64, r[10])))
        end
    end
    f = joinpath(RES, "multiperiod", "four_faces_$(lab).csv")
    open(io -> foreach(l -> println(io, l), out), f, "w")
    # table row: own-angle daily-mean sensitivity of the moving substations, reach range, measured vs predicted
    own = [mean(K[names[i]][t][i-1] for t in eachindex(K[names[i]])) for i in 2:4]
    rr(case) = [parse(Float64, r[9]) for r in meas if r[4] == case]
    dd = [100 * (parse(Float64, r[9]) / parse(Float64, r[10]) - 1) for r in meas]
    @printf("%-9s K %.0f-%.0f  idle %.2f-%.2f  dispatched %.2f-%.2f  delta %+.0f to %+.0f  (%d points; wrote %s)\n", lab,
            minimum(own), maximum(own), extrema(rr("nobatt"))..., extrema(rr("cyclic"))..., extrema(dd)..., length(meas), basename(f))
end
