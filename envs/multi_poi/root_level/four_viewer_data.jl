# Data for the interactive four-substation viewer: for each four_* stage and battery case, the predicted
# region's surface on a sphere of directions (reach from the hourly angle-to-import sensitivities, as
# predicted_reach does, with the limiting substation per point) and the measured boundary points
# (multi_windows.csv). Writes results/multiperiod/four_viewer_data.json.
#   julia --project=envs/msopf envs/multi_poi/root_level/four_viewer_data.jl
using Printf
const RES = joinpath(@__DIR__, "..", "results")
rows(f) = [split(l, ',') for l in readlines(f)[2:end]]

# (stage, system, sources, tab title)
const STAGES = [("four_ieee123a", "ieee123_5poi_1ph", "subs1:subs2:subs3:subs4", "IEEE 123-node, sources 1, 2, 3, 4"),
                ("four_ieee123b", "ieee123_5poi_1ph", "subs1:subs3:subs4:subs5", "IEEE 123-node, sources 1, 3, 4, 5"),
                ("four_ieee2522", "ieee2522C_1ph", "source:subs2:subs3:subs4", "2522-bus feeder, sources at 947, 996, 1120"),
                ("four_ieee9500", "ieee9500_1ph", "subs1:subs2:subs3:subs4", "2705-bus feeder (from IEEE 9500), its 3 + 1 added"),
                # the same with every source at 1.03 pu (rows in multi_windows_fib78_v1.03.csv)
                ("four_ieee9500_v1.03", "ieee9500_1ph", "subs1:subs2:subs3:subs4", "2705-bus feeder, sources at 1.03 pu")]

function sens(stage)
    f = first(p for (r, _, fs) in walkdir(joinpath(RES, "kamiak")) for x in fs for p in (joinpath(r, x),) if endswith(p, "multi_sens_$(stage)_T24.csv"))
    names = String[]; P = Dict{String,Vector{Float64}}(); K = Dict{String,Vector{Vector{Float64}}}()
    for r in rows(f)
        k = String(r[1]); k in names || push!(names, k)
        push!(get!(P, k, Float64[]), parse(Float64, r[3])); push!(get!(K, k, Vector{Float64}[]), parse.(Float64, r[4:end]))
    end
    return names, P, K
end

function reach(names, P, K, u; day)
    r, lim = Inf, 0
    for (i, k) in enumerate(names)
        s = [sum(K[k][t] .* u) for t in eachindex(P[k])]
        x = day ? (sum(s) < 0 ? sum(P[k]) / -sum(s) : Inf) : minimum((s[t] < 0 ? P[k][t] / -s[t] : Inf) for t in eachindex(s))
        x < r && ((r, lim) = (x, i - 1))
    end
    return r, lim
end

fmt(x) = @sprintf("%.5f", x)
arr(v) = "[" * join(v, ",") * "]"
const DSET = isempty(ARGS) ? "cube26" : ARGS[1]           # cube26 (multi_windows.csv) or fib78 (multi_windows_fib78.csv)
data = rows(joinpath(RES, "multiperiod", DSET == "cube26" ? "multi_windows.csv" : "multi_windows_$(DSET).csv"))
nth, nph = 49, 97
out = String[]
for (stage, sys, srcs, title) in STAGES
    names, P, K = sens(stage)
    tag = occursin("_v", stage) ? stage[findfirst("_v", stage)[1]:end] : ""   # e.g. _v1.03: its own rows
    src = isempty(tag) ? data : rows(joinpath(RES, "multiperiod", (DSET == "cube26" ? "multi_windows" : "multi_windows_$(DSET)") * tag * ".csv"))
    meas = [r for r in src if r[1] == sys && r[2] == srcs]
    isempty(meas) && (meas = [r for r in rows(joinpath(RES, "multiperiod", "multi_windows.csv")) if r[1] == sys && r[2] == srcs])   # not run in this set yet: the 26 directions
    rmax = maximum(parse(Float64, r[9]) for r in meas if r[12] == "false")
    cases = String[]
    for case in ("nobatt", "cyclic")
        X, Y, Z, F = String[], String[], String[], String[]
        for i in 1:nth
            th = pi * (i - 1) / (nth - 1); xs, ys, zs, fs = String[], String[], String[], String[]
            for j in 1:nph
                ph = 2pi * (j - 1) / (nph - 1)
                u = [sin(th) * cos(ph), sin(th) * sin(ph), cos(th)]
                r, lim = reach(names, P, K, u; day = case == "cyclic")
                r = min(r, 2.5rmax)
                push!(xs, fmt(r * u[1])); push!(ys, fmt(r * u[2])); push!(zs, fmt(r * u[3])); push!(fs, string(lim))
            end
            push!(X, arr(xs)); push!(Y, arr(ys)); push!(Z, arr(zs)); push!(F, arr(fs))
        end
        pts = String[]
        for r in meas
            r[4] == case || continue
            u = parse.(Float64, r[6:8]); v = parse(Float64, r[9]); pv = parse(Float64, r[10])
            lim = findfirst(==(r[11]), names) - 1
            # what bound the measured edge: the index of the substation, -1 a voltage limit, -2 unknown or open
            bl = length(r) >= 19 ? String(r[19]) : "unknown"
            bi = bl in names ? findfirst(==(bl), names) - 1 : (bl == "voltage" ? -1 : -2)
            push!(pts, "[" * join(fmt.(v .* u), ",") * ",$(lim),$(fmt(v)),$(fmt(pv)),$(r[12] == "true"),$(bi)]")
        end
        push!(cases, "\"$case\":{\"surface\":{\"x\":$(arr(X)),\"y\":$(arr(Y)),\"z\":$(arr(Z)),\"face\":$(arr(F))},\"points\":$(arr(pts))}")
    end
    lin = isempty(tag) ? [parse(Float64, r[8]) for r in rows(joinpath(RES, "multiperiod", "multi_linearity.csv")) if r[1] == sys && r[2] == srcs] : [NaN]
    # daily-mean face normals n_k = -sum_t K_k(t), for the edge-on views
    normals = [arr(fmt.(-sum(K[k]))) for k in names]; offsets = [fmt(sum(P[k])) for k in names]   # cyclic face: n.x = offset
    push!(out, "{\"id\":\"$stage\",\"title\":\"$title\",\"sources\":$(arr(["\"$n\"" for n in names])),\"normals\":$(arr(normals)),\"offsets\":$(arr(offsets)),\"linearity_max\":$(isnan(maximum(lin)) ? "null" : fmt(100maximum(lin))),\"cases\":{$(join(cases, ","))}}")
end
jf = joinpath(RES, "multiperiod", DSET == "cube26" ? "four_viewer_data.json" : "four_viewer_data_$(DSET).json")
open(jf, "w") do io
    print(io, "[", join(out, ","), "]")
end
println("wrote ", basename(jf), ", ", round(filesize(jf) / 1e6; digits = 2), " MB")
