# The three-substation window of the large decks (stages three_ieee2522, three_large10kd, three_ieee9500,
# three_ieee9500xf), T = 24, PV reactive power dispatched: measured boundary points (three_large_windows.csv)
# against the region predicted from the hourly angle-to-import sensitivities at equal angles
# (three_sens_<label>_T24.csv), batteries idle and cyclic. Prediction as in multiperiod_study.jl's
# predicted_reach: along a direction, the angle at which the first substation import reaches zero, in its
# worst hour (idle) or in its daily mean (cyclic).
#   julia --project=envs/msopf envs/multi_poi/root_level/plot_three_large.jl [label ...]
# Writes fig_three_<label>.png/.pdf and three_figure_data_<label>.csv (pv, case, kind, x, y in degrees).
using Plots, Printf
const DIR = joinpath(@__DIR__, "..", "results", "multiperiod")
const LIMIT = 20.0

rows(f) = [split(l, ',') for l in readlines(joinpath(DIR, f))[2:end]]

"Hourly imports P[k] and sensitivities K[k] = (to subs2, to subs3) per substation."
function sensitivities(label)
    P = Dict{String,Vector{Float64}}(); K = Dict{String,Vector{NTuple{2,Float64}}}()
    for r in rows("three_sens_$(label)_T24.csv")
        k = String(r[1]); push!(get!(P, k, Float64[]), parse(Float64, r[3]))
        push!(get!(K, k, NTuple{2,Float64}[]), (parse(Float64, r[4]), parse(Float64, r[5])))
    end
    return P, K
end

function reach(P, K, th; day)
    r, lim = Inf, ""
    for k in keys(P)
        slope = [a * cos(th) + b * sin(th) for (a, b) in K[k]]
        x = day ? (sum(slope) < 0 ? sum(P[k]) / -sum(slope) : Inf) :
                  minimum((slope[t] < 0 ? P[k][t] / -slope[t] : Inf) for t in eachindex(slope))
        x < r && ((r, lim) = (x, k))
    end
    return r, lim
end

function panel(label)
    P, K = sensitivities(label)
    meas = [r for r in rows("three_large_windows.csv") if r[1] == label]
    data = ["pv,case,kind,d2_minus_d1_deg,d3_minus_d1_deg"]
    rmax = maximum(parse(Float64, r[8]) for r in meas if r[11] != "true")   # zoom to the measured edges
    lims = (-1.25rmax, 1.25rmax)
    p = plot(; aspect_ratio = :equal, xlims = lims, ylims = lims, title = label, legend = :outerbottom, legend_columns = 2,
             xlabel = "δ₂ − δ₁ (deg)", ylabel = "δ₃ − δ₁ (deg)")
    ths = range(0, 2pi; length = 1441)[1:end-1]
    for (case, col, lab) in (("cyclic", colorant"#7a4fd0", "cyclic batteries"), ("nobatt", colorant"#6b6b6b", "batteries idle"))
        r = [first(reach(P, K, th; day = case == "cyclic")) for th in ths]
        r = min.(r, 3rmax)                                    # an open side (no substation limits it) drawn at the frame
        x, y = r .* cos.(ths), r .* sin.(ths)
        plot!(p, [x; x[1]], [y; y[1]]; seriestype = :shape, fillalpha = 0.12, color = col, linecolor = col, lw = 1.4, label = "predicted, " * lab)
        sel = [m for m in meas if m[5] == case]
        th = [deg2rad(parse(Float64, m[6])) for m in sel]; rm = [parse(Float64, m[8]) for m in sel]
        capped = [m[11] == "true" for m in sel]
        scatter!(p, (rm .* cos.(th))[.!capped], (rm .* sin.(th))[.!capped]; color = col, ms = 3.5, markerstrokewidth = 0, label = "measured, " * lab)
        rc = min.(rm, 1.15rmax)                               # an open direction is marked at the frame
        any(capped) && scatter!(p, (rc .* cos.(th))[capped], (rc .* sin.(th))[capped]; color = :white, markerstrokecolor = col, ms = 4,
                                label = case == "cyclic" ? "no edge within $(Int(LIMIT)) deg" : "")
        for i in 1:4:length(ths)
            push!(data, @sprintf("qfree,%s,predicted,%.6f,%.6f", case, x[i], y[i]))
        end
        for (t, v, c) in zip(th, rm, capped)
            push!(data, @sprintf("qfree,%s,%s,%.6f,%.6f", case, c ? "measured_capped" : "measured", v * cos(t), v * sin(t)))
        end
        # how far the measured points lie from the predicted boundary, by direction
        dev = [(v - first(reach(P, K, t; day = case == "cyclic"))) / first(reach(P, K, t; day = case == "cyclic")) for (t, v, c) in zip(th, rm, capped) if !c]
        @printf("%-16s %-6s measured vs predicted: median %+.1f%%, range %+.1f%% to %+.1f%% over %d directions\n",
                label, case, 100 * sort(dev)[(end + 1) ÷ 2], 100 * minimum(dev), 100 * maximum(dev), length(dev))
    end
    scatter!(p, [0], [0]; color = :black, ms = 2.5, label = "equal angles")
    open(io -> foreach(l -> println(io, l), data), joinpath(DIR, "three_figure_data_$label.csv"), "w")
    return p
end

labels = isempty(ARGS) ? sort(unique(r[1] for r in rows("three_large_windows.csv"))) : ARGS
for label in labels
    isfile(joinpath(DIR, "three_sens_$(label)_T24.csv")) || (println("skip $label: no three_sens_$(label)_T24.csv"); continue)
    p = panel(label)
    plot!(p; size = (620, 700), margin = 5Plots.mm)
    foreach(ext -> savefig(p, joinpath(DIR, "fig_three_$label.$ext")), ("png", "pdf"))
    println("wrote fig_three_$label.png/.pdf and three_figure_data_$label.csv")
end
