# The three-substation window of ieee123 (subs1 reference, subs3 and subs4 moving), T = 24, with and
# without batteries: measured boundary points against the region predicted from the angle-to-import
# sensitivities (stage "simplex"). Two panels: PV at unity power factor, where only batteries are
# dispatched and the prediction should be exact, and PV reactive power free, the paper's case.
# Also prints each predicted triangle's corners and face directions, to compare the two cases.
#   julia --project=problems/msopf problems/msopf/root_level/plot_simplex.jl
# Needs simplex_check_unity.csv (SIMPLEX_PV=unity ... simplex) and three_region.csv (stage three).
empty!(ARGS); push!(ARGS, "none")
include(joinpath(@__DIR__, "multiperiod_study.jl"))

REF, A, B = "subs1", "subs3", "subs4"
net, pvs = load_sources([REF, A, B])
rows(f) = [split(l, ',') for l in readlines(joinpath(OUT_DIR, f))[2:end]]

"Measured (theta in deg => radius) per case for PV mode `pv`."
function measured(pv)
    m = Dict("nobatt" => Dict{Float64,Float64}(), string(DAY) => Dict{Float64,Float64}())
    if pv == "qfree"
        for r in rows("three_region.csv")
            haskey(m, r[2]) && (m[r[2]][parse(Float64, r[1])] = parse(Float64, r[4]))
        end
    else
        for r in rows("simplex_check_$pv.csv")
            r[1] == "three" && haskey(m, r[4]) && (m[r[4]][parse(Float64, r[3])] = parse(Float64, r[5]))
        end
    end
    return m
end

"Predicted boundary on a fine angle grid, and the corners: where the limiting substation changes."
function predicted(sens; day)
    ths = range(0, 2pi; length = 1441)[1:end-1]
    pts = [predicted_reach(sens, Dict(A => cos(th), B => sin(th)); day) for th in ths]
    r = first.(pts); lim = last.(pts)
    corners = [(r[i] * cos(ths[i]), r[i] * sin(ths[i]), lim[mod1(i - 1, end)], lim[i]) for i in eachindex(ths) if lim[i] != lim[mod1(i - 1, end)]]
    return (; x = r .* cos.(ths), y = r .* sin.(ths), lim, corners)
end

plot_style!()
panels = []
# Figure data for the paper's own renderer: pv, case, kind (predicted | measured), x, y in degrees.
# face: the substation whose import reaches zero there (for measured points, the one predicted along that direction).
figdata = ["pv,case,kind,d3_minus_d1_deg,d4_minus_d1_deg,face"]
for (pv, title) in (("qfree", "PV reactive power dispatched"),)      # the only PV configuration in the paper (user, 2026-10-20)
    s = setup(net, pvs; T = 24, pv)
    sens = sensitivity(s, [A, B])
    meas = measured(pv)
    p = plot(; xlabel = "δ₃ − δ₁ (deg)", ylabel = "δ₄ − δ₁ (deg)", aspect_ratio = :equal, title, legend = :outerbottom, legend_columns = 2,
             xlims = (-0.45, 0.45), ylims = (-0.45, 0.55))
    say("\n", title)
    for (case, col, lab) in ((DAY, colorant"#7a4fd0", "cyclic batteries"), (:nobatt, colorant"#6b6b6b", "batteries idle"))
        pr = predicted(sens; day = case === DAY)
        plot!(p, [pr.x; pr.x[1]], [pr.y; pr.y[1]]; seriestype = :shape, fillalpha = 0.12, color = col, linecolor = col, lw = 1.4,
              label = "predicted, " * lab)
        mm = meas[string(case)]
        ths = sort(collect(keys(mm)))
        scatter!(p, [mm[t] * cosd(t) for t in ths], [mm[t] * sind(t) for t in ths]; color = col, ms = 3.5, markerstrokewidth = 0,
                 label = "measured, " * lab)
        for i in 1:4:length(pr.x)
            push!(figdata, @sprintf("%s,%s,predicted,%.6f,%.6f,%s", pv, case, pr.x[i], pr.y[i], pr.lim[i]))
        end
        for t in ths
            lt = last(predicted_reach(sens, Dict(A => cosd(t), B => sind(t)); day = case === DAY))
            push!(figdata, @sprintf("%s,%s,measured,%.6f,%.6f,%s", pv, case, mm[t] * cosd(t), mm[t] * sind(t), lt))
        end
        say("  ", rpad(string(case), 7), " corners: ", join([@sprintf("(%+.3f, %+.3f) %s|%s", c[1], c[2], c[3], c[4]) for c in pr.corners], "  "))
    end
    # Face directions: the cyclic faces use the daily-mean sensitivities, the idle ones the binding hour.
    for k in sens.names
        Kd = [sum(sens.K[(k, A)]), sum(sens.K[(k, B)])] ./ sens.T
        say(@sprintf("  face %s: daily-mean normal (%+.1f, %+.1f) kW/deg, direction %.1f deg; mean import %.1f kW, worst hour %.1f kW",
                     k, Kd[1], Kd[2], atand(Kd[2], Kd[1]), sum(sens.P[k]) / sens.T, minimum(sens.P[k])))
    end
    scatter!(p, [0], [0]; color = :black, ms = 2.5, label = "equal angles")
    push!(panels, p)
end
fig = plot(panels...; layout = (1, length(panels)), size = (560 * length(panels), 640), margin = 6Plots.mm)
for ext in ("png", "pdf")
    savefig(fig, joinpath(OUT_DIR, "fig_simplex.$ext"))
end
open(io -> foreach(l -> println(io, l), figdata), joinpath(OUT_DIR, "simplex_figure_data.csv"), "w")
say("\nwrote fig_simplex.png, fig_simplex.pdf and simplex_figure_data.csv")
