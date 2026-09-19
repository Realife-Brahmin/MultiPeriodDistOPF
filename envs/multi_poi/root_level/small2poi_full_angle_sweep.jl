# small2poi_full_angle_sweep.jl -- how the two substation powers of small2poi_1ph move
# with the source angles delta1, delta2. Full-angle AC model in JuMP (Ipopt, with
# Gurobi as a global check); every point is cross-checked against OpenDSS.
#
# Run from the repo root:
#     julia --project=envs/tadmm envs/multi_poi/root_level/small2poi_full_angle_sweep.jl
#
# Writes (gitignored) envs/multi_poi/processedData/small2poi_1ph/full_angle_sweep/:
#     sweep.csv            every sweep point: JuMP-Ipopt, JuMP-Gurobi, OpenDSS, differences
#     psubs_vs_angle.png   P_Subs and Q_Subs against delta1 - delta2
#
# Two source models, both built from the same deck:
#   deck  -- as written: each Vsource keeps OpenDSS's default series impedance, so delta
#            is the angle of the EMF behind it and the substation bus angle floats;
#   stiff -- Vsource impedance replaced in memory by a 1e-8 ohm reactance (OpenDSS's
#            stand-in for an ideal source, as ieee2522C/large10kC use), so delta is the
#            substation bus angle itself -- what the branch-flow MPOPF scripts assume.
#            The sweep also solves a truly ideal source (Zs = 0) and reports the gap.

include(joinpath(@__DIR__, "..", "full_angle_pf.jl"))
using Plots

const SYSTEM = "small2poi_1ph"
const OUT_DIR = normpath(joinpath(@__DIR__, "..", "processedData", SYSTEM, "full_angle_sweep"))
const MODELS = [("deck", false), ("stiff", true)]    # (name, stiff_sources)
const DDELTA = -4.0:0.25:4.0                         # delta1 - delta2 in degrees (delta2 = 0)
const SAME_DDELTA = [(1.0, 0.0), (0.5, -0.5), (15.0, 14.0), (-44.0, -45.0), (121.0, 120.0)]
const TOL_KW = 1e-3                                  # agreement demanded: 1 W, 1 var
const TOL_V = 1e-8                                   # agreement demanded: pu voltage phasor

# Gurobi is the global cross-check; without a valid license the sweep still runs on Ipopt
# (verified against OpenDSS) and every Gurobi column reads SKIPPED / NaN.
const USE_GUROBI, GUROBI_WHY = gurobi_usable()
USE_GUROBI || println("\nGUROBI SKIPPED -- ", GUROBI_WHY, "\nRunning Ipopt + OpenDSS only.")

"Solve one (delta1, delta2) point three ways (four with `with_ideal`: Ipopt with Zs = 0)."
function solve_all(net, d1, d2; with_ideal = false)
    delta = Dict("grid1" => d1, "grid2" => d2)
    return (; d1, d2, dd = d1 - d2,
            ip = solve_full_angle(net, delta; solver = :ipopt),
            gu = USE_GUROBI ? solve_full_angle(net, delta; solver = :gurobi) : skipped(net, :gurobi),
            ds = opendss_point(net, delta),
            id = with_ideal ? solve_full_angle(net, delta; ideal_sources = true) : skipped(net, :ipopt))
end

function print_row(label, status, x, load_bus)
    V = x.V_pu[load_bus]
    @printf "  %-8s %-17s %11.3f %11.3f %11.3f %11.3f %9.6f %9.4f\n" label status x.P_subs_kW["grid1"] x.P_subs_kW["grid2"] x.Q_subs_kvar["grid1"] x.Q_subs_kvar["grid2"] abs(V) rad2deg(angle(V))
end

header() = @printf "  %-8s %-17s %11s %11s %11s %11s %9s %9s\n" "" "status" "P_Subs1 kW" "P_Subs2 kW" "Q_Subs1 kvar" "Q_Subs2 kvar" "|V_L| pu" "∠V_L deg"

rows = []
nets = Dict{String,Any}()
for (model, stiff) in MODELS
    compile_deck(SYSTEM; stiff_sources = stiff)
    net = nets[model] = read_network()
    Set(s.name for s in net.sources) == Set(["grid1", "grid2"]) || error("expected Vsources grid1, grid2")
    load_bus = only(unique(d.bus for d in net.loads))
    P_load = sum(d.kW for d in net.loads)

    println("\n", "="^104)
    @printf "%s -- %s sources\n" SYSTEM model
    for s in net.sources
        z = s.z_series
        @printf "  Vsource.%-6s bus %-3s  %.3f pu on %.4g kV (L-N)   series Z = %.6f + j%.2e ohm%s\n" s.name s.bus s.pu net.kV_base real(z) imag(z) (stiff ? "  (stand-in, set in memory)" : "")
    end
    for l in net.lines
        z = -1 / l.Y[1, 2]
        @printf "  Line.%-9s %-3s -> %-5s  Z = %.6f + j%.6f ohm, shunt j%.3e S per end\n" l.name l.from l.to real(z) imag(z) imag(l.Y[1, 1] + l.Y[1, 2])
    end
    for d in net.loads
        @printf "  Load.%-9s bus %-5s  %.1f kW + j%.1f kvar, constant PQ for %.2f-%.2f pu of %.4g kV\n" d.name d.bus d.kW d.kvar d.vminpu d.vmaxpu d.kV
    end
    println("="^104)

    # ---- base case: delta1 = delta2 = 0 ----------------------------------------------
    b = solve_all(net, 0.0, 0.0)
    println("\nBase case, delta1 = delta2 = 0:")
    header()
    print_row("Ipopt", string(b.ip.status), b.ip, load_bus)
    USE_GUROBI ? print_row("Gurobi", string(b.gu.status), b.gu, load_bus) :
                 @printf("  %-8s %s\n", "Gurobi", "SKIPPED (see message at top)")
    print_row("OpenDSS", b.ds.converged ? "converged ($(b.ds.iterations) it)" : "NOT CONVERGED", b.ds, load_bus)
    di = max_diff(b.ip, b.ds)
    @printf "  |Ipopt - OpenDSS|:  P %.2e kW  Q %.2e kvar  V %.2e pu\n" di.dP di.dQ di.dV
    if USE_GUROBI
        dg = max_diff(b.gu, b.ds)
        @printf "  |Gurobi - OpenDSS|: P %.2e kW  Q %.2e kvar  V %.2e pu\n" dg.dP dg.dQ dg.dV
    end

    # The other power-flow root of the same equations: Ipopt started low (and, when
    # licensed, Gurobi maximizing losses, which finds it globally)
    zero_delta = Dict("grid1" => 0.0, "grid2" => 0.0)
    other = [("Ipopt, load bus started at 0.02 pu", solve_full_angle(net, zero_delta; V_load_start = 0.02))]
    USE_GUROBI && push!(other, ("Gurobi, max losses", solve_full_angle(net, zero_delta; solver = :gurobi,
                                                                       sense = MAX_SENSE)))
    for (how, x) in other
        @printf "  Other power-flow root (%s): %s, |V_L| = %.4f pu, P_Subs1 + P_Subs2 = %.0f kW\n" how string(x.status) abs(x.V_pu[load_bus]) sum(values(x.P_subs_kW))
    end

    # ---- only delta1 - delta2 should matter ---------------------------------------------
    println("\nSame delta1 - delta2 = 1 deg, different common rotation:")
    @printf "  %9s %9s %14s %14s %14s %14s\n" "delta1" "delta2" "P_Subs1 Ipopt" "P_Subs2 Ipopt" "P_Subs1 DSS" "P_Subs2 DSS"
    same = [solve_all(net, d1, d2) for (d1, d2) in SAME_DDELTA]
    for r in same
        @printf "  %9.2f %9.2f %14.4f %14.4f %14.4f %14.4f\n" r.d1 r.d2 r.ip.P_subs_kW["grid1"] r.ip.P_subs_kW["grid2"] r.ds.P_subs_kW["grid1"] r.ds.P_subs_kW["grid2"]
    end
    spread(f) = maximum(f(r) for r in same) - minimum(f(r) for r in same)
    @printf "  spread across rotations: Ipopt %.2e kW, OpenDSS %.2e kW\n" max(spread(r -> r.ip.P_subs_kW["grid1"]), spread(r -> r.ip.P_subs_kW["grid2"])) max(spread(r -> r.ds.P_subs_kW["grid1"]), spread(r -> r.ds.P_subs_kW["grid2"]))

    # ---- the sweep ------------------------------------------------------------------------
    for dd in DDELTA
        push!(rows, (; model, load_bus, P_load, solve_all(net, dd, 0.0; with_ideal = stiff)...))
    end
end

# ---- summary ----------------------------------------------------------------------------------
ok_all = true
for (model, _) in MODELS
    r = filter(x -> x.model == model, rows)
    println("\n", "="^104)
    @printf "Sweep, %s sources: delta2 = 0, delta1 = %.2f .. %.2f deg (%d points; table shows whole degrees)\n" model first(DDELTA) last(DDELTA) length(r)
    @printf "  %9s %11s %11s %12s %12s %9s %11s\n" "d1-d2 deg" "P_Subs1 kW" "P_Subs2 kW" "Q_Subs1 kvar" "Q_Subs2 kvar" "|V_L| pu" "losses kW"
    for x in r
        isinteger(x.dd) || continue
        P1, P2 = x.ip.P_subs_kW["grid1"], x.ip.P_subs_kW["grid2"]
        @printf "  %9.2f %11.1f %11.1f %12.1f %12.1f %9.5f %11.1f\n" x.dd P1 P2 x.ip.Q_subs_kvar["grid1"] x.ip.Q_subs_kvar["grid2"] abs(x.ip.V_pu[x.load_bus]) P1 + P2 - x.P_load
    end

    # Sensitivity at delta1 = delta2 and the angles where one substation stops importing
    i0 = findfirst(x -> x.dd == 0, r)
    slope = (r[i0+1].ip.P_subs_kW["grid1"] - r[i0-1].ip.P_subs_kW["grid1"]) / (r[i0+1].dd - r[i0-1].dd)
    @printf "  dP_Subs1/d(delta1 - delta2) at 0: %.0f kW/deg\n" slope
    for k in ("grid2", "grid1")
        y = [x.ip.P_subs_kW[k] for x in r]
        j = findfirst(i -> sign(y[i]) != sign(y[i+1]), 1:length(y)-1)
        if j === nothing
            @printf "  P_Subs%s does not cross zero in this range\n" last(k)
        else
            z = r[j].dd + (r[j+1].dd - r[j].dd) * y[j] / (y[j] - y[j+1])
            @printf "  P_Subs%s crosses zero at delta1 - delta2 = %+.3f deg (linear interpolation)\n" last(k) z
        end
    end

    n = length(r)
    e_ip = [max_diff(x.ip, x.ds) for x in r]
    e_gu = [max_diff(x.gu, x.ds) for x in r]
    worst(es, f) = maximum(f, es)
    n_ip, n_gu = count(x -> x.ip.ok, r), count(x -> x.gu.ok, r)
    n_conv, n_band = count(x -> x.ds.converged, r), count(x -> x.ds.pq_band, r)
    @printf "  solved: Ipopt %d/%d, OpenDSS converged %d/%d, load inside its PQ band %d/%d, Gurobi (global) %s\n" n_ip n n_conv n n_band n (USE_GUROBI ? "$n_gu/$n" : "SKIPPED")
    @printf "  worst |Ipopt - OpenDSS|:  P %.2e kW  Q %.2e kvar  V %.2e pu\n" worst(e_ip, e -> e.dP) worst(e_ip, e -> e.dQ) worst(e_ip, e -> e.dV)
    USE_GUROBI && @printf("  worst |Gurobi - OpenDSS|: P %.2e kW  Q %.2e kvar  V %.2e pu\n",
                          worst(e_gu, e -> e.dP), worst(e_gu, e -> e.dQ), worst(e_gu, e -> e.dV))
    checked = USE_GUROBI ? vcat(e_ip, e_gu) : e_ip
    passed = n_ip == n_conv == n_band == n && (!USE_GUROBI || n_gu == n) &&
             all(e -> e.dP <= TOL_KW && e.dQ <= TOL_KW && e.dV <= TOL_V, checked)
    global ok_all &= passed
    @printf "  => %s (tolerance %.0e kW/kvar, %.0e pu)\n" (passed ? "AGREE" : "DISAGREE") TOL_KW TOL_V
    ran = filter(x -> x.id.status != :SKIPPED, r)
    if !isempty(ran)
        good = filter(x -> x.id.ok, ran)
        d_ideal = isempty(good) ? NaN : maximum(max_diff(x.id, x.ip).dP for x in good)
        @printf "  ideal source (Zs = 0) vs the 1e-8 ohm stand-in: max |ΔP_Subs| %.2e kW (ideal solved %d/%d)\n" d_ideal length(good) length(ran)
    end
    d_vs = maximum(abs(x.ds.S_vsource_kVA[k] - complex(x.ds.P_subs_kW[k], x.ds.Q_subs_kvar[k]))
                   for x in r for k in ("grid1", "grid2"))
    @printf "  OpenDSS's own Vsource power reading vs line-flow reading: max %.2e kVA\n" d_vs
end

# ---- CSV --------------------------------------------------------------------------------------
mkpath(OUT_DIR)
csv = joinpath(OUT_DIR, "sweep.csv")
open(csv, "w") do io
    println(io, join(["model", "delta1_deg", "delta2_deg", "ddelta_deg",
        "P_subs1_kW", "P_subs2_kW", "Q_subs1_kvar", "Q_subs2_kvar", "Vload_pu", "Vload_deg", "network_losses_kW",
        "gurobi_P_subs1_kW", "gurobi_P_subs2_kW", "gurobi_Q_subs1_kvar", "gurobi_Q_subs2_kvar",
        "dss_P_subs1_kW", "dss_P_subs2_kW", "dss_Q_subs1_kvar", "dss_Q_subs2_kvar", "dss_Vload_pu", "dss_Vload_deg",
        "dss_converged", "dss_iterations", "dss_load_in_pq_band",
        "dss_vsource_P_subs1_kW", "dss_vsource_P_subs2_kW",
        "err_ipopt_dss_P_kW", "err_ipopt_dss_Q_kvar", "err_ipopt_dss_V_pu",
        "err_gurobi_dss_P_kW", "err_gurobi_dss_Q_kvar", "err_gurobi_dss_V_pu",
        "ideal_P_subs1_kW", "ideal_P_subs2_kW",
        "ipopt_status", "gurobi_status", "ipopt_time_s", "gurobi_time_s"], ","))
    for x in rows
        ei, eg = max_diff(x.ip, x.ds), max_diff(x.gu, x.ds)
        Vi, Vd = x.ip.V_pu[x.load_bus], x.ds.V_pu[x.load_bus]
        vals = Any[x.model, x.d1, x.d2, x.dd,
            x.ip.P_subs_kW["grid1"], x.ip.P_subs_kW["grid2"], x.ip.Q_subs_kvar["grid1"], x.ip.Q_subs_kvar["grid2"],
            abs(Vi), rad2deg(angle(Vi)), x.ip.P_subs_kW["grid1"] + x.ip.P_subs_kW["grid2"] - x.P_load,
            x.gu.P_subs_kW["grid1"], x.gu.P_subs_kW["grid2"], x.gu.Q_subs_kvar["grid1"], x.gu.Q_subs_kvar["grid2"],
            x.ds.P_subs_kW["grid1"], x.ds.P_subs_kW["grid2"], x.ds.Q_subs_kvar["grid1"], x.ds.Q_subs_kvar["grid2"],
            abs(Vd), rad2deg(angle(Vd)), x.ds.converged, x.ds.iterations, x.ds.pq_band,
            real(x.ds.S_vsource_kVA["grid1"]), real(x.ds.S_vsource_kVA["grid2"]),
            ei.dP, ei.dQ, ei.dV, eg.dP, eg.dQ, eg.dV,
            x.id.P_subs_kW["grid1"], x.id.P_subs_kW["grid2"],
            x.ip.status, x.gu.status, x.ip.time, x.gu.time]
        println(io, join((v isa AbstractFloat ? @sprintf("%.12g", v) : string(v) for v in vals), ","))
    end
end

# ---- plot ---------------------------------------------------------------------------------------
# Reference palette, light surface: categorical slots 1 and 2 (validated), chrome in ink tokens.
const SURFACE = colorant"#fcfcfb"
const INK2, MUTED = colorant"#52514e", colorant"#898781"
const GRIDC, AXISC = colorant"#e1e0d9", colorant"#c3c2b7"
const SERIES = [("grid1", "Subs 1", colorant"#2a78d6"), ("grid2", "Subs 2", colorant"#eb6834")]

"Integer tick label with thousands separators: -12500 -> \"-12,500\"."
function with_commas(v)
    s = string(round(Int, abs(v)))
    s = reverse(join((join(c) for c in Iterators.partition(reverse(s), 3)), ","))
    return (round(Int, v) < 0 ? "-" : "") * s
end

# Lines: JuMP at every sweep point. Dots: OpenDSS at whole degrees only -- every point is
# checked numerically above, and a dot on each 0.25 deg step would chop the line into dashes.
function panel(model, quantity, title; legend = false)
    r = filter(x -> x.model == model, rows)
    dots = filter(q -> isinteger(q.dd), r)
    pick(res, k) = quantity == :P ? res.P_subs_kW[k] : res.Q_subs_kvar[k]
    p = plot(; title, legend, xlabel = "δ1 − δ2 (deg)", ylabel = quantity == :P ? "kW" : "kvar",
             yformatter = with_commas, xticks = -4:1:4, xlims = (-4.3, 4.9))
    hline!(p, [0.0]; color = AXISC, linewidth = 1, label = "")
    for (k, name, c) in SERIES
        plot!(p, [q.dd for q in r], [pick(q.ip, k) for q in r]; color = c, linewidth = 2,
              label = name * " (JuMP)")
        scatter!(p, [q.dd for q in dots], [pick(q.ds, k) for q in dots]; color = c, markersize = 5,
                 markerstrokecolor = SURFACE, markerstrokewidth = 1.5, label = "")
        annotate!(p, last(r).dd + 0.15, pick(last(r).ip, k), text(name, 8, INK2, :left))
    end
    scatter!(p, [NaN], [NaN]; color = MUTED, markersize = 5, markerstrokecolor = SURFACE,
             markerstrokewidth = 1.5, label = "OpenDSS")
    return p
end

default(; fontfamily = "sans-serif", background_color = SURFACE, foreground_color_axis = AXISC,
        foreground_color_border = AXISC, foreground_color_text = INK2, foreground_color_guide = INK2,
        foreground_color_title = colorant"#0b0b0b", gridcolor = GRIDC, gridalpha = 1.0,
        gridlinewidth = 1, gridstyle = :solid, titlefontsize = 11, guidefontsize = 10,
        tickfontsize = 9, legendfontsize = 9, legend_foreground_color = GRIDC)

deck_z = nets["deck"].sources[1].z_series
zs = @sprintf("Zs = %.4f + j%.4f Ω", real(deck_z), imag(deck_z))
fig = plot(panel("deck", :P, "P_Subs — deck sources ($zs)"; legend = :top),
           panel("stiff", :P, "P_Subs — stiff sources (Zs ≈ 0)"),
           panel("deck", :Q, "Q_Subs — deck sources"),
           panel("stiff", :Q, "Q_Subs — stiff sources");
           layout = (2, 2), link = :y, size = (1400, 950), left_margin = 8Plots.mm,
           bottom_margin = 6Plots.mm, top_margin = 3Plots.mm,
           plot_title = "small2poi_1ph: substation power vs source-angle difference (δ2 = 0)",
           plot_titlefontsize = 13)
png_path = joinpath(OUT_DIR, "psubs_vs_angle.png")
savefig(fig, png_path)

println("\n", "="^104)
solvers = USE_GUROBI ? "Ipopt and Gurobi" : "Ipopt; Gurobi SKIPPED"
println(ok_all ? "ALL SWEEP POINTS AGREE with OpenDSS ($solvers)." :
                 "SOME SWEEP POINTS DISAGREE with OpenDSS ($solvers) -- see the per-model summaries above.")
println("Wrote $csv")
println("Wrote $png_path")
