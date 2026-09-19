# two_source_angle_sweep.jl -- how two substations share a network as their source angles
# move: P_Subs, Q_Subs and network losses against delta_a - delta_b, from the full-angle
# AC model in JuMP (Ipopt; Gurobi as a global check when licensed), with every point
# cross-checked against OpenDSS. Vsources other than the chosen two, and the element
# classes in the config's `disable` list (PV, storage), are switched off in memory --
# never on disk.
#
# Run from the repo root:
#     julia --project=envs/tadmm envs/multi_poi/root_level/two_source_angle_sweep.jl small2poi
#     julia --project=envs/tadmm envs/multi_poi/root_level/two_source_angle_sweep.jl ieee123
#     julia --project=envs/tadmm envs/multi_poi/root_level/two_source_angle_sweep.jl ieee123 subs1 subs4
#
# Writes (gitignored) envs/multi_poi/processedData/<system>/angle_sweep_<a>_<b>/:
#     sweep.csv   every point: JuMP-Ipopt, JuMP-Gurobi, OpenDSS, and their differences
#     sweep.png   P_Subs, Q_Subs and network losses against delta_a - delta_b
#
# ieee123 defaults to subs3 + subs4: electrically the farthest apart of the deck's five
# substations (1.48 ohm between them; every run prints the full table).
#
# Two source models per run, from the same deck:
#   stiff -- each Vsource's impedance replaced in memory by OpenDSS's 1e-8 ohm stand-in for
#            an ideal source, so delta is the substation bus angle itself, as the
#            branch-flow MPOPF formulations assume. An exactly ideal source (Zs = 0) is
#            solved alongside and the gap reported.
#   deck  -- each Vsource's impedance as written, so delta is the EMF angle behind it.
#
# Loads: OpenDSS keeps a model=1 load constant-PQ only inside [Vminpu, Vmaxpu] (0.95-1.05
# in ieee123_5poi_1ph) and makes it a constant impedance outside. The JuMP model is
# constant-PQ everywhere, like the MPOPF formulations, so the band is widened in memory for
# the cross-check; where voltages leave the deck's own band is reported separately.

include(joinpath(@__DIR__, "..", "full_angle_pf.jl"))
using Plots

# Per system: deck, default pair, what to switch off, sweep grid, plot ticks. Gurobi runs at
# the base case only on ieee123: a global solve of a 130-bus nonconvex QCQP at every sweep
# point could take hours.
const CONFIGS = Dict(
    "small2poi" => (system = "small2poi_1ph", pair = ("grid1", "grid2"), disable = String[],
                    snapshot = false, ddelta = -4.0:0.25:4.0, dot_step = 1.0, xticks = -4:1:4,
                    gurobi_everywhere = true),
    "ieee123" => (system = "ieee123_5poi_1ph", pair = ("subs3", "subs4"),
                  disable = ["PVSystem", "Storage"], snapshot = true,
                  ddelta = -2.0:0.05:2.0, dot_step = 0.25, xticks = -2:0.5:2,
                  gurobi_everywhere = false))
const MODELS = [("stiff", true), ("deck", false)]    # (name, stiff_sources)
const LOAD_BAND = (0.5, 1.5)                          # constant-PQ band used for the cross-check
const ROTATIONS = [0.0, -0.5, 14.0, -45.0, 120.0]     # common offsets for the invariance check
const TOL_KW = 1e-3                                   # agreement demanded: 1 W, 1 var
const TOL_V = 1e-8                                    # agreement demanded: pu voltage phasor

const KEY = isempty(ARGS) ? "small2poi" : lowercase(ARGS[1])
haskey(CONFIGS, KEY) || error("unknown system '$KEY'; one of: $(join(sort(collect(keys(CONFIGS))), ", "))")
const CFG = CONFIGS[KEY]
const A, B = length(ARGS) >= 3 ? (lowercase(ARGS[2]), lowercase(ARGS[3])) : CFG.pair
const OUT_DIR = normpath(joinpath(@__DIR__, "..", "processedData", CFG.system, "angle_sweep_$(A)_$(B)"))

"\"subs3\" -> \"Subs 3\", \"grid1\" -> \"Subs 1\"."
label(name) = (m = match(r"(\d+)$", name)) === nothing ? name : "Subs " * m.captures[1]
"\"subs3\" -> \"δ3\"."
dlabel(name) = (m = match(r"(\d+)$", name)) === nothing ? "δ_" * name : "δ" * m.captures[1]

const USE_GUROBI, GUROBI_WHY = gurobi_usable()
USE_GUROBI || println("\nGUROBI SKIPPED -- ", GUROBI_WHY, "\nRunning Ipopt + OpenDSS only.")

"One (delta_a, delta_b) point: Ipopt, Gurobi (if asked and licensed), OpenDSS, and Ipopt
with ideal sources (if asked)."
function solve_all(net, da, db; gurobi = true, with_ideal = false)
    delta = Dict(A => da, B => db)
    return (; da, db, dd = da - db,
            ip = solve_full_angle(net, delta; solver = :ipopt),
            gu = gurobi && USE_GUROBI ? solve_full_angle(net, delta; solver = :gurobi) : skipped(net, :gurobi),
            ds = opendss_point(net, delta),
            id = with_ideal ? solve_full_angle(net, delta; ideal_sources = true) : skipped(net, :ipopt))
end

"Voltage at every load, per unit of that load's own kV."
load_V(net, x) = [abs(x.V_pu[d.bus]) * net.kV_base / d.kV for d in net.loads]
"Does every load sit inside the (Vminpu, Vmaxpu) band given for it?"
function in_band(net, band, x)
    for (d, v) in zip(net.loads, load_V(net, x))
        lo, hi = band[d.name]
        lo <= v <= hi || return false
    end
    return true
end
losses(x, P_load) = x.P_subs_kW[A] + x.P_subs_kW[B] - P_load

"Where P_Subs of `src` crosses zero along the sweep (linear interpolation), or nothing."
function zero_crossing(r, src)
    y = [x.ip.P_subs_kW[src] for x in r]
    j = findfirst(i -> sign(y[i]) != sign(y[i+1]), 1:length(y)-1)
    j === nothing && return nothing
    return r[j].dd + (r[j+1].dd - r[j].dd) * y[j] / (y[j] - y[j+1])
end

ipopt_status(x) = string(x.status) * (x.iterations === missing ? "" : " ($(x.iterations) it)")

function print_row(label, status, x, net)
    Vl = load_V(net, x)
    @printf "  %-8s %-24s %12.3f %12.3f %12.3f %12.3f %9.5f %9.5f\n" label status x.P_subs_kW[A] x.P_subs_kW[B] x.Q_subs_kvar[A] x.Q_subs_kvar[B] minimum(Vl) maximum(Vl)
end
header() = @printf "  %-8s %-24s %12s %12s %12s %12s %9s %9s\n" "" "status" "P_$A kW" "P_$B kW" "Q_$A kvar" "Q_$B kvar" "Vmin pu" "Vmax pu"

# ---- which pair: electrical distance between every two substations in the deck -----------------
compile_deck(CFG.system; disable = CFG.disable, snapshot = CFG.snapshot)
all_sources = read_network()
if length(all_sources.sources) > 2
    snames = sort([s.name for s in all_sources.sources])
    sep = source_separation(all_sources)
    println("\n", "="^104)
    println("$(CFG.system): series impedance between substation buses (ohm, through the lines)")
    @printf "  %-8s" ""; foreach(n -> @printf("%16s", n), snames); println()
    for a in snames
        @printf "  %-8s" a
        for b in snames
            z = a == b ? nothing : sep[minmax(a, b)]
            @printf "%16s" (z === nothing ? "-" : @sprintf("%.3f+j%.3f", real(z), imag(z)))
        end
        println()
    end
    @printf "  chosen pair: %s + %s, |Z| = %.3f ohm\n" A B abs(sep[minmax(A, B)])
end

rows = []
nets, bands = Dict{String,Any}(), Dict{String,Any}()
for (model, stiff) in MODELS
    dk = compile_deck(CFG.system; sources = [A, B], disable = CFG.disable, snapshot = CFG.snapshot,
                      stiff_sources = stiff, load_band = LOAD_BAND)
    net = nets[model] = read_network()
    band = bands[model] = dk.deck_load_band
    P_load = sum(d.kW for d in net.loads)

    println("\n", "="^104)
    @printf "%s -- %s sources (%s + %s)%s\n" CFG.system model A B (isempty(CFG.disable) ? "" : "; off: " * join(CFG.disable, ", "))
    @printf "  %d buses, %d lines, %d loads: %.1f kW + j%.1f kvar%s\n" length(net.buses) length(net.lines) length(net.loads) P_load sum(d.kvar for d in net.loads) (CFG.snapshot ? " (snapshot, base load)" : "")
    for s in net.sources
        @printf "  Vsource.%-6s bus %-4s %.3f pu on %.4g kV (L-N)   Zs = %.6f + j%.2e ohm%s\n" s.name s.bus s.pu net.kV_base real(s.z_series) imag(s.z_series) (stiff ? "  (stand-in, set in memory)" : "")
    end
    println("="^104)

    # ---- base case: both angles zero ---------------------------------------------------
    b = solve_all(net, 0.0, 0.0)
    println("\nBase case, $(dlabel(A)) = $(dlabel(B)) = 0 (Vmin/Vmax over loads):")
    header()
    print_row("Ipopt", ipopt_status(b.ip), b.ip, net)
    USE_GUROBI ? print_row("Gurobi", string(b.gu.status), b.gu, net) :
                 @printf("  %-8s %s\n", "Gurobi", "SKIPPED (see message at top)")
    print_row("OpenDSS", b.ds.converged ? "converged ($(b.ds.iterations) it)" : "NOT CONVERGED", b.ds, net)
    for (who, x) in (("Ipopt", b.ip), ("Gurobi", b.gu))
        x.status == :SKIPPED && continue
        d = max_diff(x, b.ds)
        @printf "  |%s - OpenDSS|: P %.2e kW  Q %.2e kvar  V %.2e pu\n" who d.dP d.dQ d.dV
    end

    # With a single load bus the power-flow equations have exactly two roots; show the other
    # one (Ipopt started low, and Gurobi maximizing losses when licensed).
    if length(unique(d.bus for d in net.loads)) == 1
        zero_delta = Dict(A => 0.0, B => 0.0)
        other = [("Ipopt, load bus started at 0.02 pu", solve_full_angle(net, zero_delta; V_load_start = 0.02))]
        USE_GUROBI && push!(other, ("Gurobi, max losses",
                                    solve_full_angle(net, zero_delta; solver = :gurobi, sense = MAX_SENSE)))
        for (how, x) in other
            @printf "  Other power-flow root (%s): %s, |V_L| = %.4f pu, P_Subs total = %.0f kW\n" how string(x.status) minimum(load_V(net, x)) x.P_subs_kW[A] + x.P_subs_kW[B]
        end
    end

    # ---- only delta_a - delta_b should matter -----------------------------------------
    println("\nSame $(dlabel(A)) - $(dlabel(B)) = 1 deg, different common rotation:")
    @printf "  %9s %9s %15s %15s %15s %15s\n" dlabel(A) dlabel(B) "P_$A Ipopt" "P_$B Ipopt" "P_$A DSS" "P_$B DSS"
    same = [solve_all(net, 1.0 + c, c; gurobi = false) for c in ROTATIONS]
    for r in same
        @printf "  %9.2f %9.2f %15.4f %15.4f %15.4f %15.4f\n" r.da r.db r.ip.P_subs_kW[A] r.ip.P_subs_kW[B] r.ds.P_subs_kW[A] r.ds.P_subs_kW[B]
    end
    spread(f) = maximum(f(r) for r in same) - minimum(f(r) for r in same)
    @printf "  spread across rotations: Ipopt %.2e kW, OpenDSS %.2e kW\n" max(spread(r -> r.ip.P_subs_kW[A]), spread(r -> r.ip.P_subs_kW[B])) max(spread(r -> r.ds.P_subs_kW[A]), spread(r -> r.ds.P_subs_kW[B]))

    # ---- the sweep: delta_b = 0, delta_a over the grid ----------------------------------
    for dd in CFG.ddelta
        x = solve_all(net, dd, 0.0; gurobi = CFG.gurobi_everywhere, with_ideal = stiff)
        push!(rows, (; model, P_load, deck_ok = in_band(net, band, x.ds), x...))
    end
end

# ---- summary -------------------------------------------------------------------------------------
ok_all = true
for (model, _) in MODELS
    r = filter(x -> x.model == model, rows)
    net, n, P_load = nets[model], length(r), first(r).P_load
    println("\n", "="^104)
    @printf "Sweep, %s sources: %s = 0, %s = %.2f .. %.2f deg (%d points; table every %.2g deg)\n" model dlabel(B) dlabel(A) first(CFG.ddelta) last(CFG.ddelta) n CFG.dot_step
    @printf "  %9s %12s %12s %12s %12s %9s %9s %10s\n" "d_a-d_b" "P_$A kW" "P_$B kW" "Q_$A kvar" "Q_$B kvar" "Vmin pu" "Vmax pu" "losses kW"
    for x in r
        isinteger(round(x.dd / CFG.dot_step; digits = 9)) || continue
        Vl = load_V(net, x.ip)
        @printf "  %9.2f %12.2f %12.2f %12.2f %12.2f %9.5f %9.5f %10.2f\n" x.dd x.ip.P_subs_kW[A] x.ip.P_subs_kW[B] x.ip.Q_subs_kvar[A] x.ip.Q_subs_kvar[B] minimum(Vl) maximum(Vl) losses(x.ip, P_load)
    end

    i0 = findfirst(x -> x.dd == 0, r)
    slope = (r[i0+1].ip.P_subs_kW[A] - r[i0-1].ip.P_subs_kW[A]) / (r[i0+1].dd - r[i0-1].dd)
    @printf "  dP_%s/d(%s - %s) at 0: %.0f kW/deg\n" A dlabel(A) dlabel(B) slope
    zA, zB = zero_crossing(r, A), zero_crossing(r, B)
    if zA !== nothing && zB !== nothing
        w1, w2 = extrema((zA, zB))
        @printf "  both substations import for %s - %s in [%+.3f, %+.3f] deg (%.3f deg wide)\n" dlabel(A) dlabel(B) w1 w2 w2 - w1
    else
        println("  a substation does not cross zero within the sweep")
    end
    L = [losses(x.ip, P_load) for x in r]
    i = argmin(L)
    if 1 < i < n                                  # parabola through the three lowest points
        h, (y0, y1, y2) = r[i+1].dd - r[i].dd, (L[i-1], L[i], L[i+1])
        xs = r[i].dd - h * (y2 - y0) / (2 * (y2 - 2y1 + y0))
        @printf "  network losses are least, %.3f kW, at %s - %s = %+.3f deg\n" y1 - (y2 - y0)^2 / (8 * (y2 - 2y1 + y0)) dlabel(A) dlabel(B) xs
    end
    held = [x.dd for x in r if x.deck_ok]
    bs = unique(values(bands[model]))
    band_str = length(bs) == 1 ? @sprintf("%.2f-%.2f pu", bs[1][1], bs[1][2]) : "per load"
    if length(held) == n
        @printf "  every load stays inside the deck's own band (%s) across the sweep\n" band_str
    elseif isempty(held)
        @printf "  some load is outside the deck's own band (%s) at every point\n" band_str
    else
        @printf "  every load inside the deck's own band (%s) at %d/%d points, %s - %s from %+.2f to %+.2f\n" band_str length(held) n dlabel(A) dlabel(B) minimum(held) maximum(held)
    end

    e_ip = [max_diff(x.ip, x.ds) for x in r]
    ran_gu = filter(x -> x.gu.status != :SKIPPED, r)
    e_gu = [max_diff(x.gu, x.ds) for x in ran_gu]
    worst(es, f) = maximum(f, es)
    its = [x.ip.iterations for x in r if x.ip.iterations !== missing]
    n_ip, n_conv, n_pq = count(x -> x.ip.ok, r), count(x -> x.ds.converged, r), count(x -> x.ds.pq_band, r)
    @printf "  Ipopt solved %d/%d (iterations: mean %.1f, max %d; %.3f s per solve), OpenDSS converged %d/%d, constant-PQ held %d/%d\n" n_ip n (isempty(its) ? NaN : sum(its) / length(its)) (isempty(its) ? -1 : maximum(its)) sum(x.ip.time for x in r) / n n_conv n n_pq n
    @printf "  Gurobi (global): %s\n" (!USE_GUROBI ? "SKIPPED" : isempty(ran_gu) ? "base case only" : "solved $(count(x -> x.gu.ok, ran_gu))/$(length(ran_gu))")
    @printf "  worst |Ipopt - OpenDSS|:  P %.2e kW  Q %.2e kvar  V %.2e pu\n" worst(e_ip, e -> e.dP) worst(e_ip, e -> e.dQ) worst(e_ip, e -> e.dV)
    isempty(e_gu) || @printf("  worst |Gurobi - OpenDSS|: P %.2e kW  Q %.2e kvar  V %.2e pu\n",
                             worst(e_gu, e -> e.dP), worst(e_gu, e -> e.dQ), worst(e_gu, e -> e.dV))
    passed = n_ip == n_conv == n_pq == n && all(x -> x.gu.ok, ran_gu) &&
             all(e -> e.dP <= TOL_KW && e.dQ <= TOL_KW && e.dV <= TOL_V, vcat(e_ip, e_gu))
    global ok_all &= passed
    @printf "  => %s (tolerance %.0e kW/kvar, %.0e pu)\n" (passed ? "AGREE" : "DISAGREE") TOL_KW TOL_V
    ran_id = filter(x -> x.id.status != :SKIPPED, r)
    if !isempty(ran_id)
        good = filter(x -> x.id.ok, ran_id)
        d_ideal = isempty(good) ? NaN : maximum(max_diff(x.id, x.ip).dP for x in good)
        @printf "  ideal source (Zs = 0) vs the 1e-8 ohm stand-in: max |ΔP_Subs| %.2e kW (ideal solved %d/%d)\n" d_ideal length(good) length(ran_id)
    end
    d_vs = maximum(abs(x.ds.S_vsource_kVA[k] - complex(x.ds.P_subs_kW[k], x.ds.Q_subs_kvar[k])) for x in r for k in (A, B))
    @printf "  OpenDSS's own Vsource power reading vs line-flow reading: max %.2e kVA\n" d_vs
end

# ---- CSV ---------------------------------------------------------------------------------------
mkpath(OUT_DIR)
csv = joinpath(OUT_DIR, "sweep.csv")
open(csv, "w") do io
    cols = ["model", "delta_$(A)_deg", "delta_$(B)_deg", "ddelta_deg",
            "P_$(A)_kW", "P_$(B)_kW", "Q_$(A)_kvar", "Q_$(B)_kvar", "Vload_min_pu", "Vload_max_pu",
            "network_losses_kW", "deck_load_band_ok", "ipopt_status", "ipopt_iterations", "ipopt_time_s",
            "gurobi_P_$(A)_kW", "gurobi_P_$(B)_kW", "gurobi_status", "gurobi_time_s",
            "dss_P_$(A)_kW", "dss_P_$(B)_kW", "dss_Q_$(A)_kvar", "dss_Q_$(B)_kvar",
            "dss_converged", "dss_iterations", "dss_constant_pq",
            "dss_vsource_P_$(A)_kW", "dss_vsource_P_$(B)_kW",
            "err_ipopt_dss_P_kW", "err_ipopt_dss_Q_kvar", "err_ipopt_dss_V_pu",
            "err_gurobi_dss_P_kW", "err_gurobi_dss_Q_kvar", "err_gurobi_dss_V_pu",
            "ideal_P_$(A)_kW", "ideal_P_$(B)_kW"]
    println(io, join(cols, ","))
    for x in rows
        net = nets[x.model]
        ei, eg = max_diff(x.ip, x.ds), max_diff(x.gu, x.ds)
        Vl = load_V(net, x.ip)
        vals = Any[x.model, x.da, x.db, x.dd,
            x.ip.P_subs_kW[A], x.ip.P_subs_kW[B], x.ip.Q_subs_kvar[A], x.ip.Q_subs_kvar[B],
            minimum(Vl), maximum(Vl), losses(x.ip, x.P_load), x.deck_ok,
            x.ip.status, x.ip.iterations, x.ip.time,
            x.gu.P_subs_kW[A], x.gu.P_subs_kW[B], x.gu.status, x.gu.time,
            x.ds.P_subs_kW[A], x.ds.P_subs_kW[B], x.ds.Q_subs_kvar[A], x.ds.Q_subs_kvar[B],
            x.ds.converged, x.ds.iterations, x.ds.pq_band,
            real(x.ds.S_vsource_kVA[A]), real(x.ds.S_vsource_kVA[B]),
            ei.dP, ei.dQ, ei.dV, eg.dP, eg.dQ, eg.dV,
            x.id.P_subs_kW[A], x.id.P_subs_kW[B]]
        println(io, join((v isa AbstractFloat ? @sprintf("%.12g", v) : string(v) for v in vals), ","))
    end
end

# ---- plot ----------------------------------------------------------------------------------------
# Reference palette on the light surface: categorical slots 1 and 2 (validated) for the two
# substations, secondary ink for the single-series loss panels, a neutral wash for the window
# in which both substations import.
const SURFACE = colorant"#fcfcfb"
const INK2, MUTED = colorant"#52514e", colorant"#898781"
const GRIDC, AXISC, WASH = colorant"#e1e0d9", colorant"#c3c2b7", colorant"#f0efec"
const SERIES = [(A, label(A), colorant"#2a78d6"), (B, label(B), colorant"#eb6834")]

"Integer tick label with thousands separators: -12500 -> \"-12,500\"."
function with_commas(v)
    s = string(round(Int, abs(v)))
    s = reverse(join((join(c) for c in Iterators.partition(reverse(s), 3)), ","))
    return (round(Int, v) < 0 ? "-" : "") * s
end

# Lines: JuMP at every sweep point. Dots: OpenDSS every `dot_step` -- every point is checked
# numerically above, and a dot on each grid step would chop the line into dashes.
function panel(model, quantity, title; legend = false)
    r = filter(x -> x.model == model, rows)
    dots = filter(x -> isinteger(round(x.dd / CFG.dot_step; digits = 9)), r)
    span = last(CFG.ddelta) - first(CFG.ddelta)
    ylabel = quantity == :Q ? "kvar" : "kW"
    p = plot(; title, legend, xlabel = "$(dlabel(A)) − $(dlabel(B)) (deg)", ylabel,
             yformatter = with_commas, xticks = CFG.xticks,
             xlims = (first(CFG.ddelta) - 0.03span, last(CFG.ddelta) + 0.12span))
    if quantity == :P
        zA, zB = zero_crossing(r, A), zero_crossing(r, B)
        (zA === nothing || zB === nothing) ||
            vspan!(p, collect(extrema((zA, zB))); color = WASH, linecolor = WASH, label = "both import")
    end
    hline!(p, [0.0]; color = AXISC, linewidth = 1, label = "")
    if quantity == :loss
        plot!(p, [x.dd for x in r], [losses(x.ip, x.P_load) for x in r]; color = INK2, linewidth = 2, label = "")
        scatter!(p, [x.dd for x in dots], [losses(x.ds, x.P_load) for x in dots]; color = INK2,
                 markersize = 5, markerstrokecolor = SURFACE, markerstrokewidth = 1.5, label = "")
        return p
    end
    pick(res, k) = quantity == :P ? res.P_subs_kW[k] : res.Q_subs_kvar[k]
    for (k, name, c) in SERIES
        plot!(p, [x.dd for x in r], [pick(x.ip, k) for x in r]; color = c, linewidth = 2,
              label = name * " (JuMP)")
        scatter!(p, [x.dd for x in dots], [pick(x.ds, k) for x in dots]; color = c, markersize = 5,
                 markerstrokecolor = SURFACE, markerstrokewidth = 1.5, label = "")
        annotate!(p, last(r).dd + 0.02span, pick(last(r).ip, k), text(name, 8, INK2, :left))
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

zs = unique(s.z_series for s in nets["deck"].sources)
deck_title = length(zs) == 1 ? @sprintf("deck sources (Zs = %.4f %s j%.4f Ω)", real(zs[1]),
                                        imag(zs[1]) < 0 ? "−" : "+", abs(imag(zs[1]))) :
                               "deck sources (Zs as written)"
fig = plot(panel("stiff", :P, "P_Subs — stiff sources (Zs ≈ 0)"; legend = :top),
           panel("deck", :P, "P_Subs — " * deck_title),
           panel("stiff", :Q, "Q_Subs — stiff sources"),
           panel("deck", :Q, "Q_Subs — deck sources"),
           panel("stiff", :loss, "Network losses — stiff sources"),
           panel("deck", :loss, "Network losses — deck sources");
           layout = (3, 2), link = :y, size = (1400, 1350), left_margin = 8Plots.mm,
           bottom_margin = 6Plots.mm, top_margin = 3Plots.mm,
           plot_title = "$(CFG.system), $(label(A)) + $(label(B)): substation power vs source-angle difference ($(dlabel(B)) = 0)",
           plot_titlefontsize = 13)
png_path = joinpath(OUT_DIR, "sweep.png")
savefig(fig, png_path)

println("\n", "="^104)
solvers = USE_GUROBI ? "Ipopt and Gurobi" : "Ipopt; Gurobi SKIPPED"
println(ok_all ? "ALL SWEEP POINTS AGREE with OpenDSS ($solvers)." :
                 "SOME SWEEP POINTS DISAGREE with OpenDSS ($solvers) -- see the per-model summaries above.")
println("Wrote $csv")
println("Wrote $png_path")
