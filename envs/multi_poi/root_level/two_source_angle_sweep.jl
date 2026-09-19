# two_source_angle_sweep.jl -- how two substations share a network as their source angles
# move: P_Subs, Q_Subs, network losses and load voltages against delta_a - delta_b, from
# the full-angle AC model in JuMP (Ipopt; Gurobi as a global check when licensed), with
# every point cross-checked against OpenDSS. Vsources other than the chosen two, and the
# element classes in the config's `disable` list (PV, storage), are switched off in
# memory -- never on disk.
#
# Run from the repo root:
#     julia --project=envs/tadmm envs/multi_poi/root_level/two_source_angle_sweep.jl small2poi
#     julia --project=envs/tadmm envs/multi_poi/root_level/two_source_angle_sweep.jl ieee123
#     julia --project=envs/tadmm envs/multi_poi/root_level/two_source_angle_sweep.jl ieee123 subs1 subs4
#
# Writes (gitignored) envs/multi_poi/processedData/<system>/angle_sweep_<a>_<b>/:
#     sweep.csv   every point: JuMP-Ipopt, JuMP-Gurobi, OpenDSS, and their differences
#     sweep.png   P_Subs, Q_Subs, network losses and load voltages against delta_a - delta_b
#
# ieee123 defaults to subs3 + subs4: electrically the farthest apart of the deck's five
# substations (1.48 ohm between them; every run prints the full table).
#
# Sources: every substation in these decks is stiff -- the 1e-8 ohm stand-in for an ideal
# source that the single-source decks use -- so delta is the substation bus angle itself,
# as the branch-flow MPOPF formulations assume. An exactly ideal source (Zs = 0) is solved
# alongside and the gap reported.
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

# ---- the deck with only the chosen pair ----------------------------------------------------------
dk = compile_deck(CFG.system; sources = [A, B], disable = CFG.disable, snapshot = CFG.snapshot,
                  load_band = LOAD_BAND)
const NET = read_network()
const BAND = dk.deck_load_band
const P_LOAD = sum(d.kW for d in NET.loads)

println("\n", "="^104)
@printf "%s -- %s + %s%s\n" CFG.system A B (isempty(CFG.disable) ? "" : "; off: " * join(CFG.disable, ", "))
@printf "  %d buses, %d lines, %d loads: %.1f kW + j%.1f kvar%s\n" length(NET.buses) length(NET.lines) length(NET.loads) P_LOAD sum(d.kvar for d in NET.loads) (CFG.snapshot ? " (snapshot, base load)" : "")
for s in NET.sources
    @printf "  Vsource.%-6s bus %-4s %.3f pu on %.4g kV (L-N)   Zs = %.6f + j%.2e ohm\n" s.name s.bus s.pu NET.kV_base real(s.z_series) imag(s.z_series)
end
println("="^104)

# ---- base case: both angles zero -----------------------------------------------------------------
b = solve_all(NET, 0.0, 0.0)
println("\nBase case, $(dlabel(A)) = $(dlabel(B)) = 0 (Vmin/Vmax over loads):")
header()
print_row("Ipopt", ipopt_status(b.ip), b.ip, NET)
USE_GUROBI ? print_row("Gurobi", string(b.gu.status), b.gu, NET) :
             @printf("  %-8s %s\n", "Gurobi", "SKIPPED (see message at top)")
print_row("OpenDSS", b.ds.converged ? "converged ($(b.ds.iterations) it)" : "NOT CONVERGED", b.ds, NET)
for (who, x) in (("Ipopt", b.ip), ("Gurobi", b.gu))
    x.status == :SKIPPED && continue
    d = max_diff(x, b.ds)
    @printf "  |%s - OpenDSS|: P %.2e kW  Q %.2e kvar  V %.2e pu\n" who d.dP d.dQ d.dV
end

# With a single load bus the power-flow equations have exactly two roots; show the other one
# (Ipopt started low, and Gurobi maximizing losses when licensed).
if length(unique(d.bus for d in NET.loads)) == 1
    zero_delta = Dict(A => 0.0, B => 0.0)
    other = [("Ipopt, load bus started at 0.02 pu", solve_full_angle(NET, zero_delta; V_load_start = 0.02))]
    USE_GUROBI && push!(other, ("Gurobi, max losses",
                                solve_full_angle(NET, zero_delta; solver = :gurobi, sense = MAX_SENSE)))
    for (how, x) in other
        @printf "  Other power-flow root (%s): %s, |V_L| = %.4f pu, P_Subs total = %.0f kW\n" how string(x.status) minimum(load_V(NET, x)) x.P_subs_kW[A] + x.P_subs_kW[B]
    end
end

# ---- only delta_a - delta_b should matter --------------------------------------------------------
println("\nSame $(dlabel(A)) - $(dlabel(B)) = 1 deg, different common rotation:")
@printf "  %9s %9s %15s %15s %15s %15s\n" dlabel(A) dlabel(B) "P_$A Ipopt" "P_$B Ipopt" "P_$A DSS" "P_$B DSS"
same = [solve_all(NET, 1.0 + c, c; gurobi = false) for c in ROTATIONS]
for r in same
    @printf "  %9.2f %9.2f %15.4f %15.4f %15.4f %15.4f\n" r.da r.db r.ip.P_subs_kW[A] r.ip.P_subs_kW[B] r.ds.P_subs_kW[A] r.ds.P_subs_kW[B]
end
spread(f) = maximum(f(r) for r in same) - minimum(f(r) for r in same)
@printf "  spread across rotations: Ipopt %.2e kW, OpenDSS %.2e kW\n" max(spread(r -> r.ip.P_subs_kW[A]), spread(r -> r.ip.P_subs_kW[B])) max(spread(r -> r.ds.P_subs_kW[A]), spread(r -> r.ds.P_subs_kW[B]))

# ---- the sweep: delta_b = 0, delta_a over the grid -----------------------------------------------
rows = []
for dd in CFG.ddelta
    x = solve_all(NET, dd, 0.0; gurobi = CFG.gurobi_everywhere, with_ideal = true)
    push!(rows, (; deck_ok = in_band(NET, BAND, x.ds), x...))
end

# ---- summary -------------------------------------------------------------------------------------
n = length(rows)
println("\n", "="^104)
@printf "Sweep: %s = 0, %s = %.2f .. %.2f deg (%d points; table every %.2g deg)\n" dlabel(B) dlabel(A) first(CFG.ddelta) last(CFG.ddelta) n CFG.dot_step
@printf "  %9s %12s %12s %12s %12s %9s %9s %10s\n" "d_a-d_b" "P_$A kW" "P_$B kW" "Q_$A kvar" "Q_$B kvar" "Vmin pu" "Vmax pu" "losses kW"
for x in rows
    isinteger(round(x.dd / CFG.dot_step; digits = 9)) || continue
    Vl = load_V(NET, x.ip)
    @printf "  %9.2f %12.2f %12.2f %12.2f %12.2f %9.5f %9.5f %10.2f\n" x.dd x.ip.P_subs_kW[A] x.ip.P_subs_kW[B] x.ip.Q_subs_kvar[A] x.ip.Q_subs_kvar[B] minimum(Vl) maximum(Vl) losses(x.ip, P_LOAD)
end

i0 = findfirst(x -> x.dd == 0, rows)
slope = (rows[i0+1].ip.P_subs_kW[A] - rows[i0-1].ip.P_subs_kW[A]) / (rows[i0+1].dd - rows[i0-1].dd)
@printf "  dP_%s/d(%s - %s) at 0: %.0f kW/deg\n" A dlabel(A) dlabel(B) slope
zA, zB = zero_crossing(rows, A), zero_crossing(rows, B)
if zA !== nothing && zB !== nothing
    w1, w2 = extrema((zA, zB))
    @printf "  both substations import for %s - %s in [%+.3f, %+.3f] deg (%.3f deg wide)\n" dlabel(A) dlabel(B) w1 w2 w2 - w1
else
    println("  a substation does not cross zero within the sweep")
end
L = [losses(x.ip, P_LOAD) for x in rows]
i = argmin(L)
if 1 < i < n                                      # parabola through the three lowest points
    h, (y0, y1, y2) = rows[i+1].dd - rows[i].dd, (L[i-1], L[i], L[i+1])
    xs = rows[i].dd - h * (y2 - y0) / (2 * (y2 - 2y1 + y0))
    @printf "  network losses are least, %.3f kW, at %s - %s = %+.3f deg\n" y1 - (y2 - y0)^2 / (8 * (y2 - 2y1 + y0)) dlabel(A) dlabel(B) xs
end
held = [x.dd for x in rows if x.deck_ok]
bs = unique(values(BAND))
band_str = length(bs) == 1 ? @sprintf("%.2f-%.2f pu", bs[1][1], bs[1][2]) : "per load"
if length(held) == n
    @printf "  every load stays inside the deck's own band (%s) across the sweep\n" band_str
elseif isempty(held)
    @printf "  some load is outside the deck's own band (%s) at every point\n" band_str
else
    @printf "  every load inside the deck's own band (%s) at %d/%d points, %s - %s from %+.2f to %+.2f\n" band_str length(held) n dlabel(A) dlabel(B) minimum(held) maximum(held)
end

e_ip = [max_diff(x.ip, x.ds) for x in rows]
ran_gu = filter(x -> x.gu.status != :SKIPPED, rows)
e_gu = [max_diff(x.gu, x.ds) for x in ran_gu]
worst(es, f) = maximum(f, es)
its = [x.ip.iterations for x in rows if x.ip.iterations !== missing]
n_ip, n_conv, n_pq = count(x -> x.ip.ok, rows), count(x -> x.ds.converged, rows), count(x -> x.ds.pq_band, rows)
@printf "  Ipopt solved %d/%d (iterations: mean %.1f, max %d; %.3f s per solve), OpenDSS converged %d/%d, constant-PQ held %d/%d\n" n_ip n (isempty(its) ? NaN : sum(its) / length(its)) (isempty(its) ? -1 : maximum(its)) sum(x.ip.time for x in rows) / n n_conv n n_pq n
@printf "  Gurobi (global): %s\n" (!USE_GUROBI ? "SKIPPED" : isempty(ran_gu) ? "base case only" : "solved $(count(x -> x.gu.ok, ran_gu))/$(length(ran_gu))")
@printf "  worst |Ipopt - OpenDSS|:  P %.2e kW  Q %.2e kvar  V %.2e pu\n" worst(e_ip, e -> e.dP) worst(e_ip, e -> e.dQ) worst(e_ip, e -> e.dV)
isempty(e_gu) || @printf("  worst |Gurobi - OpenDSS|: P %.2e kW  Q %.2e kvar  V %.2e pu\n",
                         worst(e_gu, e -> e.dP), worst(e_gu, e -> e.dQ), worst(e_gu, e -> e.dV))
passed = n_ip == n_conv == n_pq == n && all(x -> x.gu.ok, ran_gu) &&
         all(e -> e.dP <= TOL_KW && e.dQ <= TOL_KW && e.dV <= TOL_V, vcat(e_ip, e_gu))
@printf "  => %s (tolerance %.0e kW/kvar, %.0e pu)\n" (passed ? "AGREE" : "DISAGREE") TOL_KW TOL_V
good = filter(x -> x.id.ok, rows)
d_ideal = isempty(good) ? NaN : maximum(max_diff(x.id, x.ip).dP for x in good)
@printf "  ideal source (Zs = 0) vs the decks' 1e-8 ohm stand-in: max |ΔP_Subs| %.2e kW (ideal solved %d/%d)\n" d_ideal length(good) n
d_vs = maximum(abs(x.ds.S_vsource_kVA[k] - complex(x.ds.P_subs_kW[k], x.ds.Q_subs_kvar[k])) for x in rows for k in (A, B))
@printf "  OpenDSS's own Vsource power reading vs line-flow reading: max %.2e kVA\n" d_vs

# ---- CSV ---------------------------------------------------------------------------------------
mkpath(OUT_DIR)
csv = joinpath(OUT_DIR, "sweep.csv")
open(csv, "w") do io
    cols = ["delta_$(A)_deg", "delta_$(B)_deg", "ddelta_deg",
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
        ei, eg = max_diff(x.ip, x.ds), max_diff(x.gu, x.ds)
        Vl = load_V(NET, x.ip)
        vals = Any[x.da, x.db, x.dd,
            x.ip.P_subs_kW[A], x.ip.P_subs_kW[B], x.ip.Q_subs_kvar[A], x.ip.Q_subs_kvar[B],
            minimum(Vl), maximum(Vl), losses(x.ip, P_LOAD), x.deck_ok,
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
# substations; secondary ink for the single-quantity panels; a neutral wash for the window in
# which both substations import.
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
function panel(quantity, title; legend = false)
    xs = [x.dd for x in rows]
    dots = filter(x -> isinteger(round(x.dd / CFG.dot_step; digits = 9)), rows)
    xd = [x.dd for x in dots]
    span = last(CFG.ddelta) - first(CFG.ddelta)
    endx = last(xs) + 0.02span
    p = plot(; title, legend, xlabel = "$(dlabel(A)) − $(dlabel(B)) (deg)",
             ylabel = Dict(:P => "kW", :Q => "kvar", :loss => "kW", :V => "pu")[quantity],
             yformatter = quantity == :V ? :auto : with_commas, xticks = CFG.xticks,
             xlims = (first(xs) - 0.03span, last(xs) + 0.16span))
    dot_style = (markersize = 5, markerstrokecolor = SURFACE, markerstrokewidth = 1.5, label = "")

    if quantity == :V                             # the spread of load voltages, deck band edges in view
        vmin, vmax = x -> minimum(load_V(NET, x)), x -> maximum(load_V(NET, x))
        lo, hi = [vmin(x.ip) for x in rows], [vmax(x.ip) for x in rows]
        edges = filter(e -> minimum(lo) - 0.002 <= e <= maximum(hi) + 0.002, unique(Iterators.flatten(values(BAND))))
        for edge in edges
            hline!(p, [edge]; color = MUTED, linewidth = 1, label = "")
        end
        isempty(edges) || plot!(p; title = title * " (gray line: deck's band edge, $(join(edges, ", ")) pu)")
        plot!(p, xs, hi; fillrange = lo, fillcolor = INK2, fillalpha = 0.10, color = INK2, linewidth = 2, label = "")
        plot!(p, xs, lo; color = INK2, linewidth = 2, label = "")
        scatter!(p, xd, [vmax(x.ds) for x in dots]; color = INK2, dot_style...)
        scatter!(p, xd, [vmin(x.ds) for x in dots]; color = INK2, dot_style...)
        annotate!(p, endx, last(hi), text("highest load", 8, INK2, :left))
        annotate!(p, endx, last(lo), text("lowest load", 8, INK2, :left))
        return p
    end

    if quantity == :P
        zA, zB = zero_crossing(rows, A), zero_crossing(rows, B)
        (zA === nothing || zB === nothing) ||
            vspan!(p, collect(extrema((zA, zB))); color = WASH, linecolor = WASH, label = "both import")
    end
    hline!(p, [0.0]; color = AXISC, linewidth = 1, label = "")
    if quantity == :loss
        plot!(p, xs, [losses(x.ip, P_LOAD) for x in rows]; color = INK2, linewidth = 2, label = "")
        scatter!(p, xd, [losses(x.ds, P_LOAD) for x in dots]; color = INK2, dot_style...)
        return p
    end
    pick(res, k) = quantity == :P ? res.P_subs_kW[k] : res.Q_subs_kvar[k]
    for (k, name, c) in SERIES
        plot!(p, xs, [pick(x.ip, k) for x in rows]; color = c, linewidth = 2, label = name * " (JuMP)")
        scatter!(p, xd, [pick(x.ds, k) for x in dots]; color = c, dot_style...)
        annotate!(p, endx, pick(last(rows).ip, k), text(name, 8, INK2, :left))
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

fig = plot(panel(:P, "P_Subs"; legend = :top), panel(:Q, "Q_Subs"),
           panel(:loss, "Network losses"), panel(:V, "Load voltages");
           layout = (2, 2), size = (1400, 950), left_margin = 8Plots.mm, bottom_margin = 6Plots.mm,
           top_margin = 3Plots.mm, plot_titlefontsize = 13,
           plot_title = "$(CFG.system), $(label(A)) + $(label(B)), stiff substations: power vs source-angle difference ($(dlabel(B)) = 0)")
png_path = joinpath(OUT_DIR, "sweep.png")
savefig(fig, png_path)

println("\n", "="^104)
solvers = USE_GUROBI ? "Ipopt and Gurobi" : "Ipopt; Gurobi SKIPPED"
println(passed ? "ALL SWEEP POINTS AGREE with OpenDSS ($solvers)." :
                 "SOME SWEEP POINTS DISAGREE with OpenDSS ($solvers) -- see the summary above.")
println("Wrote $csv")
println("Wrote $png_path")
