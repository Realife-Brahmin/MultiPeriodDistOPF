# battery_angle_sweep.jl -- can batteries let two substations run further apart in angle
# without either one back-feeding its grid?
#
# At each delta_a - delta_b on a grid (delta_b = 0), two solves of the full-angle model:
#   * no battery action (P_B = 0): the plain power flow of two_source_angle_sweep.jl;
#   * the OPF: every battery's output P_B free within +-kWrated, P_Subs >= 0 at both
#     substations, objective min sum(P_B^2) -- the smallest dispatch that stops back-feed.
#     Where it is infeasible, no dispatch within the ratings can.
# Both are replayed in OpenDSS (storage elements set to the dispatched kW) and compared.
# The edges of the window in which both substations import are refined by bisection.
#
# Run from the repo root:
#     julia --project=envs/tadmm envs/multi_poi/root_level/battery_angle_sweep.jl small2poi
#     julia --project=envs/tadmm envs/multi_poi/root_level/battery_angle_sweep.jl ieee123
#
# Writes (gitignored) envs/multi_poi/processedData/<system>/battery_sweep_<a>_<b>/:
#     sweep.csv, sweep.png
#
# small2poi: one test battery at the load bus rated at 100% of the load, added in memory --
# the deck on disk has none. ieee123: the deck's own 26 batteries, PV off, subs3 + subs4.
#
# A snapshot study of power only. P_B follows the repo's convention, positive discharging,
# so holding back-feed off shows up as charging, P_B < 0 -- and for how long a battery can
# keep charging is an energy question for the multi-period model.

include(joinpath(@__DIR__, "..", "full_angle_pf.jl"))
include(joinpath(@__DIR__, "..", "sweep_common.jl"))

const CONFIGS = Dict(
    "small2poi" => (system = "small2poi_1ph", pair = ("grid1", "grid2"), disable = String[],
                    snapshot = false, test_battery_bus = "1load",
                    ddelta = -8.0:0.1:8.0, dot_step = 1.0, xticks = -8:2:8),
    "ieee123" => (system = "ieee123_5poi_1ph", pair = ("subs3", "subs4"), disable = ["PVSystem"],
                  snapshot = true, test_battery_bus = nothing,
                  ddelta = -2.0:0.02:2.0, dot_step = 0.25, xticks = -2:0.5:2))
const LOAD_BAND = (0.5, 1.5)    # constant-PQ band for OpenDSS's loads and storage (see compile_deck)
const V_GUARD = 0.8             # keeps the OPF off the low-voltage power-flow root
const TOL_KW = 1e-3             # agreement demanded: 1 W, 1 var
const TOL_V = 1e-8              # agreement demanded: pu voltage phasor

const KEY = isempty(ARGS) ? "small2poi" : lowercase(ARGS[1])
haskey(CONFIGS, KEY) || error("unknown system '$KEY'; one of: $(join(sort(collect(keys(CONFIGS))), ", "))")
const CFG = CONFIGS[KEY]
const A, B = length(ARGS) >= 3 ? (lowercase(ARGS[2]), lowercase(ARGS[3])) : CFG.pair
const OUT_DIR = normpath(joinpath(@__DIR__, "..", "processedData", CFG.system, "battery_sweep_$(A)_$(B)"))

# ---- the network, with the test battery added in memory where the config asks for one -----------
setup = (; sources = [A, B], disable = CFG.disable, snapshot = CFG.snapshot, load_band = LOAD_BAND)
extra = String[]
if CFG.test_battery_bus !== nothing
    compile_deck(CFG.system; setup...)
    probe = read_network()
    rating = sum(d.kW for d in probe.loads)       # 100% of the load
    push!(extra, @sprintf("New Storage.battery1 phases=1 bus1=%s kV=%.6g kWrated=%.10g kVA=%.10g kWhrated=%.10g %%stored=50 %%reserve=0 %%IdlingkW=0 DispMode=External",
                          CFG.test_battery_bus, probe.kV_base, rating, rating, 4rating))
end
dk = compile_deck(CFG.system; extra, setup...)
const NET = read_network()
const BAND = dk.deck_load_band
const P_LOAD = sum(d.kW for d in NET.loads)
const B_RATED = sum(bt.kW_rated for bt in NET.batteries)
isempty(NET.batteries) && error("no enabled battery in $(CFG.system)")

println("\n", "="^104)
@printf "%s -- %s + %s%s\n" CFG.system A B (isempty(CFG.disable) ? "" : "; off: " * join(CFG.disable, ", "))
@printf "  %d buses, %d lines, %d loads: %.1f kW + j%.1f kvar%s\n" length(NET.buses) length(NET.lines) length(NET.loads) P_LOAD sum(d.kvar for d in NET.loads) (CFG.snapshot ? " (snapshot, base load)" : "")
@printf "  %d batter%s, %.1f kW in all (%.0f%% of the load)%s: %s\n" length(NET.batteries) (length(NET.batteries) == 1 ? "y" : "ies") B_RATED 100B_RATED / P_LOAD (isempty(extra) ? "" : ", added in memory") join(unique(bt.bus for bt in NET.batteries), ", ")
for s in NET.sources
    @printf "  Vsource.%-6s bus %-4s %.3f pu on %.4g kV (L-N)   Zs = %.6f + j%.2e ohm\n" s.name s.bus s.pu NET.kV_base real(s.z_series) imag(s.z_series)
end
println("="^104)

# ---- the sweep ---------------------------------------------------------------------------------------
opf(dd) = solve_full_angle(NET, Dict(A => dd, B => 0.0); battery = :optimize, no_backflow = true, V_guard = V_GUARD)
no_dss(ds) = merge(blank(ds), (; converged = false, pq_band = false))

rows = []
for dd in CFG.ddelta
    delta = Dict(A => dd, B => 0.0)
    pf, op = solve_full_angle(NET, delta), opf(dd)
    ds_pf = opendss_point(NET, delta)
    ds_op = op.ok ? opendss_point(NET, delta; P_B = op.P_B_kW) : no_dss(ds_pf)
    push!(rows, (; dd, pf = pf.ok ? pf : blank(pf), op = op.ok ? op : blank(op), op_status = op.status,
                 ds_pf = ds_pf.converged ? ds_pf : no_dss(ds_pf), ds_op))
end
dds = [x.dd for x in rows]
i0 = findfirst(==(0), dds)
total_B(x) = sum(values(x.P_B_kW); init = 0.0)
op_ok(x) = x.op_status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
op_cmp(x) = op_ok(x) && x.ds_op.converged && x.ds_op.pq_band     # OPF replayed, constant-PQ held
pf_cmp(x) = x.pf.ok && x.ds_pf.converged && x.ds_pf.pq_band

"Bisect between an angle where the OPF is feasible and one where it is not; returns the last
feasible angle and its solution."
function refine(inside, outside; iterations = 20)
    best = opf(inside)
    for _ in 1:iterations
        mid = (inside + outside) / 2
        trial = opf(mid)
        if trial.ok
            inside, best = mid, trial
        else
            outside = mid
        end
    end
    return inside, best
end

# ---- summary -----------------------------------------------------------------------------------------
n = length(rows)
println("\nSweep: $(dlabel(B)) = 0, $(dlabel(A)) = $(first(CFG.ddelta)) .. $(last(CFG.ddelta)) deg ($n points; table every $(CFG.dot_step) deg)")
@printf "  %8s | %11s %11s | %11s | %11s %11s | %10s %10s\n" "d_a-d_b" "P_$A" "P_$B" "battery" "P_$A" "P_$B" "losses" "losses"
@printf "  %8s | %23s | %11s | %23s | %10s %10s\n" "deg" "no battery action, kW" "P_B kW" "batteries dispatched, kW" "none, kW" "disp., kW"
for x in rows
    isinteger(round(x.dd / CFG.dot_step; digits = 9)) || continue
    if op_ok(x)
        @printf "  %8.2f | %11.1f %11.1f | %11.1f | %11.1f %11.1f | %10.2f %10.2f\n" x.dd x.pf.P_subs_kW[A] x.pf.P_subs_kW[B] total_B(x.op) x.op.P_subs_kW[A] x.op.P_subs_kW[B] losses(x.pf, P_LOAD) losses(x.op, P_LOAD)
    else
        @printf "  %8.2f | %11.1f %11.1f | %-49s | %10.2f\n" x.dd x.pf.P_subs_kW[A] x.pf.P_subs_kW[B] "  no dispatch within the ratings ($(x.op_status))" losses(x.pf, P_LOAD)
    end
end

fmt(v) = v === nothing ? "beyond the sweep" : @sprintf("%+.3f", v)
w0 = import_window(dds, [x.pf.P_subs_kW[A] for x in rows], [x.pf.P_subs_kW[B] for x in rows])
w0 === nothing && error("without battery action the substations do not both import at $(dlabel(A)) = $(dlabel(B))")
@printf "\n  no battery action:     both substations import for %s - %s from %s to %s deg\n" dlabel(A) dlabel(B) fmt(w0[1]) fmt(w0[2])
feas = [op_ok(x) for x in rows]
feas[i0] || error("the OPF is infeasible at $(dlabel(A)) = $(dlabel(B))")
l, r = stretch(feas, i0)
edges = [(l > 1 ? refine(dds[l], dds[l-1]) : nothing), (r < n ? refine(dds[r], dds[r+1]) : nothing)]
w1 = [e === nothing ? nothing : e[1] for e in edges]
@printf "  batteries dispatched: both substations import for %s - %s from %s to %s deg\n" dlabel(A) dlabel(B) fmt(w1[1]) fmt(w1[2])
for (side, e) in zip(("lower", "upper"), edges)
    e === nothing && continue
    at_limit = count(bt -> abs(abs(e[2].P_B_kW[bt.name]) - bt.kW_rated) <= 1e-3 * bt.kW_rated, NET.batteries)
    @printf "    at the %s edge the batteries charge %.1f kW (%.1f%% of their %.1f kW), %d of %d at their rating\n" side -total_B(e[2]) -100total_B(e[2]) / B_RATED B_RATED at_limit length(NET.batteries)
end
any(isnothing, (w0..., w1...)) ||
    @printf("  => the window widens from %.3f to %.3f deg (x%.1f)\n", w0[2] - w0[1], w1[2] - w1[1],
            (w1[2] - w1[1]) / (w0[2] - w0[1]))

# OPF behaviour, and the OpenDSS replays
ok_op = filter(op_ok, rows)
its = [x.op.iterations for x in ok_op if x.op.iterations !== missing]
@printf "\n  OPF solved %d/%d (iterations: mean %.1f, max %d; %.3f s per solve); the rest: %s\n" length(ok_op) n (isempty(its) ? NaN : sum(its) / length(its)) (isempty(its) ? -1 : maximum(its)) sum(x.op.time for x in ok_op; init = 0.0) / max(1, length(ok_op)) join(unique(string(x.op_status) for x in rows if !op_ok(x)), ", ")
cmp_pf = filter(pf_cmp, rows)
cmp_op = filter(op_cmp, rows)
worst(es, f) = isempty(es) ? NaN : maximum(f, es)
e_pf = [max_diff(x.pf, x.ds_pf) for x in cmp_pf]
e_op = [max_diff(x.op, x.ds_op) for x in cmp_op]
dB = worst([abs(x.ds_op.P_B_kW[k] - x.op.P_B_kW[k]) for x in cmp_op for k in keys(x.op.P_B_kW)], identity)
@printf "  OpenDSS replay, no battery action:     %d/%d compared, worst P %.2e kW  Q %.2e kvar  V %.2e pu\n" length(cmp_pf) n worst(e_pf, e -> e.dP) worst(e_pf, e -> e.dQ) worst(e_pf, e -> e.dV)
@printf "  OpenDSS replay, batteries dispatched:  %d/%d compared, worst P %.2e kW  Q %.2e kvar  V %.2e pu; battery output off by at most %.2e kW\n" length(cmp_op) length(ok_op) worst(e_op, e -> e.dP) worst(e_op, e -> e.dQ) worst(e_op, e -> e.dV) dB
passed = length(cmp_pf) == count(x -> x.pf.ok, rows) && length(cmp_op) == length(ok_op) && dB <= TOL_KW &&
         all(e -> e.dP <= TOL_KW && e.dQ <= TOL_KW && e.dV <= TOL_V, vcat(e_pf, e_op))
@printf "  => %s (tolerance %.0e kW/kvar, %.0e pu)\n" (passed ? "AGREE" : "DISAGREE") TOL_KW TOL_V
held = count(x -> in_band(NET, BAND, x.ds_op), cmp_op)
bs = unique(values(BAND))
@printf "  with the batteries dispatched, every load stays inside the deck's own band (%s) at %d/%d points\n" (length(bs) == 1 ? @sprintf("%.2f-%.2f pu", bs[1][1], bs[1][2]) : "per load") held length(cmp_op)

# ---- CSV ---------------------------------------------------------------------------------------------
mkpath(OUT_DIR)
csv = joinpath(OUT_DIR, "sweep.csv")
open(csv, "w") do io
    bn = [bt.name for bt in NET.batteries]
    println(io, join(vcat(["ddelta_deg", "P_$(A)_kW", "P_$(B)_kW", "losses_kW",
                           "opf_status", "opf_P_$(A)_kW", "opf_P_$(B)_kW", "opf_losses_kW", "opf_P_B_total_kW",
                           "opf_Vload_min_pu", "opf_Vload_max_pu",
                           "dss_P_$(A)_kW", "dss_P_$(B)_kW", "dss_opf_P_$(A)_kW", "dss_opf_P_$(B)_kW",
                           "dss_opf_P_B_total_kW", "err_pf_dss_P_kW", "err_opf_dss_P_kW", "err_opf_dss_V_pu"],
                          ["opf_P_B_$(k)_kW" for k in bn]), ","))
    for x in rows
        ep, eo = max_diff(x.pf, x.ds_pf), max_diff(x.op, x.ds_op)
        Vl = load_V(NET, x.op)
        vals = Any[x.dd, x.pf.P_subs_kW[A], x.pf.P_subs_kW[B], losses(x.pf, P_LOAD),
                   x.op_status, x.op.P_subs_kW[A], x.op.P_subs_kW[B], losses(x.op, P_LOAD), total_B(x.op),
                   minimum(Vl), maximum(Vl),
                   x.ds_pf.P_subs_kW[A], x.ds_pf.P_subs_kW[B], x.ds_op.P_subs_kW[A], x.ds_op.P_subs_kW[B],
                   total_B(x.ds_op), ep.dP, eo.dP, eo.dV]
        append!(vals, [x.op.P_B_kW[k] for k in bn])
        println(io, join((v isa AbstractFloat ? @sprintf("%.12g", v) : string(v) for v in vals), ","))
    end
end

# ---- plot --------------------------------------------------------------------------------------------
# Colour is the substation (end-labelled). The no-battery case is drawn thin and faded, the
# dispatched case bold; dots are the OpenDSS replays of the dispatched case.
plot_style!()
const SERIES = [(A, label(A), SUB_COLORS[1]), (B, label(B), SUB_COLORS[2])]
const DARK_WASH = colorant"#e1e0d9"
dots = filter(x -> isinteger(round(x.dd / CFG.dot_step; digits = 9)) && op_cmp(x), rows)
xd = [x.dd for x in dots]
span = last(dds) - first(dds)
thin = (linewidth = 1.2, linealpha = 0.45, label = "")
bold = (linewidth = 2, label = "")
dot_style = (markersize = 5, markerstrokecolor = SURFACE, markerstrokewidth = 1.5, label = "")
base_plot(title, ylabel; legend = false, yformatter = with_commas) =
    plot(; title, legend, ylabel, yformatter, xlabel = "$(dlabel(A)) − $(dlabel(B)) (deg)",
         xticks = CFG.xticks, xlims = (first(dds) - 0.03span, last(dds) + 0.14span))
function windows!(p)
    wb = [something(w1[1], first(dds)), something(w1[2], last(dds))]
    vspan!(p, wb; color = WASH, linecolor = WASH, label = "both import, with batteries")
    vspan!(p, [something(w0[1], first(dds)), something(w0[2], last(dds))]; color = DARK_WASH,
           linecolor = DARK_WASH, label = "both import, without")
    hline!(p, [0.0]; color = AXISC, linewidth = 1, label = "")
end

pP = base_plot("P_Subs", "kW"; legend = :top)
windows!(pP)
for (k, name, c) in SERIES
    plot!(pP, dds, [x.pf.P_subs_kW[k] for x in rows]; color = c, thin...)
    plot!(pP, dds, [x.op.P_subs_kW[k] for x in rows]; color = c, bold...)
    scatter!(pP, xd, [x.ds_op.P_subs_kW[k] for x in dots]; color = c, dot_style...)
    annotate!(pP, last(dds) + 0.02span, lastfinite([x.pf.P_subs_kW[k] for x in rows]), text(name, 8, INK2, :left))
end
plot!(pP, [NaN], [NaN]; color = INK2, linewidth = 1.2, linealpha = 0.45, label = "no battery")
plot!(pP, [NaN], [NaN]; color = INK2, linewidth = 2, label = "batteries dispatched")
scatter!(pP, [NaN], [NaN]; color = MUTED, markersize = 5, markerstrokecolor = SURFACE,
         markerstrokewidth = 1.5, label = "OpenDSS")

pB = base_plot("Battery output needed (P_B, discharging > 0)", "kW")
windows!(pB)
hline!(pB, [-B_RATED]; color = MUTED, linewidth = 1, label = "")
annotate!(pB, 0.0, -B_RATED, text("full charging rate, $(with_commas(B_RATED)) kW", 8, MUTED, :center, :bottom))
plot!(pB, dds, [total_B(x.op) for x in rows]; color = INK2, bold...)
scatter!(pB, xd, [total_B(x.ds_op) for x in dots]; color = INK2, dot_style...)

pL = base_plot("Network losses", "kW")
plot!(pL, dds, [losses(x.pf, P_LOAD) for x in rows]; color = INK2, thin...)
plot!(pL, dds, [losses(x.op, P_LOAD) for x in rows]; color = INK2, bold...)
scatter!(pL, xd, [losses(x.ds_op, P_LOAD) for x in dots]; color = INK2, dot_style...)

pV = base_plot("Load voltages, batteries dispatched", "pu"; yformatter = :auto)
lo, hi = [minimum(load_V(NET, x.op)) for x in rows], [maximum(load_V(NET, x.op)) for x in rows]
edges_in_view = filter(e -> minimum(filter(isfinite, lo)) - 0.002 <= e <= maximum(filter(isfinite, hi)) + 0.002,
                       unique(Iterators.flatten(values(BAND))))
foreach(e -> hline!(pV, [e]; color = MUTED, linewidth = 1, label = ""), edges_in_view)
isempty(edges_in_view) || plot!(pV; title = "Load voltages, batteries dispatched (gray: deck's band edge)")
plot!(pV, dds, hi; fillrange = lo, fillcolor = INK2, fillalpha = 0.10, color = INK2, bold...)
plot!(pV, dds, lo; color = INK2, bold...)

fig = plot(pP, pB, pL, pV; layout = (2, 2), size = (1400, 950), left_margin = 8Plots.mm,
           bottom_margin = 6Plots.mm, top_margin = 3Plots.mm, plot_titlefontsize = 13,
           plot_title = "$(CFG.system), $(label(A)) + $(label(B)): keeping both substations importing with batteries ($(dlabel(B)) = 0)")
png_path = joinpath(OUT_DIR, "sweep.png")
savefig(fig, png_path)

println("\n", "="^104)
println(passed ? "OpenDSS REPLAYS AGREE (Ipopt)." : "SOME OpenDSS REPLAYS DISAGREE -- see the summary above.")
println("Wrote $csv")
println("Wrote $png_path")
