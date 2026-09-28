# source_impedance_study.jl -- how the substation source impedance Zs shapes the no-backflow
# angle window, with and without batteries and the 0.95-1.05 pu voltage band.
#
# A point delta_a - delta_b (delta_b = 0) is FEASIBLE when, besides the AC equations, both
# substations import (P_Subs >= 0) and -- where the band is on -- every bus not set by a
# substation stays within 0.95-1.05 pu. Four windows per Zs:
#   * PF       plain power flow, batteries idle, limits checked on the solution;
#   * PF+band  the same, band included;
#   * OPF      batteries free within +-kWrated, P_Subs >= 0 imposed, min sum(P_B^2);
#   * OPF+band the same with the band imposed (V_band) in place of the low-root guard.
# Each window's edges are found by stepping outward from 0 (0.1 deg) to the first infeasible
# point and bisecting that step, so the window is the interval connected to equal angles. A
# power-flow point also has to be on the high-voltage root (every non-source |V| >= 0.8 pu):
# far out, both substations can "import" on the low root, which is not an operating point.
#
# Zs cases: stiff (Zs = 0), and |Zs| = 0.5, 0.8 (the decks' value), 2, 4 ohm at X/R = 25.
# 0.8 ohm is U^2/S_k for the typical 20 kV short-circuit level 0.50 GVA (Traupmann &
# Kienberger, Energies 2020, Table 19); 0.5-4 ohm spans their 20 kV range of 0.10-0.80 GVA
# (Table 18). The deck's own Zs is also cross-checked against OpenDSS.
#
# Run from the repo root:
#     julia --project=envs/tadmm envs/multi_poi/root_level/source_impedance_study.jl
#
# Writes (tracked) envs/multi_poi/results/source_impedance/:
#     windows.csv        every system x Zs x window
#     table_windows.tex  the same as a LaTeX table, \input by the MSOPF paper

include(joinpath(@__DIR__, "..", "full_angle_pf.jl"))

const SYSTEMS = [
    (key = "ieee123", system = "ieee123_5poi_1ph", pair = ("subs3", "subs4"),
     setup = (; disable = ["PVSystem"], snapshot = true), test_battery_bus = nothing),
    (key = "small2poi", system = "small2poi_1ph", pair = ("grid1", "grid2"),
     setup = (; disable = String[], snapshot = false), test_battery_bus = "1load")]
const ZS_OHM = [0.0, 0.5, 0.8, 2.0, 4.0]    # 0 = stiff (ideal_sources)
const X_OVER_R = 25.0
const BAND = (0.95, 1.05)
const LIMIT_DEG = 20.0                       # search limit; a window reaching it is reported as ">="
const STEP_DEG = 0.1                         # outward scan step before bisecting the last one
const BISECT = 20                            # ~1e-7 deg resolution within that step
const V_GUARD = 0.8                          # low-root guard when the band is off
const OUT_DIR = normpath(joinpath(@__DIR__, "..", "results", "source_impedance"))

"The network with every source's series impedance set to |Zs| (ohm) at X/R = X_OVER_R."
function with_zs(net, zmag)
    r = zmag / sqrt(1 + X_OVER_R^2)
    return merge(net, (; sources = [merge(s, (; z_series = complex(r, X_OVER_R * r))) for s in net.sources]))
end

"Compile a system's deck (keeping its two sources, adding the test battery if configured)."
function load_system(cfg)
    A, B = cfg.pair
    kw = (; sources = [A, B], cfg.setup..., load_band = (0.5, 1.5))
    extra = String[]
    if cfg.test_battery_bus !== nothing                  # small2poi: one battery at 100% of the load
        compile_deck(cfg.system; kw...)
        probe = read_network()
        r = sum(d.kW for d in probe.loads)
        push!(extra, @sprintf("New Storage.battery1 phases=1 bus1=%s kV=%.6g kWrated=%.10g kVA=%.10g kWhrated=%.10g %%stored=50 %%reserve=0 %%IdlingkW=0 DispMode=External",
                              cfg.test_battery_bus, probe.kV_base, r, r, 4r))
    end
    compile_deck(cfg.system; extra, kw...)
    return read_network()
end

banded_buses(net) = [b for b in net.buses if !(b in Set(s.bus for s in net.sources))]
in_band(net, sol) = all(b -> BAND[1] <= abs(sol.V_pu[b]) <= BAND[2], banded_buses(net))
"On the high-voltage power-flow root: the same V_GUARD the band-free OPF imposes."
high_root(net, sol) = all(b -> abs(sol.V_pu[b]) >= V_GUARD, banded_buses(net))

"""
Edge of the feasible interval containing 0, in direction `sgn` (+1 or -1): step outward from
0 until the first infeasible point, then bisect that last step. Stepping, rather than testing
the search limit first, matters: far out, a solver can land on the low-voltage root where both
substations import heavily -- a spurious "feasible" point disconnected from the real window.
"""
function edge(feasible, sgn)
    inside = 0.0
    while abs(inside) < LIMIT_DEG
        trial = sgn * min(abs(inside) + STEP_DEG, LIMIT_DEG)
        feasible(trial) || return (bisect(feasible, inside, trial), false)
        inside = trial
    end
    return (inside, true)                                # still feasible at the search limit
end

function bisect(feasible, inside, outside)
    for _ in 1:BISECT
        mid = (inside + outside) / 2
        feasible(mid) ? (inside = mid) : (outside = mid)
    end
    return inside
end

function window(feasible)
    feasible(0.0) || return (lo = NaN, hi = NaN, lo_capped = false, hi_capped = false, feasible = false)
    (lo, lc), (hi, hc) = edge(feasible, -1), edge(feasible, +1)
    return (; lo, hi, lo_capped = lc, hi_capped = hc, feasible = true)
end

rows = []
for cfg in SYSTEMS
    A, B = cfg.pair
    deck = load_system(cfg)
    zdeck = [s.z_series for s in deck.sources]
    println("\n", "="^96)
    @printf "%s -- %s + %s; deck Zs = %s ohm\n" cfg.system A B join(string.(round.(zdeck; sigdigits = 4)), ", ")

    # The deck's own Zs against OpenDSS
    worst_P, worst_V = 0.0, 0.0
    for dd in (-1.0, 0.0, 0.5, 2.0)
        d = Dict(A => dd, B => 0.0)
        ip, ds = solve_full_angle(deck, d), opendss_point(deck, d)
        (ip.ok && ds.converged) || continue
        worst_P = max(worst_P, maximum(abs(ip.P_subs_kW[k] - ds.P_subs_kW[k]) for k in (A, B)))
        worst_V = max(worst_V, maximum(abs(ip.V_pu[b] - ds.V_pu[b]) for b in deck.buses))
    end
    @printf "  deck Zs vs OpenDSS: worst |dP_subs| = %.2e kW, worst |dV| = %.2e pu\n" worst_P worst_V

    for zmag in ZS_OHM
        stiff = zmag == 0
        net = stiff ? deck : with_zs(deck, zmag)
        δ(dd) = Dict(A => dd, B => 0.0)
        pf(dd) = solve_full_angle(net, δ(dd); ideal_sources = stiff)
        pf_ok(dd, band) = (r = pf(dd); r.ok && high_root(net, r) && minimum(values(r.P_subs_kW)) >= 0 &&
                                         (!band || in_band(net, r)))
        opf_ok(dd, band) = solve_full_angle(net, δ(dd); ideal_sources = stiff, battery = :optimize, no_backflow = true,
                                            (band ? (; V_band = BAND) : (; V_guard = V_GUARD))...).ok
        for (case, f) in (("PF", dd -> pf_ok(dd, false)), ("PF+band", dd -> pf_ok(dd, true)),
                          ("OPF", dd -> opf_ok(dd, false)), ("OPF+band", dd -> opf_ok(dd, true)))
            w = window(f)
            push!(rows, (; system = cfg.key, zs_ohm = zmag, case, w...))
            txt = w.feasible ? @sprintf("[%s%+.3f, %s%+.3f]", w.lo_capped ? "<=" : "", w.lo, w.hi_capped ? ">=" : "", w.hi) : "infeasible at equal angles"
            @printf "  |Zs| = %-4s ohm  %-9s %s\n" (stiff ? "0" : string(zmag)) case txt
        end
    end
end

# ---- outputs -------------------------------------------------------------------------------------------
mkpath(OUT_DIR)
open(joinpath(OUT_DIR, "windows.csv"), "w") do io
    println(io, "system,zs_ohm,x_over_r,case,feasible_at_equal_angles,lo_deg,hi_deg,lo_at_search_limit,hi_at_search_limit")
    for r in rows
        @printf io "%s,%g,%g,%s,%s,%.6f,%.6f,%s,%s\n" r.system r.zs_ohm (r.zs_ohm == 0 ? NaN : X_OVER_R) r.case r.feasible r.lo r.hi r.lo_capped r.hi_capped
    end
end

cell(r) = !r.feasible ? "infeasible" :
    @sprintf("\$[%s%.2f, %s%+.2f]\$", r.lo_capped ? "\\le " : "", r.lo, r.hi_capped ? "\\ge " : "", r.hi)
open(joinpath(OUT_DIR, "table_windows.tex"), "w") do io
    println(io, "% Generated by envs/multi_poi/root_level/source_impedance_study.jl -- do not edit by hand.")
    println(io, "\\begin{table}[!t]\n\\centering\n\\caption{No-backflow window of \$\\delta_{s_a}-\\delta_{s_b}\$ (deg) versus source impedance \$|Z_s|\$ (\$X/R=25\$). ``band'': 0.95--1.05\\,pu imposed.}")
    println(io, "\\label{tab:zs-windows}\n\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{llllll}\n\\hline")
    println(io, "System & \$|Z_s|\$ (\$\\Omega\$) & PF & PF+band & Batteries & Batt.+band \\\\\n\\hline")
    for cfg in SYSTEMS, zmag in ZS_OHM
        get(c) = only(filter(r -> r.system == cfg.key && r.zs_ohm == zmag && r.case == c, rows))
        label = zmag == ZS_OHM[1] ? cfg.key : ""
        z = zmag == 0 ? "0 (stiff)" : (zmag == 0.8 ? "0.8 (deck)" : string(zmag))
        @printf io "%s & %s & %s & %s & %s & %s \\\\\n" label z cell(get("PF")) cell(get("PF+band")) cell(get("OPF")) cell(get("OPF+band"))
    end
    println(io, "\\hline\n\\end{tabular}\n\\end{table}")
end
println("\nWrote ", joinpath(OUT_DIR, "windows.csv"), " and table_windows.tex")
