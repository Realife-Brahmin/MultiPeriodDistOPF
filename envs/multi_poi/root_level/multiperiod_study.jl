# multiperiod_study.jl -- what the full-angle model says once batteries have to obey their
# state of charge over a whole day.
#
# The substation angles are held CONSTANT over the horizon. An angle configuration is feasible
# if every period satisfies the AC equations, P_Subs >= 0 at every substation and the
# 0.95-1.05 pu band, and every battery stays inside 30-95% SOC from 62.5%. For two
# substations a and b (delta_b = 0) that gives a window of delta_a - delta_b, the interval
# connected to equal angles, under four assumptions:
#   nobatt    batteries idle (PV reactive power, where PV is on, is still dispatched);
#   hourly    batteries, each hour solved alone with unlimited energy: the intersection of
#             the hourly windows -- what a snapshot study reports;
#   soft      batteries with SOC dynamics and the soft terminal SOC -- the user's standing
#             rule, and THE multi-period result reported in tables and figures;
#   cyclic    the same with a hard B[T] == B[0]; a diagnostic kept in the CSVs only.
# Containment nobatt <= cyclic <= soft <= hourly is checked and printed.
#
# PV modes: "off"; "unity" (real power only); "qfree" (each unit's reactive power a decision
# within its inverter rating -- the DDP/tADMM convention, and the default here).
#
# Stages (pass any on the command line; default: all), each writing its own CSV:
#   main      ieee123 subs3+subs4, T = 24, all three PV modes: hourly and horizon windows,
#             the cost-optimal day at two angles (replayed in OpenDSS), figures;
#   horizon   the horizon windows for T = 6, 12, 24, 48, 96;
#   zs        ... versus source impedance (T = 24);
#   duration  ... versus battery duration, 1-8 h at the same kW (T = 24);
#   pairs     ... for all ten substation pairs (T = 24);
#   three     three substations (subs1 the angle reference): the feasible REGION in the
#             plane (delta_3 - delta_1, delta_4 - delta_1), by rays from equal angles;
#   five      all five substations: one substation's angle moved against the other four;
#   bound     what limits the battery window: the circulating power per degree, the angle at
#             which a substation's daily energy import reaches zero, and the windows against
#             battery energy and power;
#   tables    rebuild the LaTeX tables from the CSVs (no solves).
#
# Run from the repo root:
#     julia --project=envs/tadmm envs/multi_poi/root_level/multiperiod_study.jl [stage ...]
#
# Writes (tracked) envs/multi_poi/results/multiperiod/.

include(joinpath(@__DIR__, "..", "multiperiod_full_angle.jl"))
include(joinpath(@__DIR__, "..", "sweep_common.jl"))
using Logging

const SYSTEM = "ieee123_5poi_1ph"
const MAIN_PAIR = ("subs3", "subs4")
const ALL_SUBS = ["subs1", "subs2", "subs3", "subs4", "subs5"]
const X_OVER_R = 25.0
const OUT_DIR = normpath(joinpath(@__DIR__, "..", "results", "multiperiod"))
const STAGES = isempty(ARGS) ? ["main", "horizon", "zs", "duration", "pairs", "three", "five", "bound", "tables"] : lowercase.(ARGS)
mkpath(OUT_DIR)
say(args...) = (println(args...); flush(stdout))

# ---- setup -------------------------------------------------------------------------------------------
function load_sources(srcs)
    kw = (; sources = collect(srcs), snapshot = true, load_band = (0.5, 1.5))
    pvs = with_logger(NullLogger()) do                  # "? PVSystem.x.Pmpp" logs a warning per unit
        read_pv(SYSTEM; kw...)
    end
    compile_deck(SYSTEM; disable = ["PVSystem"], kw...)
    return read_network(), pvs
end

with_zs(net, zmag) = (r = zmag / sqrt(1 + X_OVER_R^2);
    merge(net, (; sources = [merge(s, (; z_series = complex(r, X_OVER_R * r))) for s in net.sources])))
with_duration(net, hours) = merge(net, (; batteries = [merge(bt, (; kWh_rated = hours * bt.kW_rated)) for bt in net.batteries]))
banded(net) = [b for b in net.buses if !(b in Set(s.bus for s in net.sources))]

"Everything one configuration needs: profiles, the PV mode's inputs, and the cost weights."
function setup(net, pvs; T, pv = "qfree", ideal = false)
    lam, c, pvf = load_profile(T), price_profile(T), pv_profile(T)
    P, F = pv == "off" ? ([], zeros(T)) : (pvs, pvf)
    return (; net, T, lam, c, P, F, q = pv == "unity" ? :unity : :optimize, CB = cb_by_rule(net), GAM = gamma_by_rule(net), ideal)
end

const VTOL = 1e-4              # per unit: total limit violation below this counts as feasible
const LIMIT = 20.0             # degrees: how far a ray is followed
const RAY_OK = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)

"""
    problem(s, case; hour = nothing)

Keyword arguments of solve_multiperiod for one case of setup `s`: over the whole horizon
(:nobatt, :soft, :cyclic), or for hour `hour` alone with unlimited battery energy (:nobatt,
:batt).
"""
function problem(s, case; hour = nothing)
    common = (; pvs = s.P, pv_q = s.q, C_B = s.CB, ideal_sources = s.ideal, battery = case === :nobatt ? :idle : :optimize)
    hour === nothing && return (; common..., T = s.T, load_factor = s.lam, price = s.c, pv_factor = s.F, gamma = s.GAM,
                                terminal = case === :cyclic ? :cyclic : :soft)
    return (; common..., T = 1, dt_h = 24 / s.T, load_factor = [s.lam[hour]], price = [s.c[hour]], pv_factor = [s.F[hour]],
            soc = (-Inf, Inf), terminal = :free)
end
solve_case(s, d, case; kwargs...) = solve_multiperiod(s.net, d; problem(s, case)..., kwargs...)
equal_angles(net) = Dict(s.name => 0.0 for s in net.sources)

"""
    rays(s, case, dirs; hour = nothing, check = :inside)

How far the angles can move from equal along each direction in `dirs` (source name =>
degrees per unit of r): one feasibility solve at equal angles, then one maximize-r solve per
direction. `check = :inside` confirms each edge with the violation test at 0.99 r; `:both`
also requires a violation at r + 0.02 deg; failures are printed and counted in `bad`.
"""
function rays(s, case, dirs; hour = nothing, check = :inside)
    kw, net = problem(s, case; hour), s.net
    # A violation solve counts if Ipopt converged, including to its "acceptable" tolerance
    # (ALMOST_LOCALLY_SOLVED): feasibility is then read off the violation itself. Requiring the
    # strict status flagged one correct edge (three substations, ray 202.5 deg, 2026-10-11).
    solved(x) = x.status in RAY_OK && isfinite(x.violation)
    at0 = solve_multiperiod(net, equal_angles(net); kw..., objective = :violation)
    feasible = solved(at0) && at0.violation <= VTOL
    r, capped, bad = fill(NaN, length(dirs)), fill(false, length(dirs)), 0
    feasible || return (; feasible, r, capped, bad)
    tag = hour === nothing ? string(case) : "$case, hour $hour"
    for (i, dir) in enumerate(dirs)
        x = solve_multiperiod(net, equal_angles(net); kw..., objective = :ray, ray_dir = dir, ray_limit = LIMIT, start = at0)
        if !(x.status in RAY_OK) || !isfinite(x.r)
            bad += 1
            say("  ** ray solve failed (", tag, "): ", x.status)
            continue
        end
        r[i], capped[i] = x.r, x.r >= LIMIT - 1e-6
        (check === :none || x.r <= 1e-6) && continue
        inn = solve_multiperiod(net, Dict(k => 0.99 * x.r * v for (k, v) in dir); kw..., objective = :violation, start = x)
        good = solved(inn) && inn.violation <= VTOL
        if check === :both && !capped[i]
            out = solve_multiperiod(net, Dict(k => (x.r + 0.02) * v for (k, v) in dir); kw..., objective = :violation, start = x)
            good &= solved(out) && out.violation > VTOL
        end
        good || (bad += 1; say(@sprintf("  ** edge check failed (%s): r = %.4f deg", tag, x.r)))
    end
    return (; feasible, r, capped, bad)
end

"Window of delta_a, every other source at 0: the two rays along +a and -a."
function window(s, case, a; kwargs...)
    dir(sgn) = Dict(x.name => (x.name == a ? sgn : 0.0) for x in s.net.sources)
    y = rays(s, case, [dir(+1.0), dir(-1.0)]; kwargs...)
    return (; y.feasible, lo = -y.r[2], hi = y.r[1], lo_capped = y.capped[2], hi_capped = y.capped[1], y.bad)
end

function intersect_windows(ws)
    all(w -> w.feasible && isfinite(w.lo) && isfinite(w.hi), ws) || return (; feasible = false, lo = NaN, hi = NaN, lo_capped = false, hi_capped = false, bad = 0)
    return (; feasible = true, lo = maximum(w.lo for w in ws), hi = minimum(w.hi for w in ws),
            lo_capped = all(w -> w.lo_capped, ws), hi_capped = all(w -> w.hi_capped, ws), bad = sum(w.bad for w in ws))
end

"The horizon windows of delta_a for setup `s`; with `hourly = true` also every hour's windows."
function windows(s, a; hourly = false, check = :inside)
    w = Dict{Symbol,Any}(case => window(s, case, a; check) for case in (:nobatt, :soft, :cyclic))
    hw = nothing
    if hourly
        hw = [(; t, nobatt = window(s, :nobatt, a; hour = t), batt = window(s, :batt, a; hour = t)) for t in 1:s.T]
        w[:hourly] = intersect_windows([h.batt for h in hw])
        w[:nobatt_hourly] = intersect_windows([h.nobatt for h in hw])   # must equal :nobatt
    end
    return w, hw
end

fmt(w) = w === nothing ? "(not run)" : (w.feasible ? @sprintf("[%s%+.3f, %s%+.3f]", w.lo_capped ? "<=" : "", w.lo, w.hi_capped ? ">=" : "", w.hi) : "infeasible at equal angles")
const CASES = (:nobatt, :hourly, :soft, :cyclic)

function report(tag, w)
    g(c) = get(w, c, nothing)
    say(@sprintf("  %-36s nobatt %s | hourly %s | soft %s | cyclic %s", tag, fmt(g(:nobatt)), fmt(g(:hourly)), fmt(g(:soft)), fmt(g(:cyclic))))
    tol = 2e-3
    inside(a, b) = !a.feasible || (b.feasible && a.lo >= b.lo - tol && a.hi <= b.hi + tol)
    ok = inside(w[:nobatt], w[:cyclic]) && inside(w[:cyclic], w[:soft]) && (g(:hourly) === nothing || inside(w[:soft], w[:hourly]))
    ok || say("  ** containment nobatt <= cyclic <= soft <= hourly FAILS for ", tag)
    if haskey(w, :nobatt_hourly)
        a, b = w[:nobatt], w[:nobatt_hourly]
        same = a.feasible == b.feasible && (!a.feasible || (abs(a.lo - b.lo) < tol && abs(a.hi - b.hi) < tol))
        same || say("  ** nobatt as one T-period solve ", fmt(a), " differs from the intersection of hourly solves ", fmt(b))
    end
    nbad = sum(x.bad for x in values(w))
    nbad == 0 || say("  ** ", nbad, " failed ray solves or edge checks for ", tag)
end

"CSV rows (header: ...,case,feasible_at_equal_angles,lo_deg,hi_deg,lo_at_limit,hi_at_limit,failed_checks)."
function csv_rows(io, key, w)
    for case in CASES
        x = get(w, case, nothing)
        x === nothing && continue
        @printf io "%s,%s,%s,%.6f,%.6f,%s,%s,%d\n" key case x.feasible x.lo x.hi x.lo_capped x.hi_capped x.bad
    end
    flush(io)
end
const CSV_TAIL = ",case,feasible_at_equal_angles,lo_deg,hi_deg,lo_at_limit,hi_at_limit,failed_checks"

# ---- stage: main ---------------------------------------------------------------------------------------
if "main" in STAGES
    A, B = MAIN_PAIR
    net, pvs = load_sources([A, B])
    T = 24
    say("\n", "="^110)
    say(@sprintf("MAIN  %s, %s + %s, T = %d (dt = 1 h)", SYSTEM, A, B, T))
    say(@sprintf("  load %.0f kW; %d batteries %.0f kW / %.0f kWh; %d PV units %.0f kW / %.0f kVA; Zs = %s ohm", sum(d.kW for d in net.loads), length(net.batteries), sum(b.kW_rated for b in net.batteries), sum(b.kWh_rated for b in net.batteries), length(pvs), sum(p.kW for p in pvs), sum(p.kVA for p in pvs), string(round(net.sources[1].z_series; sigdigits = 4))))
    say(@sprintf("  C_B = %.3e \$/kW^2h, gamma = %.3e (the DDP sizing rules on this fleet)", cb_by_rule(net), gamma_by_rule(net)))
    MODES = ("off", "unity", "qfree")
    S = Dict(m => setup(net, pvs; T, pv = m) for m in MODES)
    res = Dict{String,Any}()
    for m in MODES                                     # every horizon edge checked inside AND outside
        res[m] = windows(S[m], A; hourly = true, check = :both)
        report("PV " * m, res[m][1])
    end

    open(joinpath(OUT_DIR, "main_hourly_windows.csv"), "w") do io
        println(io, "pv,hour,load_factor,pv_factor,price,nobatt_feasible,nobatt_lo_deg,nobatt_hi_deg,batt_feasible,batt_lo_deg,batt_hi_deg")
        for m in MODES, h in res[m][2]
            s = S[m]
            @printf io "%s,%d,%.6f,%.6f,%.6f,%s,%.6f,%.6f,%s,%.6f,%.6f\n" m h.t s.lam[h.t] s.F[h.t] s.c[h.t] h.nobatt.feasible h.nobatt.lo h.nobatt.hi h.batt.feasible h.batt.lo h.batt.hi
        end
    end
    open(joinpath(OUT_DIR, "main_horizon_windows.csv"), "w") do io
        println(io, "pv", CSV_TAIL)
        for m in MODES
            csv_rows(io, m, res[m][1])
        end
    end

    # The cost-optimal day (soft terminal SOC, PV reactive power free) at equal angles and at an
    # angle only batteries make feasible: halfway between the no-battery edge and the soft edge.
    s, w = S["qfree"], res["qfree"][1]
    d_star = round((w[:nobatt].hi + w[:soft].hi) / 2; digits = 2)
    days = []
    for dd in (0.0, d_star)
        d = Dict(A => dd, B => 0.0)
        x = solve_case(s, d, :soft)
        rp = replay_opendss(net, x, d; load_factor = s.lam, pvs = s.P, pv_factor = s.F)
        say(@sprintf("  day at %+.2f deg: %s, cost \$%.2f (energy %.2f, cycling %.2f); OpenDSS replay |dP| %.1e kW, |dQ| %.1e kvar, |dV| %.1e pu, converged %s, PQ %s", dd, x.status, x.cost.total, x.cost.energy, x.cost.cycling, rp.dP, rp.dQ, rp.dV, rp.converged, rp.pq_band))
        push!(days, (; dd, x, rp))
    end
    PL, PVk = sum(d.kW for d in net.loads), sum(p.kW for p in pvs)
    nb = banded(net)
    open(joinpath(OUT_DIR, "main_days.csv"), "w") do io
        println(io, "ddelta_deg,hour,load_kW,pv_kW,price,P_subs_a_kW,P_subs_b_kW,fleet_P_B_kW,fleet_energy_kWh,fleet_PV_Q_kvar,min_V_pu,max_V_pu")
        for y in days, t in 1:T
            x = y.x
            @printf io "%.2f,%d,%.3f,%.3f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.6f,%.6f\n" y.dd t s.lam[t]*PL s.F[t]*PVk s.c[t] x.P_subs_kW[A][t] x.P_subs_kW[B][t] sum(v[t] for v in values(x.P_B_kW)) sum(v[t+1] for v in values(x.B_kWh)) sum(v[t] for v in values(x.Q_PV_kvar)) minimum(abs(x.V_pu[b][t]) for b in nb) maximum(abs(x.V_pu[b][t]) for b in nb)
        end
    end
    open(joinpath(OUT_DIR, "main_opendss_check.csv"), "w") do io
        println(io, "ddelta_deg,status,opendss_converged,pq_band,max_dP_kW,max_dQ_kvar,max_dV_pu,max_dP_B_kW,cost_total,cost_energy,cost_cycling")
        for y in days
            @printf io "%.2f,%s,%s,%s,%.3e,%.3e,%.3e,%.3e,%.4f,%.4f,%.4f\n" y.dd y.x.status y.rp.converged y.rp.pq_band y.rp.dP y.rp.dQ y.rp.dV y.rp.dB y.x.cost.total y.x.cost.energy y.x.cost.cycling
        end
    end

    # ---- figures (kW and kWh, never per unit) ----
    plot_style!()
    hrs = 1:T
    BATT, SOFT, CYC = colorant"#eb6834", colorant"#7a4fd0", colorant"#1c8a5a"
    p1 = plot(hrs, s.lam .* PL; label = "Load", color = SUB_COLORS[1], lw = 2, ylabel = "kW", xlabel = "Hour", title = "Load and PV", legend = :bottomright, ylims = (0, 1.1PL))
    plot!(p1, hrs, s.F .* PVk; label = "PV", color = BATT, lw = 2)
    p2 = plot(hrs, s.c; label = "", color = INK2, lw = 2, ylabel = "\$/kWh", xlabel = "Hour", title = "Energy price (both substations)", ylims = (0, 0.22))
    savefig(plot(p1, p2; layout = (1, 2), size = (900, 300), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), joinpath(OUT_DIR, "fig_inputs.png"))

    panels = []
    for (m, ttl) in (("off", "PV off"), ("qfree", "PV on, reactive power dispatched"))
        r, hw = res[m]
        p = plot(; xlabel = "Hour", ylabel = "$(dlabel(A)) − $(dlabel(B)) (deg)", title = ttl, legend = m == "qfree" ? :outerright : false, xlims = (1, T))
        plot!(p, hrs, [h.batt.hi for h in hw]; fillrange = [h.batt.lo for h in hw], fillalpha = 0.18, color = BATT, lw = 1.2, label = "Batteries, each hour alone")
        plot!(p, hrs, [h.batt.lo for h in hw]; color = BATT, lw = 1.2, label = "")
        plot!(p, hrs, [h.nobatt.hi for h in hw]; fillrange = [h.nobatt.lo for h in hw], fillalpha = 0.30, color = SUB_COLORS[1], lw = 1.2, label = "No batteries")
        plot!(p, hrs, [h.nobatt.lo for h in hw]; color = SUB_COLORS[1], lw = 1.2, label = "")
        r[:soft].feasible && hline!(p, [r[:soft].lo, r[:soft].hi]; color = SOFT, lw = 2, label = "Batteries, whole day")
        hline!(p, [0]; color = MUTED, lw = 1, label = "")
        push!(panels, p)
    end
    savefig(plot(panels...; layout = grid(1, 2, widths = [0.40, 0.60]), size = (1150, 380), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm), joinpath(OUT_DIR, "fig_windows.png"))

    ps = []
    ER = sum(b.kWh_rated for b in net.batteries)
    for y in days
        x = y.x
        pa = plot(hrs, x.P_subs_kW[A]; label = label(A), color = SUB_COLORS[1], lw = 2, ylabel = "kW", title = @sprintf("%s − %s = %+.2f°: substation imports", dlabel(A), dlabel(B), y.dd), legend = :topright)
        plot!(pa, hrs, x.P_subs_kW[B]; label = label(B), color = SUB_COLORS[2], lw = 2)
        hline!(pa, [0]; color = MUTED, lw = 1, label = "")
        pb = plot(hrs, [sum(v[t] for v in values(x.P_B_kW)) for t in 1:T]; label = "", color = BATT, lw = 2, ylabel = "kW", title = "Battery fleet power (+ discharging)", seriestype = :steppre)
        hline!(pb, [0]; color = MUTED, lw = 1, label = "")
        pc = plot(0:T, [sum(v[t+1] for v in values(x.B_kWh)) for t in 0:T]; label = "", color = CYC, lw = 2, ylabel = "kWh", xlabel = "Hour", title = "Battery fleet energy (dashed: 30% and 95%)")
        hline!(pc, [0.30ER, 0.95ER]; color = MUTED, lw = 1, ls = :dash, label = "")
        push!(ps, pa, pb, pc)
    end
    savefig(plot(ps[1], ps[4], ps[2], ps[5], ps[3], ps[6]; layout = (3, 2), size = (1100, 760), left_margin = 6Plots.mm, bottom_margin = 4Plots.mm), joinpath(OUT_DIR, "fig_days.png"))
    say("  wrote main_*.csv, fig_inputs.png, fig_windows.png, fig_days.png")
end

# ---- stages: horizon, zs, duration ---------------------------------------------------------------------
"One sweep stage: `configs` is (key, label, setup) triples; windows of delta_a per config."
function sweep_stage(title, file, header, configs, a)
    say("\n", "="^110, "\n", title)
    open(joinpath(OUT_DIR, file), "w") do io
        println(io, header, CSV_TAIL)
        for (key, lab, s) in configs
            w, _ = windows(s, a)
            report(lab, w)
            csv_rows(io, key, w)
        end
    end
end

if "horizon" in STAGES
    A, B = MAIN_PAIR
    net, pvs = load_sources([A, B])
    sweep_stage("HORIZON  $A + $B, PV reactive power free", "horizon_windows.csv", "T",
                [(string(T), "T = $T (dt = $(24 / T) h)", setup(net, pvs; T)) for T in (6, 12, 24, 48, 96)], A)
end
if "zs" in STAGES
    A, B = MAIN_PAIR
    deck, pvs = load_sources([A, B])
    sweep_stage("ZS  $A + $B, T = 24", "zs_windows.csv", "zs_ohm",
                [(string(z), "|Zs| = $z ohm", setup(z == 0 ? deck : with_zs(deck, z), pvs; T = 24, ideal = z == 0)) for z in (0.0, 0.5, 0.8, 2.0, 4.0)], A)
end
if "duration" in STAGES
    A, B = MAIN_PAIR
    deck, pvs = load_sources([A, B])
    sweep_stage("DURATION  $A + $B, T = 24 (battery kWh = hours x kW)", "duration_windows.csv", "battery_hours",
                [(string(h), "$h h batteries", setup(with_duration(deck, h), pvs; T = 24)) for h in (1.0, 2.0, 4.0, 8.0)], A)
end

# ---- stage: pairs ------------------------------------------------------------------------------------------
if "pairs" in STAGES
    say("\n", "="^110, "\nPAIRS  T = 24, every pair of substations")
    open(joinpath(OUT_DIR, "pairs_windows.csv"), "w") do io
        println(io, "sub_a,sub_b,separation_ohm", CSV_TAIL)
        for i in 1:length(ALL_SUBS), j in i+1:length(ALL_SUBS)
            A, B = ALL_SUBS[i], ALL_SUBS[j]
            net, pvs = load_sources([A, B])
            sep = abs(only(values(source_separation(net))))
            w, _ = windows(setup(net, pvs; T = 24), A)
            report(@sprintf("%s + %s (|Z| between %.3f ohm)", A, B, sep), w)
            csv_rows(io, @sprintf("%s,%s,%.4f", A, B, sep), w)
        end
    end
end

# ---- stage: three ------------------------------------------------------------------------------------------
# Three substations: subs1 is the angle reference, and the "window" becomes a region in the
# plane (delta_3 - delta_1, delta_4 - delta_1). Mapped by rays from equal angles.
if "three" in STAGES
    REF, A, B = "subs1", "subs3", "subs4"
    net, pvs = load_sources([REF, A, B])
    s = setup(net, pvs; T = 24)
    NDIR = 16
    say("\n", "="^110, "\nTHREE  $REF (reference) + $A + $B, T = 24: feasible region by $NDIR rays")
    ths = [2pi * k / NDIR for k in 0:NDIR-1]
    dirs = [Dict(REF => 0.0, A => cos(th), B => sin(th)) for th in ths]
    reg = Dict{Symbol,Any}()
    for case in (:nobatt, :cyclic, :soft)
        reg[case] = rays(s, case, dirs)
        say(@sprintf("  %-7s feasible at equal angles %s; radius %.3f to %.3f deg over %d rays; %d failed checks", case, reg[case].feasible, minimum(reg[case].r), maximum(reg[case].r), NDIR, reg[case].bad))
    end
    open(joinpath(OUT_DIR, "three_region.csv"), "w") do io
        println(io, "ray_deg,case,feasible_at_equal_angles,radius_deg,d3_minus_d1_deg,d4_minus_d1_deg,at_limit")
        for (i, th) in enumerate(ths), case in (:nobatt, :cyclic, :soft)
            y = reg[case]
            @printf io "%.2f,%s,%s,%.6f,%.6f,%.6f,%s\n" rad2deg(th) case y.feasible y.r[i] y.r[i]*cos(th) y.r[i]*sin(th) y.capped[i]
        end
    end
    plot_style!()
    p = plot(; xlabel = "$(dlabel(A)) − $(dlabel(REF)) (deg)", ylabel = "$(dlabel(B)) − $(dlabel(REF)) (deg)", aspect_ratio = :equal, legend = :outerright,
             title = "Feasible angle region, three substations (T = 24)", size = (820, 560))
    for (case, col, lab, alpha) in ((:soft, colorant"#7a4fd0", "Batteries, whole day", 0.18), (:nobatt, SUB_COLORS[1], "No batteries", 0.35))
        reg[case].feasible || continue
        xs = reg[case].r .* cos.(ths); ys = reg[case].r .* sin.(ths)
        plot!(p, [xs; xs[1]], [ys; ys[1]]; seriestype = :shape, fillalpha = alpha, color = col, lw = 1.5, linecolor = col, label = lab)
    end
    scatter!(p, [0], [0]; color = INK, ms = 3, label = "Equal angles")
    savefig(p, joinpath(OUT_DIR, "fig_three_region.png"))
    say("  wrote three_region.csv, fig_three_region.png")
end

# ---- stage: five -------------------------------------------------------------------------------------------
# All five substations at once. The angle space is four-dimensional; reported here is one
# line through it per substation: that substation's angle moved against the other four.
if "five" in STAGES
    net, pvs = load_sources(ALL_SUBS)
    s = setup(net, pvs; T = 24)
    say("\n", "="^110, "\nFIVE  all five substations, T = 24: one substation's angle against the other four")
    open(joinpath(OUT_DIR, "five_windows.csv"), "w") do io
        println(io, "moved_substation", CSV_TAIL)
        for a in ALL_SUBS
            w, _ = windows(s, a)
            report("moving " * a, w)
            csv_rows(io, a, w)
        end
    end
end

# ---- stage: bound ------------------------------------------------------------------------------------------
# What limits the window batteries can sustain? A constant angle difference drives a roughly
# constant circulating power, and a battery can only move energy between hours. So the most
# a fleet can do is trade the worst hour's import margin for the day's average: the window
# should end where the constrained substation's DAILY ENERGY import reaches zero, and storage
# beyond what flattens the day should add nothing. Tested here:
#   (a) batteries idle, no-backflow NOT imposed: each substation's hourly import and daily
#       energy against the angle difference -- the transfer per degree and both zero crossings;
#   (b) the windows against battery energy (hours at the same kW) and power (kW at the same
#       kWh), for the soft terminal SOC and for the hard end = start diagnostic.
if "bound" in STAGES
    A, B = MAIN_PAIR
    deck, pvs = load_sources([A, B])
    s = setup(deck, pvs; T = 24)
    say("\n", "="^110, "\nBOUND  $A + $B, T = 24: circulating power, and the window against battery size")
    dds = -1.2:0.1:1.2
    sweep = []
    for dd in dds
        x = solve_multiperiod(deck, Dict(A => dd, B => 0.0); problem(s, :nobatt)..., no_backflow = false)
        x.ok || (say("  ** idle sweep failed at ", dd, ": ", x.status); continue)
        push!(sweep, (; dd, Pa = x.P_subs_kW[A], Pb = x.P_subs_kW[B]))
    end
    dt = 24 / s.T
    Ea(y) = sum(y.Pa) * dt                               # kWh over the day
    Eb(y) = sum(y.Pb) * dt
    "x where the line through the two points bracketing a sign change of f crosses zero"
    function crossing(f)
        for i in 1:length(sweep)-1
            a, b = f(sweep[i]), f(sweep[i+1])
            a * b <= 0 && a != b && return sweep[i].dd + (sweep[i+1].dd - sweep[i].dd) * a / (a - b)
        end
        return NaN
    end
    i0 = findfirst(y -> abs(y.dd) < 1e-9, sweep)
    K = [(sweep[i0+1].Pa[t] - sweep[i0-1].Pa[t]) / (sweep[i0+1].dd - sweep[i0-1].dd) for t in 1:s.T]   # kW per degree, into A
    hour_hi, hour_lo = crossing(y -> minimum(y.Pb)), crossing(y -> minimum(y.Pa))     # first hour to back-feed
    day_hi, day_lo = crossing(Eb), crossing(Ea)                                       # daily energy import to zero
    say(@sprintf("  transfer per degree: %.0f to %.0f kW/deg over the day (mean %.0f)", minimum(K), maximum(K), sum(K) / length(K)))
    say(@sprintf("  imports at equal angles: %s %.0f-%.0f kW (%.0f kWh/day), %s %.0f-%.0f kW (%.0f kWh/day)", A, extrema(sweep[i0].Pa)..., Ea(sweep[i0]), B, extrema(sweep[i0].Pb)..., Eb(sweep[i0])))
    say(@sprintf("  first hour back-feeds at   [%+.3f, %+.3f] deg   (the no-battery window)", hour_lo, hour_hi))
    say(@sprintf("  daily energy import zero at [%+.3f, %+.3f] deg   (predicted limit for batteries)", day_lo, day_hi))
    open(joinpath(OUT_DIR, "bound_sweep.csv"), "w") do io
        println(io, "ddelta_deg,daily_import_a_kWh,daily_import_b_kWh,min_hourly_import_a_kW,min_hourly_import_b_kW")
        for y in sweep
            @printf io "%.2f,%.3f,%.3f,%.4f,%.4f\n" y.dd Ea(y) Eb(y) minimum(y.Pa) minimum(y.Pb)
        end
    end

    with_power(net, k) = merge(net, (; batteries = [merge(bt, (; kW_rated = k * bt.kW_rated)) for bt in net.batteries]))
    sizes = vcat([(; kind = "hours", value = h, net = with_duration(deck, h)) for h in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)],
                 [(; kind = "power", value = k, net = with_power(deck, k)) for k in (0.25, 0.5, 2.0, 4.0)])
    rows = []
    for z in sizes
        sz = setup(z.net, pvs; T = 24)
        w = Dict(case => window(sz, case, A) for case in (:soft, :cyclic))
        say(@sprintf("  %-5s x %-5s  soft %s | cyclic %s%s", z.kind, string(z.value), fmt(w[:soft]), fmt(w[:cyclic]), w[:soft].bad + w[:cyclic].bad == 0 ? "" : "   ** failed checks"))
        push!(rows, (; z.kind, z.value, w))
    end
    open(joinpath(OUT_DIR, "bound_storage.csv"), "w") do io
        println(io, "varied,value,case,lo_deg,hi_deg,failed_checks,predicted_lo_deg,predicted_hi_deg,no_battery_lo_deg,no_battery_hi_deg")
        for r in rows, case in (:soft, :cyclic)
            @printf io "%s,%g,%s,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f\n" r.kind r.value case r.w[case].lo r.w[case].hi r.w[case].bad day_lo day_hi hour_lo hour_hi
        end
    end

    plot_style!()
    SOFT, CYC = colorant"#7a4fd0", colorant"#1c8a5a"
    x = [y.dd for y in sweep]
    p1 = plot(x, [Ea(y) / 1000 for y in sweep]; label = label(A), color = SUB_COLORS[1], lw = 2, xlabel = "$(dlabel(A)) − $(dlabel(B)) (deg)", ylabel = "MWh per day",
              title = "Daily energy imported, batteries idle", legend = :top)
    plot!(p1, x, [Eb(y) / 1000 for y in sweep]; label = label(B), color = SUB_COLORS[2], lw = 2)
    hline!(p1, [0]; color = MUTED, lw = 1, label = "")
    vline!(p1, [hour_lo, hour_hi]; color = SUB_COLORS[1], lw = 1, ls = :dash, label = "First hour back-feeds")
    vline!(p1, [day_lo, day_hi]; color = CYC, lw = 1.5, ls = :dash, label = "Daily import reaches zero")
    hrs = [r for r in rows if r.kind == "hours"]
    hv = [r.value for r in hrs]
    p2 = plot(hv, [r.w[:soft].hi for r in hrs]; label = "Soft end SOC", color = SOFT, lw = 2, marker = :circle, xscale = :log2, xticks = (hv, string.(hv)),
              xlabel = "Battery duration (h) at the same kW", ylabel = "Window edge (deg)", title = "Window against stored energy", legend = :right)
    plot!(p2, hv, [r.w[:soft].lo for r in hrs]; label = "", color = SOFT, lw = 2, marker = :circle)
    plot!(p2, hv, [r.w[:cyclic].hi for r in hrs]; label = "End = start", color = CYC, lw = 2, marker = :diamond)
    plot!(p2, hv, [r.w[:cyclic].lo for r in hrs]; label = "", color = CYC, lw = 2, marker = :diamond)
    hline!(p2, [day_lo, day_hi]; color = CYC, lw = 1.5, ls = :dash, label = "Daily import reaches zero")
    hline!(p2, [hour_lo, hour_hi]; color = SUB_COLORS[1], lw = 1, ls = :dash, label = "No batteries")
    savefig(plot(p1, p2; layout = (1, 2), size = (1150, 400), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm), joinpath(OUT_DIR, "fig_bound.png"))
    say("  wrote bound_sweep.csv, bound_storage.csv, fig_bound.png")
end

# ---- stage: tables ------------------------------------------------------------------------------------------
if "tables" in STAGES
    readcsv(f) = (p = joinpath(OUT_DIR, f); isfile(p) && filesize(p) > 0 ? [split(l, ',') for l in readlines(p)[2:end]] : nothing)
    cell(lo, hi, ok, loc, hic) = ok == "true" ? @sprintf("\$[%s%.2f, %s%+.2f]\$", loc == "true" ? "\\le " : "", parse(Float64, lo), hic == "true" ? "\\ge " : "", parse(Float64, hi)) : "infeasible"
    function pivot(rows, nkey)                           # rows of (key..., case, ok, lo, hi, lo_at_limit, hi_at_limit, bad)
        d, order = Dict(), []
        for r in rows
            k = Tuple(String.(r[1:nkey]))
            haskey(d, k) || (push!(order, k); d[k] = Dict())
            d[k][String(r[nkey+1])] = cell(r[nkey+3], r[nkey+4], r[nkey+2], r[nkey+5], r[nkey+6])
        end
        return d, order
    end
    # Tables report the soft terminal SOC only; the cyclic diagnostic stays in the CSVs.
    cols = Dict("nobatt" => "No batteries", "hourly" => "Batt., hourly", "soft" => "Batt., whole day")   # must fit one IEEE column
    three = ("nobatt", "soft")
    open(joinpath(OUT_DIR, "tables.tex"), "w") do io
        println(io, "% Generated by envs/multi_poi/root_level/multiperiod_study.jl -- do not edit by hand.")
        function table(file, nkey, caption, lab, head, keyfmt; cases = three)
            rows = readcsv(file)
            rows === nothing && return
            d, order = pivot(rows, nkey)
            println(io, "\\begin{table}[!t]\n\\centering\n\\caption{$caption}\n\\label{$lab}\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{l", "l"^length(cases), "}\n\\hline")
            println(io, head, " & ", join((cols[c] for c in cases), " & "), " \\\\\n\\hline")
            for k in order
                println(io, keyfmt(k), " & ", join((get(d[k], c, "--") for c in cases), " & "), " \\\\")
            end
            println(io, "\\hline\n\\end{tabular}\n\\end{table}\n")
        end
        pvname = Dict("off" => "Off", "unity" => "On, unity p.f.", "qfree" => "On, \$q\$ free")
        table("main_horizon_windows.csv", 1, "Feasible window of a constant \$\\delta_{s_3}-\\delta_{s_4}\$ (deg) over a 24\\,h day, ieee123. ``Hourly'': each hour solved alone, batteries without memory.", "tab:mp-main", "PV", k -> pvname[k[1]]; cases = ("nobatt", "hourly", "soft"))
        table("horizon_windows.csv", 1, "The same window versus horizon length \$T\$ (\$\\Delta t = 24/T\$\\,h).", "tab:mp-horizon", "\$T\$", k -> k[1])
        table("zs_windows.csv", 1, "The same window versus source impedance \$|Z_s|\$ (\$\\Omega\$), \$T=24\$.", "tab:mp-zs", "\$|Z_s|\$", k -> k[1] == "0.0" ? "0 (stiff)" : (k[1] == "0.8" ? "0.8 (typical)" : k[1]))
        table("duration_windows.csv", 1, "The same window versus battery duration (kWh per kW), \$T=24\$.", "tab:mp-duration", "Hours", k -> k[1])
        table("pairs_windows.csv", 3, "Windows for every substation pair, \$T=24\$; \$|Z|\$ is the feeder impedance between the pair.", "tab:mp-pairs", "Pair (\$|Z|\$, \$\\Omega\$)", k -> @sprintf("%s--%s (%.2f)", replace(k[1], "subs" => ""), replace(k[2], "subs" => ""), parse(Float64, k[3])))
        table("five_windows.csv", 1, "All five substations: window of one substation's angle against the other four (deg), \$T=24\$.", "tab:mp-five", "Moved", k -> "Subs " * replace(k[1], "subs" => ""))
    end
    say("\nwrote ", joinpath(OUT_DIR, "tables.tex"))
end
