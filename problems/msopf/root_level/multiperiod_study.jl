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
#   cyclic    batteries with SOC dynamics and a HARD terminal SOC, B[T] == B[0]: a day that
#             can repeat. THE multi-period result reported in tables and figures (user,
#             2026-10-11: feasibility is the concern here, so the terminal SOC is hard);
#   soft      the same with only the soft terminal penalty of the DDP work; a diagnostic kept
#             in the CSVs. Its window is one day's borrowing: the batteries may end fuller,
#             so it grows without bound with battery size (stage "bound").
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
#   estimate  that limit as an estimate (worst hour without batteries, daily mean with them)
#             against the measured windows of stage main, per PV mode;
#   tables    rebuild the LaTeX tables from the CSVs (no solves).
# Two more, heavy and run only when named:
#   ieee2522  the 2522-bus deck with a second substation added in memory at several buses:
#             no-battery and whole-day windows, the same limit estimate, an OpenDSS replay;
#   large10k  the same on the 10k-bus deck, a substation at the end of one feeder.
#   large10kd the 10k-bus deck with loads, PV and batteries divided by 20 (in memory).
#
# Run from the repo root:
#     julia --project=problems/msopf problems/msopf/root_level/multiperiod_study.jl [stage ...]
#
# Writes (tracked) problems/msopf/results/multiperiod/.

include(joinpath(@__DIR__, "..", "multiperiod_full_angle.jl"))
include(joinpath(@__DIR__, "..", "sweep_common.jl"))
using Logging, SparseArrays, Random

const SYSTEM = "ieee123_5poi_1ph"
const MAIN_PAIR = ("subs3", "subs4")
const ALL_SUBS = ["subs1", "subs2", "subs3", "subs4", "subs5"]
const X_OVER_R = 25.0
# MULTIPERIOD_OUT_DIR gives a run its own output folder, so concurrent jobs (Kamiak) cannot
# overwrite each other's CSVs; their rows are merged into the tracked files afterwards.
const OUT_DIR = normpath(get(ENV, "MULTIPERIOD_OUT_DIR", joinpath(@__DIR__, "..", "results", "multiperiod")))
const STAGES = isempty(ARGS) ? ["main", "horizon", "zs", "duration", "pairs", "three", "five", "bound", "estimate", "tables"] : lowercase.(ARGS)
mkpath(OUT_DIR)
say(args...) = (println(args...); flush(stdout))

# ---- setup -------------------------------------------------------------------------------------------
function load_sources(srcs; system = SYSTEM, extra = String[])
    kw = (; sources = collect(srcs), snapshot = true, load_band = (0.5, 1.5), extra)
    pvs = with_logger(NullLogger()) do                  # "? PVSystem.x.Pmpp" logs a warning per unit
        read_pv(system; kw...)
    end
    compile_deck(system; disable = ["PVSystem"], kw...)
    return read_network(), pvs
end

with_zs(net, zmag) = (r = zmag / sqrt(1 + X_OVER_R^2);
    merge(net, (; sources = [merge(s, (; z_series = complex(r, X_OVER_R * r))) for s in net.sources])))
with_duration(net, hours) = merge(net, (; batteries = [merge(bt, (; kWh_rated = hours * bt.kW_rated)) for bt in net.batteries]))
banded(net) = [b for b in net.buses if !(b in Set(s.bus for s in net.sources))]
customer_buses(net) = Set(l.to for l in net.lines if startswith(lowercase(l.name), "svc_"))
# Primary-only decks retain their historical band; service-transformer decks constrain customers only.
study_band(net) = (bs = customer_buses(net); isempty(bs) ? nothing : bs)
study_banded(s) = s.band_buses === nothing ? banded(s.net) : collect(s.band_buses)

"Everything one configuration needs: profiles, the PV mode's inputs, and the cost weights."
# `feas` (per unit) is the feasibility tolerance the user's DDP work uses, 1e-5 (2026-10-12): no
# limit violated by more than `feas`, with Ipopt's constraint tolerance at `feas` and its
# convergence tolerance a decade tighter so that a slack at zero ends well below it. `nothing`
# keeps the tight default, 1e-9, and the summed-slack test (VTOL).
function setup(net, pvs; T, pv = "qfree", ideal = false, feas = nothing, band_buses = study_band(net))
    lam, c, pvf = load_profile(T), price_profile(T), pv_profile(T)
    P, F = pv == "off" ? ([], zeros(T)) : (pvs, pvf)
    ipopt = feas === nothing ? (; tol = IPOPT_TOL, constr_viol_tol = IPOPT_TOL) : (; tol = feas / 10, constr_viol_tol = feas)
    return (; net, T, lam, c, P, F, q = pv == "unity" ? :unity : :optimize, CB = cb_by_rule(net), GAM = gamma_by_rule(net), ideal, feas, ipopt, band_buses)
end

const DAY = :cyclic            # the whole-day case that is reported: hard terminal SOC
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
    common = (; pvs = s.P, pv_q = s.q, C_B = s.CB, ideal_sources = s.ideal, band_buses = s.band_buses, battery = case === :nobatt ? :idle : :optimize, s.ipopt...,
              # MULTI_VMARGIN=k: margin solves count voltage room too (k kW per pu of |V|^2), so edges set by
              # the band are found by the same secant as those set by no-backflow (see solve_multiperiod)
              margin_band_kw = haskey(ENV, "MULTI_VMARGIN") ? parse(Float64, ENV["MULTI_VMARGIN"]) : nothing)
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
    # What the edge solution sits against: the smallest substation import over the day (kW;
    # ~0 when no-backflow binds) and the voltage extremes (pu; 0.95 / 1.05 when the band binds).
    bind = fill((; minP = NaN, vmin = NaN, vmax = NaN), length(dirs))
    feasible || return (; feasible, r, capped, bad, bind)
    tag = hour === nothing ? string(case) : "$case, hour $hour"
    nb = study_banded(s)
    for (i, dir) in enumerate(dirs)
        x = solve_multiperiod(net, equal_angles(net); kw..., objective = :ray, ray_dir = dir, ray_limit = LIMIT, start = at0)
        if !(x.status in RAY_OK) || !isfinite(x.r)
            bad += 1
            say("  ** ray solve failed (", tag, "): ", x.status)
            continue
        end
        r[i], capped[i] = x.r, x.r >= LIMIT - 1e-6
        bind[i] = (; minP = minimum(minimum(v) for v in values(x.P_subs_kW)),
                   vmin = minimum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T), vmax = maximum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T))
        (check === :none || x.r <= 1e-6) && continue
        inn = solve_multiperiod(net, Dict(k => 0.99 * x.r * v for (k, v) in dir); kw..., objective = :violation, start = x)
        good = solved(inn) && inn.violation <= VTOL
        if check === :both && !capped[i]
            out = solve_multiperiod(net, Dict(k => (x.r + 0.02) * v for (k, v) in dir); kw..., objective = :violation, start = x)
            good &= solved(out) && out.violation > VTOL
        end
        good || (bad += 1; say(@sprintf("  ** edge check failed (%s): r = %.4f deg", tag, x.r)))
    end
    return (; feasible, r, capped, bad, bind)
end

"Window of delta_a, every other source at 0: the two rays along +a and -a."
function window(s, case, a; kwargs...)
    dir(sgn) = Dict(x.name => (x.name == a ? sgn : 0.0) for x in s.net.sources)
    y = rays(s, case, [dir(+1.0), dir(-1.0)]; kwargs...)
    return (; y.feasible, lo = -y.r[2], hi = y.r[1], lo_capped = y.capped[2], hi_capped = y.capped[1], y.bad,
            lo_bind = y.bind[2], hi_bind = y.bind[1])
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

    # The cost-optimal day (hard terminal SOC, PV reactive power free) at equal angles and at an
    # angle only batteries make feasible: halfway between the no-battery edge and the day's edge.
    s, w = S["qfree"], res["qfree"][1]
    d_star = round((w[:nobatt].hi + w[DAY].hi) / 2; digits = 2)
    days = []
    for dd in (0.0, d_star)
        d = Dict(A => dd, B => 0.0)
        x = solve_case(s, d, DAY)
        rp = replay_opendss(net, x, d; load_factor = s.lam, pvs = s.P, pv_factor = s.F)
        say(@sprintf("  day at %+.2f deg: %s, cost \$%.2f (energy %.2f, cycling %.2f); OpenDSS replay |dP| %.1e kW, |dQ| %.1e kvar, |dV| %.1e pu, converged %s, PQ %s", dd, x.status, x.cost.total, x.cost.energy, x.cost.cycling, rp.dP, rp.dQ, rp.dV, rp.converged, rp.pq_band))
        push!(days, (; dd, x, rp))
    end
    PL, PVk = sum(d.kW for d in net.loads), sum(p.kW for p in pvs)
    nb = study_banded(s)
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
        r[DAY].feasible && hline!(p, [r[DAY].lo, r[DAY].hi]; color = SOFT, lw = 2, label = "Batteries, whole day")
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
    for (case, col, lab, alpha) in ((DAY, colorant"#7a4fd0", "Batteries, whole day", 0.18), (:nobatt, SUB_COLORS[1], "No batteries", 0.35))
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

# ---- stage: estimate --------------------------------------------------------------------------------------
"""
    import_estimate(s, a, b; h = 0.1)

Where the window of delta_a - delta_b should close if a constant angle difference only adds a
constant circulating power. Batteries idle and no-backflow off: the imports at equal angles
and their slope per degree (central difference over +-h). Returns the worst-hour estimate
`hour` (no batteries: the first hour to back-feed) and the daily-mean estimate `day`
(batteries that end the day where they started: the daily energy import reaches zero).
"""
function import_estimate(s, a, b; h = 0.1)
    idle(dd) = solve_multiperiod(s.net, Dict(a => dd, b => 0.0); problem(s, :nobatt)..., no_backflow = false)
    x0, xp, xm = idle(0.0), idle(+h), idle(-h)
    (x0.ok && xp.ok && xm.ok) || return nothing
    Ka = (xp.P_subs_kW[a] .- xm.P_subs_kW[a]) ./ 2h       # kW per degree into a
    Kb = (xm.P_subs_kW[b] .- xp.P_subs_kW[b]) ./ 2h       # ... and out of b
    Pa, Pb = x0.P_subs_kW[a], x0.P_subs_kW[b]
    return (; hour = (-minimum(Pa ./ Ka), minimum(Pb ./ Kb)), day = (-sum(Pa) / sum(Ka), sum(Pb) / sum(Kb)),
            K = sum(Ka) / length(Ka), Pa = sum(Pa) / length(Pa), Pb = sum(Pb) / length(Pb), minPa = minimum(Pa), minPb = minimum(Pb))
end

# The estimate against the measured windows of stage "main", per PV mode. With PV off or at
# unity power factor nothing but the batteries can be dispatched and the estimate should be
# exact; with reactive power free, the edge solve steers it and the estimate (which uses the
# cost-optimal reactive power at equal angles) does not.
function estimate_stage()
    A, B = MAIN_PAIR
    net, pvs = load_sources([A, B])
    file = joinpath(OUT_DIR, "main_horizon_windows.csv")
    isfile(file) || (say("  ** estimate needs main_horizon_windows.csv: run stage main first"); return)
    measured = Dict((String(r[1]), String(r[2])) => (parse(Float64, r[4]), parse(Float64, r[5])) for r in (split(l, ',') for l in readlines(file)[2:end]))
    say("\n", "="^110, "\nESTIMATE  $A + $B, T = 24: constant circulating power against the measured windows")
    open(joinpath(OUT_DIR, "bound_estimate.csv"), "w") do io
        println(io, "pv,case,estimate_lo_deg,estimate_hi_deg,measured_lo_deg,measured_hi_deg,kW_per_deg,mean_import_a_kW,mean_import_b_kW,min_import_a_kW,min_import_b_kW")
        for m in ("off", "unity", "qfree")
            e = import_estimate(setup(net, pvs; T = 24, pv = m), A, B)
            e === nothing && (say("  ** idle solves failed for PV ", m); continue)
            for (case, est) in (("nobatt", e.hour), (string(DAY), e.day))
                lo, hi = measured[(m, case)]
                say(@sprintf("  PV %-5s %-6s estimate [%+.3f, %+.3f], measured [%+.3f, %+.3f]", m, case, est[1], est[2], lo, hi))
                @printf io "%s,%s,%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f,%.3f,%.3f\n" m case est[1] est[2] lo hi e.K e.Pa e.Pb e.minPa e.minPb
            end
        end
    end
    say("  wrote bound_estimate.csv")
end
"estimate" in STAGES && estimate_stage()

# ---- stage: simplex ----------------------------------------------------------------------------------------
# The window with several substations, predicted from one operating point. With batteries idle each
# substation's import is close to linear in the angles, P_k(delta) ~ P_k(0) + sum_j K_kj delta_j (K: the
# angle-to-import sensitivities). Each no-backflow limit P_k >= 0 is then a half-space, so with N
# substations the window is a simplex in the N-1 angle differences: an interval, a triangle, a
# tetrahedron. Along a direction the predicted reach is the first substation to reach zero import: in
# its worst hour without batteries, in its daily mean with them (the two-substation estimate,
# generalized). The stage checks this against the measured three- and five-substation windows.

"""
    sensitivity(s, movers; h = 0.1)

Every substation's import at equal angles with the batteries idle (kW per hour) and its derivative
with respect to each moving substation's angle (kW/deg per hour), by central differences.
"""
function sensitivity(s, movers; h = 0.1, replay_file = nothing)
    names = [x.name for x in s.net.sources]
    function at(dd; start = nothing)
        angles = merge(Dict(n => 0.0 for n in names), dd)
        x = solve_multiperiod(s.net, angles; problem(s, :nobatt)..., no_backflow = false, start)
        x.ok || return x
        rp = replay_opendss(s.net, x, angles; load_factor = s.lam, pvs = s.P, pv_factor = s.F)
        good = rp.converged && rp.pq_band && rp.dP <= 1e-3 && rp.dQ <= 1e-3 && rp.dV <= 1e-7 && rp.dB <= 1e-3
        say(@sprintf("    sensitivity replay: %s, dP %.2e kW, dQ %.2e kvar, dV %.2e pu", good ? "PASS" : "FAIL", rp.dP, rp.dQ, rp.dV))
        if replay_file !== nothing
            fresh = !isfile(replay_file)
            open(replay_file, "a") do io
                fresh && println(io, "angles_deg,status,iterations,solve_seconds,converged,pq_band,dP_kW,dQ_kvar,dV_pu,dB_kW,passed")
                println(io, join([join([@sprintf("%s=%.6f", k, v) for (k, v) in sort(collect(angles))], ";"),
                    x.status, x.iterations, x.time, rp.converged, rp.pq_band, rp.dP, rp.dQ, rp.dV, rp.dB, good], ','))
            end
        end
        good || error("OpenDSS rejected sensitivity dispatch; no prediction will be written")
        return x
    end
    x0 = at(Dict{String,Float64}())
    x0.ok || return nothing
    K = Dict{Tuple{String,String},Vector{Float64}}()
    for m in movers
        xp, xm = at(Dict(m => +h); start = x0), at(Dict(m => -h); start = x0)
        (xp.ok && xm.ok) || return nothing
        for k in names
            K[(k, m)] = (xp.P_subs_kW[k] .- xm.P_subs_kW[k]) ./ 2h
        end
    end
    return (; names, P = Dict(k => x0.P_subs_kW[k] for k in names), K, T = s.T)
end

"Predicted reach along `dir` (mover => degrees per unit length) and the substation that limits it."
function predicted_reach(sens, dir; day::Bool)
    r, lim = Inf, ""
    for k in sens.names
        slope = sum(sens.K[(k, m)] .* v for (m, v) in dir)
        reach = day ? (sum(slope) < 0 ? sum(sens.P[k]) / -sum(slope) : Inf) :
                      minimum((slope[t] < 0 ? sens.P[k][t] / -slope[t] : Inf) for t in 1:sens.T)
        reach < r && ((r, lim) = (reach, k))
    end
    return r, lim
end

# SIMPLEX_PV=qfree (default) compares with the committed three_region.csv and five_windows.csv. Any other
# PV mode (unity, off) measures those windows afresh in that mode: with PV reactive power fixed only the
# batteries can be dispatched, which is where the two-substation estimate is exact.
if "simplex" in STAGES
    PVM = get(ENV, "SIMPLEX_PV", "qfree")
    committed = PVM == "qfree"
    readrows(f) = [split(l, ',') for l in readlines(joinpath(OUT_DIR, f))[2:end]]
    say("\n", "="^110, "\nSIMPLEX  the multi-substation window predicted from the angle-to-import sensitivities, T = 24, PV $PVM",
        committed ? " (measured: the committed windows)" : " (measured afresh)")
    out = ["study,pv,direction,case,measured_deg,predicted_deg,limiting_substation"]
    function compare(tag, items)
        for case in ("nobatt", string(DAY))
            sel = [x for x in items if x.case == case]
            isempty(sel) && continue
            err = [abs(x.pred - x.meas) for x in sel]
            rel = [abs(x.pred - x.meas) / abs(x.meas) for x in sel]
            say(@sprintf("  %-6s %-7s %2d directions: |predicted - measured| max %.4f deg, mean %.4f deg; relative max %.1f%%, mean %.1f%%",
                         tag, case, length(sel), maximum(err), sum(err) / length(err), 100maximum(rel), 100sum(rel) / length(rel)))
            for x in sel
                push!(out, @sprintf("%s,%s,%s,%s,%.6f,%.6f,%s", tag, PVM, x.dir, case, x.meas, x.pred, x.lim))
            end
        end
    end
    # three substations, the 16 directions of stage "three"
    REF, A, B = "subs1", "subs3", "subs4"
    net, pvs = load_sources([REF, A, B])
    s3 = setup(net, pvs; T = 24, pv = PVM)
    sens = sensitivity(s3, [A, B])
    if sens === nothing
        say("  ** idle solves failed (three substations)")
    else
        for k in sens.names
            say(@sprintf("    %s: mean import %.1f kW; sensitivity to %s %+.1f, to %s %+.1f kW/deg (daily mean)", k, sum(sens.P[k]) / sens.T,
                         A, sum(sens.K[(k, A)]) / sens.T, B, sum(sens.K[(k, B)]) / sens.T))
        end
        ths = [2pi * k / 16 for k in 0:15]                # the 16 directions of stage "three"
        measured = Dict{Tuple{String,Float64},Float64}()
        if committed
            for r in readrows("three_region.csv")
                measured[(String(r[2]), parse(Float64, r[1]))] = parse(Float64, r[4])
            end
        else
            for case in (:nobatt, DAY)
                y = rays(s3, case, [Dict(REF => 0.0, A => cos(th), B => sin(th)) for th in ths])
                y.bad == 0 || say("  ** ", y.bad, " failed ray solves or edge checks (three, ", case, ")")
                for (i, th) in enumerate(ths)
                    measured[(string(case), round(rad2deg(th); digits = 2))] = y.r[i]
                end
            end
        end
        items = []
        for case in ("nobatt", string(DAY)), th in ths
            key = (case, round(rad2deg(th); digits = 2))
            haskey(measured, key) || continue
            p, lim = predicted_reach(sens, Dict(A => cos(th), B => sin(th)); day = case == string(DAY))
            push!(items, (; dir = @sprintf("%.2f", key[2]), case, meas = measured[key], pred = p, lim))
        end
        compare("three", items)
    end
    # five substations, each moved alone against the other four (stage "five")
    net5, pvs5 = load_sources(ALL_SUBS)
    s5 = setup(net5, pvs5; T = 24, pv = PVM)
    sens5 = sensitivity(s5, ALL_SUBS)
    if sens5 === nothing
        say("  ** idle solves failed (five substations)")
    else
        measured5 = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
        if committed
            for r in readrows("five_windows.csv")
                measured5[(String(r[1]), String(r[2]))] = (parse(Float64, r[4]), parse(Float64, r[5]))
            end
        else
            for m in ALL_SUBS, case in (:nobatt, DAY)
                w = window(s5, case, m)
                w.bad == 0 || say("  ** ", w.bad, " failed ray solves or edge checks (five, ", m, ", ", case, ")")
                measured5[(m, string(case))] = (w.lo, w.hi)
            end
        end
        items = []
        for m in ALL_SUBS, case in ("nobatt", string(DAY))
            haskey(measured5, (m, case)) || continue
            lo, hi = measured5[(m, case)]
            day = case == string(DAY)
            up, limu = predicted_reach(sens5, Dict(m => 1.0); day)
            dn, limd = predicted_reach(sens5, Dict(m => -1.0); day)
            push!(items, (; dir = m * " up", case, meas = hi, pred = up, lim = limu))
            push!(items, (; dir = m * " down", case, meas = -lo, pred = dn, lim = limd))
        end
        compare("five", items)
    end
    file = committed ? "simplex_check.csv" : "simplex_check_$PVM.csv"
    open(io -> foreach(l -> println(io, l), out), joinpath(OUT_DIR, file), "w")
    say("  wrote ", file)
end

# ---- stages: ieee2522, large10k ---------------------------------------------------------------------------
# The two big decks have one substation. A second one is added in memory (the decks on disk are
# untouched): a Vsource at 1.05 pu behind the source impedance, tied to a chosen bus by a
# 1e-6 ohm jumper. Where it goes is a dimension of the study, so each deck lists several buses.
#   ieee2522: Zs = 0.311 ohm, the 0.5 GVA short-circuit level of the 20 kV decks applied at
#             12.47 kV (an assumption: the literature value is for 20 kV). Buses 947, 996, 1120
#             and 1216 lie along the main trunk (3.4 to 10.1 ohm from the substation); 546 is
#             the electrically farthest bus, the tip of a thin lateral.
#   large10k: 102 identical feeders in a star, 1032 MW in total; the head bus is an aggregate,
#             so the sources stay stiff. Bus 1120 is the end of one feeder.
#   large10kd: large10kC with every load, PV unit and battery (power and energy) divided by 20,
#             in memory: 51.6 MW in total, about 0.5 MW per feeder (user, 2026-10-12: the C deck's
#             assets suit solver tests more than a realistic feeder). Rows carry "large10kD_1ph".
#             Its line sections also get R/X = 1 (user, 2026-10-13): the C deck has R/X = 7 on
#             every ordinary section (0.0175 + j0.0025 ohm), a synthetic number no 12.47 kV
#             trunk has; the review in resources/rx_ratio_report.md finds 0.3-3.4 at medium
#             voltage and recommends 1.0 as a single value. R is kept, X raised to 0.0175 ohm.
# Not in the default stage list (heavy): pass the stage name. LARGE_AT=bus[,bus] and LARGE_T=n
# restrict a run; rows are merged into the CSV by (system, bus, T), so runs can be split.
const LARGE = Dict(
    "ieee2522" => (; system = "ieee2522C_1ph", label = "ieee2522C_1ph", scale = 1.0, rx = nothing, ray = true, feas = nothing, kV = 7.2, zs = 0.311, T = 12, at = ["947", "996", "1120", "1216", "546"]),
    "large10k" => (; system = "large10kC_1ph", label = "large10kC_1ph", scale = 1.0, rx = nothing, ray = false, feas = 1e-5, kV = 12.47, zs = 0.0, T = 24, at = ["1120"]),
    "large10kd" => (; system = "large10kC_1ph", label = "large10kD_1ph", scale = 20.0, rx = 1.0, ray = false, feas = 1e-5, kV = 12.47, zs = 0.0, T = 24, at = ["1120"]))
const LARGE_HEADER = "system,bus,T,case,feasible_at_equal_angles,lo_deg,hi_deg,lo_at_limit,hi_at_limit,failed_checks," *
                     "predicted_lo_deg,predicted_hi_deg,lo_min_import_kW,lo_vmin_pu,lo_vmax_pu,hi_min_import_kW,hi_vmin_pu,hi_vmax_pu," *
                     "z_between_ohm,kW_per_deg,mean_import_new_kW,mean_import_old_kW,load_kW,seconds"

"Impedance between the two sources' buses by one sparse solve (source_separation is dense: 1.7 GB at 10k buses)."
function separation_sparse(net)
    idx, n = Dict(b => i for (i, b) in enumerate(net.buses)), length(net.buses)
    I, J, V = Int[], Int[], ComplexF64[]
    for l in net.lines
        y, i, j = -l.Y[1, 2], idx[l.from], idx[l.to]
        append!(I, (i, j, i, j)); append!(J, (i, j, j, i)); append!(V, (y, y, -y, -y))
    end
    L = sparse(I, J, V, n, n)[1:n-1, 1:n-1]            # last bus as the reference node
    a, b = idx[net.sources[1].bus], idx[net.sources[2].bus]
    e = zeros(ComplexF64, n - 1)
    a < n && (e[a] += 1); b < n && (e[b] -= 1)
    x = L \ e
    return (a < n ? x[a] : 0.0im) - (b < n ? x[b] : 0.0im)
end

"Replace the rows of `path` whose first `nkey` fields match any of `rows`, keep the rest."
function merge_csv(path, header, nkey, rows)
    key(l) = join(split(l, ',')[1:nkey], ',')
    old = isfile(path) ? readlines(path)[2:end] : String[]
    new = Set(key.(rows))
    open(path, "w") do io
        println(io, header)
        foreach(l -> println(io, l), filter(l -> !(key(l) in new), old))
        foreach(l -> println(io, l), rows)
    end
end

"""
    window_bracketed(s, case, a; ray = true)

`window` with each edge tested from outside. A maximize-angle solve can stop early at Ipopt's
"acceptable" level, run out of iterations or end locally infeasible: on the large10k decks all
three happened (2026-10-12; one edge came back as +4.9 deg while every angle up to +45 is
feasible). So the violation test is run 0.02 deg beyond each edge. If that point is feasible
the edge was not one: step outward (0.5 deg, doubling) to the first infeasible point or the
search limit, then bisect six times (to 0.5 / 64 deg or better near zero).
`ray = false` skips the maximize-angle solves and brackets both edges from equal angles: the
choice for the large10k decks, where those solves cost up to 25 minutes each and then fail.
Assumes one feasible interval around zero; a doubling step can jump a narrow infeasible gap.
"""
function window_bracketed(s, case, a; ray = true)
    kw, net, nb = problem(s, case), s.net, study_banded(s)
    function test(sgn, r, prev)
        x = solve_multiperiod(net, Dict(k.name => (k.name == a ? sgn * r : 0.0) for k in net.sources); kw..., objective = :violation, start = prev)
        within = s.feas === nothing ? x.violation <= VTOL : x.max_violation <= s.feas
        return x.status in RAY_OK && isfinite(x.violation) && within, x
    end
    bind_of(x) = (; minP = minimum(minimum(v) for v in values(x.P_subs_kW)),
                  vmin = minimum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T), vmax = maximum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T))
    nobind = (; minP = NaN, vmin = NaN, vmax = NaN)
    if ray
        w = window(s, case, a)
        w.feasible || return w
        x0 = nothing
    else
        ok0, x0 = test(+1.0, 0.0, nothing)
        ok0 || return (; feasible = false, lo = NaN, hi = NaN, lo_capped = false, hi_capped = false, bad = 0, lo_bind = nobind, hi_bind = nobind)
        b0 = bind_of(x0)
        w = (; hi = 0.0, lo = 0.0, hi_capped = false, lo_capped = false, hi_bind = b0, lo_bind = b0, bad = 0)
    end
    function settle(sgn, r0, capped0, bind0)
        r = isfinite(r0) ? r0 : 0.0
        (capped0 || r >= LIMIT - 1e-6) && return (; r, capped = true, bind = bind0, moved = false)
        ok, x = test(sgn, r + 0.02, x0)
        ok || return (; r, capped = false, bind = bind0, moved = !isfinite(r0))
        r, last, step = r + 0.02, x, 0.5
        while r < LIMIT - 1e-6
            trial = min(r + step, LIMIT)
            ok, x = test(sgn, trial, last)
            if !ok
                rin, rout = r, trial                      # last feasible, first infeasible
                for _ in 1:6
                    mid = (rin + rout) / 2
                    ok, x = test(sgn, mid, last)
                    ok ? ((rin, last) = (mid, x)) : (rout = mid)
                end
                return (; r = rin, capped = false, bind = bind_of(last), moved = true)
            end
            r, last, step = trial, x, 2 * step
        end
        return (; r = LIMIT, capped = true, bind = bind_of(last), moved = true)
    end
    hi = settle(+1.0, w.hi, w.hi_capped, w.hi_bind)
    lo = settle(-1.0, -w.lo, w.lo_capped, w.lo_bind)
    for (side, e, old) in (("upper", hi, w.hi), ("lower", lo, -w.lo))
        ray && e.moved && say(@sprintf("  ** %s: the %s edge from the maximize-angle solve (%s deg) was not an edge; the outside test moved it to %s%.3f deg",
                                case, side, isfinite(old) ? @sprintf("%.3f", old) : "failed", e.capped ? ">= " : "", e.r))
    end
    # A failed ray solve is superseded once its edge has been bracketed.
    bad = hi.moved && lo.moved ? 0 : (hi.moved || lo.moved ? max(w.bad - 1, 0) : w.bad)
    return (; feasible = true, lo = -lo.r, hi = hi.r, lo_capped = lo.capped, hi_capped = hi.capped, bad, lo_bind = lo.bind, hi_bind = hi.bind)
end

const MTOL_KW = 0.05           # kW: an import margin this small is the edge (2e-4 deg at 239 kW/deg)

"""
    window_secant(s, case, a, K, pred_hi)

The window of delta_a from the import margin instead of a search. `objective = :margin` gives
the largest import every substation can be guaranteed at a fixed angle; it is zero at an edge
set by no-backflow and falls almost linearly with the angle, at the angle-to-import
sensitivity `K` (kW per degree, from `import_estimate`). So:
  1. the margin at equal angles, m0 (negative: infeasible there);
  2. lower edge: first guess m0 / K, then secant steps on the margin until |m| <= MTOL_KW;
  3. upper edge: the estimate `pred_hi`; beyond the search limit, one feasibility test there;
  4. each finite edge is confirmed by the violation test 0.02 deg inside and outside.
A margin solve that fails, a margin that does not fall, or a failed confirmation (an edge set
by the voltage band, say) hands the whole window to `window_bracketed`.
"""
function window_secant(s, case, a, K, pred_hi)
    kw, net, nb = problem(s, case), s.net, study_banded(s)
    at(sgn, r) = Dict(k.name => (k.name == a ? sgn * r : 0.0) for k in net.sources)
    vtol = something(s.feas, 1e-5)
    nm, nt = Ref(0), Ref(0)
    function margin(sgn, r, prev)
        x = solve_multiperiod(net, at(sgn, r); kw..., objective = :margin, start = prev); nm[] += 1
        return x.status in RAY_OK && isfinite(x.margin_kW) && x.max_violation <= vtol, x
    end
    function test(sgn, r, prev)
        x = solve_multiperiod(net, at(sgn, r); kw..., objective = :violation, start = prev); nt[] += 1
        within = s.feas === nothing ? x.violation <= VTOL : x.max_violation <= s.feas
        return x.status in RAY_OK && isfinite(x.violation) && within, x
    end
    bind_of(x) = (; minP = minimum(minimum(v) for v in values(x.P_subs_kW)),
                  vmin = minimum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T), vmax = maximum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T))
    giveup(why) = (say("  ** ", case, ": secant on the import margin gave up (", why, "); bracketing instead"); window_bracketed(s, case, a; ray = false))
    ok0, x0 = margin(+1.0, 0.0, nothing)
    ok0 || return giveup("margin solve at equal angles: " * string(x0.status))
    nobind = (; minP = NaN, vmin = NaN, vmax = NaN)
    x0.margin_kW >= -MTOL_KW || return (; feasible = false, lo = NaN, hi = NaN, lo_capped = false, hi_capped = false, bad = 0, lo_bind = nobind, hi_bind = nobind)
    function edge(sgn, guess)
        if !isfinite(guess) || guess >= LIMIT             # nothing binds before the search limit, if the estimate is right
            okc, xc = test(sgn, LIMIT, x0)
            okc && return (; r = LIMIT, capped = true, bind = bind_of(xc))
            guess = LIMIT / 2
        end
        ra, ma, rb, xb, mb = 0.0, x0.margin_kW, clamp(guess, 0.02, LIMIT), x0, NaN
        for _ in 1:6
            okb, xb = margin(sgn, rb, xb)
            okb || return "margin solve at $(round(sgn * rb; digits = 3)) deg: $(xb.status), max violation $(xb.max_violation)"
            mb = xb.margin_kW
            abs(mb) <= MTOL_KW && break
            slope = (mb - ma) / (rb - ra)
            slope < -1e-6 || (slope = -abs(K))            # the margin must fall with distance: use the sensitivity
            rn = clamp(rb - mb / slope, 0.0, LIMIT)
            ra, ma, rb = rb, mb, rn
        end
        abs(mb) <= MTOL_KW || return "no convergence in 6 margin solves"
        okin, _ = test(sgn, max(rb - 0.02, 0.0), xb)
        okout, _ = test(sgn, rb + 0.02, xb)
        (okin && !okout) || return @sprintf("edge %.4f deg not confirmed: inside feasible %s, outside feasible %s", sgn * rb, okin, okout)
        return (; r = rb, capped = false, bind = bind_of(xb))
    end
    lo = edge(-1.0, x0.margin_kW / abs(K))
    lo isa String && return giveup("lower edge, " * lo)
    hi = edge(+1.0, pred_hi)
    hi isa String && return giveup("upper edge, " * hi)
    say(@sprintf("           %s: import margin at equal angles %.1f kW; edges from %d margin solves and %d feasibility tests", case, x0.margin_kW, nm[], nt[]))
    return (; feasible = true, lo = -lo.r, hi = hi.r, lo_capped = lo.capped, hi_capped = hi.capped, bad = 0, lo_bind = lo.bind, hi_bind = hi.bind)
end

function large_stage(name)
    cfg = LARGE[name]
    T = parse(Int, get(ENV, "LARGE_T", string(cfg.T)))
    at = haskey(ENV, "LARGE_AT") ? String.(split(ENV["LARGE_AT"], ',')) : cfg.at
    A, B = "subs2", "source"                              # the added substation and the deck's own
    rz = cfg.zs / sqrt(1 + X_OVER_R^2)
    imp = cfg.zs == 0 ? "R1=0 X1=1e-8 R0=0 X0=1e-8" : @sprintf("R1=%.8g X1=%.8g R0=%.8g X0=%.8g", rz, X_OVER_R * rz, rz, X_OVER_R * rz)
    say("\n", "="^110, "\n", uppercase(name), "  $(cfg.label) with a second substation, T = $T (dt = $(24 / T) h), Zs = $(cfg.zs) ohm",
        cfg.scale == 1 ? "" : "; loads, PV and batteries of $(cfg.system) divided by $(cfg.scale)",
        cfg.rx === nothing ? "" : "; line sections at R/X = $(cfg.rx)")
    # Sections whose R/X exceeds the target get X = R / target; R and every other section stay.
    reactance = String[]
    if cfg.rx !== nothing
        for l in eachline(joinpath(REPO_ROOT, "test_networks", cfg.system, "BranchData.dss"))
            m = match(r"^New Line\.(\S+) .*\br1=(\S+) x1=(\S+)"i, l)
            m === nothing && continue
            r1, x1 = parse(Float64, m[2]), parse(Float64, m[3])
            r1 / x1 > cfg.rx + 1e-9 && push!(reactance, @sprintf("Edit Line.%s x1=%.10g", m[1], r1 / cfg.rx))
        end
        say("  ", length(reactance), " line sections edited to R/X = ", cfg.rx)
    end
    for far in at
        t0 = time()
        extra = [reactance;
                 "New Line.poi2 Phases=1 Bus1=2s Bus2=$far r1=1e-06 x1=1e-06 c1=0 c0=0";
                 "New Vsource.subs2 Phases=1 Bus1=2s basekv=$(cfg.kV) pu=1.05 angle=0 $imp";
                 "Edit Vsource.source $imp"]
        net, pvs = load_sources([B, A]; system = cfg.system, extra)
        if cfg.scale != 1                                 # in memory; the OpenDSS replay scales its loads by LoadMult
            k = cfg.scale
            net = merge(net, (; loads = [merge(d, (; kW = d.kW / k, kvar = d.kvar / k)) for d in net.loads],
                                batteries = [merge(b, (; kW_rated = b.kW_rated / k, kVA = b.kVA / k, kWh_rated = b.kWh_rated / k)) for b in net.batteries]))
            pvs = [merge(p, (; kW = p.kW / k, kVA = p.kVA / k)) for p in pvs]
        end
        s = setup(net, pvs; T, cfg.feas)
        PL = sum(d.kW for d in net.loads)
        zc = separation_sparse(net)
        zsep = abs(zc)
        zloop = abs(zc + sum(x.z_series for x in net.sources))
        say(@sprintf("\n  substation at bus %s: %d buses, load %.0f kW, %d batteries %.0f kW / %.0f kWh, %d PV %.0f kW; |Z| between the two %.3f ohm (R/X %.2f; loop %.3f)",
                     far, length(net.buses), PL, length(net.batteries), sum(b.kW_rated for b in net.batteries), sum(b.kWh_rated for b in net.batteries), length(pvs), sum(p.kW for p in pvs), zsep, real(zc) / imag(zc), zloop))

        e = import_estimate(s, A, B)
        pred = (; nobatt = (NaN, NaN), cyclic = (NaN, NaN))
        K = Pa = Pb = NaN
        if e !== nothing
            K, Pa, Pb, pred = e.K, e.Pa, e.Pb, (; nobatt = e.hour, cyclic = e.day)
            say(@sprintf("    idle: imports at equal angles %s %.0f kW mean (%.0f min), %s %.0f kW mean (%.0f min); %.0f kW/deg (V^2/|Z_loop| = %.0f)",
                         A, Pa, e.minPa, B, Pb, e.minPb, K, (1.05 * cfg.kV)^2 * 1e3 / zloop * pi / 180))
        else
            say("    ** idle solves failed")
        end

        rows = String[]
        for case in (:nobatt, DAY)
            tc = time()
            p = case === :nobatt ? pred.nobatt : pred.cyclic
            edges = Symbol(get(ENV, "LARGE_EDGES", cfg.ray ? "ray" : "secant"))     # ray | bracket | secant
            w = edges === :secant && e !== nothing ? window_secant(s, case, A, K, p[2]) : window_bracketed(s, case, A; ray = edges === :ray)
            say(@sprintf("    %-6s window %s deg, predicted [%+.3f, %+.3f]%s   (%.0f s, peak RSS %.1f GB)", case, fmt(w), p[1], p[2], w.bad == 0 ? "" : "   ** failed checks", time() - tc, Sys.maxrss() / 2^30))
            w.feasible && say(@sprintf("           lower edge: smallest import %.2f kW, V %.4f-%.4f pu; upper edge: smallest import %.2f kW, V %.4f-%.4f pu",
                                       w.lo_bind.minP, w.lo_bind.vmin, w.lo_bind.vmax, w.hi_bind.minP, w.hi_bind.vmin, w.hi_bind.vmax))
            push!(rows, @sprintf("%s,%s,%d,%s,%s,%.6f,%.6f,%s,%s,%d,%.6f,%.6f,%.4f,%.6f,%.6f,%.4f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f,%.3f,%.0f",
                                 cfg.label, far, T, case, w.feasible, w.lo, w.hi, w.lo_capped, w.hi_capped, w.bad, p[1], p[2],
                                 w.lo_bind.minP, w.lo_bind.vmin, w.lo_bind.vmax, w.hi_bind.minP, w.hi_bind.vmin, w.hi_bind.vmax, zsep, K, Pa, Pb, PL, time() - tc))
        end
        merge_csv(joinpath(OUT_DIR, "large_windows.csv"), LARGE_HEADER, 3, rows)

        # Standing practice: check a dispatch against OpenDSS before trusting the instance.
        d0 = equal_angles(net)
        x = solve_case(s, d0, DAY)
        if x.ok
            rp = replay_opendss(net, x, d0; load_factor = s.lam ./ cfg.scale, pvs = s.P, pv_factor = s.F)
            say(@sprintf("    cost-optimal day at equal angles: %s, \$%.2f; OpenDSS replay |dP| %.1e kW, |dQ| %.1e kvar, |dV| %.1e pu, converged %s, PQ %s",
                         x.status, x.cost.total, rp.dP, rp.dQ, rp.dV, rp.converged, rp.pq_band))
        else
            say("    ** cost-optimal day at equal angles failed: ", x.status)
        end
        say(@sprintf("    bus %s done in %.0f s", far, time() - t0))
    end
    say("  wrote large_windows.csv")
end
foreach(large_stage, filter(in(STAGES), ["ieee2522", "large10k", "large10kd"]))

# ---- stages: three_ieee2522, three_large10kd -----------------------------------------------------------------
# Three substations on the big decks: the deck's own plus two added in memory (as in the large stages), and
# the window in the two independent angle differences along 16 directions, against the simplex predicted
# from the angle-to-import sensitivities (stage "simplex"). User, 2026-10-19: up to three substations on every
# test system. ieee2522: buses 947 and 1120 on the main trunk (both backflow-limited edges at T = 12).
# large10kD: buses 1120 (end of feeder 1) and 1060 (its middle), on one feeder so the two interact; feeders
# meet only at the stiff head bus. Edges: the maximize-angle solve with both checks on ieee2522; on
# large10kD, where that solve is unreliable, feasibility tests bracketing the predicted reach.
# THREE_DIRS=i:j:... (0-15; colons, since sbatch --export splits on commas), THREE_CASES=nobatt|cyclic and
# THREE_METHOD=ray|bracket split or override a run; LARGE_T sets T.
const THREE = Dict(
    "three_ieee2522" => (; base = "ieee2522", at = ["947", "1120"], method = :ray),
    "three_large10kd" => (; base = "large10kd", at = ["1120", "1060"], method = :bracket),
    # ieee9500_1ph has its three substations in the deck; the five normally-open switches that join
    # different feeders are closed (the four inside a feeder stay open).
    "three_ieee9500" => (; base = nothing, system = "ieee9500_1ph", at = ["sub2", "sub3"], method = :secant, feas = 1e-5,
                           ties = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]),
    # The same deck with its service transformers modeled (loads on their secondary buses, DERs unchanged): does
    # leaving the secondaries out change the window?
    "three_ieee9500xf" => (; base = nothing, system = "ieee9500xf_1ph", at = ["sub2", "sub3"], method = :secant, feas = 1e-5,
                             ties = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]))
const THREE_HEADER = "system,bus_a,bus_b,T,case,ray_deg,feasible_at_equal_angles,measured_deg,predicted_deg,limiting_substation,at_limit,failed_checks,min_import_kW,vmin_pu,vmax_pu,seconds"

"""
The margin at equal angles for `case`, solved without a cap; every later solve is then capped at 8x
its time (at least 120 s), and a capped solve is parked in parked_solves.csv with what it takes to
rerun it, instead of grinding on (a test on the infeasible side took 10x a feasible one on ieee2522).
MULTIPERIOD_CAP=0 turns the cap off.
"""
function reference_solve!(x0s, s, case, names, stage)
    haskey(x0s, case) && return x0s[case]
    delete!(ENV, "MULTIPERIOD_MAX_WALL")
    x = solve_multiperiod(s.net, Dict(k => 0.0 for k in names); problem(s, case)..., objective = :margin)
    PARK[] = (; file = joinpath(OUT_DIR, "parked_solves.csv"), stage, case = string(case))
    if get(ENV, "MULTIPERIOD_CAP", "1") == "1"
        cap = max(120.0, 8 * x.time)
        ENV["MULTIPERIOD_MAX_WALL"] = string(round(cap))
        say(@sprintf("    %s: reference solve %.0f s; later solves capped at %.0f s, capped ones parked in parked_solves.csv", case, x.time, cap))
    end
    return x0s[case] = x
end

"""
A margin solve that, if it hits the time cap, is retried once with MULTIPERIOD_RETRY (default 4) times
the cap from the same start: a capped solve means "not resolved yet", not "infeasible". A solve still
capped after that (or at once, with MULTIPERIOD_RETRY=1: park without retrying) counts in UNRESOLVED (the
search then flags its direction) and stays in parked_solves.csv for a later uncapped replay.
"""
const UNRESOLVED = Ref(0)
const CAPPED = (MOI.TIME_LIMIT, MOI.ITERATION_LIMIT)
function margin_solve(net, angles, kw; start = nothing)
    x = solve_multiperiod(net, angles; kw..., objective = :margin, start)
    x.status in CAPPED || return x
    f = parse(Float64, get(ENV, "MULTIPERIOD_RETRY", "4"))
    f <= 1 && (UNRESOLVED[] += 1; say("      capped solve: parked, direction flagged"); return x)
    cap = get(ENV, "MULTIPERIOD_MAX_WALL", "")
    isempty(cap) || (ENV["MULTIPERIOD_MAX_WALL"] = string(f * parse(Float64, cap)))
    say("      capped solve: retrying once with $(f)x the time")
    try
        x = solve_multiperiod(net, angles; kw..., objective = :margin, start)
    finally
        isempty(cap) || (ENV["MULTIPERIOD_MAX_WALL"] = cap)
    end
    x.status in CAPPED && (UNRESOLVED[] += 1; say("      still capped: unresolved"))
    return x
end

"""
The largest current over the day (A, in the deck's base: S_base / kV_base) in each class of branch, with the
element carrying it: `source` (substations, and the poi lines joining added ones), `tie` (the switches closed
for the study, `s.tie_names`, and the converter's `tie_*` lines), `service` (service transformers, `svc_*`)
and `primary` (every other line). No limit is imposed; this only keeps tabs, since the classes will need
different ratings (2026-10-26).
"""
function current_tabs(s, x)
    Ib = s.net.S_base_kVA / s.net.kV_base
    ties = get(s, :tie_names, Set{String}())
    cls(n) = startswith(n, "svc_") ? "service" : (n in ties || startswith(n, "tie_")) ? "tie" : startswith(n, "poi") ? "source" : "primary"
    best = Dict(c => (0.0, "") for c in CURRENT_CLASSES)
    for (n, r) in x.raw.lr
        m = Ib * maximum(hypot.(r, x.raw.li[n]); init = 0.0)
        c = cls(lowercase(n)); m > best[c][1] && (best[c] = (m, n))
    end
    for (n, r) in x.raw.ir
        m = Ib * maximum(hypot.(r, x.raw.ii[n]); init = 0.0)
        m > best["source"][1] && (best["source"] = (m, n))
    end
    return best
end
const CURRENT_CLASSES = ("source", "tie", "primary", "service")
const CURRENT_COLUMNS = join(["max_A_$c,at_$c" for c in CURRENT_CLASSES], ",")
current_cells(bind) = haskey(bind, :amps) ? join([@sprintf("%.2f,%s", bind.amps[c]...) for c in CURRENT_CLASSES], ",") :
                      join(fill("NaN,", length(CURRENT_CLASSES)), ",")

"""
Reach along `dir` by a secant on the import margin (`objective = :margin`, the largest import every
substation can be guaranteed): it is zero at an edge set by no-backflow and nearly linear in the angle
(the imports are), so from the margin at equal angles `x0` and one solve at the predicted reach `guess`,
secant steps reach |margin| <= MTOL_KW in a few solves (MTOL_KW / K ~ 1e-4 deg). Each solve starts from
the nearest solution: `prev` (the neighbouring direction's edge) for the first, then the last iterate.
The edge solution is itself feasible; one violation test just outside confirms that the edge is real
(a voltage limit binding first, or a margin solve failing, hands the direction to `reach_bracket`).
Returns the edge solution as `x`, to start the next direction from.
"""
function reach_secant(s, case, dir, guess, x0; prev = nothing)
    kw, nb = problem(s, case), study_banded(s)
    vtol = something(s.feas, 1e-5)
    at(r) = Dict(k => r * v for (k, v) in dir)
    n = Ref(0)
    function margin(r, st)
        x = margin_solve(s.net, at(r), kw; start = st); n[] += 1
        return x.status in RAY_OK && isfinite(x.margin_kW) && x.max_violation <= vtol, x
    end
    function test(r, st)                                   # feasible iff the margin is >= 0 within the limits
        x = solve_multiperiod(s.net, at(r); kw..., objective = :margin, start = st); n[] += 1
        return x.status in RAY_OK && isfinite(x.margin_kW) && x.margin_kW >= -MTOL_KW && x.max_violation <= vtol, x
    end
    bind_of(x) = (; minP = minimum(minimum(v) for v in values(x.P_subs_kW)),
                  vmin = minimum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T), vmax = maximum(abs(x.V_pu[b][t]) for b in nb, t in 1:x.T),
                  imports = Dict(k => minimum(v) for (k, v) in x.P_subs_kW), amps = current_tabs(s, x))
    fallback(why, st) = (say("      secant gave up (", why, "); bracketing"); b = reach_bracket(s, case, dir, guess);
                         (; b..., solves = n[], method = :bracket, x = nothing))
    x0.margin_kW >= -MTOL_KW || return (; feasible = false, r = NaN, capped = false, bind = bind_of(x0), solves = n[], method = :secant, x = nothing)
    g = isfinite(guess) ? guess : LIMIT
    if g >= LIMIT                                           # nothing binds before the search limit, if the prediction is right
        ok, x = test(LIMIT, something(prev, x0))
        ok && return (; feasible = true, r = LIMIT, capped = true, bind = bind_of(x), solves = n[], method = :secant, x)
        g = LIMIT / 2
    end
    ra, ma, rb, xb, mb = 0.0, x0.margin_kW, g, something(prev, x0), NaN
    rfeas, xfeas = 0.0, x0                                  # the farthest point known feasible, for an unresolved direction
    for _ in 1:8
        ok, xn = margin(rb, xb)
        # still capped after the retry: keep the farthest known-feasible reach, flagged, rather than bracketing on
        xn.status in CAPPED && return (; feasible = true, r = rfeas, capped = false, bind = bind_of(xfeas), solves = n[], method = :unresolved, x = nothing)
        ok || return fallback(@sprintf("margin solve at %.4f deg: %s, max violation %.1e", rb, xn.status, xn.max_violation), xb)
        xn.margin_kW >= -MTOL_KW && rb > rfeas && ((rfeas, xfeas) = (rb, xn))
        xb, mb = xn, xn.margin_kW
        abs(mb) <= MTOL_KW && break
        slope = (mb - ma) / (rb - ra)
        slope < -1e-6 || return fallback("margin does not fall with the angle", xb)
        ra, ma, rb = rb, mb, clamp(rb - mb / slope, 0.0, LIMIT)
    end
    abs(mb) <= MTOL_KW || return fallback("no convergence in 8 margin solves", xb)
    # Just outside: a margin solve (always feasible as an optimization, fast) must give a negative margin.
    _, xo = margin(rb * 1.005 + 1e-3, xb)
    xo.status in RAY_OK || return fallback("confirmation solve: " * string(xo.status) * " (parked)", xb)
    okout = xo.margin_kW >= -MTOL_KW && xo.max_violation <= vtol
    okout && return fallback(@sprintf("still feasible just outside %.4f deg", rb), xb)
    return (; feasible = true, r = rb, capped = false, bind = bind_of(xb), solves = n[], method = :secant, x = xb)
end

"Reach along `dir` by feasibility tests: start at the predicted reach, step out (x1.15) or in (x0.85) to a bracket, bisect."
function reach_bracket(s, case, dir, guess; bisections = 7)
    kw, nb = problem(s, case), study_banded(s)
    last = Ref{Any}(nothing)
    function test(r)
        # a margin solve (always feasible as an optimization): feasible here iff the margin is >= 0 within the limits
        x = margin_solve(s.net, Dict(k => r * v for (k, v) in dir), kw; start = last[])
        ok = x.status in RAY_OK && isfinite(x.margin_kW) && x.margin_kW >= -MTOL_KW && x.max_violation <= something(s.feas, 1e-5)
        ok && (last[] = x)
        return ok
    end
    bind() = last[] === nothing ? (; minP = NaN, vmin = NaN, vmax = NaN, imports = Dict{String,Float64}()) :
             (; minP = minimum(minimum(v) for v in values(last[].P_subs_kW)),
                vmin = minimum(abs(last[].V_pu[b][t]) for b in nb, t in 1:last[].T), vmax = maximum(abs(last[].V_pu[b][t]) for b in nb, t in 1:last[].T),
                imports = Dict(k => minimum(v) for (k, v) in last[].P_subs_kW), amps = current_tabs(s, last[]))
    test(0.0) || return (; feasible = false, r = NaN, capped = false, bind = bind())
    g = isfinite(guess) ? min(guess, LIMIT) : LIMIT
    if g >= LIMIT - 1e-9
        test(LIMIT) && return (; feasible = true, r = LIMIT, capped = true, bind = bind())
        g = LIMIT / 2
    end
    lo, hi = 0.0, LIMIT
    if test(g)
        lo = g
        while true
            r = min(lo * 1.15, LIMIT)
            r <= lo + 1e-9 && return (; feasible = true, r = LIMIT, capped = true, bind = bind())
            test(r) ? (lo = r) : (hi = r; break)
        end
    else
        hi = g
        while true
            r = hi * 0.85
            r < 1e-3 && (lo = 0.0; break)
            test(r) ? (lo = r; break) : (hi = r)
        end
    end
    for _ in 1:bisections
        mid = (lo + hi) / 2
        test(mid) ? (lo = mid) : (hi = mid)
    end
    test(lo)                                              # leave the edge solution as the last one, for bind()
    return (; feasible = true, r = lo, capped = false, bind = bind())
end

function three_stage(name)
    tc = THREE[name]
    tc.base === nothing && return three_deck_stage(name, tc)
    cfg = LARGE[tc.base]
    T = parse(Int, get(ENV, "LARGE_T", "24"))
    dirs_i = haskey(ENV, "THREE_DIRS") ? parse.(Int, split(ENV["THREE_DIRS"], r"[,: ]"; keepempty = false)) : collect(0:15)
    method = Symbol(get(ENV, "THREE_METHOD", string(tc.method)))     # ray | bracket
    cases = haskey(ENV, "THREE_CASES") ? Symbol.(split(ENV["THREE_CASES"], ',')) : [:nobatt, DAY]
    S0, A1, A2 = "source", "subs2", "subs3"
    rz = cfg.zs / sqrt(1 + X_OVER_R^2)
    imp = cfg.zs == 0 ? "R1=0 X1=1e-8 R0=0 X0=1e-8" : @sprintf("R1=%.8g X1=%.8g R0=%.8g X0=%.8g", rz, X_OVER_R * rz, rz, X_OVER_R * rz)
    reactance = String[]
    if cfg.rx !== nothing
        for l in eachline(joinpath(REPO_ROOT, "test_networks", cfg.system, "BranchData.dss"))
            m = match(r"^New Line\.(\S+) .*\br1=(\S+) x1=(\S+)"i, l)
            m === nothing && continue
            r1, x1 = parse(Float64, m[2]), parse(Float64, m[3])
            r1 / x1 > cfg.rx + 1e-9 && push!(reactance, @sprintf("Edit Line.%s x1=%.10g", m[1], r1 / cfg.rx))
        end
    end
    extra = [reactance;
             "New Line.poi2 Phases=1 Bus1=2s Bus2=$(tc.at[1]) r1=1e-06 x1=1e-06 c1=0 c0=0";
             "New Vsource.subs2 Phases=1 Bus1=2s basekv=$(cfg.kV) pu=1.05 angle=0 $imp";
             "New Line.poi3 Phases=1 Bus1=3s Bus2=$(tc.at[2]) r1=1e-06 x1=1e-06 c1=0 c0=0";
             "New Vsource.subs3 Phases=1 Bus1=3s basekv=$(cfg.kV) pu=1.05 angle=0 $imp";
             "Edit Vsource.source $imp"]
    net, pvs = load_sources([S0, A1, A2]; system = cfg.system, extra)
    if cfg.scale != 1
        k = cfg.scale
        net = merge(net, (; loads = [merge(d, (; kW = d.kW / k, kvar = d.kvar / k)) for d in net.loads],
                            batteries = [merge(b, (; kW_rated = b.kW_rated / k, kVA = b.kVA / k, kWh_rated = b.kWh_rated / k)) for b in net.batteries]))
        pvs = [merge(p, (; kW = p.kW / k, kVA = p.kVA / k)) for p in pvs]
    end
    s = setup(net, pvs; T, cfg.feas)
    say("\n", "="^110, "\n", uppercase(name), "  $(cfg.label), substations at buses $(tc.at[1]) (subs2) and $(tc.at[2]) (subs3) besides the deck's own, T = $T; directions $(join(dirs_i, ',')), cases $(join(cases, ','))")
    three_measure(name, cfg.label, tc.at, s, T, dirs_i, method, cases, S0, A1, A2)
end

"A deck that has its three substations already (ieee9500_1ph): close the ties, then measure as above."
function three_deck_stage(name, tc)
    T = parse(Int, get(ENV, "LARGE_T", "24"))
    dirs_i = haskey(ENV, "THREE_DIRS") ? parse.(Int, split(ENV["THREE_DIRS"], r"[,: ]"; keepempty = false)) : collect(0:15)
    method = Symbol(get(ENV, "THREE_METHOD", string(tc.method)))
    cases = haskey(ENV, "THREE_CASES") ? Symbol.(split(ENV["THREE_CASES"], ',')) : [:nobatt, DAY]
    net, pvs = load_sources(["subs1", "subs2", "subs3"]; system = tc.system, extra = ["Line.$t.enabled=true" for t in tc.ties])
    s = merge(setup(net, pvs; T, tc.feas), (; tie_names = Set(lowercase.(tc.ties))))
    sep =Dict((a.name, b.name) => abs(separation_sparse(merge(net, (; sources = [a, b])))) for a in net.sources, b in net.sources if a.name < b.name)
    say("\n", "="^110, "\n", uppercase(name), "  $(tc.system), its own three substations, ties ", join(tc.ties, ","), " closed, T = $T; directions $(join(dirs_i, ',')), cases $(join(cases, ','))")
    say("    |Z| between substations: ", join([@sprintf("%s-%s %.3f ohm", k[1], k[2], z) for (k, z) in sort(collect(sep))], ", "))
    three_measure(name, tc.system, tc.at, s, T, dirs_i, method, cases, "subs1", "subs2", "subs3")
end

function three_measure(name, label, at, s, T, dirs_i, method, cases, S0, A1, A2)
    t0 = time()
    sens = sensitivity(s, [A1, A2])
    sens === nothing && (say("  ** idle solves failed"); return)
    for k in sens.names
        say(@sprintf("    %-6s mean import %8.1f kW (worst hour %8.1f); sensitivity to subs2 %+8.1f, to subs3 %+8.1f kW/deg (daily mean)", k,
                     sum(sens.P[k]) / sens.T, minimum(sens.P[k]), sum(sens.K[(k, A1)]) / sens.T, sum(sens.K[(k, A2)]) / sens.T))
    end
    say(@sprintf("    sensitivities: %.0f s", time() - t0))
    # Hourly imports and sensitivities at equal angles, from which the predicted region is drawn
    # (THREE_DIRS= with no directions runs only this).
    open(joinpath(OUT_DIR, "three_sens_$(label)_T$(T).csv"), "w") do io
        println(io, "substation,hour,import_kW,K_subs2_kW_per_deg,K_subs3_kW_per_deg")
        for k in sens.names, t in 1:sens.T
            @printf(io, "%s,%d,%.6f,%.6f,%.6f\n", k == S0 ? "ref" : k, t, sens.P[k][t], sens.K[(k, A1)][t], sens.K[(k, A2)][t])
        end
    end
    rows = String[]
    x0s = Dict{Any,Any}(); prev = Dict{Any,Any}()             # per case: the margin at equal angles; the last edge
    for case in cases, i in dirs_i
        th = 2pi * i / 16
        dir = Dict(S0 => 0.0, A1 => cos(th), A2 => sin(th))
        pred, lim = predicted_reach(sens, Dict(A1 => cos(th), A2 => sin(th)); day = case === DAY)
        tr = time()
        if method === :ray
            y = rays(s, case, [dir]; check = :both)
            w = (; feasible = y.feasible, r = y.r[1], capped = y.capped[1], bad = y.bad, bind = y.bind[1])
        elseif method === :secant
            reference_solve!(x0s, s, case, sens.names, name)
            b = reach_secant(s, case, dir, pred, x0s[case]; prev = get(prev, case, nothing))
            b.x === nothing || (prev[case] = b.x)
            say(@sprintf("      %s: %d solves (%s)", case, b.solves, b.method))
            w = (; b.feasible, b.r, b.capped, bad = 0, b.bind)
        else
            b = reach_bracket(s, case, dir, pred)
            w = (; b.feasible, b.r, b.capped, bad = 0, b.bind)
        end
        say(@sprintf("    %-6s ray %6.1f deg: measured %s%.4f, predicted %.4f (%s)%s; min import %.2f kW, V %.4f-%.4f  (%.0f s, peak RSS %.1f GB)",
                     case, rad2deg(th), w.capped ? ">=" : "", w.r, pred, lim, w.bad == 0 ? "" : "  ** $(w.bad) failed checks",
                     w.bind.minP, w.bind.vmin, w.bind.vmax, time() - tr, Sys.maxrss() / 2^30))
        push!(rows, @sprintf("%s,%s,%s,%d,%s,%.2f,%s,%.6f,%.6f,%s,%s,%d,%.4f,%.6f,%.6f,%.0f", label, at[1], at[2], T, case, rad2deg(th),
                             w.feasible, w.r, pred, lim, w.capped, w.bad, w.bind.minP, w.bind.vmin, w.bind.vmax, time() - tr))
        merge_csv(joinpath(OUT_DIR, "three_large_currents.csv"), "system,bus_a,bus_b,T,case,ray_deg,measured_deg,at_limit,$CURRENT_COLUMNS", 6,
                  [@sprintf("%s,%s,%s,%d,%s,%.2f,%.6f,%s,%s", label, at[1], at[2], T, case, rad2deg(th), w.r, w.capped, current_cells(w.bind))])
    end
    isempty(rows) || merge_csv(joinpath(OUT_DIR, "three_large_windows.csv"), THREE_HEADER, 6, rows)
    say("  wrote three_large_windows.csv; ", name, " done in ", round(Int, time() - t0), " s")
end
foreach(three_stage, filter(in(STAGES), collect(keys(THREE))))

# ---------------------------------------------------------------------------------------------
# N substations (four here): the region in the N-1 angle differences should be a simplex with one
# planar face per substation (a tetrahedron for four), predicted from the angle-to-import
# sensitivities at equal angles. Measured along the 26 directions of a cube's faces, edges and
# corners (normalized) by the secant search, idle and cyclic; plus a direct test of the affine
# import model at random angle points inside the predicted idle region.
#   MULTI_DIRS=i:j:... (1-26), MULTI_CASES=nobatt|cyclic, MULTI_SAMPLES=n (idle linearity points, 30;
#   0 skips), LARGE_T. Rows merge into multi_windows.csv, samples into multi_linearity.csv.
# A variant master contains only switching commands after Redirect Master.dss.
# Keep the converter's topology authoritative; do not duplicate its tie list in the study.
function variant_switches(system, master)
    commands = String[]
    for line in eachline(joinpath(REPO_ROOT, "test_networks", system, master))
        text = strip(line)
        occursin(r"(?i)^(Line|Vsource)\.\S+\.enabled=(true|false)$", text) && push!(commands, String(text))
    end
    isempty(commands) && error("No switching commands in $master")
    return commands
end
const TIES_9500 = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]
const MULTI = Dict(
    "four_ieee123a" => (; system = SYSTEM, label = "ieee123_5poi_1ph", sources = ["subs1", "subs2", "subs3", "subs4"], add = nothing, feas = nothing),
    "four_ieee123b" => (; system = SYSTEM, label = "ieee123_5poi_1ph", sources = ["subs1", "subs3", "subs4", "subs5"], add = nothing, feas = nothing),
    # ieee2522 with three added substations (0.311 ohm each, as in stage ieee2522) at buses where no-backflow sets the edges.
    "four_ieee2522" => (; system = "ieee2522C_1ph", label = "ieee2522C_1ph", sources = ["source", "subs2", "subs3", "subs4"],
                         add = ["947", "996", "1120"], feas = 1e-5),
    # ieee9500 with the five ties between feeders closed and a fourth substation at l2804256, a three-phase bus
    # 6.4 / 6.8 / 5.0 ohm from the three real ones (choose_fourth_9500.jl), behind the same transformer impedance.
    "four_ieee9500" => (; system = "ieee9500_1ph", label = "ieee9500_1ph", sources = ["subs1", "subs2", "subs3", "subs4"],
                         add = ["l2804256"], feas = 1e-5, zs = complex(0.104186, 0.617339), kv = 7.1996, edit_source = false,
                         ties = TIES_9500),
    # The service-transformer edition retains the primary graph, including l2804256;
    # adding pendant secondary branches does not change its series effective distance.
    "four_ieee9500xf" => (; system = "ieee9500xf_1ph", label = "ieee9500xf_1ph", sources = ["subs1", "subs2", "subs3", "subs4"],
                         add = ["l2804256"], feas = 1e-5, zs = complex(0.104186, 0.617339), kv = 7.1996, edit_source = false, ties = TIES_9500),
    "four_smartds_p1u_typeiii" => (; system = "smartds_p1u_xf_1ph", label = "smartds_p1u_xf_typeIII", sources = ["subs1", "subs2", "subs3", "subs4"],
                         add = nothing, feas = 1e-5, variant = "Master_4sub_typeIII.dss"),
    "four_smartds_p1u_typeii" => (; system = "smartds_p1u_xf_1ph", label = "smartds_p1u_xf_typeII", sources = ["subs1", "subs2", "subs4", "subs5"],
                         add = nothing, feas = 1e-5, variant = "Master_4sub_typeII.dss"))
const MULTI_HEADER = "system,sources,T,case,direction,u1,u2,u3,measured_deg,predicted_deg,limiting_substation,at_limit,min_import_kW,vmin_pu,vmax_pu,solves,method,seconds,binding,min_imports_kW"
const CUBE_DIRS = [v ./ sqrt(sum(abs2, v)) for v in ([a, b, c] for a in -1:1, b in -1:1, c in -1:1) if any(!=(0), v)]
# n directions spread evenly over the sphere (Fibonacci lattice): MULTI_DIRSET=fib78 triples the cube's 26.
fib_dirs(n) = [(z = 1 - 2(i + 0.5) / n; r = sqrt(1 - z^2); ph = i * pi * (3 - sqrt(5)); [r * cos(ph), r * sin(ph), z]) for i in 0:n-1]
const DIRSETS = Dict("cube26" => CUBE_DIRS, "fib78" => fib_dirs(78))

"What bound at an edge: the substation whose smallest import is ~0 (< 1 kW), else `voltage` at the band edge, else `none`."
function binding_label(bind, capped)
    capped && return "open"
    isempty(bind.imports) && return "unknown"
    k = argmin(bind.imports)
    bind.imports[k] < 1.0 && return k
    (bind.vmin <= 0.9501 || bind.vmax >= 1.0499) && return "voltage"
    return "none"
end

function multi_stage(name)
    mc = MULTI[name]
    T = parse(Int, get(ENV, "LARGE_T", "24"))
    dset = get(ENV, "MULTI_DIRSET", "cube26"); DIRS = DIRSETS[dset]
    di = haskey(ENV, "MULTI_DIRS") ? parse.(Int, split(ENV["MULTI_DIRS"], r"[,: ]"; keepempty = false)) : collect(1:length(DIRS))
    cases = haskey(ENV, "MULTI_CASES") ? Symbol.(split(ENV["MULTI_CASES"], ',')) : [:nobatt, DAY]
    nsamp = parse(Int, get(ENV, "MULTI_SAMPLES", "30"))
    extra = haskey(mc, :variant) ? variant_switches(mc.system, mc.variant) : String[]
    # MULTI_VSRC=v: every source (the deck's own and the added ones) at v pu instead of 1.05; outputs
    # carry the tag _v<v>, so they never mix with the 1.05 pu rows
    vs = get(ENV, "MULTI_VSRC", "1.05"); vtag = haskey(ENV, "MULTI_VSRC") ? "_v$vs" : ""
    if mc.add !== nothing
        # added substations take the names after the deck's own (subs2.. on ieee2522, subs4 on ieee9500)
        zs = get(mc, :zs, complex(0.311 / sqrt(1 + X_OVER_R^2), X_OVER_R * 0.311 / sqrt(1 + X_OVER_R^2)))
        imp = @sprintf("R1=%.8g X1=%.8g R0=%.8g X0=%.8g", real(zs), imag(zs), real(zs), imag(zs))
        off = length(mc.sources) - length(mc.add)
        for (i, b) in enumerate(mc.add)
            push!(extra, "New Line.poi$(i+off) Phases=1 Bus1=$(i+off)s Bus2=$b r1=1e-06 x1=1e-06 c1=0 c0=0",
                         "New Vsource.subs$(i+off) Phases=1 Bus1=$(i+off)s basekv=$(get(mc, :kv, 7.2)) pu=$vs angle=0 $imp")
        end
        get(mc, :edit_source, true) && push!(extra, "Edit Vsource.source $imp")
    end
    haskey(mc, :ties) && append!(extra, ["Line.$t.enabled=true" for t in mc.ties])
    own = mc.add === nothing ? mc.sources : mc.sources[1:end-length(mc.add)]
    isempty(vtag) || append!(extra, ["Edit Vsource.$k pu=$vs" for k in own])
    net, pvs = load_sources(mc.sources; system = mc.system, extra)
    isempty(vtag) || say("    sources at ", join(sort(unique([@sprintf("%.3f", x.pu) for x in net.sources])), ", "), " pu (MULTI_VSRC)")
    s = setup(net, pvs; T, feas = mc.feas)
    # the switches closed for the study (ties, and a variant master's bus tie), for the current tabs
    s = merge(s, (; tie_names = Set(lowercase(m[1]) for c in extra for m in (match(r"(?i)^Line\.(\S+)\.enabled=true$", c),) if m !== nothing)))
    S0, movers =mc.sources[1], mc.sources[2:end]
    say("\n", "="^110, "\n", uppercase(name), "  $(mc.label), sources ", join(mc.sources, ","), mc.add === nothing ? "" : " (added at buses $(join(mc.add, ",")))",
        ", reference $S0, T = $T; directions $(join(di, ',')), cases $(join(cases, ','))")
    t0 = time()
    say("    voltage band on ", length(study_banded(s)), " buses", s.band_buses === nothing ? " (primary-only legacy deck)" : " (customers only)")
    sens = sensitivity(s, movers; h = parse(Float64, get(ENV, "MULTI_SENS_H", "0.1")),
                       replay_file = joinpath(OUT_DIR, "multi_replay_$(name)$(vtag)_T$(T).csv"))
    sens === nothing && error("Idle sensitivity solves failed for $name; no prediction available")
    for k in sens.names
        say(@sprintf("    %-6s mean import %8.1f kW (worst hour %8.1f); sensitivity %s kW/deg (daily mean)", k, sum(sens.P[k]) / T, minimum(sens.P[k]),
                     join([@sprintf("%s %+.1f", m, sum(sens.K[(k, m)]) / T) for m in movers], ", ")))
    end
    colsum = [sum(sum(sens.K[(k, m)]) / T for k in sens.names) for m in movers]
    say("    sensitivities summed over substations (0 if raising an angle only shifts power): ", join([@sprintf("%+.1f", x) for x in colsum], ", "), " kW/deg")
    open(joinpath(OUT_DIR, "multi_sens_$(name)$(vtag)_T$(T).csv"), "w") do io
        println(io, "substation,hour,import_kW," * join(["K_$(m)_kW_per_deg" for m in movers], ","))
        for k in sens.names, t in 1:T
            println(io, k, ",", t, ",", @sprintf("%.6f", sens.P[k][t]), ",", join([@sprintf("%.6f", sens.K[(k, m)][t]) for m in movers], ","))
        end
    end
    # Predict all directions even when running only the sensitivity stage. These are
    # backflow-based estimates; voltage limits and redispatch can move the actual edge.
    open(joinpath(OUT_DIR, "multi_predicted_$(name)$(vtag)_$(dset)_T$(T).csv"), "w") do io
        println(io, "case,direction,u1,u2,u3,predicted_deg,limiting_substation")
        for case in cases, (i, u) in enumerate(DIRS)
            r, lim = predicted_reach(sens, Dict(zip(movers, u)); day = case === DAY)
            println(io, join([case, i, u..., r, lim], ','))
        end
    end
    get(ENV, "MULTI_SENS_ONLY", "0") == "1" && (say("  sensitivity and predicted reaches complete; edge search skipped (MULTI_SENS_ONLY)"); return)
    # the affine import model at random points inside the predicted idle region
    if nsamp > 0 && :nobatt in cases
        rng = Random.MersenneTwister(2027); lin = String[]; errs = Float64[]
        for j in 1:nsamp
            u = randn(rng, length(movers)); u ./= sqrt(sum(abs2, u))
            r = first(predicted_reach(sens, Dict(zip(movers, u)); day = false)) * (0.2 + 0.7 * rand(rng))
            d = Dict(zip(movers, r .* u))
            x = solve_multiperiod(s.net, merge(Dict(S0 => 0.0), d); problem(s, :nobatt)..., no_backflow = false)
            x.ok || continue
            e = maximum(abs(x.P_subs_kW[k][t] - (sens.P[k][t] + sum(sens.K[(k, m)][t] * d[m] for m in movers))) for k in sens.names, t in 1:T)
            sc = maximum(maximum(sens.P[k]) for k in sens.names)
            push!(errs, e / sc)
            push!(lin, @sprintf("%s,%s,%d,%d,%s,%.6f,%.4f,%.6f", mc.label, join(mc.sources, ":"), T, j, join([@sprintf("%.6f", d[m]) for m in movers], ":"), r, e, e / sc))
        end
        if !isempty(errs)
            say(@sprintf("    affine import model at %d random idle points: largest hourly error %.2f%% of the largest import (median %.2f%%)",
                         length(errs), 100maximum(errs), 100sort(errs)[(end + 1) ÷ 2]))
            merge_csv(joinpath(OUT_DIR, "multi_linearity$(vtag).csv"), "system,sources,T,sample,angles_deg,radius_deg,max_abs_error_kW,max_rel_error", 4, lin)
        end
    end
    rows = String[]; x0s = Dict{Any,Any}(); prev = Dict{Any,Any}()
    fname = (dset == "cube26" ? "multi_windows" : "multi_windows_$dset") * vtag * ".csv"
    for case in cases, i in di
        u = DIRS[i]
        pred, lim = predicted_reach(sens, Dict(zip(movers, u)); day = case === DAY)
        dir = merge(Dict(S0 => 0.0), Dict(zip(movers, u)))
        tr = time()
        reference_solve!(x0s, s, case, sens.names, name)
        u0 = UNRESOLVED[]
        b = reach_secant(s, case, dir, pred, x0s[case]; prev = get(prev, case, nothing))
        b.x === nothing || (prev[case] = b.x)
        # a direction with a solve still capped after its retry: its edge is a lower bound, flagged in the method column
        nu = UNRESOLVED[] - u0
        nu > 0 && (b = merge(b, (; method = Symbol(string(b.method) * "+unresolved$nu"))))
        bl = binding_label(b.bind, b.capped)
        say(@sprintf("    %-6s dir %2d (%+.2f,%+.2f,%+.2f): measured %s%.4f, predicted %.4f (%s), bound by %s; min import %.2f kW, V %.4f-%.4f; %d solves (%s, %.0f s)",
                     case, i, u..., b.capped ? ">=" : "", b.r, pred, lim, bl, b.bind.minP, b.bind.vmin, b.bind.vmax, b.solves, b.method, time() - tr))
        push!(rows, @sprintf("%s,%s,%d,%s,%d,%.4f,%.4f,%.4f,%.6f,%.6f,%s,%s,%.4f,%.6f,%.6f,%d,%s,%.0f,%s,%s", mc.label, join(mc.sources, ":"), T, case, i, u...,
                             b.r, pred, lim, b.capped, b.bind.minP, b.bind.vmin, b.bind.vmax, b.solves, b.method, time() - tr, bl,
                             join([@sprintf("%s=%.3f", k, v) for (k, v) in sort(collect(b.bind.imports))], ";")))
        merge_csv(joinpath(OUT_DIR, fname), MULTI_HEADER, 5, rows[end:end])     # each row as soon as it is measured
        # largest currents at the edge solution, by branch class (kept apart so multi_windows keeps its columns)
        merge_csv(joinpath(OUT_DIR, replace(fname, "multi_windows" => "multi_currents")), "system,sources,T,case,direction,measured_deg,at_limit,$CURRENT_COLUMNS", 5,
                  [@sprintf("%s,%s,%d,%s,%d,%.6f,%s,%s", mc.label, join(mc.sources, ":"), T, case, i, b.r, b.capped, current_cells(b.bind))])
    end
    # (rows were written one by one above)
    say("  wrote $fname; ", name, " done in ", round(Int, time() - t0), " s")
end
foreach(multi_stage, filter(in(STAGES), sort(collect(keys(MULTI)))))

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
    # Tables report the hard terminal SOC (DAY) only; the soft diagnostic stays in the CSVs.
    day = string(DAY)
    cols = Dict("nobatt" => "No batteries", "hourly" => "Batt., hourly", day => "Batt., whole day")   # must fit one IEEE column
    three = ("nobatt", day)
    open(joinpath(OUT_DIR, "tables.tex"), "w") do io
        println(io, "% Generated by problems/msopf/root_level/multiperiod_study.jl -- do not edit by hand.")
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
        table("main_horizon_windows.csv", 1, "Feasible window of a constant \$\\delta_{s_3}-\\delta_{s_4}\$ (deg) over a 24\\,h day, ieee123. ``Hourly'': each hour solved alone, batteries without memory.", "tab:mp-main", "PV", k -> pvname[k[1]]; cases = ("nobatt", "hourly", day))
        table("horizon_windows.csv", 1, "The same window versus horizon length \$T\$ (\$\\Delta t = 24/T\$\\,h).", "tab:mp-horizon", "\$T\$", k -> k[1])
        table("zs_windows.csv", 1, "The same window versus source impedance \$|Z_s|\$ (\$\\Omega\$), \$T=24\$.", "tab:mp-zs", "\$|Z_s|\$", k -> k[1] == "0.0" ? "0 (stiff)" : (k[1] == "0.8" ? "0.8 (typical)" : k[1]))
        table("duration_windows.csv", 1, "The same window versus battery duration (kWh per kW), \$T=24\$.", "tab:mp-duration", "Hours", k -> k[1])
        table("pairs_windows.csv", 3, "Windows for every substation pair, \$T=24\$; \$|Z|\$ is the feeder impedance between the pair.", "tab:mp-pairs", "Pair (\$|Z|\$, \$\\Omega\$)", k -> @sprintf("%s--%s (%.2f)", replace(k[1], "subs" => ""), replace(k[2], "subs" => ""), parse(Float64, k[3])))
        table("five_windows.csv", 1, "All five substations: window of one substation's angle against the other four (deg), \$T=24\$.", "tab:mp-five", "Moved", k -> "Subs " * replace(k[1], "subs" => ""))
        est = readcsv("bound_estimate.csv")
        if est !== nothing
            iv(lo, hi) = @sprintf("\$[%.3f, %+.3f]\$", parse(Float64, lo), parse(Float64, hi))
            short = Dict("off" => "Off", "unity" => "Unity p.f.", "qfree" => "\$q\$ free")
            println(io, "\\begin{table}[!t]\n\\centering\n\\caption{Constant circulating power as an estimate of the window (deg), \$T=24\$: the worst hour without batteries, the daily mean with them.}\n\\label{tab:mp-estimate}\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{llll}\n\\hline")
            println(io, "PV & Batteries & Estimate & Measured \\\\\n\\hline")
            for r in est
                println(io, short[r[1]], " & ", r[2] == "nobatt" ? "No" : "Whole day", " & ", iv(r[3], r[4]), " & ", iv(r[5], r[6]), " \\\\")
            end
            println(io, "\\hline\n\\end{tabular}\n\\end{table}\n")
        end
        big = readcsv("large_windows.csv")              # columns: see LARGE_HEADER
        for name in ("ieee2522", "large10k", "large10kd")
            big === nothing && break
            cfg = LARGE[name]
            rows = [r for r in big if r[1] == cfg.label && r[3] == string(cfg.T)]
            isempty(rows) && continue
            title = cfg.scale == 1 ? name : @sprintf("%s (loads, PV and batteries divided by %g%s)", replace(cfg.label, "_1ph" => ""), cfg.scale, cfg.rx === nothing ? "" : @sprintf("; line \$R/X=%g\$", cfg.rx))
            println(io, "\\begin{table}[!t]\n\\centering\n\\caption{$title with a second substation at the given bus (\$|Z|\$ to the first, \$\\Omega\$): window of the new substation's angle against the first (deg), \$T=$(cfg.T)\$. Estimate: the daily-mean rule.}\n\\label{tab:mp-$name}\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{llll}\n\\hline")
            println(io, "Bus (\$|Z|\$) & No batteries & Batt., whole day & Estimate \\\\\n\\hline")
            for bus in unique(String(r[2]) for r in rows)
                pick(case) = (m = [r for r in rows if r[2] == bus && r[4] == case]; isempty(m) ? nothing : m[end])
                nb, dy = pick("nobatt"), pick(day)
                c(r) = r === nothing ? "--" : cell(r[6], r[7], r[5], r[8], r[9])
                est = dy === nothing ? "--" : @sprintf("\$[%.2f, %+.2f]\$", parse(Float64, dy[11]), parse(Float64, dy[12]))
                println(io, @sprintf("%s (%.1f)", bus, parse(Float64, something(nb, dy)[19])), " & ", c(nb), " & ", c(dy), " & ", est, " \\\\")
            end
            println(io, "\\hline\n\\end{tabular}\n\\end{table}\n")
        end
    end
    say("\nwrote ", joinpath(OUT_DIR, "tables.tex"))
end
