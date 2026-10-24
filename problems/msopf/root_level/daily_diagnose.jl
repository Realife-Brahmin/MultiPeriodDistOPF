# Where the daily certificate loses its margin (2026-10-22). DIAGNOSTICS, NOT CERTIFICATES where marked.
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/daily_diagnose.jl
# The unresolved daily queries of daily_pilot.jl (1.02 r*, T = 24, batteries dispatched, cyclic SOC). Four
# relaxations of the restoration problem at the same angle vector, all with both nodal strengthenings and scaled rows:
#   V  valid:      derived bounds for dispatched batteries, batteries free (the pilot's run);
#   D1 diagnostic: the idle-battery bounds of each hour imposed although batteries move (NOT valid for this model):
#                  what the wider device-current sets cost;
#   D2 diagnostic: batteries fixed to the validated exact restoration schedule, dispatched-battery bounds:
#                  what battery flexibility costs inside the relaxation;
#   D3 diagnostic: batteries fixed to that schedule, bounds recomputed for that fixed schedule (valid for the
#                  fixed-schedule model): the per-hour network relaxation alone.
# The fixed-schedule model's exact optimum equals the free model's (the schedule is the free model's optimum).
ENV["SCB_LIB"] = "1"
include(joinpath(@__DIR__, "source_current_bounds.jl"))
include(joinpath(@__DIR__, "..", "current_bounds.jl"))

const DD = joinpath(OUT_DIR, "daily_diagnose")

function relax_case(ang, kw, bnd, rub; fixed = nothing)
    S = net.S_base_kVA; t0 = time()
    kwb = (; kw..., source_ampacity_pu = bnd)
    R = relax_margin(net, ang, kwb; k = KV, cuts = CUTS, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, solver_attributes = CLARABEL_DEFAULTS, src_scale = true)
    if fixed !== nothing
        for kk in axes(R.P_B, 1), t in 1:kw.T; fix(R.P_B[kk, t], fixed[kk][t] / S; force = true); end
    end
    optimize!(R.model); ts = time() - t0
    tighten_boxes!(R, net, kwb, ang)
    st = termination_status(R.model)
    db = has_duals(R.model) ? dual_bound(R.model; box = R.box) : (; bound = -Inf, box_term = NaN)
    rho = has_values(R.model) ? -value(R.m) * S : NaN
    # fictitious losses: line losses booked beyond what the relaxation's own currents imply
    fict = NaN
    if has_values(R.model)
        fict = 0.0
        for l in net.lines, t in 1:kw.T
            Rl = real(-net.Y_base / l.Y[1, 2]); g(sy) = value(R.Lv[sy][l.name, t])
            fict += Rl * (g(:lrlr) + g(:lili) - g(:lr)^2 - g(:li)^2)
        end
        fict *= S
    end
    return (; st, rho, lb = db.bound * S, box = db.box_term * S, fict, ts)
end

function dd_main()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw24 = problem(s24, :cyclic); T = kw24.T
    q = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "daily_pilot_queries.csv"))[2:end]]
    RM = response_matrix()
    open(io -> println(io, "theta_deg,point,case,valid,relax_status,rho_relax_kW,rho_lb_kW,box_term_kW,fictitious_losses_kWh,solve_s,rho_exact_kW"), DD * ".csv", "w")
    for row in q
        (row[2] == "1.02" && occursin("unresolved", row[end])) || continue
        th = parse(Float64, row[1]); r = parse(Float64, row[3]); ang = at(r, th)
        xr = solve_multiperiod(net, ang; kw24..., objective = :margin, margin_band_hard = true, margin_reg = 0.0)
        vr = validate_dispatch(net, ang, kw24, xr; k = KV); rho_ex = max(0.0, -vr.min_import_kW); rub = rho_ex + 1.0
        say(@sprintf("\n== theta %.1f, 1.02 r* (r = %.6f deg): exact restoration rho %.4f kW (equations %.1e, SOC %.1e)", th, r, rho_ex, vr.equations, vr.viol["soc"]))
        bd = Dict(s.name => [current_bounds(RM, ang, device_rectangles(net, kw24, t), t; ndir_rho = 720).out[s.name].sect for t in 1:T] for s in net.sources)
        kwI = (; kw24..., battery = :idle)
        bi = Dict(s.name => [current_bounds(RM, ang, device_rectangles(net, kwI, t), t; ndir_rho = 720).out[s.name].sect for t in 1:T] for s in net.sources)
        # fixed-schedule bounds: each battery bus's power is its fixed schedule value (a degenerate rectangle)
        Sb = net.S_base_kVA
        bf = Dict(s.name => Float64[] for s in net.sources)
        for t in 1:T
            rects = device_rectangles(net, kwI, t)
            add = Dict(b => 0.0 for b in net.buses)
            for bt in net.batteries; add[bt.bus] += xr.P_B_kW[bt.name][t] / Sb; end
            seen = Set(rr[1] for rr in rects)
            rects2 = [(rr[1], rr[2] + add[rr[1]], rr[3] + add[rr[1]], rr[4], rr[5]) for rr in rects]
            for b in net.buses
                (add[b] != 0 && !(b in seen)) && push!(rects2, (b, add[b], add[b], 0.0, 0.0))
            end
            cb = current_bounds(RM, ang, rects2, t; ndir_rho = 720)
            for s in net.sources; push!(bf[s.name], cb.out[s.name].sect); end
        end
        for (name, valid, bnd, fx) in (("V: dispatched bounds, batteries free", "yes", bd, nothing), ("D1: idle-battery bounds, batteries free", "NO (diagnostic)", bi, nothing),
                                       ("D2: dispatched bounds, batteries fixed to the witness schedule", "fixed-schedule model", bd, xr.P_B_kW),
                                       ("D3: fixed-schedule bounds, batteries fixed", "fixed-schedule model", bf, xr.P_B_kW))
            o = relax_case(ang, kw24, bnd, rub; fixed = fx)
            say(@sprintf("  %-62s valid: %-20s %s, rho relax %8.4f kW, checked lower bound %8.4f kW (box term %.2e), fictitious losses %.2f kWh, %.1f s",
                         name, valid, o.st, o.rho, o.lb, o.box, o.fict, o.ts))
            open(io -> println(io, join([th, "1.02", "\"" * name * "\"", valid, o.st, @sprintf("%.4f", o.rho), @sprintf("%.4f", o.lb), @sprintf("%.3e", o.box), @sprintf("%.3f", o.fict),
                                         @sprintf("%.1f", o.ts), @sprintf("%.4f", rho_ex)], ',')), DD * ".csv", "a")
        end
    end
end
dd_main()
