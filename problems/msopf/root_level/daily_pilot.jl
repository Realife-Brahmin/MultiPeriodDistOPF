# Daily certificate pilot: T = 24, batteries dispatched, hard cyclic SOC (2026-10-22).
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/daily_pilot.jl
# IEEE 123-node, sources 1, 3, 4, PV reactive power dispatched, batteries dispatched (lossless, |P_B| <= rating,
# SOC in [0.30, 0.95], B[0] = 0.625, B[24] = B[0]); no current caps or ampacity anywhere. Per direction
# (292.5 and 157.5 deg):
#   * the locally measured daily reach r*: secant on the validated exact margin, best of a cold and a warm start
#     at every iterate;
#   * source-current bounds per source and hour from current_bounds.jl, with the battery power range in the device
#     rectangles (and, for comparison only, the idle-battery bounds of the same hour);
#   * the full coupled 24-hour restoration relaxation (band, device limits, SOC dynamics and cyclic SOC hard;
#     min rho with P_s,t >= -rho for every source and hour; rho <= rho_exact + 1 kW, which keeps every rho = 0
#     candidate), both nodal strengthenings, the derived bounds as lifted source-current constraints;
#   * 0.98 and 1.02 r*, and 1.10 r* only if 1.02 r* stays unresolved.
# Writes daily_pilot_{reach,bounds,queries}.csv, certificates/daily_*.txt and daily_pilot.log.
ENV["SCB_LIB"] = "1"
include(joinpath(@__DIR__, "source_current_bounds.jl"))
include(joinpath(@__DIR__, "..", "current_bounds.jl"))
include(joinpath(@__DIR__, "verify_certificate.jl"))

const DP = joinpath(OUT_DIR, "daily_pilot")
const CERT_DIR = joinpath(OUT_DIR, "certificates")
dp_csv(file, row) = open(io -> println(io, join(row, ',')), DP * "_" * file * ".csv", "a")

function daily_margin(kw, r, th, warm)
    ang = at(r, th); best = (-Inf, nothing, NaN, NaN)
    vals = Float64[]
    for st in (warm === nothing ? (nothing,) : (nothing, warm))
        x = solve_multiperiod(net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0, start = st)
        v = validate_dispatch(net, ang, kw, x; k = KV)
        L = x.status in RAY_OK && v.equations <= 1e-6 && v.viol["soc"] <= 1e-6 ? v.margin_kW : -Inf
        push!(vals, L)
        L > best[1] && (best = (L, x, v.equations, v.viol["soc"]))
    end
    return best..., vals
end

cyclic_res(x) = maximum(abs(x.B_kWh[k][end] - x.B_kWh[k][1]) for k in keys(x.B_kWh))

function dp_main()
    mkpath(CERT_DIR)
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw24 = problem(s24, :cyclic)
    kw24.battery === :optimize && kw24.terminal === :cyclic || error("expected dispatched batteries and cyclic SOC")
    T = kw24.T; S = net.S_base_kVA
    for f in ("reach", "bounds", "queries"); rm(DP * "_" * f * ".csv"; force = true); end
    dp_csv("reach", ["theta_deg", "r_deg", "margin_cold_kW", "margin_warm_kW"])
    dp_csv("bounds", ["theta_deg", "point", "hour", "source", "bound_dispatched_A", "bound_idle_A", "termwise_A", "no_load_A", "full_sectors", "device_buses", "witness_A"])
    dp_csv("queries", ["theta_deg", "point", "r_deg", "d3_deg", "d4_deg", "r_star_deg", "L_kW", "equations", "soc_residual", "cyclic_residual_kWh", "rho_exact_kW", "rho_max_kW",
                       "bound_s", "build_s", "solve_s", "verify_s", "relax_status", "dual_status", "rho_relax_kW", "rho_lb_kW", "box_term_kW", "rounding_kW", "primal_residual",
                       "dual_psd_violation", "dual_sign_violation", "farkas_bound", "lifted_worst", "lifted_at", "lifted_box_worst", "exact_lb_kW", "exact_verify_s", "class"])
    tp = time(); RM = response_matrix(); tpre = time() - tp
    say(@sprintf("reusable preprocessing: network response (BigFloat 256) %.1f s", tpre))
    # SOC coupling: can it shrink any hour's battery power range below the rating?
    soc = get(kw24, :soc, (0.30, 0.95)); soc0 = get(kw24, :soc0, 0.625)
    hrs = minimum(min(soc0 - soc[1], soc[2] - soc0) * bt.kWh_rated / bt.kW_rated for bt in net.batteries)
    say(@sprintf("SOC coupling: every battery can run at full rating for at least %.2f h in either direction from B[0] (dt = %.1f h), and the SOC window is %.2f-%.2f; no hour's power range is tightened by SOC, so the per-hour ranges [-P_rated, P_rated] are used",
                 hrs, get(kw24, :dt_h, 24 / T), soc...))
    cw = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "certified_window.csv"))[2:end]]
    for th in (292.5, 157.5)
        rh = parse(Float64, only(r for r in cw if parse(Float64, r[1]) == th && r[2] == "sources+lines 1.2x" && r[4] == "")[5])
        # daily reach: secant on the best validated margin
        r0, r1 = rh, 2rh
        m0, x0, _, _, v0 = daily_margin(kw24, r0, th, nothing); dp_csv("reach", [th, r0, v0[1], length(v0) > 1 ? v0[2] : NaN])
        m1, x1, _, _, v1 = daily_margin(kw24, r1, th, x0); dp_csv("reach", [th, r1, v1[1], length(v1) > 1 ? v1[2] : NaN])
        for _ in 1:15
            (isfinite(m0) && isfinite(m1) && m1 != m0) || break
            r2 = max(r1 - m1 * (r1 - r0) / (m1 - m0), 1e-4)
            r0, m0, x0 = r1, m1, x1
            r1 = r2; m1, x1, _, _, v1 = daily_margin(kw24, r1, th, x0); dp_csv("reach", [th, r1, v1[1], length(v1) > 1 ? v1[2] : NaN])
            abs(m1) < 1e-4 && break
        end
        rs = r1
        say(@sprintf("\n== theta %.1f deg: daily reach r* = %.6f deg (validated margin %+.1e kW; hour-12 reach %.6f deg)", th, rs, m1, rh))
        facs = [0.98, 1.02]; k = 0
        while k < length(facs)
            k += 1; fac = facs[k]; r = fac * rs; ang = at(r, th); pname = @sprintf("%.2f", fac)
            L, xm, eqm, socm, _ = daily_margin(kw24, r, th, x1)
            xr = solve_multiperiod(net, ang; kw24..., objective = :margin, margin_band_hard = true, margin_reg = 0.0, start = xm)
            vr = validate_dispatch(net, ang, kw24, xr; k = KV)
            okr = xr.status in RAY_OK && vr.equations <= 1e-6 && vr.viol["band"] <= 1e-6 && vr.viol["soc"] <= 1e-6
            rho_ex = okr ? max(0.0, -vr.min_import_kW) : NaN; rub = (okr ? rho_ex : 50.0) + 1.0
            say(@sprintf("  -- %s r* (r = %.6f deg; delta3 %+.6f, delta4 %+.6f): exact margin L %+.4f kW (equations %.1e, SOC %.1e, cyclic %.1e kWh); exact restoration rho %.4f kW (validated %s; SOC %.1e, cyclic %.1e kWh); rho <= %.4f kW",
                         fac, r, ang["subs3"], ang["subs4"], L, eqm, socm, xm === nothing ? NaN : cyclic_res(xm), rho_ex, okr, vr.viol["soc"], cyclic_res(xr), rub))
            # derived bounds, per hour
            tb = time()
            CBs = [current_bounds(RM, ang, device_rectangles(net, kw24, t), t; ndir = 3600, ndir_rho = 720) for t in 1:T]
            tbound = time() - tb
            kwI = (; kw24..., battery = :idle)
            CIs = [current_bounds(RM, ang, device_rectangles(net, kwI, t), t; ndir = 3600, ndir_rho = 720) for t in 1:T]
            bnd = Dict(s.name => [CBs[t].out[s.name].sect for t in 1:T] for s in net.sources)
            wit = xm === nothing ? nothing : xm
            for t in 1:T, s in net.sources
                w = wit === nothing ? NaN : abs(complex(wit.raw.ir[s.name][t], wit.raw.ii[s.name][t]))
                dp_csv("bounds", [th, pname, t, s.name, @sprintf("%.3f", bnd[s.name][t] * IbA), @sprintf("%.3f", CIs[t].out[s.name].sect * IbA), @sprintf("%.3f", CBs[t].out[s.name].term * IbA),
                                  @sprintf("%.3f", CBs[t].out[s.name].a0 * IbA), CBs[t].full_sectors, CBs[t].nd, @sprintf("%.3f", w * IbA)])
                isfinite(w) && w > bnd[s.name][t] && say("  !! witness current above its derived bound: ", s.name, " hour ", t)
            end
            for s in net.sources
                b = bnd[s.name] .* IbA; bi = [CIs[t].out[s.name].sect * IbA for t in 1:T]
                say(@sprintf("     %s bounds by hour (A): %s | idle-battery bounds %.1f-%.1f A | dispatched/idle up to %.2fx (hour %d)", s.name, join([@sprintf("%.0f", x) for x in b], " "),
                             minimum(bi), maximum(bi), maximum(b ./ bi), argmax(b ./ bi)))
            end
            say(@sprintf("     bound construction %.1f s (24 hours); sectors unrestricted in phase (device set contains zero power): %s of %d device buses by hour",
                         tbound, join(unique([CBs[t].full_sectors for t in 1:T]), "/"), CBs[1].nd))
            # the coupled restoration relaxation
            kwb = (; kw24..., source_ampacity_pu = bnd)
            tw = time()
            c = certify(net, ang, kwb; k = KV, cuts = CUTS, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, tighten = true, src_scale = true)   # rows scaled (certificate_audit)
            infeas = c.status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
            lb = infeas ? NaN : (isfinite(c.db.bound) ? c.db.bound * S : -Inf)
            cls = infeas ? (c.farkas.bound > 0 ? "no-backflow day impossible (checked Farkas)" : "solver-only infeasibility (ray not validated)") :
                  lb > 0 ? "no-backflow day impossible (checked lower bound on rho > 0)" : "unresolved (rho lower bound <= 0)"
            lv = (; worst = NaN, at = "", box_worst = NaN, box_at = "")
            if isfinite(L) && L >= 0 && xm !== nothing
                lv = lifted_violation(c.R, net, xm, kwb; margin_kW = 0.0)
                say(@sprintf("     lifted %s r* daily dispatch (rho = 0): largest constraint violation %.2e (%s); largest box violation %.2e (%s)", pname, lv.worst, lv.at, lv.box_worst, lv.box_at))
            end
            exact_lb, tex = NaN, NaN
            if occursin("impossible", cls)
                fname = joinpath(CERT_DIR, @sprintf("daily_restoration_theta%.1f_%s.txt", th, pname))
                export_certificate(c.R.model, c.R.box, fname; c = infeas ? zero(AffExpr) : objective_function(c.R.model),
                                   note = @sprintf("IEEE123 T=24, batteries dispatched, cyclic SOC, theta %.1f deg, r = %.6f deg (%s r*), restoration", th, r, pname))
                ve = verify_file(fname); exact_lb = Float64(BigFloat(ve.bound)) * S; tex = ve.seconds
                say(@sprintf("     exact re-check of the exported certificate: %.6f kW (widened boxes %.6f kW), %d of %d PSD duals shifted, %.1f s",
                             exact_lb, Float64(BigFloat(ve.bound_widened)) * S, ve.nshift, ve.npsd, ve.seconds))
            end
            say(@sprintf("     relaxation (%d PSD blocks, %d variables): %s / %s, rho relax %.4f kW, checked lower bound %.4f kW%s (box term %.2e kW, rounding %.1e kW, residual %.1e, dual psd/sign %.1e/%.1e); build %.1f, solve %.1f, verify %.1f s -> %s",
                         c.nblocks, c.nvars, c.status, c.dual_status, infeas ? NaN : -c.primal, lb, infeas ? @sprintf(", Farkas %.3e", c.farkas.bound) : "", c.box_term_kW, c.rounding_kW,
                         c.primal_residual, c.psd_dual_violation, c.dual_sign_violation, c.build_s, c.solve_s, c.verify_s, cls))
            dp_csv("queries", [th, pname, @sprintf("%.6f", r), @sprintf("%.6f", ang["subs3"]), @sprintf("%.6f", ang["subs4"]), @sprintf("%.6f", rs), @sprintf("%.4f", L), @sprintf("%.1e", eqm),
                               @sprintf("%.1e", socm), @sprintf("%.1e", xm === nothing ? NaN : cyclic_res(xm)), @sprintf("%.4f", rho_ex), @sprintf("%.4f", rub), @sprintf("%.1f", tbound),
                               @sprintf("%.1f", c.build_s), @sprintf("%.1f", c.solve_s), @sprintf("%.1f", c.verify_s), string(c.status), string(c.dual_status), @sprintf("%.4f", infeas ? NaN : -c.primal),
                               @sprintf("%.4f", lb), @sprintf("%.3e", c.box_term_kW), @sprintf("%.1e", c.rounding_kW), @sprintf("%.1e", c.primal_residual), @sprintf("%.1e", c.psd_dual_violation),
                               @sprintf("%.1e", c.dual_sign_violation), @sprintf("%.3e", c.farkas.bound), @sprintf("%.2e", lv.worst), replace(string(lv.at), ',' => ';'), @sprintf("%.2e", lv.box_worst),
                               @sprintf("%.6f", exact_lb), @sprintf("%.1f", tex), cls])
            fac == 1.02 && !occursin("impossible", cls) && push!(facs, 1.10)
        end
    end
end
dp_main()
