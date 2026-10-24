# Daily pilot with the refined battery-bus enclosure (2026-10-22).
#   julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_refined.jl
# Same queries and settings as daily_pilot.jl (T = 24, batteries dispatched, cyclic SOC, no caps), one change:
# current_bounds(...; corners = true). For a fixed V the device current J = conj(S) V / |V|^2 is a rotated, scaled copy
# of the device power rectangle, so its support function is a max over the rectangle's corners (corner_support! in
# current_bounds.jl), instead of a sector or, when the rectangle contains zero power (load + PV + battery buses), the
# full disc of radius Smax / Vmin the pilot used. Still an outer enclosure of every allowed state. Checks: the hour-12
# idle-battery corner bounds against every current the falsification search found (source_current_bounds_falsify.csv);
# the 0.98 r* daily witnesses against the refined bounds, and lifted into the relaxation.
ENV["SCB_LIB"] = "1"
include(joinpath(@__DIR__, "source_current_bounds.jl"))
include(joinpath(@__DIR__, "..", "current_bounds.jl"))
include(joinpath(@__DIR__, "verify_certificate.jl"))

const DR = joinpath(OUT_DIR, "daily_refined")
const CERT_DIR = joinpath(OUT_DIR, "certificates")

function dr_main()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw24 = problem(s24, :cyclic); T = kw24.T; S = net.S_base_kVA
    q = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "daily_pilot_queries.csv"))[2:end]]
    RM = response_matrix()
    # validity cross-check on hour 12 (idle batteries): corner bounds vs the largest validated currents found earlier
    kw1 = hour_kw(s24, 12)
    fz = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "source_current_bounds_falsify.csv"))[2:end]]
    worst = 0.0
    for key in unique([(f[1], f[2], f[3]) for f in fz])
        ang1 = at(parse(Float64, key[3]), parse(Float64, key[1]))
        cb = current_bounds(RM, ang1, device_rectangles(net, kw1, 1), 1; corners = true)
        for s in net.sources
            found = maximum(parse(Float64, f[10]) for f in fz if (f[1], f[2], f[3]) == key && f[4] == s.name)
            sect = parse(Float64, first(f[11] for f in fz if (f[1], f[2], f[3]) == key && f[4] == s.name))
            worst = max(worst, found / (cb.out[s.name].sect * IbA))
            say(@sprintf("  hour 12, theta %s, %s r*, %s: corner bound %.3f A (sector bound %.3f A); largest validated current found %.3f A", key[1], key[2], s.name, cb.out[s.name].sect * IbA, sect, found))
        end
    end
    say(@sprintf("hour-12 cross-check: largest found current / corner bound = %.4f (%s)", worst, worst <= 1 ? "no corner bound violated" : "A CORNER BOUND IS VIOLATED"))
    worst <= 1 || error("corner enclosure violated by a validated dispatch")
    open(io -> println(io, "theta_deg,point,r_deg,rho_exact_kW,bound_s,max_bound_A_pilot,max_bound_A_refined,double_sector_buses,full_sector_buses,witness_over_bound,relax_status,rho_relax_kW,rho_lb_kW,box_term_kW,rounding_kW,primal_residual,dual_psd_violation,build_s,solve_s,verify_s,lifted_worst,lifted_box_worst,exact_lb_kW,class"), DR * ".csv", "w")
    for row in q
        th = parse(Float64, row[1]); pname = row[2]; r = parse(Float64, row[3]); ang = at(r, th); rho_ex = parse(Float64, row[11]); rub = parse(Float64, row[12])
        xm = solve_multiperiod(net, ang; kw24..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        vm = validate_dispatch(net, ang, kw24, xm; k = KV)
        L = xm.status in RAY_OK && vm.equations <= 1e-6 && vm.viol["soc"] <= 1e-6 ? vm.margin_kW : NaN
        tb = time()
        CBs = [current_bounds(RM, ang, device_rectangles(net, kw24, t), t; ndir = 3600, ndir_rho = 720, corners = true) for t in 1:T]
        tbound = time() - tb
        bnd = Dict(s.name => [CBs[t].out[s.name].sect for t in 1:T] for s in net.sources)
        old = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "daily_pilot_bounds.csv"))[2:end]]
        oldmax = maximum(parse(Float64, o[5]) for o in old if parse(Float64, o[1]) == th && o[2] == pname)
        wr = isfinite(L) ? maximum(abs(complex(xm.raw.ir[s.name][t], xm.raw.ii[s.name][t])) / bnd[s.name][t] for s in net.sources, t in 1:T) : NaN
        say(@sprintf("\n== theta %.1f, %s r* (r = %.6f deg): exact margin L %+.4f kW, exact restoration rho %.4f kW; refined bounds %.1f-%.1f A (pilot up to %.1f A), %d double-sector and %d full-disc device buses; %.1f s%s",
                     th, pname, r, L, rho_ex, minimum(minimum(v) for v in values(bnd)) * IbA, maximum(maximum(v) for v in values(bnd)) * IbA, oldmax,
                     CBs[1].double_sectors, CBs[1].full_sectors, tbound, isfinite(wr) ? @sprintf("; witness currents at most %.4f of the refined bounds", wr) : ""))
        for s in net.sources
            say(@sprintf("     %s refined bounds by hour (A): %s", s.name, join([@sprintf("%.0f", x * IbA) for x in bnd[s.name]], " ")))
        end
        kwb = (; kw24..., source_ampacity_pu = bnd)
        c = certify(net, ang, kwb; k = KV, cuts = CUTS, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, tighten = true, src_scale = true)
        infeas = c.status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
        lb = infeas ? NaN : (isfinite(c.db.bound) ? c.db.bound * S : -Inf)
        cls = infeas ? (c.farkas.bound > 0 ? "no-backflow day impossible (checked Farkas)" : "solver-only infeasibility (ray not validated)") :
              lb > 0 ? "no-backflow day impossible (checked lower bound on rho > 0)" : "unresolved (rho lower bound <= 0)"
        lv = (; worst = NaN, box_worst = NaN, at = "", box_at = "")
        if isfinite(L) && L >= 0
            lv = lifted_violation(c.R, net, xm, kwb; margin_kW = 0.0)
            say(@sprintf("     lifted %s r* daily dispatch: largest constraint violation %.2e (%s); largest box violation %.2e (%s)", pname, lv.worst, lv.at, lv.box_worst, lv.box_at))
        end
        exact_lb = NaN
        if occursin("impossible", cls)
            fname = joinpath(CERT_DIR, @sprintf("daily_refined_restoration_theta%.1f_%s.txt", th, pname))
            export_certificate(c.R.model, c.R.box, fname; c = infeas ? zero(AffExpr) : objective_function(c.R.model),
                               note = @sprintf("IEEE123 T=24, batteries dispatched, cyclic SOC, refined battery enclosure, theta %.1f deg, r = %.6f deg (%s r*), restoration", th, r, pname))
            ve = verify_file(fname); exact_lb = Float64(BigFloat(ve.bound)) * S
            say(@sprintf("     exact re-check: %.6f kW (widened boxes %.6f kW), %d of %d PSD duals shifted, %.1f s", exact_lb, Float64(BigFloat(ve.bound_widened)) * S, ve.nshift, ve.npsd, ve.seconds))
        end
        say(@sprintf("     relaxation: %s / %s, rho relax %.4f kW, checked lower bound %.4f kW (box term %.2e kW, rounding %.1e kW, residual %.1e, dual psd %.1e); build %.1f, solve %.1f, verify %.1f s -> %s",
                     c.status, c.dual_status, infeas ? NaN : -c.primal, lb, c.box_term_kW, c.rounding_kW, c.primal_residual, c.psd_dual_violation, c.build_s, c.solve_s, c.verify_s, cls))
        open(io -> println(io, join([th, pname, @sprintf("%.6f", r), @sprintf("%.4f", rho_ex), @sprintf("%.1f", tbound), @sprintf("%.2f", oldmax), @sprintf("%.2f", maximum(maximum(v) for v in values(bnd)) * IbA),
                                     CBs[1].double_sectors, CBs[1].full_sectors, @sprintf("%.4f", wr), c.status, @sprintf("%.4f", infeas ? NaN : -c.primal), @sprintf("%.4f", lb),
                                     @sprintf("%.3e", c.box_term_kW), @sprintf("%.1e", c.rounding_kW), @sprintf("%.1e", c.primal_residual), @sprintf("%.1e", c.psd_dual_violation),
                                     @sprintf("%.1f", c.build_s), @sprintf("%.1f", c.solve_s), @sprintf("%.1f", c.verify_s), @sprintf("%.2e", lv.worst), @sprintf("%.2e", lv.box_worst),
                                     @sprintf("%.6f", exact_lb), cls], ',')), DR * ".csv", "a")
    end
end
dr_main()
