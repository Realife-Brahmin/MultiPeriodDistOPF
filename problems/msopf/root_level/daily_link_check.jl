# Validity and robustness of the dispatch link (2026-10-22).
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/daily_link_check.jl
# (a) For the validated exact dispatches at every daily query (margin witness; restoration witness, which lies in the
#     restoration domain rho <= rho_max), evaluate |I_s,t - a_s - sum_b k_sb conj(S_b,t)| / eps_s,t from the
#     dispatch itself: a real state above 1 would show eps invalid.
# (b) At 157.5 deg, 1.02 r*: the link relaxation with eps x 1.1, 1.25, 1.5 (all valid), and eps x 1 under two other
#     Clarabel regularizations, to see whether the certificate is robust or a numerical accident.
ENV["DL_LIB"] = "1"
include(joinpath(@__DIR__, "daily_link.jl"))
function link_residual(x, ang, kw, CBs)
    S = net.S_base_kVA; T = kw.T; worst = 0.0; at_ = ""
    for t in 1:T, s in net.sources
        cb = CBs[t]; a = cb.out[s.name].a; cs = cb.out[s.name].cs
        I = complex(x.raw.ir[s.name][t], x.raw.ii[s.name][t])
        A = a
        for (q, b) in enumerate(cb.D)
            Sb = sum((-kw.load_factor[t] * complex(d.kW, d.kvar) / S for d in net.loads if d.bus == b); init = 0.0im) +
                 sum((kw.pv_factor[t] * p.kW / S + im * x.Q_PV_kvar[p.name][t] / S for p in kw.pvs if p.bus == b); init = 0.0im) +
                 sum((x.P_B_kW[bt.name][t] / S for bt in net.batteries if bt.bus == b); init = 0.0im)
            A += cs[q] / conj(cb.V0D[q]) * conj(Sb)
        end
        ep = sum(abs(cs[q]) * cb.Smax[q] * cb.rho[q] / (cb.vmn[q] * abs(cb.V0D[q])) for q in eachindex(cb.D))
        rr = abs(I - A) / ep
        rr > worst && (worst = rr; at_ = "$(s.name) hour $t")
    end
    return worst, at_
end
function chk_main()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw24 = problem(s24, :cyclic); T = kw24.T
    q = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "daily_pilot_queries.csv"))[2:end]]
    RM = response_matrix()
    for row in q
        th = parse(Float64, row[1]); r = parse(Float64, row[3]); ang = at(r, th)
        CBs = [current_bounds(RM, ang, device_rectangles(net, kw24, t), t; ndir = 3600, ndir_rho = 720, corners = true) for t in 1:T]
        xm = solve_multiperiod(net, ang; kw24..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        xr = solve_multiperiod(net, ang; kw24..., objective = :margin, margin_band_hard = true, margin_reg = 0.0)
        vm = validate_dispatch(net, ang, kw24, xm; k = KV); vr = validate_dispatch(net, ang, kw24, xr; k = KV)
        wm, am = link_residual(xm, ang, kw24, CBs); wr, ar = link_residual(xr, ang, kw24, CBs)
        say(@sprintf("theta %s, %s r*: margin witness (equations %.1e) uses at most %.4f of eps (%s); restoration witness (rho %.4f kW, equations %.1e, band %.1e) uses at most %.4f of eps (%s)",
                     row[1], row[2], vm.equations, wm, am, max(0.0, -vr.min_import_kW), vr.equations, vr.viol["band"], wr, ar))
    end
    row = only(r for r in q if r[1] == "157.5" && r[2] == "1.02")
    th = 157.5; r = parse(Float64, row[3]); ang = at(r, th); rub = parse(Float64, row[12])
    CBs = [current_bounds(RM, ang, device_rectangles(net, kw24, t), t; ndir = 3600, ndir_rho = 720, corners = true) for t in 1:T]
    bnd = Dict(s.name => [CBs[t].out[s.name].sect for t in 1:T] for s in net.sources)
    regs = Dict("reg 1e-6" => [("max_iter", 400), ("static_regularization_constant", 1e-6), ("tol_gap_abs", 1e-6), ("tol_gap_rel", 1e-6), ("tol_feas", 1e-6)],
                "reg 1e-8" => [("max_iter", 400), ("static_regularization_constant", 1e-8), ("tol_gap_abs", 1e-6), ("tol_gap_rel", 1e-6), ("tol_feas", 1e-6)])
    open(io -> println(io, "theta_deg,point,eps_scale,solver,relax_status,rho_relax_kW,rho_lb_kW,fictitious_losses_kWh,max_source_moment_gap_pu2,solve_s"), DL * "_robust.csv", "w")
    for (sc, rn) in ((1.1, "default"), (1.25, "default"), (1.5, "default"), (1.0, "reg 1e-6"), (1.0, "reg 1e-8"), (2.0, "reg 1e-6"))
        o = run_relax(ang, kw24, bnd, rub, CBs; link = true, eps_scale = sc, attrs = rn == "default" ? CLARABEL_DEFAULTS : regs[rn])
        say(@sprintf("  157.5, 1.02 r*, eps x %.2f, %s: %s, rho relax %.4f kW, checked lower bound %.4f kW, fictitious losses %.1f kWh, source moment gap %.3f pu^2, solve %.1f s",
                     sc, rn, o.st, o.rho, o.lb, o.fict, o.srcgap, o.ts))
        dl_csv("_robust.csv", [157.5, "1.02", sc, rn, o.st, @sprintf("%.4f", o.rho), @sprintf("%.4f", o.lb), @sprintf("%.2f", o.fict), @sprintf("%.3e", o.srcgap), @sprintf("%.1f", o.ts)])
    end
end
chk_main()
