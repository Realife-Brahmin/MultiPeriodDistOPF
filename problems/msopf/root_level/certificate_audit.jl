# Audit of the two near-boundary certificates (2026-10-22).
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/certificate_audit.jl
# IEEE 123-node, sources 1, 3, 4, hour 12 alone, batteries idle, PV reactive power dispatched; directions 292.5 and
# 157.5 deg; the points 1.02 r* (certified with derived source-current bounds and both nodal strengthenings) and
# 0.98 r* (feasible controls). For each:
#   (a) the derived source-current bounds recomputed: the committed Float64 code, the generic rectangle code in
#       Float64, and the same code in BigFloat (256-bit), same direction grid -- the Float64 values used must be
#       at least the BigFloat ones;
#   (b) the restoration certificate re-solved in four numerically different but equivalent forms (as committed;
#       the source-current rows scaled by |z_s|^2; tighter Clarabel tolerances; both), each exported to a text
#       file and re-checked in exact rational arithmetic by verify_certificate.jl (no solver);
#   (c) the 0.98 r* validated dispatch lifted into the strengthened restoration relaxation: every constraint and
#       every tightened verifier box checked.
# Writes certificate_audit.csv, certificate_audit_bounds.csv, certificates/*.txt and certificate_audit.log.
ENV["SCB_LIB"] = "1"
include(joinpath(@__DIR__, "source_current_bounds.jl"))
include(joinpath(@__DIR__, "..", "current_bounds.jl"))
include(joinpath(@__DIR__, "verify_certificate.jl"))

const CA = joinpath(OUT_DIR, "certificate_audit")
const CERT_DIR = joinpath(OUT_DIR, "certificates")
const TIGHT = [("max_iter", 1000), ("static_regularization_constant", 1e-8), ("tol_gap_abs", 1e-9), ("tol_gap_rel", 1e-9), ("tol_feas", 1e-9), ("tol_ktratio", 1e-8)]
ca_csv(file, row) = open(io -> println(io, join(row, ',')), file, "a")

function audit()
    mkpath(CERT_DIR)
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    cw = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "certified_window.csv"))[2:end]]
    F = CA * ".csv"; FB = CA * "_bounds.csv"
    open(io -> println(io, "theta_deg,point,variant,relax_status,dual_status,rho_relax_kW,lb_float_kW,rounding_allowance_kW,lb_exact_kW,exact_minus_float_kW,lb_exact_widened_kW,box_term_kW,raw_psd_dual_min_eig,psd_blocks_shifted,max_shift_rel,sign_violation,primal_residual,build_s,solve_s,verify_float_s,verify_exact_s,file"), F, "w")
    open(io -> println(io, "theta_deg,point,source,committed_float_A,rect_float_A,bigfloat_A,committed_minus_bigfloat_rel,termwise_A,no_load_A"), FB, "w")
    t0 = time(); RM = response_matrix(); say(@sprintf("network response (BigFloat 256): %.1f s; Re(Y) min pivot %.3e; |YZ - I| %.1e", time() - t0, Float64(RM.minpiv), Float64(RM.resid)))
    for th in (292.5, 157.5)
        rs = parse(Float64, only(r for r in cw if parse(Float64, r[1]) == th && r[2] == "sources+lines 1.2x" && r[4] == "")[5])
        for (pname, fac) in (("1.02", 1.02), ("0.98", 0.98))
            r = fac * rs; ang = at(r, th)
            say(@sprintf("\n== theta %.1f, %s r* (r = %.6f deg; delta3 %+.6f, delta4 %+.6f deg)", th, pname, r, ang["subs3"], ang["subs4"]))
            # (a) bounds three ways
            AB = analytic_bounds(RM, ang, kw1)
            rects = device_rectangles(net, kw1, 1)
            tf = time(); CF = current_bounds(RM, ang, rects, 1; T = Float64); tf = time() - tf
            tb = time(); CB = setprecision(BigFloat, 256) do; current_bounds(RM, ang, rects, 1; T = BigFloat); end; tb = time() - tb
            bnd = Dict{String,Float64}()
            for s in net.sources
                old = AB.out[s.name].sect; nf = CF.out[s.name].sect; bf = CB.out[s.name].sect
                bnd[s.name] = old
                say(@sprintf("  %s: committed %.9f A | rectangle code, Float64 %.9f A | BigFloat %.9f A | committed - BigFloat = %+.2e relative (%s)",
                             s.name, old * IbA, nf * IbA, Float64(bf) * IbA, (old - bf) / bf, old >= bf ? "committed value covers the high-precision bound" : "COMMITTED VALUE BELOW THE HIGH-PRECISION BOUND"))
                ca_csv(FB, [@sprintf("%.1f", th), pname, s.name, @sprintf("%.9f", old * IbA), @sprintf("%.9f", nf * IbA), @sprintf("%.9f", Float64(bf) * IbA),
                            @sprintf("%.3e", Float64((old - bf) / bf)), @sprintf("%.6f", CF.out[s.name].term * IbA), @sprintf("%.6f", CF.out[s.name].a0 * IbA)])
            end
            say(@sprintf("  bound construction: Float64 %.1f s, BigFloat %.1f s; disc radius history (BigFloat) %s pu", tf, tb, join([@sprintf("%.6f", Float64(x)) for x in CB.rho_hist], " -> ")))
            # exact witnesses
            xm = solve_multiperiod(net, ang; kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
            vm = validate_dispatch(net, ang, kw1, xm; k = KV)
            xr = solve_multiperiod(net, ang; kw1..., objective = :margin, margin_band_hard = true, margin_reg = 0.0)
            vr = validate_dispatch(net, ang, kw1, xr; k = KV)
            rho_ex = max(0.0, -vr.min_import_kW); rub = rho_ex + 1.0
            say(@sprintf("  witnesses: margin L %+.4f kW (equations %.1e, band %.1e); restoration rho %.4f kW (equations %.1e, band %.1e); rho <= %.4f kW in the relaxation",
                         vm.margin_kW, vm.equations, vm.viol["band"], rho_ex, vr.equations, vr.viol["band"], rub))
            kwb = (; kw1..., source_ampacity_pu = copy(bnd))
            variants = pname == "1.02" ? [("committed", CLARABEL_DEFAULTS, false), ("rows scaled", CLARABEL_DEFAULTS, true), ("tight tolerances", TIGHT, false), ("scaled + tight", TIGHT, true)] :
                                         [("committed", CLARABEL_DEFAULTS, false)]
            for (vn, attrs, sc) in variants
                tw = time()
                c = certify(net, ang, kwb; k = KV, cuts = CUTS, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, tighten = true, solver_attributes = attrs, src_scale = sc)
                S = net.S_base_kVA
                infeas = c.status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
                fname = joinpath(CERT_DIR, @sprintf("restoration_theta%.1f_%s_%s.txt", th, pname, replace(vn, ' ' => '_', '+' => "and")))
                export_certificate(c.R.model, c.R.box, fname; c = infeas ? zero(AffExpr) : objective_function(c.R.model),
                                   note = @sprintf("IEEE123 hour 12, idle batteries, theta %.1f deg, r = %.6f deg (%s r*), restoration (min rho = -m), %s; %s", th, r, pname, vn, infeas ? "Farkas (c = 0)" : "objective"))
                ve = verify_file(fname)
                lbf = infeas ? c.farkas.bound : c.db.bound
                say(@sprintf("  [%s] %s / %s, rho relax %.4f kW | float check %.6f kW (allowance %.1e kW) | EXACT %.9f kW (exact - float %+.2e kW; boxes widened 1e-12: %.9f kW) | box term %.3e kW | raw PSD duals: min eig %.1e, %d of %d blocks shifted (largest %.1e) | sign viol %.1e | residual %.1e | build %.1f, solve %.1f, float verify %.2f, exact verify %.1f s",
                             vn, c.status, c.dual_status, infeas ? NaN : -c.primal, lbf * S, c.rounding_kW, Float64(BigFloat(ve.bound)) * S, (Float64(BigFloat(ve.bound)) - lbf) * S,
                             Float64(BigFloat(ve.bound_widened)) * S, Float64(ve.boxterm) * S, ve.psd_min_eig, ve.nshift, ve.npsd, ve.maxtau, ve.sign_viol, c.primal_residual,
                             c.build_s, c.solve_s, c.verify_s, ve.seconds))
                ca_csv(F, [@sprintf("%.1f", th), pname, vn, string(c.status), string(c.dual_status), @sprintf("%.6f", infeas ? NaN : -c.primal), @sprintf("%.9f", lbf * S),
                           @sprintf("%.2e", c.rounding_kW), @sprintf("%.9f", Float64(BigFloat(ve.bound)) * S), @sprintf("%.2e", (Float64(BigFloat(ve.bound)) - lbf) * S),
                           @sprintf("%.9f", Float64(BigFloat(ve.bound_widened)) * S), @sprintf("%.4e", Float64(ve.boxterm) * S), @sprintf("%.2e", ve.psd_min_eig), string(ve.nshift),
                           @sprintf("%.1e", ve.maxtau), @sprintf("%.1e", ve.sign_viol), @sprintf("%.1e", c.primal_residual), @sprintf("%.2f", c.build_s), @sprintf("%.2f", c.solve_s),
                           @sprintf("%.2f", c.verify_s), @sprintf("%.1f", ve.seconds), relpath(fname, OUT_DIR)])
                if pname == "0.98"
                    lv = lifted_violation(c.R, net, xm, kwb; margin_kW = 0.0)
                    ws = maximum(abs(complex(xm.raw.ir[s.name][1], xm.raw.ii[s.name][1])) / bnd[s.name] for s in net.sources)
                    say(@sprintf("  lifted 0.98 r* dispatch (rho = 0): largest constraint violation %.2e (%s); largest box violation %.2e (%s); source currents at most %.4f of their derived bounds",
                                 lv.worst, lv.at, lv.box_worst, lv.box_at, ws))
                end
            end
        end
    end
end
audit()
