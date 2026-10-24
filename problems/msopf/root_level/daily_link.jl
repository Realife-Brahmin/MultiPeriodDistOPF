# Dispatch-dependent network-response constraints in the daily relaxation (2026-10-22).
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/daily_link.jl
# IEEE 123, sources 1, 3, 4, T = 24, PV reactive power dispatched, batteries dispatched, hard cyclic SOC; the daily
# queries 0.98 and 1.02 r* at 292.5 and 157.5 deg (daily_pilot_queries.csv). Baseline = the corner-enclosure run of
# daily_refined.jl (derived per-hour source-current bounds, both nodal strengthenings, scaled rows). Added: for every
# source s and hour t the second-order cone
#     | I_s,t - a_s - sum_b k_sb conj(S_b,t) | <= eps_s,t,    k_sb = c_sb / conj(V0_b),
#     eps_s,t = sum_b |c_sb| Smax_b,t rho_b,t / (vmin_b,t |V0_b|) + delta,
# with S_b,t the model's own device power (fixed load and PV real power, the PV reactive-power and battery-power
# variables), I_s,t = (E_s - V_s,t) / z_s from the terminal voltage moments, rho, vmin, Smax from current_bounds (corner
# enclosure) and delta an allowance for rounding the coefficients to Float64. Derivation: daily_link.md.
# Also: the thin 292.5 deg 1.02 r* certificate under three modest solver variations; the 1.10 r* certificates as a
# regression check of the changed formulation. Writes daily_link{.csv,_eps.csv,.log} and certificates/daily_link_*.txt.
ENV["SCB_LIB"] = "1"
include(joinpath(@__DIR__, "source_current_bounds.jl"))
include(joinpath(@__DIR__, "..", "current_bounds.jl"))
include(joinpath(@__DIR__, "verify_certificate.jl"))

const DL = joinpath(OUT_DIR, "daily_link")
const CERT_DIR = joinpath(OUT_DIR, "certificates")
dl_csv(f, row) = open(io -> println(io, join(row, ',')), DL * f, "a")

"Add the dispatch-link cones to relaxation R; returns eps (pu) by (source, hour) and the cone references."
function add_dispatch_links!(R, ang, kw, CBs; eps_scale = 1.0)
    S = net.S_base_kVA; T = kw.T; lam = kw.load_factor; pvf = get(kw, :pv_factor, zeros(T))
    model = R.model; eps = Dict{Tuple{String,Int},Float64}(); cons = []
    pvat = Dict(b => [p.name for p in kw.pvs if p.bus == b] for b in net.buses)
    btat = Dict(b => [bt.name for bt in net.batteries if bt.bus == b] for b in net.buses)
    for t in 1:T
        cb = CBs[t]; D = cb.D
        minV0 = minimum(abs.(cb.V0D)); minV0 > 0.5 || error("no-load voltage too small at a device bus (hour $t): the link needs V0 != 0")
        for s in net.sources
            a = cb.out[s.name].a; cs = cb.out[s.name].cs
            z = s.z_series * net.Y_base; y = 1 / z; E = s.pu * cis(deg2rad(angle_at(ang, s.name, t))); bt = s.bus
            # I_s = y E - y V_s
            Dre = AffExpr(real(y * E) - real(a)); Dim = AffExpr(imag(y * E) - imag(a))
            add_to_expression!(Dre, -real(y), R.e[bt, t]); add_to_expression!(Dre, imag(y), R.f[bt, t])
            add_to_expression!(Dim, -imag(y), R.e[bt, t]); add_to_expression!(Dim, -real(y), R.f[bt, t])
            ep = 0.0; mag = abs(a)
            for (q, b) in enumerate(D)
                k = cs[q] / conj(cb.V0D[q]); kr, ki = real(k), imag(k)
                # conj(S_b) = P_b - j Q_b;  k conj(S) = (kr P + ki Q) + j (ki P - kr Q)
                P0 = sum((-lam[t] * d.kW / S for d in net.loads if d.bus == b); init = 0.0) + sum((pvf[t] * p.kW / S for p in kw.pvs if p.bus == b); init = 0.0)
                Q0 = sum((-lam[t] * d.kvar / S for d in net.loads if d.bus == b); init = 0.0)
                add_to_expression!(Dre, -(kr * P0 + ki * Q0)); add_to_expression!(Dim, -(ki * P0 - kr * Q0))
                for kk in btat[b]; add_to_expression!(Dre, -kr, R.P_B[kk, t]); add_to_expression!(Dim, -ki, R.P_B[kk, t]); end
                for kk in pvat[b]; add_to_expression!(Dre, -ki, R.Q_PV[kk, t]); add_to_expression!(Dim, kr, R.Q_PV[kk, t]); end
                ep += abs(cs[q]) * cb.Smax[q] * cb.rho[q] / (cb.vmn[q] * abs(cb.V0D[q]))
                mag += abs(k) * cb.Smax[q]
            end
            delta = 1e-12 * (mag + abs(y) * 3)                 # Float64 rounding of a, k, y and the products (generous)
            e_st = (ep * (1 + 1e-9) + delta) * eps_scale       # eps_scale < 1 is a DIAGNOSTIC only (not valid)
            eps[(s.name, t)] = e_st
            sc = abs(z)                                        # scale the cone by |z_s| (unit-size voltage coefficients)
            push!(cons, @constraint(model, [sc * e_st; sc * Dre; sc * Dim] in SecondOrderCone()))
        end
    end
    return eps, cons
end

"Restoration relaxation with or without the links; returns status, checked bound (kW), times, diagnostics."
function run_relax(ang, kw, bnd, rub, CBs; link, attrs = CLARABEL_DEFAULTS, export_to = nothing, eps_scale = 1.0)
    S = net.S_base_kVA; T = kw.T
    kwb = (; kw..., source_ampacity_pu = bnd)
    t0 = time()
    R = relax_margin(net, ang, kwb; k = KV, cuts = CUTS, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, solver_attributes = attrs, src_scale = true)
    eps = nothing; link && ((eps, _) = add_dispatch_links!(R, ang, kw, CBs; eps_scale))
    tb = time() - t0
    optimize!(R.model); ts = time() - t0 - tb
    tv = time()
    tighten_boxes!(R, net, kwb, ang)
    st = termination_status(R.model); infeas = st in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
    db = has_duals(R.model) ? dual_bound(R.model; box = R.box) : (; bound = -Inf, box_term = NaN, rounding = NaN, psd_viol = NaN, sign_viol = NaN)
    fb = infeas && dual_status(R.model) in (MOI.INFEASIBILITY_CERTIFICATE, MOI.NEARLY_INFEASIBILITY_CERTIFICATE) ? dual_bound(R.model; c = zero(AffExpr), box = R.box).bound : NaN
    tver = time() - tv
    lb = infeas ? NaN : db.bound * S
    certified = infeas ? fb > 0 : lb > 0
    # where the relaxation's freedom sits: fictitious line losses, and at the sources the lifted |I_s|^2 beyond the
    # square of its own first moment (pu^2), and the first-moment source current against the link's affine centre
    fict = NaN; srcgap = NaN; rho = NaN; pres = NaN
    if has_values(R.model)
        rho = -value(R.m) * S
        pres = maximum(values(primal_feasibility_report(R.model; atol = 0.0)); init = 0.0)
        fict = 0.0
        for l in net.lines, tt in 1:T
            Rl = real(-net.Y_base / l.Y[1, 2]); g(sy) = value(R.Lv[sy][l.name, tt])
            fict += Rl * (g(:lrlr) + g(:lili) - g(:lr)^2 - g(:li)^2)
        end
        fict *= S
        srcgap = 0.0
        for s in net.sources, tt in 1:T
            z = s.z_series * net.Y_base; E = s.pu * cis(deg2rad(angle_at(ang, s.name, tt))); b = s.bus
            ev, fv = value(R.e[b, tt]), value(R.f[b, tt]); W = value(R.Wee[b, tt]) + value(R.Wff[b, tt])
            lifted = (abs2(E) - 2 * (real(E) * ev + imag(E) * fv) + W) / abs2(z); first = abs2((E - complex(ev, fv)) / z)
            srcgap = max(srcgap, lifted - first)
        end
    end
    exact = NaN; exs = NaN; nsoc = 0; nfix = 0
    if certified && export_to !== nothing
        export_certificate(R.model, R.box, export_to; c = infeas ? zero(AffExpr) : objective_function(R.model), note = "daily restoration relaxation, dispatch link $(link)")
        ve = verify_file(export_to); exact = Float64(BigFloat(ve.bound)) * S; exs = ve.seconds; nsoc = ve.nsoc; nfix = ve.nsocfix
    end
    return (; st, dst = dual_status(R.model), lb, fb, certified, rho, fict, srcgap, pres, box = db.box_term * S, rounding = db.rounding * S,
            psd = db.psd_viol, sign = db.sign_viol, tb, ts, tver, eps, exact, exs, nsoc, nfix, R)
end

function dl_main()
    mkpath(CERT_DIR)
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw24 = problem(s24, :cyclic); T = kw24.T
    q = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "daily_pilot_queries.csv"))[2:end]]
    tp = time(); RM = response_matrix(); say(@sprintf("network response (BigFloat 256): %.1f s", time() - tp))
    open(io -> println(io, "theta_deg,point,case,link,relax_status,dual_status,rho_relax_kW,rho_lb_kW,exact_lb_kW,box_term_kW,rounding_kW,primal_residual,dual_psd_violation,dual_sign_violation,fictitious_losses_kWh,max_source_moment_gap_pu2,eps_min_A,eps_max_A,eps_mean_A,bounds_s,build_s,solve_s,verify_s,exact_verify_s,soc_multipliers_adjusted,lifted_worst,lifted_at,lifted_box_worst,class"), DL * ".csv", "w")
    open(io -> println(io, "theta_deg,point,source,hour,eps_A,eps_pu,current_bound_A"), DL * "_eps.csv", "w")
    queries = [(th, pt) for th in ("292.5", "157.5") for pt in ("0.98", "1.02")]
    push!(queries, ("292.5", "1.10"), ("157.5", "1.10"))                       # regression: the formulation changed
    for (th_s, pt) in queries
        row = only(r for r in q if r[1] == th_s && r[2] == pt)
        th = parse(Float64, th_s); r = parse(Float64, row[3]); ang = at(r, th); rho_ex = parse(Float64, row[11]); rub = parse(Float64, row[12])
        xm = solve_multiperiod(net, ang; kw24..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        vm = validate_dispatch(net, ang, kw24, xm; k = KV)
        L = xm.status in RAY_OK && vm.equations <= 1e-6 && vm.viol["soc"] <= 1e-6 ? vm.margin_kW : NaN
        tb0 = time(); CBs = [current_bounds(RM, ang, device_rectangles(net, kw24, t), t; ndir = 3600, ndir_rho = 720, corners = true) for t in 1:T]; tbound = time() - tb0
        bnd = Dict(s.name => [CBs[t].out[s.name].sect for t in 1:T] for s in net.sources)
        say(@sprintf("\n== theta %s, %s r* (r = %.6f deg; delta3 %+.6f, delta4 %+.6f): exact margin L %+.4f kW, exact restoration rho %.4f kW, rho <= %.4f kW; corner bounds %.1f s; min |V0| at device buses %.4f pu",
                     th_s, pt, r, ang["subs3"], ang["subs4"], L, rho_ex, rub, tbound, minimum(minimum(abs.(cb.V0D)) for cb in CBs)))
        cases = Any[("baseline (corner bounds)", false, CLARABEL_DEFAULTS), ("dispatch link", true, CLARABEL_DEFAULTS)]
        if th_s == "292.5" && pt == "1.02"
            push!(cases, ("baseline, static reg 1e-8", false, [("max_iter", 400), ("static_regularization_constant", 1e-8), ("tol_gap_abs", 1e-6), ("tol_gap_rel", 1e-6), ("tol_feas", 1e-6)]),
                         ("baseline, static reg 1e-6", false, [("max_iter", 400), ("static_regularization_constant", 1e-6), ("tol_gap_abs", 1e-6), ("tol_gap_rel", 1e-6), ("tol_feas", 1e-6)]),
                         ("baseline, no equilibration", false, [CLARABEL_DEFAULTS; ("equilibrate_enable", false)]))
        end
        pt == "1.10" && (cases = cases[2:2])
        for (cn, link, attrs) in cases
            fn = joinpath(CERT_DIR, @sprintf("daily_link_theta%s_%s_%s.txt", th_s, pt, replace(cn, r"[ ,()]+" => "_")))
            o = run_relax(ang, kw24, bnd, rub, CBs; link, attrs, export_to = fn)
            cls = o.certified ? "no-backflow day impossible (checked)" : "unresolved"
            lv = (; worst = NaN, at = "", box_worst = NaN)
            if isfinite(L) && L >= 0
                lv = lifted_violation(o.R, net, xm, (; kw24..., source_ampacity_pu = bnd); margin_kW = 0.0)
            end
            ev = o.eps === nothing ? Float64[] : collect(values(o.eps)) .* IbA
            if o.eps !== nothing && cn == "dispatch link"
                for ((sn, t), e) in sort(collect(o.eps))
                    dl_csv("_eps.csv", [th_s, pt, sn, t, @sprintf("%.4f", e * IbA), @sprintf("%.6e", e), @sprintf("%.3f", bnd[sn][t] * IbA)])
                end
            end
            say(@sprintf("  [%s] %s / %s: rho relax %.4f kW, checked lower bound %.4f kW%s%s | box term %.2e, rounding %.1e kW | residual %.1e, dual psd/sign %.1e/%.1e | fictitious losses %.1f kWh, largest source |I|^2 moment gap %.2e pu^2%s | build %.1f, solve %.1f, verify %.1f s%s -> %s",
                         cn, o.st, o.dst, o.rho, o.lb, isfinite(o.fb) ? @sprintf(", Farkas %.3e", o.fb) : "", isfinite(o.exact) ? @sprintf(" (exact re-check %.6f kW, %.1f s, %d SOC multipliers, %d adjusted)", o.exact, o.exs, o.nsoc, o.nfix) : "",
                         o.box, o.rounding, o.pres, o.psd, o.sign, o.fict, o.srcgap,
                         isempty(ev) ? "" : @sprintf(" | eps %.3f-%.3f A (mean %.3f A)", minimum(ev), maximum(ev), sum(ev) / length(ev)),
                         o.tb, o.ts, o.tver, isfinite(lv.worst) ? @sprintf(" | lifted 0.98 witness: constraints %.2e (%s), boxes %.2e", lv.worst, lv.at, lv.box_worst) : "", cls))
            dl_csv(".csv", [th_s, pt, "\"" * cn * "\"", link, o.st, o.dst, @sprintf("%.4f", o.rho), @sprintf("%.4f", o.lb), @sprintf("%.6f", o.exact), @sprintf("%.3e", o.box), @sprintf("%.1e", o.rounding),
                            @sprintf("%.1e", o.pres), @sprintf("%.1e", o.psd), @sprintf("%.1e", o.sign), @sprintf("%.2f", o.fict), @sprintf("%.3e", o.srcgap),
                            isempty(ev) ? "" : @sprintf("%.4f", minimum(ev)), isempty(ev) ? "" : @sprintf("%.4f", maximum(ev)), isempty(ev) ? "" : @sprintf("%.4f", sum(ev) / length(ev)),
                            @sprintf("%.1f", tbound), @sprintf("%.1f", o.tb), @sprintf("%.1f", o.ts), @sprintf("%.1f", o.tver), @sprintf("%.1f", o.exs), o.nfix,
                            @sprintf("%.2e", lv.worst), replace(string(lv.at), ',' => ';'), @sprintf("%.2e", lv.box_worst), cls])
        end
    end
end
get(ENV, "DL_LIB", "0") == "1" || dl_main()
