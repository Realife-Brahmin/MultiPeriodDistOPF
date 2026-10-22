# Source-current bounds derived from the original model, and whether they keep the certificates (2026-10-22).
#   julia --project=envs/msopf_sdp envs/multi_poi/root_level/source_current_bounds.jl
# IEEE 123-node, sources 1, 3, 4, hour 12, batteries idle, PV reactive power dispatched; directions
# 292.5 and 157.5 deg. No ampacity or empirical cap anywhere in the model or in the bounds. Derivation:
# results/multiperiod/source_current_bounds.md. Three bounds per source and query:
#   naive    |I_s| <= (|E_s| + max|V_s|) / |z_s|, the model's box |e|, |f| <= 1.5 at the terminal (KVL and the
#            triangle inequality alone; the bound the relaxation already had);
#   termwise the network eliminated, I_s = a_s + sum_b c_sb J_b, then |I_s| <= |a_s| + sum_b |c_sb| Jmax_b with
#            Jmax_b = Smax_b / Vmin_b (the :jbound bound);
#   sector   the same linear relation kept whole: each J_b lies in a circular sector (its phase is arg V_b - arg S_b,
#            with V_b in a disc around the no-load voltage and S_b on its PV reactive segment); the bound is the
#            support function of the Minkowski sum, maximized over a direction grid and divided by cos(grid/2);
# then optimization-based tightening (OBBT): maximize the lifted |I_s|^2 over the restoration relaxation (band
# hard, P_s >= -rho_max, both nodal additions) with the previous round's bounds, and take the checked upper
# bound from dual_bound. Finally the restoration test with the final bounds.
ENV["CV_EXTRA"] = "lib"
include(joinpath(@__DIR__, "certificate_validation.jl"))

const SCB = joinpath(OUT_DIR, "source_current_bounds")
const NDIR = 3600                                   # direction grid, 0.1 deg
const SAFE = 1 + 1e-9                               # outward rounding allowance on the Float64 support sums
const CUTS = [:kcl0, :jbound]
const IbA = net.S_base_kVA / net.kV_base            # A per pu
setprecision(BigFloat, 256)
scb_csv(file, row) = open(io -> println(io, join(row, ',')), SCB * "_" * file * ".csv", "a")

"Bus admittance matrix of the model's current-form network equations, with the sources' series admittances, in BigFloat."
function response_matrix()
    buses = net.buses; n = length(buses); id = Dict(b => i for (i, b) in enumerate(buses)); Yb = big(net.Y_base)
    Y = zeros(Complex{BigFloat}, n, n)
    for l in net.lines
        Yl = Complex{BigFloat}.(l.Y) ./ Yb; ys = -Yl[1, 2]; i, j = id[l.from], id[l.to]
        Y[i, i] += ys + Yl[1, 1] + Yl[1, 2]; Y[j, j] += ys + Yl[2, 2] + Yl[2, 1]; Y[i, j] -= ys; Y[j, i] -= ys
    end
    for s in net.sources
        Y[id[s.bus], id[s.bus]] += 1 / (big(s.z_series) * Yb)
    end
    ch = cholesky(Symmetric(real.(Y)))              # Re(Y) positive definite => Y nonsingular
    minpiv = minimum(diag(ch.U))^2
    Z = inv(Y)
    resid = maximum(abs.(Y * Z - LinearAlgebra.I))
    return (; Y, Z, id, minpiv, resid)
end

"Device power at each bus: S_b = P_b + j(Q0_b + q), q in [-qd_b, qd_b] (loads fixed, PV real fixed, PV reactive free, batteries idle)."
function device_segments(kw)
    S = net.S_base_kVA; lam = kw.load_factor[1]; pvf = kw.pv_factor[1]
    kw.battery === :idle || error("batteries must be idle here")
    P = Dict(b => 0.0 for b in net.buses); Q0 = Dict(b => 0.0 for b in net.buses); qd = Dict(b => 0.0 for b in net.buses)
    for d in net.loads; P[d.bus] -= lam * d.kW / S; Q0[d.bus] -= lam * d.kvar / S; end
    for p in kw.pvs; P[p.bus] += pvf * p.kW / S; qd[p.bus] += sqrt(max(p.kVA^2 - (pvf * p.kW)^2, 0.0)) / S; end
    dev = [b for b in net.buses if P[b] != 0 || Q0[b] != 0 || qd[b] != 0]
    src = Set(s.bus for s in net.sources); any(b in src for b in dev) && error("device at a source terminal")
    return (; P, Q0, qd, dev)
end

"Add the support function of the hull of the sector {r e^{j phi}: r <= R, |phi - phic| <= w} to acc over the grid."
function support!(acc, psis, phic, w, R)
    w >= pi && (acc .+= R; return acc)
    @inbounds for k in eachindex(psis)
        d = abs(rem2pi(psis[k] - phic, RoundNearest))
        acc[k] += d <= w ? R : (d - w < pi / 2 ? R * cos(d - w) : 0.0)
    end
    return acc
end

"The three analytical bounds (pu) at angle vector ang, plus the sectors used."
function analytic_bounds(RM, ang, kw)
    lo = get(kw, :V_band, (0.95, 1.05))[1]
    dv = device_segments(kw); id = RM.id; Yb = big(net.Y_base)
    n = length(net.buses); b0 = zeros(Complex{BigFloat}, n)
    a = Dict{String,ComplexF64}(); c = Dict{String,Vector{ComplexF64}}(); ysd = Dict{String,Complex{BigFloat}}(); Ed = Dict{String,Complex{BigFloat}}()
    for s in net.sources
        E = big(s.pu) * cis(deg2rad(big(angle_at(ang, s.name, 1)))); y = 1 / (big(s.z_series) * Yb)
        b0[id[s.bus]] += y * E; ysd[s.name] = y; Ed[s.name] = E
    end
    V0b = RM.Z * b0
    for s in net.sources
        t = id[s.bus]; y = ysd[s.name]
        a[s.name] = ComplexF64(y * (Ed[s.name] - V0b[t]))
        c[s.name] = ComplexF64.(-y .* RM.Z[t, :])
    end
    V0 = ComplexF64.(V0b); Zf = ComplexF64.(RM.Z)
    D = dv.dev; k_ = [id[b] for b in D]
    Slo = [complex(dv.P[b], dv.Q0[b] - dv.qd[b]) for b in D]; Shi = [complex(dv.P[b], dv.Q0[b] + dv.qd[b]) for b in D]
    Smax = max.(abs.(Slo), abs.(Shi))
    # arc of arg S along the segment (it misses the origin unless P = 0; then the arc [-pi/2, pi/2] covers both ends)
    sarc = map(eachindex(D)) do q
        if dv.P[D[q]] == 0 && real(Slo[q]) == 0 && imag(Slo[q]) <= 0 <= imag(Shi[q]); (0.0, pi / 2)
        else
            p1, p2 = angle(Slo[q]), angle(Shi[q]); sp = rem2pi(p2 - p1, RoundNearest); (p1 + sp / 2, abs(sp) / 2)
        end
    end
    psis = collect(range(0, 2pi, length = NDIR + 1))[1:NDIR]; cw = cos(pi / NDIR)
    Jmax0 = Smax ./ lo                                        # the existing |J| <= Smax / Vmin, Vmin = band
    rho = [sum(abs(Zf[id[b], k_[q]]) * Jmax0[q] for q in eachindex(D)) for b in D]       # termwise disc radius
    hist = [maximum(rho)]
    phiJ = zeros(length(D)); wJ = fill(Float64(pi), length(D)); RJ = copy(Jmax0)
    for it in 1:8
        for q in eachindex(D)
            v0 = abs(V0[k_[q]]); al = rho[q] < v0 ? asin(rho[q] / v0) : pi
            RJ[q] = Smax[q] / max(lo, v0 - rho[q])
            phiJ[q] = angle(V0[k_[q]]) - sarc[q][1]; wJ[q] = min(al + sarc[q][2], pi)
        end
        new = similar(rho)
        acc = zeros(NDIR)
        for (p, b) in enumerate(D)
            fill!(acc, 0.0)
            for q in eachindex(D)
                zz = Zf[id[b], k_[q]]; support!(acc, psis, phiJ[q] + angle(zz), wJ[q], abs(zz) * RJ[q])
            end
            new[p] = maximum(acc) / cw * SAFE
        end
        rho = min.(rho, new); push!(hist, maximum(rho))
        hist[end] > hist[end-1] * (1 - 1e-4) && break
    end
    out = Dict{String,Any}()
    for s in net.sources
        cs = c[s.name][k_]
        naive = (s.pu + sqrt(2) * 1.5) / abs(s.z_series * net.Y_base)
        term = abs(a[s.name]) + sum(abs.(cs) .* Jmax0)
        acc = abs(a[s.name]) .* cos.(psis .- angle(a[s.name]))
        for q in eachindex(D); support!(acc, psis, phiJ[q] + angle(cs[q]), wJ[q], abs(cs[q]) * RJ[q]); end
        sect = max(maximum(acc), 0.0) / cw * SAFE
        # the same sectors with the device terms bounded separately (triangle inequality per device)
        sep = abs(a[s.name]) + sum(abs.(cs) .* RJ)
        out[s.name] = (; naive, term = term * SAFE, sep = sep * SAFE, sect, a0 = abs(a[s.name]))
    end
    return (; out, rho_hist = hist, V0, D, wJ)
end

"Lifted |I_s|^2 of source s at period 1 in relaxation R."
function lifted_isq(R, s, ang)
    E = s.pu * cis(deg2rad(angle_at(ang, s.name, 1))); z = s.z_series * net.Y_base; b = s.bus
    return (abs2(E) - 2 * (real(E) * R.e[b, 1] + imag(E) * R.f[b, 1]) + R.Wee[b, 1] + R.Wff[b, 1]) / abs2(z)
end

"OBBT: maximize each lifted |I_s|^2 over the restoration relaxation with the current bounds; checked upper bounds."
function obbt(ang, kw, rub, bnd; rounds = 4)
    hist = Any[]; bnd = copy(bnd)
    for rd in 1:rounds
        old = copy(bnd)
        for s in net.sources
            kwb = (; kw..., source_ampacity_pu = copy(bnd))
            t0 = time()
            R = relax_margin(net, ang, kwb; k = KV, cuts = CUTS, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, solver_attributes = CLARABEL_DEFAULTS)
            isq = lifted_isq(R, s, ang)
            @objective(R.model, Min, -isq)
            optimize!(R.model); ts = time() - t0
            tighten_boxes!(R, net, kwb, ang)
            st = termination_status(R.model)
            if st in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
                fb = dual_status(R.model) in (MOI.INFEASIBILITY_CERTIFICATE, MOI.NEARLY_INFEASIBILITY_CERTIFICATE) ? dual_bound(R.model; c = zero(AffExpr), box = R.box).bound : NaN
                push!(hist, (; round = rd, source = s.name, status = st, solver = NaN, ub = NaN, farkas = fb, solve_s = ts))
                fb > 0 && return (; bnd, hist, certificate = (; source = s.name, round = rd, farkas = fb))
                continue
            end
            db = has_duals(R.model) ? dual_bound(R.model; box = R.box) : (; bound = -Inf, box_term = NaN, rounding = NaN)
            ubsq = -db.bound
            solver = has_values(R.model) ? sqrt(max(value(isq), 0.0)) : NaN
            push!(hist, (; round = rd, source = s.name, status = st, solver, ub = ubsq >= 0 ? sqrt(ubsq) * (1 + 1e-12) : -1.0, farkas = NaN, solve_s = ts,
                         box_term = db.box_term, rounding = db.rounding))
            ubsq < 0 && return (; bnd, hist, certificate = (; source = s.name, round = rd, farkas = -ubsq))   # max |I|^2 < 0: empty domain
            isfinite(ubsq) && (bnd[s.name] = min(bnd[s.name], sqrt(ubsq) * (1 + 1e-12)))
        end
        maximum(old[k] / bnd[k] - 1 for k in keys(bnd)) < 0.005 && break
    end
    return (; bnd, hist, certificate = nothing)
end

src_currents(x) = Dict(s.name => abs(complex(x.raw.ir[s.name][1], x.raw.ii[s.name][1])) for s in net.sources)

function restoration_test(ang, kw, rub, bnd)
    kwb = (; kw..., source_ampacity_pu = copy(bnd))
    c = certify(net, ang, kwb; k = KV, cuts = CUTS, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, tighten = true)
    S = net.S_base_kVA
    infeas = c.status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
    lb = infeas ? NaN : (isfinite(c.db.bound) ? c.db.bound * S : -Inf)
    cls = infeas ? (c.farkas.bound > 0 ? "no-backflow impossible (checked Farkas)" : "solver-only infeasibility (ray not validated)") :
          lb > 0 ? "no-backflow impossible (checked lower bound on rho > 0)" : "unresolved (rho lower bound <= 0)"
    return (; c, lb, cls, infeas, rho_relax = infeas ? NaN : -c.primal)
end

function scb_main()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    env2 = Dict(String(r[1]) => parse(Float64, r[4]) for r in (split(l, ',') for l in readlines(PFX * "_envelope_v2.csv")[2:end]))
    cw = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "certified_window.csv"))[2:end]]
    open(io -> println(io, "theta_deg,factor,r_deg,source,envelope_v2_A,cap_1.5x_A,naive_A,termwise_A,sector_separate_A,sector_combined_A,no_load_term_A,obbt_A,final_A,final_pu,obbt_rounds,witness_feasible_A,witness_restoration_A"), SCB * "_bounds.csv", "w")
    open(io -> println(io, "theta_deg,factor,r_deg,L_exact_kW,exact_equations,rho_exact_kW,rho_max_kW,stage,relax_status,rho_relax_kW,rho_lb_kW,box_term_kW,rounding_kW,primal_residual,dual_psd_violation,dual_sign_violation,farkas_bound,build_s,solve_s,bound_s,class"), SCB * "_queries.csv", "w")
    t0 = time(); RM = response_matrix()
    say(@sprintf("network response: Re(Y) Cholesky min pivot %.3e (positive definite => Y nonsingular); |YZ - I| = %.1e (BigFloat 256); %.1f s", Float64(RM.minpiv), Float64(RM.resid), time() - t0))
    for th in (292.5, 157.5)
        rs = parse(Float64, only(r for r in cw if parse(Float64, r[1]) == th && r[2] == "sources+lines 1.2x" && r[4] == "")[5])
        fc = parse(Float64, only(r for r in cw if parse(Float64, r[1]) == th && r[2] == "sources only 1.5x" && r[4] == "kcl0+jbound")[11])
        say(@sprintf("\n== theta %.1f deg: r* = %.6f deg (locally measured); source-only 1.5x certified point %.3f r*", th, rs, fc))
        for fac in (0.98, 1.02, fc, 1.10, 1.25)
            r = fac * rs; ang = at(r, th)
            xm = solve_multiperiod(net, ang; kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
            vm = validate_dispatch(net, ang, kw1, xm; k = KV)
            L = xm.status in RAY_OK && vm.equations <= 1e-6 ? vm.margin_kW : NaN
            xr = solve_multiperiod(net, ang; kw1..., objective = :margin, margin_band_hard = true, margin_reg = 0.0)
            vr = validate_dispatch(net, ang, kw1, xr; k = KV)
            okr = xr.status in RAY_OK && vr.equations <= 1e-6 && vr.viol["band"] <= 1e-6
            rho_ex = okr ? max(0.0, -vr.min_import_kW) : NaN
            rub = (okr ? rho_ex : 10.0) + 1.0
            say(@sprintf("  -- %.3f r* (r = %.6f deg): exact margin L %+.4f kW (%s, equations %.1e); exact restoration rho %.4f kW (%s, validated %s); rho_max %.4f kW",
                         fac, r, L, xm.status, vm.equations, rho_ex, xr.status, okr, rub))
            ta = time(); AB = analytic_bounds(RM, ang, kw1); tab = time() - ta
            bnd0 = Dict(s => AB.out[s].sect for s in keys(AB.out))
            say(@sprintf("     analytic bounds (%.1f s; disc radius around the no-load voltages %s pu):", tab, join([@sprintf("%.4f", v) for v in AB.rho_hist], " -> ")))
            wf = isfinite(L) ? src_currents(xm) : Dict(); wr = okr ? src_currents(xr) : Dict()
            for s in net.sources
                o = AB.out[s.name]
                say(@sprintf("       %s: naive %.1f A | termwise %.2f A | sectors, devices separately %.2f A | sectors, combined %.2f A (no-load term %.2f A) | envelope %.2f A, 1.5x cap %.2f A | witnesses %.2f / %.2f A",
                             s.name, o.naive * IbA, o.term * IbA, o.sep * IbA, o.sect * IbA, o.a0 * IbA, env2["source:"*s.name] * IbA, 1.5env2["source:"*s.name] * IbA,
                             get(wf, s.name, NaN) * IbA, get(wr, s.name, NaN) * IbA))
            end
            for (stage, bnd, ob) in (("analytic", bnd0, nothing), ("obbt", nothing, true))
                tb = 0.0; obr = nothing
                if ob !== nothing
                    tq = time(); obr = obbt(ang, kw1, rub, bnd0); tb = time() - tq; bnd = obr.bnd
                    for h in obr.hist
                        say(@sprintf("     OBBT round %d %s: %s, solver max %.3f A, checked upper bound %.3f A (box term %.2e, rounding %.1e), %.1f s",
                                     h.round, h.source, h.status, h.solver * IbA, h.ub * IbA, get(h, :box_term, NaN), get(h, :rounding, NaN), h.solve_s))
                    end
                    if obr.certificate !== nothing
                        say(@sprintf("     OBBT problem itself certified empty (source %s, round %d, checked %.3e > 0): no-backflow impossible", obr.certificate.source, obr.certificate.round, obr.certificate.farkas))
                    end
                end
                rt = restoration_test(ang, kw1, rub, bnd); c = rt.c
                say(@sprintf("     restoration [%s bounds %s A]: %s, rho relax %.4f kW, checked lower bound %.4f kW (box term %.2e kW, rounding %.1e kW, residual %.1e, dual psd/sign %.1e/%.1e%s), %.1f + %.1f s -> %s",
                             stage, join([@sprintf("%.2f", bnd[s.name] * IbA) for s in net.sources], "/"), c.status, rt.rho_relax, rt.lb, c.box_term_kW, c.rounding_kW,
                             c.primal_residual, c.psd_dual_violation, c.dual_sign_violation, rt.infeas ? @sprintf(", Farkas %.3e", c.farkas.bound) : "", c.build_s, c.solve_s, rt.cls))
                scb_csv("queries", [@sprintf("%.1f", th), @sprintf("%.3f", fac), @sprintf("%.6f", r), @sprintf("%.4f", L), @sprintf("%.1e", vm.equations), @sprintf("%.4f", rho_ex), @sprintf("%.4f", rub),
                                    stage, string(c.status), @sprintf("%.4f", rt.rho_relax), @sprintf("%.4f", rt.lb), @sprintf("%.3e", c.box_term_kW), @sprintf("%.1e", c.rounding_kW),
                                    @sprintf("%.1e", c.primal_residual), @sprintf("%.1e", c.psd_dual_violation), @sprintf("%.1e", c.dual_sign_violation), @sprintf("%.3e", c.farkas.bound),
                                    @sprintf("%.2f", c.build_s), @sprintf("%.2f", c.solve_s), @sprintf("%.1f", stage == "obbt" ? tb : tab), rt.cls])
                if stage == "obbt"
                    for s in net.sources
                        o = AB.out[s.name]; nr = obr === nothing ? 0 : maximum((h.round for h in obr.hist); init = 0)
                        scb_csv("bounds", [@sprintf("%.1f", th), @sprintf("%.3f", fac), @sprintf("%.6f", r), s.name, @sprintf("%.3f", env2["source:"*s.name] * IbA), @sprintf("%.3f", 1.5env2["source:"*s.name] * IbA),
                                           @sprintf("%.1f", o.naive * IbA), @sprintf("%.3f", o.term * IbA), @sprintf("%.3f", o.sep * IbA), @sprintf("%.3f", o.sect * IbA), @sprintf("%.3f", o.a0 * IbA),
                                           @sprintf("%.3f", bnd[s.name] * IbA), @sprintf("%.3f", bnd[s.name] * IbA), @sprintf("%.6f", bnd[s.name]), string(nr),
                                           @sprintf("%.3f", get(wf, s.name, NaN) * IbA), @sprintf("%.3f", get(wr, s.name, NaN) * IbA)])
                        isfinite(get(wf, s.name, NaN)) && wf[s.name] > bnd[s.name] && say("     !! feasible witness exceeds the derived bound at ", s.name)
                        isfinite(get(wr, s.name, NaN)) && wr[s.name] > bnd[s.name] && say("     !! restoration witness exceeds the derived bound at ", s.name)
                    end
                end
            end
        end
    end
end
scb_main()
