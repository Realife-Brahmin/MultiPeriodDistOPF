# Attainable currents in the exact AC model (2026-10-22): how large can source and line currents get at
# feasible operating points, and do cushioned empirical limits drawn from them help the SDP relaxation?
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/attainable_currents.jl
# Case: IEEE 123-node, sources 1, 3, 4, PV reactive power dispatched, batteries idle, the limiting hour of
# the recent experiment (hour 12 of the 292.5-deg idle edge), alone. Every original constraint is kept
# (network equations, source EMFs behind their impedances, the 0.95-1.05 band, PV reactive limits,
# no-backflow); no current limit while the maxima are discovered.
# Study domain: substation 1 is the reference (0 deg); substations 3 and 4 move freely in [-5, +5] deg.
# Each maximization (|I|^2 of one source or one line) is nonconvex: the results are the largest
# validated currents found from the starts tried, not upper bounds.
# Writes results/multiperiod/attainable_currents_{runs,lines,limits}.csv and a log.
empty!(ARGS); push!(ARGS, "none")
include(joinpath(@__DIR__, "feasibility_poc.jl"))     # its setup helpers and the hour-1 conventions (mode "none")
banded(net) = [b for b in net.buses if !(b in Set(s.bus for s in net.sources))]

const DOMAIN = (-5.0, 5.0)
const REFNAME = "subs1"

function currents_study()
    th0 = 292.5
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw24 = problem(s24, :nobatt)
    reach = measured_reach()
    meas, _ = reach[("nobatt", th0)]
    at(r, th) = Dict("subs1" => 0.0, "subs3" => r * cosd(th), "subs4" => r * sind(th))
    xd = solve_multiperiod(net, at(0.98meas, th0); kw24..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
    tstar = argmin([minimum(xd.P_subs_kW[n][t] for n in keys(xd.P_subs_kW)) for t in 1:24])
    kw1 = (; pvs = s24.P, pv_q = :optimize, C_B = s24.CB, ideal_sources = false, battery = :idle, tol = IPOPT_TOL, constr_viol_tol = IPOPT_TOL,
           T = 1, dt_h = 1.0, load_factor = [s24.lam[tstar]], price = [s24.c[tstar]], pv_factor = [s24.F[tstar]], gamma = 0.0, terminal = :free)
    S = net.S_base_kVA; Vb = net.kV_base; Ib = S / Vb; Zb = Vb^2 / (S / 1000)       # A per pu, ohm per pu
    Sload = abs(sum(complex(d.kW, d.kvar) for d in net.loads)) * s24.lam[tstar]
    Irep = Sload / Vb
    say(@sprintf("bases: S_base %.0f kVA per phase, V_base %.4f kV line-to-neutral (%.2f kV line-to-line), I_base %.3f A, Z_base %.2f ohm", S, Vb, sqrt(3) * Vb, Ib, Zb))
    say(@sprintf("hour %d: load factor %.4f, PV factor %.4f; feeder apparent demand %.1f kVA per phase (%.1f kVA three-phase): a demand base would give I_rep = %.2f A = %.4f pu", tstar,
                 s24.lam[tstar], s24.F[tstar], Sload, 3Sload, Irep, Irep / Ib))
    lines = [l.name for l in net.lines]
    srcs = [s.name for s in net.sources]
    # validated schedules already in hand: the hour alone at the inside/outside points of hour1, and the idle
    # daily feasible schedules of the T = 24 experiment (equal angles, 98% of both reaches, beyond prediction)
    env = Dict{String,Tuple{Float64,String}}()                   # element => (largest |I| pu, where)
    function record!(x, tag; T = 1)
        for l in lines, t in 1:T
            v = abs(complex(x.raw.lr[l][t], x.raw.li[l][t])); v > get(env, l, (0.0, ""))[1] && (env[l] = (v, tag))
        end
        for n in srcs, t in 1:T
            v = abs(complex(x.raw.ir[n][t], x.raw.ii[n][t])); v > get(env, "source:" * n, (0.0, ""))[1] && (env["source:" * n] = (v, tag))
        end
    end
    hourstar = 0.118008                                          # the hour's local boundary along 292.5 deg (hour1)
    starts = Any[("equal angles (cold)", nothing, Dict("subs1" => 0.0, "subs3" => 0.0, "subs4" => 0.0))]
    for (nm, r, th, kwx, T) in (("hour inside 292.5", 0.98hourstar, th0, kw1, 1), ("hour outside 292.5", 1.02hourstar, th0, kw1, 1))
        x = solve_multiperiod(net, at(r, th); kwx..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        v = validate_dispatch(net, at(r, th), kwx, x; k = KV)
        say(@sprintf("schedule %-20s margin %+.3f kW, equations %.1e", nm, v.margin_kW, v.equations))
        v.margin_kW >= 0 && v.equations <= 1e-6 && (record!(x, nm); push!(starts, (nm, x, at(r, th))))
    end
    for (q, r, th) in (("daily equal", 0.0, 0.0), ("daily 98% 292.5", 0.98meas, th0), ("daily 98% 157.5", 0.98reach[("nobatt", 157.5)][1], 157.5),
                       ("daily beyond pred 292.5", (meas + reach[("nobatt", th0)][2]) / 2, th0))
        x = solve_multiperiod(net, at(r, th); kw24..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        v = validate_dispatch(net, at(r, th), kw24, x; k = KV)
        v.margin_kW >= 0 && v.equations <= 1e-6 && record!(x, q; T = 24)
        say(@sprintf("schedule %-20s margin %+.3f kW, equations %.1e (24 h)", q, v.margin_kW, v.equations))
    end
    for th in (0.0, 90.0, 225.0)                                 # more starting points: 95% of the daily reach in other directions
        r = 0.95reach[("nobatt", th)][1]
        x = solve_multiperiod(net, at(r, th); kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        v = validate_dispatch(net, at(r, th), kw1, x; k = KV)
        v.margin_kW >= 0 && v.equations <= 1e-6 && (record!(x, "hour 95% $(th)"); push!(starts, ("hour 95% $(th)", x, at(r, th))))
    end
    runs = joinpath(OUT_DIR, "attainable_currents_runs.csv")
    open(io -> println(io, "target,start,status,iterations,solve_s,build_s,I_pu,I_A,I_rep_pu,d3_deg,d4_deg,domain_active,equations,max_violation,min_import_kW," *
                       "P1_kW,Q1_kvar,S1_kVA,I1_A,P3_kW,Q3_kvar,S3_kVA,I3_A,P4_kW,Q4_kvar,S4_kVA,I4_A,losses_kW,vmin_band,vmax_band,V1s,V3s,V4s,max_line,max_line_A"), runs, "w")
    best = Dict{String,Any}()
    function maxrun(kind, nm)
        for (sname, sx, sang) in starts
            free = Dict("subs3" => DOMAIN, "subs4" => DOMAIN)
            t0 = time()
            x = solve_multiperiod(net, sang; kw1..., objective = :maxcurrent, current_target = (kind, nm), free_angles = free, start = sx)
            tt = time() - t0
            ok = x.status in RAY_OK
            ang = x.angles_deg
            v = ok ? validate_dispatch(net, ang, kw1, x; k = KV) : nothing
            valid = ok && v.equations <= 1e-6 && v.viol["band"] <= 1e-6 && v.viol["backflow"] <= 1e-6
            Ipu = !ok ? NaN : kind === :source ? abs(complex(x.raw.ir[nm][1], x.raw.ii[nm][1])) : abs(complex(x.raw.lr[nm][1], x.raw.li[nm][1]))
            valid && record!(x, "max $(kind) $(nm) from $(sname)")
            if ok
                sp(n) = complex(x.P_subs_kW[n][1], x.Q_subs_kvar[n][1])
                Isrc(n) = abs(complex(x.raw.ir[n][1], x.raw.ii[n][1])) * Ib
                loss = sum(real(sp(n)) for n in srcs) - sum(s24.lam[tstar] * d.kW for d in net.loads) + sum(s24.F[tstar] * p.kW for p in s24.P)
                nb = banded(net); vb = [abs(x.V_pu[b][1]) for b in nb]
                ml = argmax(l -> abs(complex(x.raw.lr[l][1], x.raw.li[l][1])), lines)
                dact = any(abs(ang[n] - lim) < 1e-6 for n in ("subs3", "subs4") for lim in DOMAIN)
                row = [string(kind, ":", nm), replace(sname, "," => ";"), string(x.status), string(something(x.iterations, -1)), @sprintf("%.2f", something(x.time, NaN)), @sprintf("%.2f", tt - something(x.time, 0.0)),
                       @sprintf("%.4f", Ipu), @sprintf("%.1f", Ipu * Ib), @sprintf("%.4f", Ipu * Ib / Irep), @sprintf("%.5f", ang["subs3"]), @sprintf("%.5f", ang["subs4"]), string(dact),
                       @sprintf("%.1e", v.equations), @sprintf("%.1e", v.max_violation), @sprintf("%.3f", v.min_import_kW),
                       (@sprintf("%.2f,%.2f,%.2f,%.1f", real(sp(n)), imag(sp(n)), abs(sp(n)), Isrc(n)) for n in ("subs1", "subs3", "subs4"))...,
                       @sprintf("%.2f", loss), @sprintf("%.5f", minimum(vb)), @sprintf("%.5f", maximum(vb)),
                       (@sprintf("%.5f", abs(x.V_pu[s.bus][1])) for s in sort(net.sources, by = s -> s.name))...,
                       ml, @sprintf("%.1f", abs(complex(x.raw.lr[ml][1], x.raw.li[ml][1])) * Ib)]
                open(io -> println(io, join(row, ',')), runs, "a")
                say(@sprintf("  max %-6s %-6s from %-22s %-22s |I| %.4f pu = %7.1f A (%.3f demand-base pu) at d3 %+.4f d4 %+.4f, eq %.1e, min import %.2f kW, losses %.1f kW, %.1f s%s",
                             kind, nm, sname, x.status, Ipu, Ipu * Ib, Ipu * Ib / Irep, ang["subs3"], ang["subs4"], v.equations, v.min_import_kW, loss, tt, dact ? " (domain bound active)" : ""))
                if valid && Ipu > get(best, string(kind, ":", nm), (; I = -Inf)).I
                    best[string(kind, ":", nm)] = (; I = Ipu, x, ang, start = sname)
                end
            else
                say(@sprintf("  max %-6s %-6s from %-22s %s, %.1f s", kind, nm, sname, x.status, tt))
            end
        end
    end
    for n in ("subs1", "subs3", "subs4"); maxrun(:source, n); end
    # lines: the fictitious-loss lines of the SDP experiment and the largest internal currents found so far
    internal = sort([(v[1], l) for (l, v) in env if !startswith(l, "source:") && !startswith(l, "poi")]; rev = true)
    targets = unique([["l115", "l116", "l46", "l110", "l111"]; [l for (_, l) in internal[1:3]]])
    say("line targets: ", join(targets, ", "))
    for l in targets; maxrun(:line, l); end
    # amperes check: replay the largest source-current dispatch in OpenDSS and compare one line's current there
    k = argmax(kk -> best[kk].I, collect(keys(best)))
    bx = best[k]
    rp = replay_opendss(net, bx.x, bx.ang; load_factor = kw1.load_factor, pvs = kw1.pvs, pv_factor = kw1.pv_factor)
    for l in ("l115", "l116")
        ODD.Circuit.SetActiveElement("Line." * l)
        ia = ODD.CktElement.CurrentsMagAng()[1]
        say(@sprintf("amperes check (%s dispatch, OpenDSS replay dP %.1e kW): line %s OpenDSS %.2f A, model %.4f pu x %.3f A = %.2f A", k, rp.dP, l, ia,
                     abs(complex(bx.x.raw.lr[l][1], bx.x.raw.li[l][1])), Ib, abs(complex(bx.x.raw.lr[l][1], bx.x.raw.li[l][1])) * Ib))
    end
    open(joinpath(OUT_DIR, "attainable_currents_lines.csv"), "w") do io
        println(io, "element,max_I_pu,max_I_A,max_I_demand_pu,where")
        for (el, (v, w)) in sort(collect(env), by = p -> -p[2][1])
            println(io, @sprintf("%s,%.5f,%.2f,%.5f,%s", el, v, v * Ib, v * Ib / Irep, replace(w, ',' => ';')))
        end
    end
    # exploratory limits from the envelope, cushioned; applied to the exact model and to the relaxation's moments
    floorpu = 0.05                                               # 4.3 A: no line limited below this
    scen = [("src x1.5", 1.5, nothing), ("src+lines x1.5", 1.5, 1.5), ("src+lines x1.2", 1.2, 1.2)]
    lim = joinpath(OUT_DIR, "attainable_currents_limits.csv")
    open(io -> println(io, "scenario,point,L_kW,nlp_status,orig_equations,limits_violated_by_validated,variant,relax_status,primal_kW,gap_kW,primal_residual,psd_primal_min_eig,fictitious_loss_kW,vsq_gap,farkas_bound,lifted_violation,build_s,solve_s,class"), lim, "w")
    for (sn, cs, cl) in scen
        samp = Dict(n => cs * env["source:" * n][1] for n in srcs)
        lamp = cl === nothing ? nothing : Dict(l => max(cl * get(env, l, (0.0, ""))[1], floorpu) for l in lines)
        kwl = (; kw1..., source_ampacity_pu = samp, ampacity_pu = lamp)
        say(@sprintf("\nscenario %s: source limits %s A; %s", sn, join([@sprintf("%s %.1f", n, samp[n] * Ib) for n in srcs], ", "),
                     lamp === nothing ? "no line limits" : @sprintf("line limits %.1f-%.1f A (floor %.1f A)", minimum(values(lamp)) * Ib, maximum(values(lamp)) * Ib, floorpu * Ib)))
        for (pname, r) in (("inside", 0.98hourstar), ("outside", 1.02hourstar))
            x = solve_multiperiod(net, at(r, th0); kwl..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
            v = validate_dispatch(net, at(r, th0), kwl, x; k = KV)
            L = x.status in RAY_OK && v.equations <= 1e-6 ? v.margin_kW : NaN
            for (vname, cv) in (("base", Symbol[]), ("both", [:kcl0, :jbound]))
                c = certify(net, at(r, th0), kwl; k = KV, cuts = cv, m_lo_kW = 0.0)
                lv = isfinite(L) && L >= 0 ? lifted_violation(c.R, net, x, kwl; margin_kW = L) : (; worst = NaN, at = "")
                infeas = c.status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
                cls = isfinite(L) && L >= 0 ? "feasible (validated L)" : (infeas && c.farkas.bound > 0 ? "infeasible (validated Farkas)" : (infeas ? "solver infeasible, certificate not validated" : "unresolved"))
                row = [sn, pname, @sprintf("%.4f", L), string(x.status), @sprintf("%.1e", v.equations), @sprintf("%.1e", v.viol["ampacity"]), vname, string(c.status),
                       @sprintf("%.3f", c.primal), @sprintf("%.2e", c.gap), @sprintf("%.2e", c.primal_residual), @sprintf("%.2e", c.psd_primal_min_eig),
                       @sprintf("%.2f", c.fictitious_loss_kW), @sprintf("%.3e", c.vsq_gap), @sprintf("%.3e", c.farkas.bound), @sprintf("%.2e", lv.worst),
                       @sprintf("%.2f", c.build_s), @sprintf("%.2f", c.solve_s), cls]
                open(io -> println(io, join(row, ',')), lim, "a")
                say(@sprintf("  %-7s %-4s L %+.3f kW | relaxation %-22s primal %8.2f kW, fictitious %8.2f kW, |V|^2 gap %.2e, farkas %.2e, lifted %.1e, %.1f+%.1f s -> %s",
                             pname, vname, L, string(c.status), c.primal, c.fictitious_loss_kW, c.vsq_gap, c.farkas.bound, lv.worst, c.build_s, c.solve_s, cls))
            end
        end
    end
end
# Farkas check for the "solver infeasible" outside point: the ray's certificate value against its box term,
# with Clarabel's default and tighter infeasibility tolerances (the scenario src+lines x1.2, no additions).
function farkas_check()
    env = Dict(r[1] => parse(Float64, r[2]) for r in (split(l, ',') for l in readlines(joinpath(OUT_DIR, "attainable_currents_lines.csv"))[2:end]))
    s24 = setup(net, pvs; T = 24, pv = "qfree")
    kw1 = (; pvs = s24.P, pv_q = :optimize, C_B = s24.CB, ideal_sources = false, battery = :idle, tol = IPOPT_TOL, constr_viol_tol = IPOPT_TOL,
           T = 1, dt_h = 1.0, load_factor = [s24.lam[12]], price = [s24.c[12]], pv_factor = [s24.F[12]], gamma = 0.0, terminal = :free)
    kwl = (; kw1..., source_ampacity_pu = Dict(s.name => 1.2 * env["source:" * s.name] for s in net.sources),
           ampacity_pu = Dict(l.name => max(1.2 * get(env, l.name, 0.0), 0.05) for l in net.lines))
    ang = Dict("subs1" => 0.0, "subs3" => 1.02 * 0.118008 * cosd(292.5), "subs4" => 1.02 * 0.118008 * sind(292.5))
    for (nm, opts) in (("default", []), ("tol_infeas 1e-10", [("tol_infeas_abs", 1e-10), ("tol_infeas_rel", 1e-10)]),
                       ("tol_infeas 1e-12, max_iter 1000", [("tol_infeas_abs", 1e-12), ("tol_infeas_rel", 1e-12), ("max_iter", 1000)]))
        c = certify(net, ang, kwl; k = KV, cuts = Symbol[], m_lo_kW = 0.0, solver_attributes = [CLARABEL_DEFAULTS; opts])
        say(@sprintf("farkas check (%s): %s / %s; ray certificate value %.3e, box term %.3e, bound %.3e (> 0 certifies), dual PSD violation %.1e",
                     nm, c.status, c.dual_status, c.farkas.constant, c.farkas.box_term, c.farkas.bound, c.farkas.psd_viol))
    end
end
currents_study()
farkas_check()
