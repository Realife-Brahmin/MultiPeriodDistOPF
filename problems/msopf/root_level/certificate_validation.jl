# Certificate validation under exploratory current limits, and a challenge to those limits (2026-10-22).
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/certificate_validation.jl
# Frozen comparison cases: IEEE 123-node, sources 1, 3, 4, hour 12 alone, batteries idle, PV reactive power
# dispatched; the hour-1 points 98% / 102% of the hour's boundary along 292.5 deg (r* = 0.118008 deg);
# limits = cushion x the frozen envelope in attainable_currents_lines.csv, no line below 0.05 pu.
#   A: sources and lines 1.2x, no nodal additions          B: sources and lines 1.5x, both nodal additions
# Parts: (0) sign and PSD-dual conventions on toy problems; (1-2) A and B with the loose and the tightened
# verifier boxes, and what dominates the Farkas correction; (3) the restoration problem (band and limits
# hard, min rho with P_s >= -rho) as a finite-objective certificate; (4) a challenge to the envelope:
# cross-seeded current maximizations with no limit; (5) cushion 1.2 / 1.5 / 2 x {sources, sources and
# lines}; (6) 157.5 deg, only if some certificate validates. Writes certificate_validation_*.csv and a log.
empty!(ARGS); push!(ARGS, "none")
include(joinpath(@__DIR__, "feasibility_poc.jl"))
banded(net) = [b for b in net.buses if !(b in Set(s.bus for s in net.sources))]
using Random

const RSTAR = 0.118008                    # hour 12's local boundary along 292.5 deg (hour1)
const TH0 = 292.5
const PFX = joinpath(OUT_DIR, "certificate_validation")
const FLOOR = 0.05
const DOMAIN = (-5.0, 5.0)                # deg, substations 3 and 4 (substation 1 is the reference)
csvw(file, header) = open(io -> println(io, header), PFX * "_" * file * ".csv", "w")
csva(file, row) = open(io -> println(io, join(row, ',')), PFX * "_" * file * ".csv", "a")
cl(x) = replace(string(x), ',' => ';')

hour_kw(s24, t) = (; pvs = s24.P, pv_q = :optimize, C_B = s24.CB, ideal_sources = false, battery = :idle, tol = IPOPT_TOL, constr_viol_tol = IPOPT_TOL,
                   T = 1, dt_h = 1.0, load_factor = [s24.lam[t]], price = [s24.c[t]], pv_factor = [s24.F[t]], gamma = 0.0, terminal = :free)
at(r, th) = Dict("subs1" => 0.0, "subs3" => r * cosd(th), "subs4" => r * sind(th))
read_env(f) = Dict(String(r[1]) => parse(Float64, r[2]) for r in (split(l, ',') for l in readlines(f)[2:end]))
limits(env, cs, cl_) = (; source_ampacity_pu = Dict(s.name => cs * env["source:" * s.name] for s in net.sources),
                         ampacity_pu = cl_ === nothing ? nothing : Dict(l.name => max(cl_ * get(env, l.name, 0.0), FLOOR) for l in net.lines))

function conventions()
    say("== (0) conventions on toy problems (Clarabel, the same dual_bound)")
    m = Model(Clarabel.Optimizer); set_silent(m); @variable(m, x); c = @constraint(m, Symmetric([1 x; x 0.25]) in PSDCone()); @objective(m, Min, x); optimize!(m)
    db = dual_bound(m; box = Dict(x => (-10.0, 10.0)))
    raw = MOI.get(backend(m), MOI.ConstraintDual(), index(c)); jd = dual(c)
    say(@sprintf("  min x s.t. [1 x; x 0.25] PSD: optimum %.6f, dual_bound %.6f (must be <= -0.5 and close); raw triangle dual %s, JuMP dual matrix %s",
                 objective_value(m), db.bound, string(round.(raw, digits = 4)), string(round.(collect(jd), digits = 4))))
    m2 = Model(Clarabel.Optimizer); set_silent(m2); @variable(m2, y); @constraint(m2, Symmetric([1 y; y 0.25]) in PSDCone()); @constraint(m2, y >= 1); optimize!(m2)
    f2 = dual_bound(m2; c = zero(AffExpr), box = Dict(y => (-10.0, 10.0)))
    m3 = Model(Clarabel.Optimizer); set_silent(m3); @variable(m3, z); @constraint(m3, Symmetric([1 z; z 0.25]) in PSDCone()); @constraint(m3, z <= -1); optimize!(m3)
    f3 = dual_bound(m3; c = zero(AffExpr), box = Dict(z => (-10.0, 10.0)))
    say(@sprintf("  infeasible (y >= 1): %s / %s, Farkas bound %.4e (> 0 certifies); infeasible (z <= -1): %s / %s, Farkas bound %.4e",
                 termination_status(m2), dual_status(m2), f2.bound, termination_status(m3), dual_status(m3), f3.bound))
    return (; ok = db.bound <= objective_value(m) + 1e-6 && f2.bound > 0 && f3.bound > 0)
end

"One case: exact margin (validated L), the relaxation with tightened boxes (and the loose ones), lifted and box checks."
function run_case(tag, kw, pname, r, th, cuts; detail = 12)
    ang = at(r, th)
    x = solve_multiperiod(net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
    v = validate_dispatch(net, ang, kw, x; k = KV)
    L = x.status in RAY_OK && v.equations <= 1e-6 ? v.margin_kW : NaN
    c = certify(net, ang, kw; k = KV, cuts, m_lo_kW = 0.0, m_hi_kW = 2000.0, tighten = true, detail)
    lv = isfinite(L) && L >= 0 ? lifted_violation(c.R, net, x, kw; margin_kW = L) : (; worst = NaN, at = "", box_worst = NaN, box_at = "")
    infeas = c.status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
    dloose = has_duals(c.R.model) && !infeas ? dual_bound(c.R.model; box = c.box_loose) : (; bound = NaN, box_term = NaN)
    fk = c.farkas
    cls = isfinite(L) && L >= 0 ? "feasible (validated exact dispatch)" :
          infeas ? (fk.bound > 0 ? "INFEASIBLE under the limits (validated Farkas, tightened boxes)" : "solver-only infeasibility (ray not validated)") :
          (isfinite(c.db.bound) && -c.db.bound * net.S_base_kVA < -UTOL ? "INFEASIBLE under the limits (validated dual bound)" : "unresolved (relaxation feasible)")
    row = [tag, pname, @sprintf("%.6f", r), @sprintf("%.1f", th), join(string.(cuts), "+"), @sprintf("%.4f", L), string(x.status), @sprintf("%.1e", v.equations),
           @sprintf("%.1e", v.viol["ampacity"]), string(c.status), string(c.dual_status),
           infeas ? "n/a" : @sprintf("%.4f", c.primal), infeas ? "n/a" : @sprintf("%.4f", -c.db.bound * net.S_base_kVA),
           infeas ? "n/a" : @sprintf("%.4f", -dloose.bound * net.S_base_kVA), infeas ? "n/a" : @sprintf("%.2f", c.fictitious_loss_kW), infeas ? "n/a" : @sprintf("%.2e", c.vsq_gap),
           infeas ? "n/a" : @sprintf("%.2e", c.primal_residual),
           @sprintf("%.3e", fk.constant), @sprintf("%.3e", fk.box_term_loose), @sprintf("%.3e", fk.bound_loose), @sprintf("%.3e", fk.box_term), @sprintf("%.3e", fk.bound),
           @sprintf("%.1e", fk.psd_viol), @sprintf("%.2e", lv.worst), @sprintf("%.2e", lv.box_worst), cl(lv.box_at), @sprintf("%.2f", c.build_s), @sprintf("%.2f", c.solve_s), cls]
    csva("cases", row)
    say(@sprintf("  %-26s %-7s %-12s L %+8.3f | %-18s %s", tag, pname, join(string.(cuts), "+"), L, string(c.status),
                 infeas ? @sprintf("Farkas: value %.3e, box term loose %.3e -> bound %.3e; tightened %.3e -> bound %.3e", fk.constant, fk.box_term_loose, fk.bound_loose, fk.box_term, fk.bound) :
                          @sprintf("primal %.3f kW, dual bound (tight) %.3f kW, (loose) %.3e kW, fictitious %.2f kW", c.primal, -c.db.bound * net.S_base_kVA, -dloose.bound * net.S_base_kVA, c.fictitious_loss_kW)))
    isfinite(L) && L >= 0 && say(@sprintf("      lifted inside dispatch: constraints %.1e, outside a tightened box by %.1e (%s)", lv.worst, lv.box_worst, lv.box_at))
    say("      -> ", cls)
    if infeas && detail > 0
        for (k, (val, nm, a, l, u)) in enumerate(fk.top)
            csva("contributors", [tag, pname, join(string.(cuts), "+"), "tightened", string(k), cl(nm), @sprintf("%.4e", a), @sprintf("%.6g", l), @sprintf("%.6g", u), @sprintf("%.4e", val)])
            k <= 6 && say(@sprintf("      box term #%d: %-34s coefficient %+.3e, box [%.5g, %.5g], contributes %.3e", k, nm, a, l, u, val))
        end
        say("      certificate value by constraint kind: ", join([@sprintf("%s %+.3e", k, v) for (k, v) in sort(collect(fk.csrc), by = p -> -abs(p[2]))], ", "))
    end
    return (; c, L, x, cls, infeas)
end

"The restoration problem: band and current limits hard, min rho >= 0 with P_s >= -rho."
function restoration(tag, kw, pname, r, th, cuts)
    ang = at(r, th)
    x = solve_multiperiod(net, ang; kw..., objective = :margin, margin_band_hard = true, margin_reg = 0.0)
    v = validate_dispatch(net, ang, kw, x; k = KV)
    ok = x.status in RAY_OK && v.equations <= 1e-6 && v.viol["band"] <= 1e-6 && v.viol["ampacity"] <= 1e-6
    rho_exact = ok ? max(0.0, -v.min_import_kW) : NaN
    if !ok
        say(@sprintf("  restoration %s %s: exact solve %s, equations %.1e, band %.1e, limits %.1e -- no validated restoration dispatch", tag, pname, x.status, v.equations, v.viol["band"], v.viol["ampacity"]))
        csva("restoration", [tag, pname, join(string.(cuts), "+"), string(x.status), "NaN", "n/a", "n/a", "n/a", "n/a", "n/a", "n/a", "no validated exact restoration dispatch"])
        return nothing
    end
    rub = rho_exact + 1.0                                       # kW: the domain rho <= rho_exact + 1 contains the optimum
    c = certify(net, ang, kw; k = KV, cuts, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, tighten = true, detail = 8)
    S = net.S_base_kVA
    lbt = isfinite(c.db.bound) ? c.db.bound * S : -Inf
    lbl = has_duals(c.R.model) ? dual_bound(c.R.model; box = c.box_loose).bound * S : -Inf
    cls = lbt > 0 ? "no-backflow IMPOSSIBLE under the limits (validated lower bound on rho > 0)" : "not certified (rho lower bound <= 0)"
    csva("restoration", [tag, pname, join(string.(cuts), "+"), string(x.status), @sprintf("%.4f", rho_exact), string(c.status), @sprintf("%.4f", -c.primal),
                         @sprintf("%.4f", lbt), @sprintf("%.4f", lbl), @sprintf("%.3e", c.db.box_term * S), @sprintf("%.1e", c.db.rounding * S), cls])
    say(@sprintf("  restoration %-26s %-7s %-12s exact rho %.4f kW (validated) | relaxation %s rho %.4f kW; validated lower bound: tightened %.4f kW, loose %.3e kW (box term %.3e kW) -> %s",
                 tag, pname, join(string.(cuts), "+"), rho_exact, c.status, -c.primal, lbt, lbl, c.db.box_term * S, cls))
    if lbt <= 0
        for (k, (val, nm, a, l, u)) in enumerate(c.db.top[1:min(4, end)])
            say(@sprintf("      box term #%d: %-34s coefficient %+.3e, box [%.5g, %.5g], contributes %.3e kW", k, nm, a, l, u, val * S))
        end
    end
    return (; rho_exact, lbt, c)
end

function main()
    conv = conventions()
    say("conventions consistent: ", conv.ok)
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    env0 = read_env(joinpath(OUT_DIR, "attainable_currents_lines.csv"))
    csvw("cases", "case,point,r_deg,theta_deg,additions,L_kW,nlp_status,orig_equations,limit_violation,relax_status,dual_status,primal_kW,U_tight_kW,U_loose_kW,fictitious_kW,vsq_gap,primal_residual,farkas_value,farkas_boxterm_loose,farkas_bound_loose,farkas_boxterm_tight,farkas_bound_tight,dual_psd_violation,lifted_violation,lifted_box_violation,lifted_box_at,build_s,solve_s,class")
    csvw("contributors", "case,point,additions,boxes,rank,variable,residual_coefficient,lower,upper,contribution")
    csvw("restoration", "case,point,additions,exact_status,rho_exact_kW,relax_status,rho_relax_kW,rho_lb_tight_kW,rho_lb_loose_kW,boxterm_kW,rounding_kW,class")
    say("\n== (1-3) frozen cases A and B, hour 12, 292.5 deg")
    cfg = [("A: src+lines 1.2x", 1.2, 1.2, Symbol[]), ("B: src+lines 1.5x", 1.5, 1.5, [:kcl0, :jbound])]
    certified = false
    for (tag, cs, cl_, cuts) in cfg
        kwl = (; kw1..., limits(env0, cs, cl_)...)
        for (pname, r) in (("inside", 0.98RSTAR), ("outside", 1.02RSTAR))
            o = run_case(tag, kwl, pname, r, TH0, cuts)
            certified |= occursin("INFEASIBLE under", o.cls)
            pname == "outside" && (z = restoration(tag, kwl, pname, r, TH0, cuts); z !== nothing && z.lbt > 0 && (certified = true))
            pname == "inside" && restoration(tag, kwl, pname, r, TH0, cuts)
        end
    end
    say("\n== (4) challenge the envelope: no current limit, cross-seeded maximizations")
    csvw("envelope_runs", "target,start,status,I_pu,I_A,d3_deg,d4_deg,equations,band_violation,min_import_kW,max_line,max_line_A,frozen_limit_1.2x_A,exceeds_frozen_envelope,exceeds_1.2x")
    Ib = net.S_base_kVA / net.kV_base
    lines = [l.name for l in net.lines]; srcs = [s.name for s in net.sources]
    newenv = Dict{String,Tuple{Float64,String}}()
    function rec!(x, tag)
        for l in lines; v = abs(complex(x.raw.lr[l][1], x.raw.li[l][1])); v > get(newenv, l, (0.0, ""))[1] && (newenv[l] = (v, tag)); end
        for n in srcs; v = abs(complex(x.raw.ir[n][1], x.raw.ii[n][1])); v > get(newenv, "source:" * n, (0.0, ""))[1] && (newenv["source:" * n] = (v, tag)); end
    end
    starts = Any[("equal (cold)", nothing, at(0.0, 0.0))]
    xi = solve_multiperiod(net, at(0.98RSTAR, TH0); kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0); push!(starts, ("hour inside 292.5", xi, at(0.98RSTAR, TH0))); rec!(xi, "hour inside 292.5")
    rng = MersenneTwister(7)
    for k in 1:10                                               # random feasible angle vectors
        a = Dict("subs1" => 0.0, "subs3" => 0.6 * rand(rng) - 0.3, "subs4" => 0.6 * rand(rng) - 0.3)
        x = solve_multiperiod(net, a; kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        v = validate_dispatch(net, a, kw1, x; k = KV)
        (x.status in RAY_OK && v.equations <= 1e-6 && v.margin_kW >= 0) && (push!(starts, ("random $k", x, a)); rec!(x, "random $k"))
    end
    say("  starts: ", length(starts), " (equal, the inside dispatch, ", length(starts) - 2, " random feasible angle vectors)")
    targets = [[(:source, n) for n in srcs]; [(:line, l) for l in ("l1", "l2", "l119", "l54", "l58", "l5", "l9", "l118", "l123", "l115", "l116", "l110", "l46", "l111")]]
    maxers = Any[]
    function sweep(starts_)
        for (kind, nm) in targets, (sname, sx, sang) in starts_
            x = solve_multiperiod(net, sang; kw1..., objective = :maxcurrent, current_target = (kind, nm), free_angles = Dict("subs3" => DOMAIN, "subs4" => DOMAIN), start = sx)
            x.status in RAY_OK || continue
            v = validate_dispatch(net, x.angles_deg, kw1, x; k = KV)
            (v.equations <= 1e-6 && v.viol["band"] <= 1e-6 && v.viol["backflow"] <= 1e-6) || continue
            rec!(x, "max $(kind) $(nm) from $(sname)")
            push!(maxers, ("max $(nm) from $(sname)", x, x.angles_deg))
            I = kind === :source ? abs(complex(x.raw.ir[nm][1], x.raw.ii[nm][1])) : abs(complex(x.raw.lr[nm][1], x.raw.li[nm][1]))
            key = kind === :source ? "source:" * nm : nm
            ml = argmax(l -> abs(complex(x.raw.lr[l][1], x.raw.li[l][1])), lines)
            csva("envelope_runs", [string(kind, ":", nm), cl(sname), string(x.status), @sprintf("%.5f", I), @sprintf("%.2f", I * Ib), @sprintf("%.5f", x.angles_deg["subs3"]), @sprintf("%.5f", x.angles_deg["subs4"]),
                                   @sprintf("%.1e", v.equations), @sprintf("%.1e", v.viol["band"]), @sprintf("%.3f", v.min_import_kW), ml, @sprintf("%.2f", abs(complex(x.raw.lr[ml][1], x.raw.li[ml][1])) * Ib),
                                   @sprintf("%.2f", max(1.2 * get(env0, key, 0.0), FLOOR) * Ib), string(I > get(env0, key, 0.0) + 1e-6), string(I > max(1.2 * get(env0, key, 0.0), FLOOR) + 1e-6)])
        end
    end
    sweep(starts); n1 = length(maxers)
    sweep(copy(maxers))                                         # cross-seeding: every maximizer as a start for every target
    say("  validated maximizations: ", n1, " from the first starts, ", length(maxers) - n1, " cross-seeded")
    exc = [(el, v, w) for (el, (v, w)) in newenv if v > get(env0, el, 0.0) + 1e-6]
    exc12 = [(el, v, w) for (el, v, w) in exc if v > max(1.2 * get(env0, el, 0.0), FLOOR) + 1e-6]
    say("  elements above the frozen envelope: ", length(exc), "; above 1.2x of it: ", length(exc12))
    for (el, v, w) in sort(exc, by = z -> -(z[2] / max(get(env0, z[1], 0.0), 1e-9)))[1:min(8, end)]
        say(@sprintf("    %-14s %.2f A (frozen %.2f A, x%.3f) at %s", el, v * Ib, get(env0, el, 0.0) * Ib, v / max(get(env0, el, 0.0), 1e-9), w))
    end
    open(PFX * "_envelope.csv", "w") do io
        println(io, "element,frozen_I_A,challenge_I_A,updated_I_A,updated_I_pu,where_challenge")
        for el in sort(unique([collect(keys(env0)); collect(keys(newenv))]))
            f0 = get(env0, el, 0.0); n0 = get(newenv, el, (0.0, ""))
            println(io, @sprintf("%s,%.3f,%.3f,%.3f,%.6f,%s", el, f0 * Ib, n0[1] * Ib, max(f0, n0[1]) * Ib, max(f0, n0[1]), cl(n0[2])))
        end
    end
    env = Dict(el => max(get(env0, el, 0.0), get(newenv, el, (0.0, ""))[1]) for el in unique([collect(keys(env0)); collect(keys(newenv))]))
    if !isempty(exc12)
        say("  envelope updated transparently; rerunning A and B with it")
        for (tag, cs, cl_, cuts) in cfg
            kwl = (; kw1..., limits(env, cs, cl_)...)
            for (pname, r) in (("inside", 0.98RSTAR), ("outside", 1.02RSTAR))
                o = run_case(tag * " (updated envelope)", kwl, pname, r, TH0, cuts)
                pname == "outside" && restoration(tag * " (updated envelope)", kwl, pname, r, TH0, cuts)
            end
        end
    end
    say("\n== (5) cushion sensitivity (envelope as used by the rest of the run)")
    for cs in (1.2, 1.5, 2.0), (kindl, cl_) in (("sources", nothing), ("sources+lines", cs))
        kwl = (; kw1..., limits(env, cs, cl_)...)
        variants = kindl == "sources" ? (Symbol[], [:kcl0, :jbound]) : (Symbol[],)
        for cuts in variants, (pname, r) in (("inside", 0.98RSTAR), ("outside", 1.02RSTAR))
            tag = @sprintf("%s %.1fx", kindl, cs)
            o = run_case(tag, kwl, pname, r, TH0, cuts; detail = 0)
            certified |= occursin("INFEASIBLE under", o.cls)
            pname == "outside" && (z = restoration(tag, kwl, pname, r, TH0, cuts); z !== nothing && z.lbt > 0 && (certified = true))
        end
    end
    say("\nany validated certificate: ", certified)
    if certified
        say("\n== (6) 157.5 deg, same envelope and cushions")
        marg(r) = (x = solve_multiperiod(net, at(r, 157.5); kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0); validate_dispatch(net, at(r, 157.5), kw1, x; k = KV).margin_kW)
        r0, r1 = 0.09, 0.11; m0, m1 = marg(r0), marg(r1)
        for _ in 1:10; r2 = r1 - m1 * (r1 - r0) / (m1 - m0); r0, m0 = r1, m1; r1 = r2; m1 = marg(r1); abs(m1) < 1e-4 && break; end
        say(@sprintf("  hour 12 boundary along 157.5 deg: r* = %.6f deg (margin %+.1e kW)", r1, m1))
        for (tag, cs, cl_, cuts) in cfg
            kwl = (; kw1..., limits(env, cs, cl_)...)
            for (pname, r) in (("inside", 0.98r1), ("outside", 1.02r1))
                run_case(tag * " 157.5", kwl, pname, r, 157.5, cuts)
                pname == "outside" && restoration(tag * " 157.5", kwl, pname, r, 157.5, cuts)
            end
        end
    end
end
"Extra cell (run after main, appends to its CSVs): 2x sources and lines with both nodal additions, updated envelope."
function extra()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    env = Dict(String(r[1]) => parse(Float64, r[5]) for r in (split(l, ',') for l in readlines(PFX * "_envelope.csv")[2:end]))
    say("== (5, extra) sources+lines 2.0x with both nodal additions (updated envelope)")
    kwl = (; kw1..., limits(env, 2.0, 2.0)...)
    for (pname, r) in (("inside", 0.98RSTAR), ("outside", 1.02RSTAR))
        run_case("sources+lines 2.0x", kwl, pname, r, TH0, [:kcl0, :jbound]; detail = 6)
        pname == "outside" && restoration("sources+lines 2.0x", kwl, pname, r, TH0, [:kcl0, :jbound])
    end
end
"Extra 2: the current caps that carry the validated Farkas rays (A, B outside), and a challenge to those lines."
function extra2()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    env = Dict(String(r[1]) => parse(Float64, r[5]) for r in (split(l, ',') for l in readlines(PFX * "_envelope.csv")[2:end]))
    Ib = net.S_base_kVA / net.kV_base
    infl = Dict{String,Float64}()
    for (tag, cs, cl_, cuts) in (("A", 1.2, 1.2, Symbol[]), ("B", 1.5, 1.5, [:kcl0, :jbound]))
        kwl = (; kw1..., limits(env, cs, cl_)...)
        c = certify(net, at(1.02RSTAR, TH0), kwl; k = KV, cuts, m_lo_kW = 0.0, m_hi_kW = 2000.0, tighten = true)
        w = Tuple{Float64,String}[]
        for cref in all_constraints(c.R.model, AffExpr, MOI.LessThan{Float64})
            fn = constraint_object(cref).func
            nm = [name(v) for v in keys(fn.terms)]
            any(startswith(x, "lrlr[") for x in nm) || continue                   # a line cap: lrlr + lili <= I_max^2
            ln = split(split(first(x for x in nm if startswith(x, "lrlr[")), '[')[2], ',')[1]
            push!(w, (abs(dual(cref)) * constraint_object(cref).set.upper, String(ln)))
        end
        sort!(w; rev = true); tot = sum(first.(w))
        say(@sprintf("  %s outside, %s: line caps carrying the ray (share of the caps' total weight): ", tag, c.status), join([@sprintf("%s %.0f%%", l, 100v / tot) for (v, l) in w[1:8]], ", "))
        for (v, l) in w[1:8]; infl[l] = max(get(infl, l, 0.0), v / tot); end
    end
    done = Set(["l1", "l2", "l119", "l54", "l58", "l5", "l9", "l118", "l123", "l115", "l116", "l110", "l46", "l111"])
    newt = [l for l in keys(infl) if !(l in done)]
    say("  influential caps not challenged before: ", isempty(newt) ? "none" : join(newt, ", "))
    isempty(newt) && return
    starts = Any[("equal (cold)", nothing, at(0.0, 0.0))]
    xi = solve_multiperiod(net, at(0.98RSTAR, TH0); kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0); push!(starts, ("hour inside 292.5", xi, at(0.98RSTAR, TH0)))
    for n in ("subs1", "subs3", "subs4"), (sn, sx, sa) in copy(starts)[1:2]
        x = solve_multiperiod(net, sa; kw1..., objective = :maxcurrent, current_target = (:source, n), free_angles = Dict("subs3" => DOMAIN, "subs4" => DOMAIN), start = sx)
        x.status in RAY_OK && push!(starts, ("max $(n) from $(sn)", x, x.angles_deg))
    end
    for l in newt, (sn, sx, sa) in starts
        x = solve_multiperiod(net, sa; kw1..., objective = :maxcurrent, current_target = (:line, l), free_angles = Dict("subs3" => DOMAIN, "subs4" => DOMAIN), start = sx)
        x.status in RAY_OK || continue
        v = validate_dispatch(net, x.angles_deg, kw1, x; k = KV)
        (v.equations <= 1e-6 && v.viol["band"] <= 1e-6 && v.viol["backflow"] <= 1e-6) || continue
        I = abs(complex(x.raw.lr[l][1], x.raw.li[l][1]))
        say(@sprintf("    max %-6s from %-28s %.2f A (envelope %.2f A, 1.2x %.2f A)%s", l, sn, I * Ib, get(env, l, 0.0) * Ib, max(1.2get(env, l, 0.0), FLOOR) * Ib,
                     I > max(1.2get(env, l, 0.0), FLOOR) + 1e-6 ? "  EXCEEDS 1.2x" : (I > get(env, l, 0.0) + 1e-6 ? "  above the envelope" : "")))
    end
end
"Extra 3: fold the second challenge's maxima into the envelope (v2) and rerun A and B."
function extra3()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    rows = [split(l, ',') for l in readlines(PFX * "_envelope.csv")[2:end]]
    env = Dict(String(r[1]) => parse(Float64, r[5]) for r in rows)
    Ib = net.S_base_kVA / net.kV_base
    starts = Any[("equal (cold)", nothing, at(0.0, 0.0))]
    xi = solve_multiperiod(net, at(0.98RSTAR, TH0); kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0); push!(starts, ("hour inside 292.5", xi, at(0.98RSTAR, TH0)))
    for n in ("subs1", "subs3", "subs4"), (sn, sx, sa) in copy(starts)[1:2]
        x = solve_multiperiod(net, sa; kw1..., objective = :maxcurrent, current_target = (:source, n), free_angles = Dict("subs3" => DOMAIN, "subs4" => DOMAIN), start = sx)
        x.status in RAY_OK && push!(starts, ("max $(n) from $(sn)", x, x.angles_deg))
    end
    upd = Dict{String,Tuple{Float64,String}}()
    for l in ("l107", "l100", "l71", "l61", "l120", "l76", "l125", "l67", "l12", "l114", "l112"), (sn, sx, sa) in starts
        x = solve_multiperiod(net, sa; kw1..., objective = :maxcurrent, current_target = (:line, l), free_angles = Dict("subs3" => DOMAIN, "subs4" => DOMAIN), start = sx)
        x.status in RAY_OK || continue
        v = validate_dispatch(net, x.angles_deg, kw1, x; k = KV)
        (v.equations <= 1e-6 && v.viol["band"] <= 1e-6 && v.viol["backflow"] <= 1e-6) || continue
        for el in [[ln.name for ln in net.lines]; ["source:" * s.name for s in net.sources]]
            I = startswith(el, "source:") ? abs(complex(x.raw.ir[el[8:end]][1], x.raw.ii[el[8:end]][1])) : abs(complex(x.raw.lr[el][1], x.raw.li[el][1]))
            I > max(get(env, el, 0.0), get(upd, el, (0.0, ""))[1]) + 1e-9 && (upd[el] = (I, "max $(l) from $(sn)"))
        end
    end
    open(PFX * "_envelope_v2.csv", "w") do io
        println(io, "element,v1_I_A,v2_I_A,v2_I_pu,raised_by")
        for el in sort(collect(keys(env)))
            v2 = max(env[el], get(upd, el, (0.0, ""))[1])
            println(io, @sprintf("%s,%.3f,%.3f,%.6f,%s", el, env[el] * Ib, v2 * Ib, v2, cl(get(upd, el, (0.0, ""))[2])))
        end
    end
    say("== envelope v2: ", length(upd), " elements raised: ", join([@sprintf("%s %.2f->%.2f A", el, env[el] * Ib, v * Ib) for (el, (v, _)) in sort(collect(upd), by = p -> -p[2][1] / max(env[p[1]], 1e-9))][1:min(10, end)], ", "))
    env2 = Dict(el => max(env[el], get(upd, el, (0.0, ""))[1]) for el in keys(env))
    for (tag, cs, cl_, cuts) in (("A: src+lines 1.2x (envelope v2)", 1.2, 1.2, Symbol[]), ("B: src+lines 1.5x (envelope v2)", 1.5, 1.5, [:kcl0, :jbound]))
        kwl = (; kw1..., limits(env2, cs, cl_)...)
        for (pname, r) in (("inside", 0.98RSTAR), ("outside", 1.02RSTAR))
            run_case(tag, kwl, pname, r, TH0, cuts; detail = 4)
            pname == "outside" && restoration(tag, kwl, pname, r, TH0, cuts)
        end
    end
end
"Extra 4: 157.5 deg with envelope v2 (the boundary r* found by the main run's secant)."
function extra4()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    env2 = Dict(String(r[1]) => parse(Float64, r[4]) for r in (split(l, ',') for l in readlines(PFX * "_envelope_v2.csv")[2:end]))
    r1 = 0.100075
    for (tag, cs, cl_, cuts) in (("A: src+lines 1.2x (envelope v2) 157.5", 1.2, 1.2, Symbol[]), ("B: src+lines 1.5x (envelope v2) 157.5", 1.5, 1.5, [:kcl0, :jbound]))
        kwl = (; kw1..., limits(env2, cs, cl_)...)
        for (pname, r) in (("inside", 0.98r1), ("outside", 1.02r1))
            run_case(tag, kwl, pname, r, 157.5, cuts; detail = 4)
            pname == "outside" && restoration(tag, kwl, pname, r, 157.5, cuts)
        end
    end
end
get(ENV, "CV_EXTRA", "0") == "lib" ? nothing : get(ENV, "CV_EXTRA", "0") == "4" ? extra4() : get(ENV, "CV_EXTRA", "0") == "3" ? extra3() : get(ENV, "CV_EXTRA", "0") == "2" ? extra2() : get(ENV, "CV_EXTRA", "0") == "1" ? extra() : main()
