# Certificates on a deck with service transformers (2026-10-26): the hour-12-style pilot of
# certificate_validation.jl (IEEE 123-node) moved to ieee9500xf_1ph, three substations, five ties closed,
# the 0.95-1.05 band at customer buses only, batteries idle, PV reactive power dispatched.
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/certificate_secondary_pilot.jl
# Steps:
#  1. The limiting hour t* of the T = 24 idle edge along the pilot direction (CERT_THETA, default 135 deg):
#     the hour with the smallest import at 98% of the campaign's measured reach (regions_oct25).
#  2. At t* alone (T = 1), the hour's own boundary r*_t along each of the 16 directions (bisection on the
#     sign of the exact margin), with no current limit.
#  3. The current envelope: each line's and source's largest current over a frozen set of validated
#     dispatches, those at equal angles and at 98% of r*_t in all 16 directions. The limits are the
#     envelope itself (1.0x, the tightest that keeps every one of those dispatches feasible), with a floor
#     of CERT_FLOOR pu (default 0.01) for lines that carry almost nothing. These are empirical limits, not
#     ratings: a certificate proves the model WITH these limits infeasible, nothing about the model without them.
#  4. Inside (98%) and outside (102%) points along the pilot directions: exact margin L (Ipopt, re-validated)
#     and the relaxation's bound U (dual_bound / Farkas), classified feasible / infeasible / unresolved.
# Writes certificate_secondary_{hours,envelope,cases}.csv and a log in results/multiperiod.
include(joinpath(@__DIR__, "..", "multiperiod_full_angle.jl"))
include(joinpath(@__DIR__, "..", "feasibility_certificate.jl"))
using Logging

const OUT = normpath(get(ENV, "MULTIPERIOD_OUT_DIR", joinpath(@__DIR__, "..", "results", "multiperiod")))
mkpath(OUT)
const PFX = joinpath(OUT, "certificate_secondary")
say(args...) = (println(args...); flush(stdout))
const RAY_OK = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
const KV = 20000.0                                   # kW of margin per pu of |V|^2, as MULTI_VMARGIN
const UTOL = 0.01                                    # kW: U must be below -UTOL to call a point infeasible
const SYSTEM = "ieee9500xf_1ph"
const TIES = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]
const S0, A1, A2 = "subs1", "subs2", "subs3"
const THETAS = [22.5 * i for i in 0:15]
const PILOT = parse.(Float64, split(get(ENV, "CERT_THETA", "135,315"), ','))
const FLOOR = parse(Float64, get(ENV, "CERT_FLOOR", "0.01"))
# infeasible: a margin bound U < -UTOL, or a validated Farkas bound > 0 (both from dual_bound, not from the solver status)
classify(L, U, fk) = isfinite(L) && L >= 0 ? "feasible" : ((isfinite(U) && U < -UTOL) || (isfinite(fk) && fk > 0) ? "infeasible" : "unresolved")
at(r, th) = Dict(S0 => 0.0, A1 => r * cosd(th), A2 => r * sind(th))

# ---- deck, profiles, the campaign's measured idle edges -----------------------------------------------
kw0 = (; sources = [S0, A1, A2], snapshot = true, load_band = (0.5, 1.5), extra = ["Line.$t.enabled=true" for t in TIES])
pvs = with_logger(NullLogger()) do; read_pv(SYSTEM; kw0...); end
compile_deck(SYSTEM; disable = ["PVSystem"], kw0...)
net = read_network()
cust = Set(l.to for l in net.lines if startswith(lowercase(l.name), "svc_"))
lam, price, pvf = load_profile(24), price_profile(24), pv_profile(24)
reach = Dict{Float64,Float64}()
for (r, _, fs) in walkdir(joinpath(@__DIR__, "..", "results", "kamiak", "regions_oct25")), f in fs
    f == "three_large_windows.csv" || continue
    for l in readlines(joinpath(r, f))[2:end]
        c = split(l, ','); (c[1] == SYSTEM && c[5] == "nobatt" && c[11] == "false") && (reach[parse(Float64, c[6])] = parse(Float64, c[8]))
    end
end
say("measured idle reaches (T = 24): ", join([@sprintf("%.1f: %.3f", k, reach[k]) for k in sort(collect(keys(reach)))], ", "))
tol = (; tol = 1e-8, constr_viol_tol = 1e-8)
kw24 = (; pvs, pv_q = :optimize, C_B = cb_by_rule(net), ideal_sources = false, battery = :idle, band_buses = cust, tol...,
        T = 24, load_factor = lam, price, pv_factor = pvf, gamma = gamma_by_rule(net), terminal = :soft)
margin(ang, kw; start = nothing) = solve_multiperiod(net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0, start)

# ---- 1. the limiting hour -----------------------------------------------------------------------------
# CERT_REUSE=1: take t*, the hour boundaries and the envelope from the files of an earlier full run, and only redo
# step 4 (with CERT_FLOOR, CERT_CUTS=kcl0,jbound,... and CERT_TAG naming its cases file).
const REUSE = get(ENV, "CERT_REUSE", "0") == "1"
const CUTS = Symbol.(split(get(ENV, "CERT_CUTS", ""), ','; keepempty = false))
const TAG = get(ENV, "CERT_TAG", "")
if !REUSE
th0 = PILOT[1]
x24 = margin(at(0.98reach[th0], th0), kw24)
x24.status in RAY_OK || error("T = 24 solve at 98% of the measured reach failed: $(x24.status)")
hourmin = [minimum(x24.P_subs_kW[n][t] for n in keys(x24.P_subs_kW)) for t in 1:24]
tstar = argmin(hourmin)
say(@sprintf("limiting hour along %.1f deg: t* = %d (smallest import %.1f kW at 98%% of the measured reach)", th0, tstar, hourmin[tstar]))
open(PFX * "_hours.csv", "w") do io
    println(io, "hour,min_import_kW"); for t in 1:24; @printf(io, "%d,%.4f\n", t, hourmin[t]); end
end
else
    hm = [parse(Float64, split(l, ',')[2]) for l in readlines(PFX * "_hours.csv")[2:end]]; tstar = argmin(hm)
    say("reusing t* = ", tstar, ", the hour boundaries and the envelope of the earlier run")
end
kw1 = merge(kw24, (; T = 1, dt_h = 1.0, load_factor = [lam[tstar]], price = [price[tstar]], pv_factor = [pvf[tstar]], gamma = 0.0, terminal = :free))

# ---- 2. the hour's own boundary along every direction (no current limit) ------------------------------
"""
Validated margin of an exact solve, and the solution. The margin problem is always feasible as an optimization, so a
failed solve is numerical: it is retried cold and from `alt` before the point counts as failed (NaN, logged).
"""
const NFAIL = Ref(0)
function vmargin(ang, kw; start = nothing, alt = nothing)
    x = nothing
    for st in unique([start, nothing, alt])
        x = margin(ang, kw; start = st)
        x.status in RAY_OK || continue
        v = validate_dispatch(net, ang, kw, x; k = KV)
        v.equations <= 1e-6 && return v.margin_kW, x
    end
    NFAIL[] += 1; say("      solve failed from every start (", x.status, ")")
    return NaN, x
end
rstar = Dict{Float64,Float64}(); inside = Dict{Float64,Any}()
env = Dict{String,Float64}(); senv = Dict{String,Float64}(); who = Dict{String,String}()
if !REUSE
for th in THETAS
    # batteries idle: every hour is feasible wherever the day is, so the day's edge is a lower bound here (checked)
    lo, hi = 0.0, 1.2reach[th]; xs = nothing
    m, x = vmargin(at(0.999reach[th], th), kw1)
    isfinite(m) && m >= 0 ? (lo = 0.999reach[th]; xs = x) : say(@sprintf("  %5.1f deg: the day's edge does not validate at this hour; bracketing from zero", th))
    while true                                                       # an infeasible point beyond the edge
        m, _ = vmargin(at(hi, th), kw1); (isfinite(m) && m < 0) && break
        lo = hi; hi *= 1.5; hi > 40 && error("no edge within 40 deg along $th")
    end
    for _ in 1:16                                                    # a point failed from every start counts as outside (logged)
        mid = (lo + hi) / 2; m, x = vmargin(at(mid, th), kw1; start = xs, alt = xs)
        isfinite(m) && m >= 0 ? (lo = mid; xs = x) : (hi = mid)
    end
    rstar[th] = lo
    m, x = vmargin(at(0.98lo, th), kw1; start = xs, alt = xs)
    isfinite(m) && m >= 0 || (say(@sprintf("  %5.1f deg: 98%% point does not validate; left out of the envelope", th)); continue)
    inside[th] = x
    say(@sprintf("  %5.1f deg: hour boundary r*_t = %.5f deg (T = 24 reach %.4f); 98%% point margin %.3f kW", th, lo, reach[th], m))
end
m0, x0 = vmargin(at(0.0, 0.0), kw1)
say(@sprintf("  equal angles: margin %.3f kW", m0))

open(PFX * "_rstar.csv", "w") do io
    println(io, "theta_deg,hour_boundary_deg,day_reach_deg"); for th in THETAS; haskey(rstar, th) && @printf(io, "%.1f,%.6f,%.6f\n", th, rstar[th], reach[th]); end
end

# ---- 3. the envelope ----------------------------------------------------------------------------------
for (tag, x) in [("equal", x0); [(@sprintf("%.1f", th), inside[th]) for th in THETAS if haskey(inside, th)]]
    for (n, r) in x.raw.lr
        I = hypot(r[1], x.raw.li[n][1]); I > get(env, n, -1.0) && (env[n] = I; who[n] = tag)
    end
    for (n, r) in x.raw.ir
        I = hypot(r[1], x.raw.ii[n][1]); senv[n] = max(get(senv, n, 0.0), I)
    end
end
Ib = net.S_base_kVA / net.kV_base
open(PFX * "_envelope.csv", "w") do io
    println(io, "element,kind,envelope_pu,envelope_A,limit_pu,from_dispatch")
    for n in sort(collect(keys(env))); @printf(io, "%s,line,%.6f,%.3f,%.6f,%s\n", n, env[n], env[n] * Ib, max(env[n], FLOOR), who[n]); end
    for n in sort(collect(keys(senv))); @printf(io, "%s,source,%.6f,%.3f,%.6f,\n", n, senv[n], senv[n] * Ib, senv[n]); end
end
else
    for l in readlines(PFX * "_rstar.csv")[2:end]; c = parse.(Float64, split(l, ',')); rstar[c[1]] = c[2]; end
    for l in readlines(PFX * "_envelope.csv")[2:end]
        c = split(l, ','); c[2] == "line" ? (env[String(c[1])] = parse(Float64, c[3])) : (senv[String(c[1])] = parse(Float64, c[3]))
    end
end
Ib = net.S_base_kVA / net.kV_base
amp =Dict(n => max(v, FLOOR) for (n, v) in env)
kwL = merge(kw1, (; ampacity_pu = amp, source_ampacity_pu = senv))
say(@sprintf("envelope: %d lines (%d at the %.3f pu floor), %d sources; base %.2f A per pu", length(amp), count(<(FLOOR), values(env)), FLOOR, length(senv), Ib))

# ---- 4. inside and outside points under the envelope ---------------------------------------------------
open(PFX * "_cases" * TAG * ".csv", "w") do io
    println(io, "theta_deg,point,r_deg,L_kW,nlp_status,max_violation,U_kW,relax_status,farkas,box_term_kW,class,nlp_s,relax_build_s,relax_solve_s,opendss")
end
for th in PILOT, (name, f) in (("inside 98%", 0.98), ("outside 102%", 1.02))
    r = f * rstar[th]; ang = at(r, th)
    t0 = time(); x = margin(ang, kwL; start = get(inside, th, nothing)); tn = time() - t0
    ok = x.status in RAY_OK
    v = ok ? validate_dispatch(net, ang, kwL, x; k = KV) : nothing
    L = ok && v.equations <= 1e-6 ? v.margin_kW : NaN
    c = certify(net, ang, kwL; k = KV, cuts = CUTS)
    cls = classify(L, c.U, c.farkas.bound)
    rp = ""
    if ok && isfinite(L) && L >= 0
        rr = replay_opendss(net, x, ang; load_factor = [lam[tstar]], pvs, pv_factor = [pvf[tstar]])
        rp = @sprintf("dP %.1e kW; dV %.1e pu; converged %s", rr.dP, rr.dV, rr.converged)
    end
    fk = isfinite(c.farkas.bound) ? @sprintf("%.4g", c.farkas.bound) : ""
    open(PFX * "_cases" * TAG * ".csv", "a") do io
        println(io, join([@sprintf("%.1f", th), name, @sprintf("%.6f", r), @sprintf("%.4f", L), string(x.status), ok ? @sprintf("%.2e", v.max_violation) : "NaN",
                          @sprintf("%.4f", c.U), string(c.status), fk, @sprintf("%.3e", c.box_term_kW), cls, @sprintf("%.1f", tn),
                          @sprintf("%.1f", c.build_s), @sprintf("%.1f", c.solve_s), rp], ','))
    end
    say(@sprintf("  %5.1f deg %-12s r %.5f | L %+9.3f kW (%s) | U %+9.3f kW (%s%s) | %s | nlp %.0f s, relax %.0f+%.0f s %s",
                 th, name, r, L, x.status, c.U, c.status, isempty(fk) ? "" : ", Farkas " * fk, cls, tn, c.build_s, c.solve_s, rp))
end
say("failed solves (every start): ", NFAIL[], "; wrote ", PFX, "_{hours,envelope,cases}.csv")
