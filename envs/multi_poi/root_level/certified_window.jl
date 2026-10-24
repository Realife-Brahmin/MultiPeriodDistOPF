# Certified window versus current-limit cushion (2026-10-22).
#   julia --project=envs/msopf_sdp envs/multi_poi/root_level/certified_window.jl
# IEEE 123-node, sources 1, 3, 4, hour 12 alone, batteries idle, PV reactive power dispatched.
# Limits = cushion x envelope v2 (certificate_validation_envelope_v2.csv), no line below 0.05 pu.
# For each direction theta (delta3 = r cos theta, delta4 = r sin theta) and each limit configuration:
#   (a) the exact window edge r*_lim with the limits on, against r* without limits, and how close any
#       limit comes to binding at the edge dispatch (largest |I| / I_max);
#   (b) the certified edge r_cert: the smallest r found at which the relaxation proves no-backflow
#       impossible (restoration problem: band and limits hard, validated lower bound on rho > 0, or a
#       validated Farkas ray). Grid 1.02 .. 3 x r*, then bisection. Points beyond r_cert are proven
#       infeasible under the limits; the band [r*_lim, r_cert) is the relaxation gap.
# Writes certified_window.csv and certified_window.log.
ENV["CV_EXTRA"] = "lib"
include(joinpath(@__DIR__, "certificate_validation.jl"))

const CW = joinpath(OUT_DIR, "certified_window")
const DIRS = [22.5 + 45k for k in 0:7]
const GRID = [1.02, 1.05, 1.1, 1.25, 1.5, 2.0, 3.0]

"Exact margin at (r, theta) under kw: validated L in kW (NaN if not validated) and the dispatch."
function exact_margin(kw, r, th)
    ang = at(r, th)
    x = solve_multiperiod(net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
    v = validate_dispatch(net, ang, kw, x; k = KV)
    return (x.status in RAY_OK && v.equations <= 1e-6 && v.viol["ampacity"] <= 1e-6 ? v.margin_kW : NaN), x
end

"Secant on the exact margin along theta."
function edge(kw, th; r0 = 0.05, r1 = 0.15)
    m0, _ = exact_margin(kw, r0, th); m1, x1 = exact_margin(kw, r1, th)
    for _ in 1:20
        (isfinite(m0) && isfinite(m1) && m1 != m0) || return (NaN, nothing)
        r2 = max(r1 - m1 * (r1 - r0) / (m1 - m0), 1e-4)
        r0, m0 = r1, m1; r1 = r2; m1, x1 = exact_margin(kw, r1, th)
        abs(m1) < 1e-5 && return (r1, x1)
    end
    return (NaN, nothing)
end

"Largest |I| / I_max over the limited elements at dispatch x."
function utilization(kw, x)
    best, at_ = 0.0, ""
    sa, la = get(kw, :source_ampacity_pu, nothing), get(kw, :ampacity_pu, nothing)
    if sa !== nothing
        for (n, a) in sa
            u = abs(complex(x.raw.ir[n][1], x.raw.ii[n][1])) / a; u > best && ((best, at_) = (u, "source:" * n))
        end
    end
    if la !== nothing
        for (n, a) in la
            u = abs(complex(x.raw.lr[n][1], x.raw.li[n][1])) / a; u > best && ((best, at_) = (u, n))
        end
    end
    return best, at_
end

"Is (r, theta) certified infeasible under kw (relaxation with cuts)? Returns (certified, how, lower bound on rho in kW)."
function certified_at(kw, r, th, cuts)
    ang = at(r, th); S = net.S_base_kVA
    x = solve_multiperiod(net, ang; kw..., objective = :margin, margin_band_hard = true, margin_reg = 0.0)
    v = validate_dispatch(net, ang, kw, x; k = KV)
    ok = x.status in RAY_OK && v.equations <= 1e-6 && v.viol["band"] <= 1e-6 && v.viol["ampacity"] <= 1e-6
    rub = ok ? max(0.0, -v.min_import_kW) + 1.0 : 2000.0       # any rho cap >= 0 keeps the test valid
    c = certify(net, ang, kw; k = KV, cuts, m_lo_kW = -rub, m_hi_kW = 0.0, band_in_margin = false, tighten = true)
    infeas = c.status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
    infeas && return (c.farkas.bound > 0, "Farkas", NaN)
    lb = isfinite(c.db.bound) ? c.db.bound * S : -Inf
    return (lb > 0, "rho bound", lb)
end

function cw_main()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    env2 = Dict(String(r[1]) => parse(Float64, r[4]) for r in (split(l, ',') for l in readlines(PFX * "_envelope_v2.csv")[2:end]))
    cfgs = Any[]
    for (cs, lines_too) in ((1.2, true), (1.5, true), (2.0, true), (3.0, true), (1.5, false)), cuts in (Symbol[], [:kcl0, :jbound])
        push!(cfgs, (@sprintf("%s %.1fx", lines_too ? "sources+lines" : "sources only", cs), cs, lines_too, cuts))
    end
    push!(cfgs, ("no limits", NaN, false, Symbol[]), ("no limits", NaN, false, [:kcl0, :jbound]))
    open(io -> println(io, "theta_deg,limits,cushion,additions,r_star_nolim_deg,r_star_lim_deg,edge_shift_pct,max_util_at_edge,util_element,r_cert_deg,r_cert_over_r_star,how,rho_lb_kW_at_r_cert,note"), CW * ".csv", "w")
    for th in DIRS
        rn, _ = edge(kw1, th)
        say(@sprintf("\n== theta %.1f deg: exact edge without limits r* = %.6f deg", th, rn))
        isfinite(rn) || continue
        for (lab, cs, lines_too, cuts) in cfgs
            kwl = isnan(cs) ? kw1 : (; kw1..., limits(env2, cs, lines_too ? cs : nothing)...)
            rl, xl = isnan(cs) ? (rn, nothing) : edge(kwl, th; r0 = 0.98rn, r1 = 1.02rn)
            u, ue = xl === nothing ? (NaN, "") : utilization(kwl, xl)
            rc, how, lb, note = NaN, "", NaN, ""
            prev = 1.0
            for g in GRID
                ok, h, l = certified_at(kwl, g * rn, th, cuts)
                if ok
                    lo, hi = prev, g; how, lb = h, l
                    for _ in 1:6
                        mid = (lo + hi) / 2; okm, hm, lm = certified_at(kwl, mid * rn, th, cuts)
                        okm ? ((hi, how, lb) = (mid, hm, lm)) : (lo = mid)
                    end
                    rc = hi * rn; note = @sprintf("first certified grid point %.2fx; bisection to %.1f%%", g, 100 * (hi - lo))
                    break
                end
                prev = g
            end
            isnan(rc) && (note = "not certified up to 3 x r*")
            row = [@sprintf("%.1f", th), lab, isnan(cs) ? "none" : @sprintf("%.1f", cs), join(string.(cuts), "+"), @sprintf("%.6f", rn), @sprintf("%.6f", rl),
                   @sprintf("%.3f", 100 * (rl / rn - 1)), @sprintf("%.3f", u), ue, @sprintf("%.6f", rc), @sprintf("%.3f", rc / rn), how, @sprintf("%.4f", lb), note]
            open(io -> println(io, join(row, ',')), CW * ".csv", "a")
            say(@sprintf("  %-20s %-12s exact edge %.6f (%+.3f%%; max |I|/Imax %.2f at %s) | certified from %s", lab, join(string.(cuts), "+"), rl, 100 * (rl / rn - 1), u, ue,
                         isnan(rc) ? "-- (none up to 3 r*)" : @sprintf("%.4f deg = %.3f r* (%s)", rc, rc / rn, how)))
        end
    end
end
cw_main()
