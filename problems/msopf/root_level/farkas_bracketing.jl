# Farkas bracketing (2026-10-26): the certified feasibility map of a system, batteries idle, T = 24 day.
#   FB_SYSTEM=<ieee123|med3980|sfp1u>_<two|three|four> FB_PHASE=<A|B|C> FB_DIR=<i> \
#     julia --project=problems/msopf/environments/sdp problems/msopf/root_level/farkas_bracketing.jl
# Along each measured direction of the region campaigns (16 for three sources, the fib78 set for four), the
# boundary is bracketed from both sides: inside by a validated AC dispatch, outside by the relaxation proving the
# point infeasible (a margin bound U < -UTOL or a validated Farkas bound > 0; feasibility_certificate.jl).
#  A (per direction): the limiting hour t* (the hour with the smallest import of a T = 24 solve at 98% of the
#    measured reach), that hour's own boundary r*_t by bisection on the sign of the exact, re-validated margin
#    (starting from the day's reach, a lower bound with idle batteries), and the currents of the validated
#    dispatches at 98% of r*_t and at equal angles. Writes fb_A_<i>.csv (summary) and fb_Acur_<i>.csv (currents).
#  B (once): the current envelope, each line's and source's largest current over all phase-A dispatches. The
#    limits are this envelope itself (1.0x), with a floor of FB_FLOOR pu (1e-4) for lines that carry almost nothing.
#    Writes fb_envelope.csv. These are empirical limits, not ratings: a certificate proves the model WITH them
#    infeasible, nothing about the model without them.
#  C (per direction): at t*, under the envelope, the certified edge: GRID x r*_t until one point certifies, then
#    five bisections. Writes fb_C_<i>.csv. An hour infeasible alone makes the idle-battery day infeasible, so the
#    certified edge holds for the day.
# Output folder: MULTIPERIOD_OUT_DIR (all phases of one system share it).
include(joinpath(@__DIR__, "..", "multiperiod_full_angle.jl"))
include(joinpath(@__DIR__, "..", "feasibility_certificate.jl"))
using Logging

const SYS = get(ENV, "FB_SYSTEM", "med3980_four")
const PHASE = get(ENV, "FB_PHASE", "A")
const OUT = normpath(get(ENV, "MULTIPERIOD_OUT_DIR", joinpath(@__DIR__, "..", "results", "multiperiod", "farkas_" * SYS)))
mkpath(OUT)
say(args...) = (println(args...); flush(stdout))
const RES = joinpath(@__DIR__, "..", "results", "multiperiod")
const RAY_OK = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
const KV = 20000.0                                   # kW of margin per pu of |V|^2, as MULTI_VMARGIN
const UTOL = 0.01                                    # kW: U must be below -UTOL to call a point infeasible
const FLOOR = parse(Float64, get(ENV, "FB_FLOOR", "1e-4"))
const CUTS = Symbol.(split(get(ENV, "FB_CUTS", "kcl0,jbound"), ','; keepempty = false))
# Largest accepted equation residual (pu). 1e-6 rejected med3980 hour-13 solves at a balance residual of 1.1e-6
# (all other residuals zero); 1e-5 is the repository's feasibility rule for the medium and large systems.
const EQTOL = parse(Float64, get(ENV, "FB_EQTOL", "1e-5"))
const GRID = [1.02, 1.05, 1.1, 1.25, 1.5, 2.0]
const LIMIT = 20.0                                   # deg: the search limit of the region campaigns
classify(L, U, fk) = isfinite(L) && L >= 0 ? "feasible" : ((isfinite(U) && U < -UTOL) || (isfinite(fk) && fk > 0) ? "infeasible" : "unresolved")

# ---- the systems ---------------------------------------------------------------------------------------
"Switching commands of a variant master (as multiperiod_study.jl's variant_switches)."
variant_switches(system, master) = [String(strip(l)) for l in eachline(joinpath(REPO_ROOT, "test_networks", system, master))
                                    if occursin(r"(?i)^(Line|Vsource)\.\S+\.enabled=(true|false)$", strip(l))]
const TIES_9500 = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]
"Deck, sources (reference first), switching/added-source commands, whether the band is at customer buses only."
function config(name)
    if name == "ieee123_three"          # the paper's three-source 123-node case, substations 1, 3, 4
        return (; deck = "ieee123_5poi_1ph", sources = ["subs1", "subs3", "subs4"], extra = String[], customers = false)
    elseif name == "ieee123_four"       # the paper's four-source 123-node case (stage four_ieee123a), substations 1-4
        return (; deck = "ieee123_5poi_1ph", sources = ["subs1", "subs2", "subs3", "subs4"], extra = String[], customers = false)
    elseif name == "med3980_four"       # stage four_ieee9500xf: three substations, a fourth at l2804256, five ties closed
        zs = complex(0.104186, 0.617339)
        imp = @sprintf("R1=%.8g X1=%.8g R0=%.8g X0=%.8g", real(zs), imag(zs), real(zs), imag(zs))
        extra = vcat(["New Line.poi4 Phases=1 Bus1=4s Bus2=l2804256 r1=1e-06 x1=1e-06 c1=0 c0=0",
                      "New Vsource.subs4 Phases=1 Bus1=4s basekv=7.1996 pu=1.05 angle=0 $imp"], ["Line.$t.enabled=true" for t in TIES_9500])
        return (; deck = "ieee9500xf_1ph", sources = ["subs1", "subs2", "subs3", "subs4"], extra, customers = true)
    elseif name == "sfp1u_four"         # stage four_smartds_p1u_typeiii
        return (; deck = "smartds_p1u_xf_1ph", sources = ["subs1", "subs2", "subs3", "subs4"],
                extra = variant_switches("smartds_p1u_xf_1ph", "Master_4sub_typeIII.dss"), customers = true)
    elseif (m = match(r"^(ieee123|med3980|sfp1u)_s([1-4]{2,3})$", name)) !== nothing
        # Every subset of a four-source case (user, 2026-10-29): e.g. med3980_s124 keeps substations 1, 2 and 4 of
        # med3980_four, the lowest-numbered one the reference; the dropped ones' sections are fed through the closed ties.
        full = config(Dict("ieee123" => "ieee123_four", "med3980" => "med3980_four", "sfp1u" => "sfp1u_four")[m[1]])
        keep = ["subs$c" for c in m[2]]
        extra = '4' in m[2] ? full.extra : filter(e -> !occursin("subs4", e) && !occursin("poi4", e), full.extra)
        return (; deck = full.deck, sources = keep, extra, customers = full.customers)
    # Two and three sources (2026-10-29), each a subset of a four-source case; the dropped substations' sections
    # are fed through the same closed ties. ieee123: the paper's two-source pair. med3980: its own substations.
    # sf6937: neighbours in the tie chain hs0-hs2-hs1-hs3, i.e. subs4-subs1-subs2-subs3 (subs1 = hs2 the reference).
    elseif name == "ieee123_two"
        return (; deck = "ieee123_5poi_1ph", sources = ["subs3", "subs4"], extra = String[], customers = false)
    elseif name in ("med3980_two", "med3980_three")
        src = name == "med3980_two" ? ["subs1", "subs2"] : ["subs1", "subs2", "subs3"]
        return (; deck = "ieee9500xf_1ph", sources = src, extra = ["Line.$t.enabled=true" for t in TIES_9500], customers = true)
    elseif name in ("sfp1u_two", "sfp1u_three")
        src = name == "sfp1u_two" ? ["subs1", "subs4"] : ["subs1", "subs2", "subs4"]
        return (; deck = "smartds_p1u_xf_1ph", sources = src,
                extra = variant_switches("smartds_p1u_xf_1ph", "Master_4sub_typeIII.dss"), customers = true)
    end
    error("unknown FB_SYSTEM $name")
end
"Measured idle directions of the T = 24 campaigns: (id, movers' unit vector, reach), at_limit ones flagged."
function directions(name, movers)
    if name == "ieee123_three"
        seen = Dict{Float64,Float64}()
        for l in readlines(joinpath(RES, "simplex_check.csv"))[2:end]
            r = split(l, ','); (r[1] == "three" && r[2] == "qfree" && r[4] == "nobatt") && (seen[parse(Float64, r[3])] = parse(Float64, r[5]))
        end
        return [(; id = i, u = [cosd(th), sind(th)], reach = seen[th], at_limit = false, label = @sprintf("%.1f deg", th))
                for (i, th) in enumerate(sort(collect(keys(seen))))]
    end
    # Cases without a measured campaign under these settings: phase A finds the day's reach itself (reach NaN).
    if length(movers) == 1
        return [(; id = i, u = [s], reach = NaN, at_limit = false, label = s > 0 ? "+$(movers[1])" : "-$(movers[1])")
                for (i, s) in enumerate((1.0, -1.0))]
    elseif endswith(name, "_three") || length(movers) == 2
        return [(; id = i, u = [cosd(th), sind(th)], reach = NaN, at_limit = false, label = @sprintf("%.1f deg", th))
                for (i, th) in enumerate(0:22.5:337.5)]
    end
    label = Dict("ieee123_four" => "ieee123_5poi_1ph", "med3980_four" => "ieee9500xf_1ph", "sfp1u_four" => "smartds_p1u_xf_typeIII")[name]
    ds = NamedTuple[]
    for l in readlines(joinpath(RES, "multi_windows_fib78.csv"))[2:end]
        r = split(l, ','); (r[1] == label && r[2] == "subs1:subs2:subs3:subs4" && r[4] == "nobatt") || continue
        push!(ds, (; id = parse(Int, r[5]), u = parse.(Float64, r[6:8]), reach = parse(Float64, r[9]), at_limit = r[12] == "true",
                    label = "direction " * r[5]))
    end
    return sort(ds; by = d -> d.id)
end

cfg = config(SYS)
kw0 = (; sources = cfg.sources, snapshot = true, load_band = (0.5, 1.5), extra = cfg.extra)
pvs = with_logger(NullLogger()) do; read_pv(cfg.deck; kw0...); end
compile_deck(cfg.deck; disable = ["PVSystem"], kw0...)
net = read_network()
S0, movers = cfg.sources[1], cfg.sources[2:end]
band = cfg.customers ? Set(l.to for l in net.lines if startswith(lowercase(l.name), "svc_")) : nothing
lam, price, pvf = load_profile(24), price_profile(24), pv_profile(24)
kw24 = (; pvs, pv_q = :optimize, C_B = cb_by_rule(net), ideal_sources = false, battery = :idle, band_buses = band, tol = 1e-8, constr_viol_tol = 1e-8,
        T = 24, load_factor = lam, price, pv_factor = pvf, gamma = gamma_by_rule(net), terminal = :soft)
hour_kw(t) = merge(kw24, (; T = 1, dt_h = 1.0, load_factor = [lam[t]], price = [price[t]], pv_factor = [pvf[t]], gamma = 0.0, terminal = :free))
DIRS = directions(SYS, movers)
angles(d, r) = merge(Dict(S0 => 0.0), Dict(m => r * d.u[k] for (k, m) in enumerate(movers)))
margin(ang, kw; start = nothing) = solve_multiperiod(net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0, start)
"Validated margin of an exact solve (retried cold and from `alt`; NaN if no start validates), and the solution."
function vmargin(ang, kw; start = nothing, alt = nothing)
    x = nothing; v = nothing
    for st in unique([start, nothing, alt])
        x = margin(ang, kw; start = st)
        x.status in RAY_OK || continue
        v = validate_dispatch(net, ang, kw, x; k = KV)
        v.equations <= EQTOL && return v.margin_kW, x
    end
    say("      solve failed from every start (", x.status, ")", v === nothing ? "" :
        "; residuals " * join((@sprintf("%s %.1e", c, r) for (c, r) in sort(collect(v.viol))), ", "))
    return NaN, x
end
say("== Farkas bracketing ", SYS, " phase ", PHASE, ": ", length(net.buses), " buses, sources ", join(cfg.sources, ","),
    ", band on ", band === nothing ? "every non-source bus" : "$(length(band)) customer buses", "; ", length(DIRS), " directions")

"The phase, in a function so that its loops update their own variables."
function main()
    if PHASE == "A"
        i = parse(Int, ENV["FB_DIR"]); d = DIRS[i]
        fa = joinpath(OUT, "fb_A_$(d.id).csv")
        open(io -> println(io, "id,label,day_reach_deg,at_limit,tstar,hour_boundary_deg,margin98_kW,u1,u2,u3,seconds"), fa, "w")
        t0 = time()
        if isnan(d.reach)
            # The day's reach by bisection on the validated T = 24 margin (batteries idle, so it is the worst hour's).
            # A point that does not validate (after the retries of probe) counts as beyond the edge: only validated
            # points move the lower end, and phase C certifies the other side.
            # A solve that does not validate at one angle often does a fraction of a percent away (at 1.0 and 0.5 deg on
            # med3980 the residuals stayed at 1e-3 while every nearby angle validated), so a failed point is retried at
            # r(1 -+ 0.3%) before it counts as beyond.
            function probe(r, xs)
                m, x = vmargin(angles(d, r), kw24; start = xs, alt = xs)
                for f in (0.997, 1.003)
                    isfinite(m) && return m, x, r
                    say(@sprintf("  %s: no validated solve at %.4f deg; trying %.4f", d.label, r, f * r))
                    r2 = f * r; m, x = vmargin(angles(d, r2), kw24; start = xs, alt = xs)
                    isfinite(m) && return m, x, r2
                end
                return m, x, r
            end
            lo, hi, xs = 0.0, 1.0, nothing
            while true
                m, x, r = probe(hi, xs)
                (isfinite(m) && m >= 0) || (hi = r; break)
                lo = r; xs = x
                hi >= LIMIT && (d = merge(d, (; reach = LIMIT, at_limit = true)); break)
                hi = min(2hi, LIMIT)
            end
            if !d.at_limit
                for _ in 1:14
                    mid = (lo + hi) / 2; m, x, r = probe(mid, xs)
                    isfinite(m) && m >= 0 ? (lo = r; xs = x) : (hi = r)
                end
                d = merge(d, (; reach = lo))
            end
            say(@sprintf("  %s: day's reach %.5f deg%s (%.0f s)", d.label, d.reach, d.at_limit ? " (search limit)" : "", time() - t0))
        end
        if d.at_limit
            open(io -> @printf(io, "%d,%s,%.6f,true,0,NaN,NaN,%s,0\n", d.id, d.label, d.reach, join(vcat(d.u, fill(0.0, 3 - length(d.u))), ',')), fa, "a")
            say("  ", d.label, ": feasible to the 20 deg search limit; nothing to bracket"); exit(0)
        end
        x24 = margin(angles(d, 0.98d.reach), kw24)
        x24.status in RAY_OK || error("T = 24 solve at 98% of the measured reach failed: $(x24.status)")
        hourmin = [minimum(x24.P_subs_kW[n][t] for n in keys(x24.P_subs_kW)) for t in 1:24]
        tstar = argmin(hourmin); kw1 = hour_kw(tstar)
        say(@sprintf("  %s: limiting hour %d (smallest import %.1f kW at 98%% of the day's reach %.5f deg)", d.label, tstar, hourmin[tstar], d.reach))
        lo, hi = 0.0, 1.2d.reach; xs = nothing
        m, x = vmargin(angles(d, 0.999d.reach), kw1)
        isfinite(m) && m >= 0 ? (lo = 0.999d.reach; xs = x) : say("  the day's edge does not validate at this hour; bracketing from zero")
        while true
            m, _ = vmargin(angles(d, hi), kw1); (isfinite(m) && m < 0) && break
            lo = hi; hi *= 1.5; hi > 60 && error("no edge within 60 deg")
        end
        for _ in 1:16
            mid = (lo + hi) / 2; m, x = vmargin(angles(d, mid), kw1; start = xs, alt = xs)
            isfinite(m) && m >= 0 ? (lo = mid; xs = x) : (hi = mid)
        end
        m98, x98 = vmargin(angles(d, 0.98lo), kw1; start = xs, alt = xs)
        m0, x0 = vmargin(angles(d, 0.0), kw1)
        (isfinite(m98) && m98 >= 0 && isfinite(m0) && m0 >= 0) || error("the 98% point or equal angles does not validate at hour $tstar")
        rp = replay_opendss(net, x98, angles(d, 0.98lo); load_factor = [lam[tstar]], pvs, pv_factor = [pvf[tstar]])
        open(io -> @printf(io, "%d,%s,%.6f,false,%d,%.6f,%.4f,%s,%.0f\n", d.id, d.label, d.reach, tstar, lo, m98,
                           join(vcat(d.u, fill(0.0, 3 - length(d.u))), ','), time() - t0), fa, "a")
        open(joinpath(OUT, "fb_Acur_$(d.id).csv"), "w") do io
            println(io, "element,kind,I_pu")
            for x in (x98, x0)
                for (n, r) in x.raw.lr; @printf(io, "%s,line,%.8f\n", n, hypot(r[1], x.raw.li[n][1])); end
                for (n, r) in x.raw.ir; @printf(io, "%s,source,%.8f\n", n, hypot(r[1], x.raw.ii[n][1])); end
            end
        end
        say(@sprintf("  hour boundary %.5f deg (day %.5f); 98%% point margin %.3f kW, OpenDSS replay dP %.1e kW dV %.1e pu (%s); %.0f s",
                     lo, d.reach, m98, rp.dP, rp.dV, rp.converged ? "converged" : "NOT converged", time() - t0))

    elseif PHASE == "B"
        env = Dict{Tuple{String,String},Float64}(); nfiles = 0
        for f in readdir(OUT); startswith(f, "fb_Acur_") || continue; nfiles += 1
            for l in readlines(joinpath(OUT, f))[2:end]
                c = split(l, ','); k = (String(c[1]), String(c[2])); env[k] = max(get(env, k, 0.0), parse(Float64, c[3]))
            end
        end
        open(joinpath(OUT, "fb_envelope.csv"), "w") do io
            println(io, "element,kind,envelope_pu,limit_pu")
            for k in sort(collect(keys(env))); @printf(io, "%s,%s,%.8f,%.8f\n", k[1], k[2], env[k], k[2] == "line" ? max(env[k], FLOOR) : env[k]); end
        end
        say(@sprintf("  envelope from %d directions' dispatches: %d lines (%d at the %.0e pu floor), %d sources", nfiles,
                     count(k -> k[2] == "line", keys(env)), count(k -> k[2] == "line" && env[k] < FLOOR, keys(env)), FLOOR, count(k -> k[2] == "source", keys(env))))

    elseif PHASE == "C"
        i = parse(Int, ENV["FB_DIR"]); d = DIRS[i]
        fa = joinpath(OUT, "fb_A_$(d.id).csv"); a = split(readlines(fa)[2], ',')
        fc = joinpath(OUT, "fb_C_$(d.id).csv")
        open(io -> println(io, "id,label,tstar,hour_boundary_deg,r_cert_deg,r_unresolved_deg,ratio_cert,U_at_cert_kW,farkas_at_cert,relaxations,seconds"), fc, "w")
        a[4] == "true" && (say("  ", d.label, ": at the search limit; skipped"); exit(0))
        tstar = parse(Int, a[5]); rs = parse(Float64, a[6])
        amp = Dict{String,Float64}(); senv = Dict{String,Float64}()
        for l in readlines(joinpath(OUT, "fb_envelope.csv"))[2:end]
            c = split(l, ','); c[2] == "line" ? (amp[String(c[1])] = parse(Float64, c[4])) : (senv[String(c[1])] = parse(Float64, c[4]))
        end
        kwL = merge(hour_kw(tstar), (; ampacity_pu = amp, source_ampacity_pu = senv))
        certified(r) = (c = certify(net, angles(d, r), kwL; k = KV, cuts = CUTS); (classify(NaN, c.U, c.farkas.bound) == "infeasible", c))
        t0 = time(); lo = rs; hi = NaN; Uhi = NaN; Fhi = NaN; n = 0
        for g in GRID
            ok, c = certified(g * rs); n += 1
            say(@sprintf("  %s: %.3f x r* = %.5f deg -> U %+.4g kW, Farkas %.4g (%s)%s", d.label, g, g * rs, c.U, c.farkas.bound, c.status, ok ? "  certified" : ""))
            ok ? (hi = g * rs; Uhi = c.U; Fhi = c.farkas.bound; break) : (lo = g * rs)
        end
        if isfinite(hi)
            for _ in 1:5
                mid = (lo + hi) / 2; ok, c = certified(mid); n += 1
                ok ? (hi = mid; Uhi = c.U; Fhi = c.farkas.bound) : (lo = mid)
            end
        end
        open(io -> @printf(io, "%d,%s,%d,%.6f,%.6f,%.6f,%.5f,%.6g,%.6g,%d,%.0f\n", d.id, d.label, tstar, rs, hi, lo, hi / rs, Uhi, Fhi, n, time() - t0), fc, "a")
        say(@sprintf("  %s: r* %.5f, certified from %.5f (x%.4f), unresolved up to %.5f; %d relaxations, %.0f s", d.label, rs, hi, hi / rs, lo, n, time() - t0))
    end
end
main()
