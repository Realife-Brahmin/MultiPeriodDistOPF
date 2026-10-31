# Farkas bracketing with the batteries DISPATCHED (2026-10-30): the certified feasibility map of the real problem,
# T = 24 day, hard cyclic terminal B[T] == B[0], PV reactive power dispatched. farkas_bracketing.jl maps the region
# with the batteries idle, an inner bound (an idle day meets the cyclic terminal); this one maps the region itself.
# With the batteries dispatched the hours no longer decouple, so every solve and every relaxation spans the day.
#   AD  (FB_DIR = i) the day's reach along direction i: starts at the idle reach (feasible, fb_A_<i>.csv of the idle
#       campaign in FB_IDLE, from the stitched idle hours), steps out by 1.25 until a validated margin is negative or
#       nothing validates, then regula falsi (Illinois) on the validated day margin to 0.5% of the reach. Writes
#       fb_AD_<i>.csv, the currents of the dispatch at 98% of the reach (fb_ADcur_<i>.csv, largest over the day) and
#       the OpenDSS replay of that day.
#   BD  the current envelope over every fb_ADcur and the idle campaign's fb_Acur files (as phase B).
#   CD  (FB_DIR = i) certifies infeasibility of the T = 24 relaxation outward from the dispatched reach (grid, then
#       five bisections, as phase C), with the envelope as current limits. Writes fb_CD_<i>.csv.
#   FB_SYSTEM=sfp1u_four FB_PHASE=AD FB_DIR=1 FB_IDLE=<idle campaign folder> julia --project=... farkas_dispatched.jl
const DPHASE = get(ENV, "FB_PHASE", "AD")
ENV["FB_PHASE"] = "X"                                # farkas_bracketing.jl: load its setup, run none of its phases
include(joinpath(@__DIR__, "farkas_bracketing.jl"))
const IDLE = get(ENV, "FB_IDLE", OUT)
# Ipopt to 1e-6 (tol and the primal infeasibility constr_viol_tol; the idle campaign used 1e-8): every point is
# validated against the original constraints at EQTOL = 1e-5 anyway, and the edge is sought to 0.5% (user, 2026-10-30).
kwb = merge(kw24, (; battery = :optimize, terminal = :cyclic, gamma = 0.0, tol = 1e-6, constr_viol_tol = 1e-6))

"A T = 24 start stitched from the hourly solutions `xs` (batteries at their initial charge, no power)."
function stitch(xs)
    cat24(f) = Dict(k => [f(xs[t])[k][1] for t in 1:24] for k in keys(f(xs[1])))
    return (; V_pu = cat24(x -> x.V_pu), P_subs_kW = cat24(x -> x.P_subs_kW), Q_subs_kvar = cat24(x -> x.Q_subs_kvar),
            Q_PV_kvar = cat24(x -> x.Q_PV_kvar), P_B_kW = Dict(k => zeros(24) for k in keys(xs[1].P_B_kW)),
            B_kWh = Dict(k => fill(v[1], 25) for (k, v) in xs[1].B_kWh),
            raw = (; lr = cat24(x -> x.raw.lr), li = cat24(x -> x.raw.li), ir = cat24(x -> x.raw.ir), ii = cat24(x -> x.raw.ii)))
end

"Validated dispatched day margin at angle r along d, from `start`, cold only if that does not validate (cold retries of
negative margins agreed to 2 kW and doubled the cost, 2026-10-30); NaN if neither validates."
function dmargin(d, r, start)
    for st in unique([start, nothing])
        t1 = time(); x = margin(angles(d, r), kwb; start = st)
        if x.status in RAY_OK
            v = validate_dispatch(net, angles(d, r), kwb, x; k = KV)
            say(@sprintf("    %.5f deg (%s): margin %.1f kW, residual %.1e (%.0f s)", r, st === nothing ? "cold" : "warm", v.margin_kW, v.equations, time() - t1))
            if v.equations <= EQTOL
                return v.margin_kW, x
            end
        else
            say(@sprintf("    %.5f deg (%s): %s (%.0f s)", r, st === nothing ? "cold" : "warm", x.status, time() - t1))
        end
    end
    return NaN, nothing
end

function dmain()
    if DPHASE == "AD"
        i = parse(Int, ENV["FB_DIR"]); d = DIRS[i]; t0 = time()
        a = split(readlines(joinpath(IDLE, "fb_A_$(d.id).csv"))[2], ',')
        ridle = parse(Float64, a[3]); ucols = join(vcat(d.u, fill(0.0, 3 - length(d.u))), ',')
        fa = joinpath(OUT, "fb_AD_$(d.id).csv")
        open(io -> println(io, "id,label,idle_reach_deg,reach_deg,beyond_deg,at_limit,solves,u1,u2,u3,seconds"), fa, "w")
        if a[4] == "true"                            # idle feasible to the limit: so is the dispatched day
            open(io -> @printf(io, "%d,%s,%.6f,%.6f,NaN,true,0,%s,0\n", d.id, d.label, ridle, LIMIT, ucols), fa, "a")
            say("  ", d.label, ": idle feasible to the search limit; so is the dispatched day"); exit(0)
        end
        isfinite(ridle) || error("direction $(d.id) has no idle reach (equal angles infeasible when idle): not handled yet")
        _, _, _, xs = daymargin(angles(d, ridle); all = true)
        any(x -> x === nothing, xs) && error("the idle hours at the idle reach do not validate")
        n = 0; lo, mlo, xlo = ridle, NaN, stitch(xs)
        m, x = dmargin(d, lo, xlo); n += 1
        (isfinite(m) && m >= 0) || error(@sprintf("the dispatched day does not validate as feasible at the idle reach (margin %.1f kW)", m))
        mlo, xlo = m, x
        hi, mhi = NaN, NaN; r = lo
        while true                                   # step out until the margin is negative or nothing validates
            r = min(1.25r, LIMIT); m, x = dmargin(d, r, xlo); n += 1
            if isfinite(m) && m >= 0
                lo, mlo, xlo = r, m, x
                r >= LIMIT && break
            else
                hi, mhi = r, m; break
            end
        end
        if isfinite(hi)
            side = 0
            while (hi - lo) > 5e-3 * lo
                # Illinois regula falsi when the far side has a validated margin; bisection when it does not
                r = isfinite(mhi) ? lo + (hi - lo) * mlo / (mlo - mhi) : (lo + hi) / 2
                r = clamp(r, lo + 0.05(hi - lo), hi - 0.05(hi - lo))
                m, x = dmargin(d, r, xlo); n += 1
                if isfinite(m) && m >= 0
                    lo, mlo, xlo = r, m, x; side == 1 && isfinite(mhi) && (mhi /= 2); side = 1
                else
                    hi, mhi = r, m; side == -1 && (mlo /= 2); side = -1
                end
            end
        end
        atlim = !isfinite(hi)
        say(@sprintf("  %s: dispatched reach %.5f deg (idle %.5f, x%.3f)%s; %d solves, %.0f s", d.label, lo, ridle, lo / ridle,
                     atlim ? " (search limit)" : "", n, time() - t0))
        m98, x98 = dmargin(d, 0.98lo, xlo)
        (isfinite(m98) && m98 >= 0) || error("the 98% point does not validate")
        rp = replay_opendss(net, x98, angles(d, 0.98lo); load_factor = lam, pvs, pv_factor = pvf)
        open(io -> @printf(io, "%d,%s,%.6f,%.6f,%.6f,%s,%d,%s,%.0f\n", d.id, d.label, ridle, lo, hi, atlim, n, ucols, time() - t0), fa, "a")
        open(joinpath(OUT, "fb_ADcur_$(d.id).csv"), "w") do io
            println(io, "element,kind,I_pu")
            for (nm, v) in x98.raw.lr; @printf(io, "%s,line,%.8f\n", nm, maximum(hypot.(v, x98.raw.li[nm]))); end
            for (nm, v) in x98.raw.ir; @printf(io, "%s,source,%.8f\n", nm, maximum(hypot.(v, x98.raw.ii[nm]))); end
        end
        say(@sprintf("  98%% point margin %.1f kW; battery swing %.0f kWh; OpenDSS replay of the day dP %.1e kW dV %.1e pu (%s); %.0f s",
                     m98, sum(maximum(b) - minimum(b) for b in values(x98.B_kWh)), rp.dP, rp.dV, rp.converged ? "converged" : "NOT converged", time() - t0))

    elseif DPHASE == "BD"
        env = Dict{Tuple{String,String},Float64}(); nfiles = 0
        for (dir, pre) in ((OUT, "fb_ADcur_"), (IDLE, "fb_Acur_")), f in readdir(dir)
            startswith(f, pre) || continue; nfiles += 1
            for l in readlines(joinpath(dir, f))[2:end]
                c = split(l, ','); k = (String(c[1]), String(c[2])); env[k] = max(get(env, k, 0.0), parse(Float64, c[3]))
            end
        end
        open(joinpath(OUT, "fb_envelope_disp.csv"), "w") do io
            println(io, "element,kind,envelope_pu,limit_pu")
            for k in sort(collect(keys(env))); @printf(io, "%s,%s,%.8f,%.8f\n", k[1], k[2], env[k], k[2] == "line" ? max(env[k], FLOOR) : env[k]); end
        end
        say(@sprintf("  envelope from %d files: %d lines, %d sources", nfiles, count(k -> k[2] == "line", keys(env)), count(k -> k[2] == "source", keys(env))))

    elseif DPHASE == "CD"
        i = parse(Int, ENV["FB_DIR"]); d = DIRS[i]
        a = split(readlines(joinpath(OUT, "fb_AD_$(d.id).csv"))[2], ',')
        fc = joinpath(OUT, "fb_CD_$(d.id).csv")
        open(io -> println(io, "id,label,reach_deg,r_cert_deg,r_unresolved_deg,ratio_cert,U_at_cert_kW,farkas_at_cert,relaxations,seconds,maxrss_GB"), fc, "w")
        a[6] == "true" && (say("  ", d.label, ": at the search limit; skipped"); exit(0))
        rs = parse(Float64, a[4])
        envf = get(ENV, "FB_ENVELOPE", joinpath(OUT, "fb_envelope_disp.csv"))
        amp = Dict{String,Float64}(); senv = Dict{String,Float64}()
        for l in readlines(envf)[2:end]
            c = split(l, ','); c[2] == "line" ? (amp[String(c[1])] = parse(Float64, c[4])) : (senv[String(c[1])] = parse(Float64, c[4]))
        end
        kwL = merge(kwb, (; ampacity_pu = amp, source_ampacity_pu = senv))
        function certified(r)
            t1 = time(); c = certify(net, angles(d, r), kwL; k = KV, cuts = CUTS)
            ok = classify(NaN, c.U, c.farkas.bound) == "infeasible"
            say(@sprintf("  %s: %.4f x r = %.5f deg -> U %+.4g kW, Farkas %.4g (%s)%s; build %.0f s, solve %.0f s, %d blocks, max RSS %.1f GB",
                         d.label, r / rs, r, c.U, c.farkas.bound, c.status, ok ? "  certified" : "", c.build_s, c.solve_s, c.nblocks, Sys.maxrss() / 2^30))
            return ok, c
        end
        # FB_CD_GRID (comma list) replaces the grid, FB_CD_BISECT (default 5) the number of bisections.
        grid = haskey(ENV, "FB_CD_GRID") ? parse.(Float64, split(ENV["FB_CD_GRID"], ',')) : GRID
        nbis = parse(Int, get(ENV, "FB_CD_BISECT", "5"))
        t0 = time(); lo = rs; hi = NaN; Uhi = NaN; Fhi = NaN; n = 0
        for g in grid
            ok, c = certified(g * rs); n += 1
            ok ? (hi = g * rs; Uhi = c.U; Fhi = c.farkas.bound; break) : (lo = g * rs)
            get(ENV, "FB_CD_ONCE", "0") == "1" && break
        end
        if isfinite(hi) && get(ENV, "FB_CD_ONCE", "0") != "1"
            for _ in 1:nbis
                mid = (lo + hi) / 2; ok, c = certified(mid); n += 1
                ok ? (hi = mid; Uhi = c.U; Fhi = c.farkas.bound) : (lo = mid)
            end
        end
        open(io -> @printf(io, "%d,%s,%.6f,%.6f,%.6f,%.5f,%.6g,%.6g,%d,%.0f,%.1f\n", d.id, d.label, rs, hi, lo, hi / rs, Uhi, Fhi, n, time() - t0, Sys.maxrss() / 2^30), fc, "a")
        say(@sprintf("  %s: reach %.5f, certified from %.5f (x%.4f), unresolved up to %.5f; %d relaxations, %.0f s", d.label, rs, hi, hi / rs, lo, n, time() - t0))
    end
end
dmain()
