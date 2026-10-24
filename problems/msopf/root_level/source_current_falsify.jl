# Falsification test for the derived source-current bounds (2026-10-22).
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/source_current_falsify.jl
# At each query angle of source_current_bounds.jl, maximize every source current in the exact model at those
# FIXED angles, over the domain the analytical bound claims (network equations, band, device limits; no-backflow
# dropped, so the search covers every state the bound must hold for), from several starts. A validated dispatch
# above a derived bound would show the bound wrong. Local maxima: a pass is evidence, not a proof.
ENV["CV_EXTRA"] = "lib"
include(joinpath(@__DIR__, "certificate_validation.jl"))
const SCB = joinpath(OUT_DIR, "source_current_bounds")
const IbA = net.S_base_kVA / net.kV_base

function falsify()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw1 = hour_kw(s24, 12)
    rows = [split(l, ',') for l in readlines(SCB * "_bounds.csv")[2:end]]
    out = SCB * "_falsify.csv"
    open(io -> println(io, "theta_deg,factor,r_deg,source,no_backflow,start,status,equations,band_violation,I_A,derived_bound_A,ratio"), out, "w")
    worst = 0.0
    for key in unique([(r[1], r[2], r[3]) for r in rows])
        th, fac, r = parse(Float64, key[1]), key[2], parse(Float64, key[3]); ang = at(r, th)
        bnd = Dict(String(q[4]) => parse(Float64, q[13]) for q in rows if (q[1], q[2], q[3]) == key)       # A
        starts = Any[("cold", nothing)]
        xm = solve_multiperiod(net, ang; kw1..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0); xm.status in RAY_OK && push!(starts, ("margin", xm))
        xr = solve_multiperiod(net, ang; kw1..., objective = :margin, margin_band_hard = true, margin_reg = 0.0); xr.status in RAY_OK && push!(starts, ("restoration", xr))
        best = Dict(s.name => 0.0 for s in net.sources)
        for nb in (false, true), pass in 1:2, s in net.sources, (sn, sx) in copy(starts)
            x = solve_multiperiod(net, ang; kw1..., objective = :maxcurrent, current_target = (:source, s.name), start = sx, no_backflow = nb)
            x.status in RAY_OK || continue
            v = validate_dispatch(net, ang, kw1, x; k = KV)
            (v.equations <= 1e-6 && v.viol["band"] <= 1e-6 && (!nb || v.viol["backflow"] <= 1e-6)) || continue
            pass == 1 && push!(starts, ("max $(s.name) ($(nb ? "nb" : "free"))", x))
            for t in net.sources
                I = abs(complex(x.raw.ir[t.name][1], x.raw.ii[t.name][1])) * IbA
                if t.name == s.name
                    open(io -> println(io, join([key..., s.name, nb, replace(sn, ',' => ';'), x.status, @sprintf("%.1e", v.equations), @sprintf("%.1e", v.viol["band"]),
                                                 @sprintf("%.3f", I), @sprintf("%.3f", bnd[s.name]), @sprintf("%.4f", I / bnd[s.name])], ',')), out, "a")
                end
                best[t.name] = max(best[t.name], I); worst = max(worst, I / bnd[t.name])
            end
        end
        say(@sprintf("theta %.1f, %s r*: largest validated source currents %s (derived bounds %s A)", th, fac,
                     join([@sprintf("%s %.2f", n, best[n]) for n in sort(collect(keys(best)))], ", "), join([@sprintf("%.2f", bnd[n]) for n in sort(collect(keys(bnd)))], "/")))
    end
    say(@sprintf("largest ratio of a validated source current to its derived bound: %.4f (%s)", worst, worst <= 1 ? "no bound violated" : "A BOUND IS VIOLATED"))
end
falsify()
