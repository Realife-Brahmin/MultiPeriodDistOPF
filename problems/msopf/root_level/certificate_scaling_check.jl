# Which data make the relaxation badly scaled on ieee9500xf_1ph (2026-10-26)? Builds the hour-13 relaxation of
# certificate_secondary_pilot.jl at equal angles, with the envelope limits it wrote, and lists the
# coefficient range of the constraints by kind, the line and source impedances at both extremes, and
# how many lines have near-zero impedance.
#   julia --project=problems/msopf/environments/sdp problems/msopf/root_level/certificate_scaling_check.jl
include(joinpath(@__DIR__, "..", "multiperiod_full_angle.jl"))
include(joinpath(@__DIR__, "..", "feasibility_certificate.jl"))
using Logging
const SYSTEM = "ieee9500xf_1ph"
const TIES = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]
kw0 = (; sources = ["subs1", "subs2", "subs3"], snapshot = true, load_band = (0.5, 1.5), extra = ["Line.$t.enabled=true" for t in TIES])
pvs = with_logger(NullLogger()) do; read_pv(SYSTEM; kw0...); end
compile_deck(SYSTEM; disable = ["PVSystem"], kw0...)
net = read_network()
cust = Set(l.to for l in net.lines if startswith(lowercase(l.name), "svc_"))
Yb = net.Y_base
z = Dict(l.name => abs(-Yb / l.Y[1, 2]) for l in net.lines)
zs = sort(collect(z); by = last)
println("lines: ", length(zs), "; |z| pu: min ", zs[1], ", median ", zs[end ÷ 2][2], ", max ", zs[end])
for thr in (1e-6, 1e-4, 1e-3)
    println(@sprintf("  |z| < %.0e pu: %d lines", thr, count(p -> p[2] < thr, zs)))
end
println("  smallest: ", join([@sprintf("%s %.2e", n, v) for (n, v) in zs[1:8]], ", "))
for s in net.sources; println(@sprintf("source %s |z_s| %.3e pu (1/|z| %.3e)", s.name, abs(s.z_series * Yb), 1 / abs(s.z_series * Yb))); end
env = Dict{String,Float64}(); senv = Dict{String,Float64}()
for l in readlines(joinpath(@__DIR__, "..", "results", "multiperiod", "certificate_secondary_envelope.csv"))[2:end]
    c = split(l, ','); c[2] == "line" ? (env[String(c[1])] = parse(Float64, c[5])) : (senv[String(c[1])] = parse(Float64, c[4]))
end
lam, price, pvf = load_profile(24), price_profile(24), pv_profile(24); t = 13
kw = (; pvs, pv_q = :optimize, C_B = cb_by_rule(net), ideal_sources = false, battery = :idle, band_buses = cust,
      T = 1, dt_h = 1.0, load_factor = [lam[t]], price = [price[t]], pv_factor = [pvf[t]], gamma = 0.0, terminal = :free,
      ampacity_pu = env, source_ampacity_pu = senv)
R = relax_margin(net, Dict("subs1" => 0.0, "subs2" => 0.0, "subs3" => 0.0), kw; k = 20000.0)
lo, hi = Dict{String,Float64}(), Dict{String,Float64}(); wlo, whi = Dict{String,String}(), Dict{String,String}()
for (F, S) in list_of_constraint_types(R.model)
    F == VariableRef && continue
    kind = string(nameof(S))
    for c in all_constraints(R.model, F, S)
        f = constraint_object(c).func
        terms = f isa AffExpr ? collect(values(f.terms)) : f isa AbstractVector ? [a for g in f for a in values(convert(AffExpr, g).terms)] :
                [t.scalar_term.coefficient for t in f.terms]
        for a in terms
            a == 0 && continue; v = abs(a)
            v < get(lo, kind, Inf) && (lo[kind] = v; wlo[kind] = first(string(c), 90))
            v > get(hi, kind, 0.0) && (hi[kind] = v; whi[kind] = first(string(c), 90))
        end
    end
end
for k in sort(collect(keys(lo)))
    println(@sprintf("%-28s |coef| %.2e .. %.2e (ratio %.1e)", k, lo[k], hi[k], hi[k] / lo[k]))
    println("    largest in: ", whi[k])
end
