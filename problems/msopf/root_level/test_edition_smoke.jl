# Smoke test of the test-edition decks (2026-10-24): do they solve, and how do their voltages and imports sit?
#   julia --project=problems/msopf problems/msopf/root_level/test_edition_smoke.jl [case ...]
# Each deck with secondaries, all its sources at equal angles, T = 24, PV reactive power dispatched, batteries
# dispatched with the hard cyclic SOC, no-backflow at every source and hour:
#   (a) cost-optimal day with no voltage band (where do the voltages sit on their own?);
#   (b) cost-optimal day with the 0.95-1.05 band at the customer (secondary) buses only (the 2026-10-24 rule);
#   (c) if (b) does not solve: the violation day (least total backflow and band violation; zero = feasible).
# Every solved day is replayed in OpenDSS. Reports status, Ipopt iterations, time, the customer and primary voltage
# ranges, the smallest import, and the replay errors. Writes results/test_networks/smoke.csv and smoke.log.
include(joinpath(@__DIR__, "..", "multiperiod_full_angle.jl"))
using Logging, Printf

const OUT = normpath(joinpath(@__DIR__, "..", "results", "test_networks"))
mkpath(OUT)
say(args...) = (println(args...); flush(stdout))
const OKST = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)

"The variant master's switching: (sources kept, extra commands)."
function variant(deck, file)
    cmds = [String(strip(l)) for l in readlines(joinpath(REPO_ROOT, "test_networks", deck, file)) if occursin(r"(?i)^(Line|Vsource)\.\S+\.enabled=", strip(l))]
    off = [lowercase(match(r"(?i)^Vsource\.(\S+)\.enabled=false", c)[1]) for c in cmds if occursin(r"(?i)^Vsource\..*enabled=false", c)]
    master = read(joinpath(REPO_ROOT, "test_networks", deck, "Master.dss"), String)
    srcs = [lowercase(m[1]) for m in eachmatch(r"(?i)New\s+Vsource\.(\S+)", master)]
    return [s for s in srcs if !(s in off)], cmds
end

const P1U3 = variant("smartds_p1u_xf_1ph", "Master_4sub_typeIII.dss")
const CASES = [
    (; label = "K1", deck = "epri_k1xf_1ph", sources = ["subs1"], extra = String[]),
    (; label = "Ckt24", deck = "epri_ckt24xf_1ph", sources = ["subs1"], extra = String[]),
    (; label = "8500", deck = "ieee8500xf_1ph", sources = ["subs1"], extra = String[]),
    (; label = "9500 (3 substations, 5 ties closed)", deck = "ieee9500xf_1ph", sources = ["subs1", "subs2", "subs3"],
       extra = ["Line.$t.enabled=true" for t in ("a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw")]),
    (; label = "P1U Type III (4 substations)", deck = "smartds_p1u_xf_1ph", sources = P1U3[1], extra = P1U3[2]),
]

function run_case(c)
    say("\n== ", c.label, " (", c.deck, ")")
    t0 = time()
    kw0 = (; sources = c.sources, snapshot = true, load_band = (0.5, 1.5), extra = c.extra)
    pvs = with_logger(NullLogger()) do; read_pv(c.deck; kw0...); end
    compile_deck(c.deck; disable = ["PVSystem"], kw0...)
    net = read_network()
    tload = time() - t0
    T = 24
    lam, price, pvf = load_profile(T), price_profile(T), pv_profile(T)
    cust = Set(l.to for l in net.lines if startswith(lowercase(l.name), "svc_"))
    srcb = Set(s.bus for s in net.sources)
    prim = [b for b in net.buses if !(b in cust) && !(b in srcb)]
    say(@sprintf("  %d buses (%d customer, %d primary, %d source), %d lines, %d PV, %d batteries; loaded in %.1f s",
                 length(net.buses), length(cust), length(prim), length(srcb), length(net.lines), length(pvs), length(net.batteries), tload))
    kw = (; T, load_factor = lam, price, pvs, pv_factor = pvf, pv_q = :optimize, battery = :optimize, terminal = :cyclic,
          C_B = cb_by_rule(net), gamma = 0.0)
    ang = Dict(s.name => 0.0 for s in net.sources)
    rows = Any[]
    function report(tag, x)
        okk = x.status in OKST
        vr(bs) = isempty(bs) || !okk ? (NaN, NaN) : extrema(abs(x.V_pu[b][t]) for b in bs, t in 1:T)
        vc, vp = vr(collect(cust)), vr(prim)
        minimp = okk ? minimum(minimum(v) for v in values(x.P_subs_kW)) : NaN
        rp = okk ? replay_opendss(net, x, ang; load_factor = lam, pvs, pv_factor = pvf) : nothing
        say(@sprintf("  (%s) %s, %s iterations, %.1f s | customer V %.4f-%.4f, primary V %.4f-%.4f pu | smallest import %.1f kW%s%s",
                     tag, x.status, string(x.iterations), x.time, vc..., vp..., minimp,
                     hasproperty(x, :violation) && isfinite(something(x.violation, NaN)) ? @sprintf(" | violation %.3e", x.violation) : "",
                     rp === nothing ? "" : @sprintf(" | OpenDSS replay: converged %s, dP %.1e kW, dV %.1e pu", rp.converged, rp.dP, rp.dV)))
        push!(rows, ["\"" * c.label * "\"", c.deck, tag, string(x.status), string(x.iterations), @sprintf("%.1f", x.time), @sprintf("%.4f", vc[1]), @sprintf("%.4f", vc[2]),
                     @sprintf("%.4f", vp[1]), @sprintf("%.4f", vp[2]), @sprintf("%.2f", minimp), rp === nothing ? "" : string(rp.converged),
                     rp === nothing ? "" : @sprintf("%.2e", rp.dP), rp === nothing ? "" : @sprintf("%.2e", rp.dV), length(net.buses), length(cust)])
        return okk
    end
    xa = solve_multiperiod(net, ang; kw..., V_band = nothing, objective = :cost)
    report("a: no band", xa)
    xb = solve_multiperiod(net, ang; kw..., V_band = (0.95, 1.05), band_buses = cust, objective = :cost, start = xa.status in OKST ? xa : nothing)
    okb = report("b: band at customer buses", xb)
    if !okb
        xc = solve_multiperiod(net, ang; kw..., V_band = (0.95, 1.05), band_buses = cust, objective = :violation, start = xa.status in OKST ? xa : nothing)
        report("c: violation", xc)
    end
    open(joinpath(OUT, "smoke.csv"), "a") do io
        foreach(r -> println(io, join(r, ',')), rows)
    end
end

function main(sel)
    f = joinpath(OUT, "smoke.csv")
    isfile(f) || open(io -> println(io, "case,deck,solve,status,iterations,seconds,cust_vmin,cust_vmax,prim_vmin,prim_vmax,min_import_kW,replay_converged,replay_dP_kW,replay_dV_pu,buses,customer_buses"), f, "w")
    for c in CASES
        (isempty(sel) || any(s -> occursin(lowercase(s), lowercase(c.label * c.deck)), sel)) && run_case(c)
    end
end
main(ARGS)
