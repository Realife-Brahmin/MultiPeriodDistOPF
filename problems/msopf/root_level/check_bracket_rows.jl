# Two rows of multi_windows_fib78_v1.03.csv found by the bracket fallback (2705-bus feeder, four
# substations, sources at 1.03 pu, batteries dispatched, directions 47 and 65; binding = none, about
# 1.15 kW of import left). Are they slightly interior endpoints, a labelling-threshold effect, or
# search failures? Margin solves (MULTI_VMARGIN = 20000, as in the run) at 100, 100.1, 100.2 and
# 100.5% of the recorded reach, warm-started along the direction, and at 100 and 100.1% from a cold
# start. The recorded rows are not touched; this writes check_bracket_rows.csv.
#   MULTI_VSRC=1.03 julia --project=problems/msopf problems/msopf/root_level/check_bracket_rows.jl
empty!(ARGS); push!(ARGS, "none")
include(joinpath(@__DIR__, "multiperiod_study.jl"))
function check_rows()
    mc = MULTI["four_ieee9500"]; vs = get(ENV, "MULTI_VSRC", "1.03"); T = 24
    zs = mc.zs; imp = @sprintf("R1=%.8g X1=%.8g R0=%.8g X0=%.8g", real(zs), imag(zs), real(zs), imag(zs))
    extra = ["New Line.poi4 Phases=1 Bus1=4s Bus2=$(mc.add[1]) r1=1e-06 x1=1e-06 c1=0 c0=0",
             "New Vsource.subs4 Phases=1 Bus1=4s basekv=$(mc.kv) pu=$vs angle=0 $imp"]
    append!(extra, ["Line.$t.enabled=true" for t in mc.ties]); append!(extra, ["Edit Vsource.$k pu=$vs" for k in mc.sources[1:3]])
    net, pvs = load_sources(mc.sources; system = mc.system, extra)
    s = setup(net, pvs; T, feas = mc.feas)
    kw = problem(s, :cyclic)
    rows = [split(l, ',') for l in readlines(joinpath(@__DIR__, "..", "results", "multiperiod", "multi_windows_fib78_v1.03.csv"))[2:end]]
    out = joinpath(OUT_DIR, "check_bracket_rows.csv")
    open(io -> println(io, "direction,recorded_deg,fraction,r_deg,start,status,iterations,margin_kW,min_import_kW,max_violation,seconds"), out, "w")
    for d in (47, 65)
        r = only([x for x in rows if x[4] == "cyclic" && x[5] == string(d)])
        u = parse.(Float64, r[6:8]); r0 = parse(Float64, r[9])
        at(f) = Dict("subs1" => 0.0, "subs2" => f * r0 * u[1], "subs3" => f * r0 * u[2], "subs4" => f * r0 * u[3])
        prev = nothing
        runs = [(1.0, :warm), (1.001, :warm), (1.002, :warm), (1.005, :warm), (1.0, :cold), (1.001, :cold)]
        for (fr, how) in runs
            how === :warm && prev === nothing && (prev = solve_multiperiod(net, at(1.0); kw..., objective = :margin, margin_band_kw = 20000.0))
            t0 = time()
            x = solve_multiperiod(net, at(fr); kw..., objective = :margin, margin_band_kw = 20000.0, start = how === :warm ? prev : nothing)
            how === :warm && x.status in RAY_OK && (prev = x)
            mi = minimum(minimum(v) for v in values(x.P_subs_kW))
            line = @sprintf("%d,%.6f,%.4f,%.6f,%s,%s,%s,%.4f,%.4f,%.2e,%.1f", d, r0, fr, fr * r0, how, x.status, string(something(x.iterations, -1)), x.margin_kW, mi, x.max_violation, time() - t0)
            open(io -> println(io, line), out, "a"); say("  ", line)
        end
    end
end
check_rows()
