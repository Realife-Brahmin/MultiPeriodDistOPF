# Where a relaxation beyond the edge finds its slack (2026-11-03): at FB_RATIO x the limiting hour's edge of direction
# FB_DIR (phase A's file), solve phase C's relaxation and list the lines with the largest fictitious loss
# R (|I|^2 moment - |first moment|^2), their relaxed currents and ratings, and the largest |V|^2 gaps.
#   FB_SYSTEM=ieee123_s13 FB_RATINGS=1 FB_DIR=1 FB_RATIO=1.10 MULTIPERIOD_OUT_DIR=<phase A folder> \
#     julia --project=<sdp env> problems/msopf/root_level/relaxation_slack.jl
ENV["FB_PHASE"] = "X"
include(joinpath(@__DIR__, "farkas_bracketing.jl"))
let i = parse(Int, ENV["FB_DIR"]), d = DIRS[i], ratio = parse(Float64, get(ENV, "FB_RATIO", "1.10"))
    a = split(readlines(joinpath(OUT, "fb_A_$(d.id).csv"))[2], ','); tstar = parse(Int, a[5]); rs = parse(Float64, a[6])
    kwL = merge(hour_kw(tstar), (; source_ampacity_pu = SRCAMP))
    x = certify(net, angles(d, ratio * rs), kwL; k = KV, cuts = CUTS)
    S = net.S_base_kVA; M = x.R.model
    g(sy, n) = value(variable_by_name(M, string(sy) * "[" * n * ",1]"))
    rows = []
    for l in net.lines
        Rl = real(-net.Y_base / l.Y[1, 2]); m2 = g(:lrlr, l.name) + g(:lili, l.name); m1 = g(:lr, l.name)^2 + g(:li, l.name)^2
        push!(rows, (Rl * (m2 - m1) * S, l.name, sqrt(max(m2, 0)), sqrt(max(m1, 0)), get(kwL.ampacity_pu, l.name, NaN)))
    end
    sort!(rows; rev = true)
    @printf("%s %.3f x r* = %.4f deg, hour %d: %s, U %+.3f kW, fictitious loss %.3f kW total\n", d.label, ratio, ratio * rs, tstar, x.status, x.U, sum(r[1] for r in rows))
    println("  line      loss kW   sqrt(moment) pu   |first moment| pu   rating pu")
    for r in rows[1:12]; @printf("  %-8s %9.3f %12.4f %18.4f %14.4f\n", r[2], r[1], r[3], r[4], r[5]); end
    @printf("  sources: %s\n", join(["$(s.name) bound $(round(SRCAMP[s.name]; digits = 3)) pu" for s in net.sources], ", "))
end
