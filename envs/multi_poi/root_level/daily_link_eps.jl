# Is the link's error allowance what limits the bound? (2026-10-22). DIAGNOSTIC where marked.
#   julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_link_eps.jl
# At theta 157.5, 1.02 r*: the dispatch-link relaxation with eps scaled by 2 (still valid: a larger allowance),
# 0.1 and 0 (NOT valid: they would exclude real states whose voltages differ from the no-load voltages). If the bound
# barely moves, the allowance is not the limit and the remaining freedom is elsewhere (the second moments).
ENV["DL_LIB"] = "1"
include(joinpath(@__DIR__, "daily_link.jl"))
function eps_main()
    s24 = setup(net, pvs; T = 24, pv = "qfree"); kw24 = problem(s24, :cyclic); T = kw24.T
    q = [split(l, ',') for l in readlines(joinpath(OUT_DIR, "daily_pilot_queries.csv"))[2:end]]
    RM = response_matrix()
    row = only(r for r in q if r[1] == "157.5" && r[2] == "1.02")
    th = 157.5; r = parse(Float64, row[3]); ang = at(r, th); rub = parse(Float64, row[12])
    CBs = [current_bounds(RM, ang, device_rectangles(net, kw24, t), t; ndir = 3600, ndir_rho = 720, corners = true) for t in 1:T]
    bnd = Dict(s.name => [CBs[t].out[s.name].sect for t in 1:T] for s in net.sources)
    open(io -> println(io, "theta_deg,point,eps_scale,valid,relax_status,rho_relax_kW,rho_lb_kW,fictitious_losses_kWh,max_source_moment_gap_pu2,solve_s"), DL * "_eps_scan.csv", "w")
    for (sc, valid) in ((2.0, "yes"), (1.0, "yes"), (0.1, "NO (diagnostic)"), (0.0, "NO (diagnostic)"))
        o = run_relax(ang, kw24, bnd, rub, CBs; link = true, eps_scale = sc)
        say(@sprintf("  eps x %.1f (valid: %s): %s, rho relax %.4f kW, checked lower bound %.4f kW, fictitious losses %.1f kWh, source moment gap %.3f pu^2, solve %.1f s",
                     sc, valid, o.st, o.rho, o.lb, o.fict, o.srcgap, o.ts))
        dl_csv("_eps_scan.csv", [th, "1.02", sc, valid, o.st, @sprintf("%.4f", o.rho), @sprintf("%.4f", o.lb), @sprintf("%.2f", o.fict), @sprintf("%.3e", o.srcgap), @sprintf("%.1f", o.ts)])
    end
end
eps_main()
