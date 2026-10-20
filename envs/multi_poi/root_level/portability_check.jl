# Portability check for a new machine (Kamiak, 2026-10-19): recompute two committed IEEE123
# results with the stage-main settings and compare.
#   - the T = 24 windows of delta_s3 - delta_s4 with PV reactive power free, batteries idle and
#     cyclic (main_horizon_windows.csv, rows "qfree,nobatt" and "qfree,cyclic");
#   - the cost-optimal cyclic day at equal angles, hour by hour (main_days.csv, ddelta 0.00).
# Run from the repo root:
#   julia --project=envs/tadmm envs/multi_poi/root_level/portability_check.jl
# The solver follows IPOPT_EXTRA_OPTIONS, so the same script checks MUMPS and HSL builds.
empty!(ARGS); push!(ARGS, "none")
include(joinpath(@__DIR__, "multiperiod_study.jl"))

ref_dir = normpath(joinpath(@__DIR__, "..", "results", "multiperiod"))
rows(f) = [split(l, ',') for l in readlines(joinpath(ref_dir, f))[2:end]]
A, B = MAIN_PAIR
net, pvs = load_sources([A, B])
s = setup(net, pvs; T = 24, pv = "qfree")
say("PORTABILITY  ieee123 ", A, " + ", B, ", T = 24, PV q free; IPOPT_EXTRA_OPTIONS = \"", get(ENV, "IPOPT_EXTRA_OPTIONS", ""), "\"")
worst = 0.0
for case in (:nobatt, :cyclic)
    t0 = time()
    w = window(s, case, A; check = :both)
    r = only(filter(r -> r[1] == "qfree" && r[2] == string(case), rows("main_horizon_windows.csv")))
    lo, hi = parse(Float64, r[4]), parse(Float64, r[5])
    d = max(abs(w.lo - lo), abs(w.hi - hi))
    global worst = max(worst, d)
    say(@sprintf("  %-6s window [%+.6f, %+.6f] deg, committed [%+.6f, %+.6f], max difference %.2e deg, failed checks %d, %.0f s",
                 case, w.lo, w.hi, lo, hi, d, w.bad, time() - t0))
end
t0 = time()
x = solve_case(s, equal_angles(net), :cyclic)
ref = filter(r -> r[1] == "0.00", rows("main_days.csv"))
dP = maximum(max(abs(x.P_subs_kW[A][t] - parse(Float64, ref[t][6])), abs(x.P_subs_kW[B][t] - parse(Float64, ref[t][7]))) for t in 1:24)
say(@sprintf("  cost-optimal day at equal angles: %s, \$%.4f; substation imports against main_days.csv: max |dP| %.2e kW, %.0f s",
             x.status, x.cost.total, dP, time() - t0))
say(@sprintf("PORTABILITY RESULT: windows within %.2e deg, imports within %.2e kW", worst, dP))
