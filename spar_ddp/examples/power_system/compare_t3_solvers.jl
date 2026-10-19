# Identical adjusted T=3 instance solved by centralized Ipopt and FilterDDP.

include("convergence_trace.jl")
include("mp_optimality.jl")

cost, demand = tadmm_tutorial3()
central = solve_explicit(3; bounded=true, cost=cost, pL=demand)
filtered = solve_ddp(3; ps_lo=0.0, pb_lo=PB_MIN_KKT, pb_hi=PB_MAX_KKT,
                     B_lo=B_MIN_KKT, B_hi=B_MAX_KKT,
                     c=cost, pL=demand, tol=1e-10, verbose=false)

J_filter = sum(cost .* filtered.psub .* DT_KKT) +
           C_B_KKT * sum(filtered.pb .^ 2) * DT_KKT +
           W_KKT * (filtered.B[4] - B_INIT_KKT)^2
gap = abs(J_filter - central.objective)
pb_error = maximum(abs.(filtered.pb .- central.P_B))
B_error = maximum(abs.(filtered.B .- central.B))
data = filtered.solver.data

println("centralized Ipopt objective = ", central.objective)
println("FilterDDP objective         = ", J_filter)
println("objective gap              = ", gap)
println("maximum P_B error          = ", pb_error)
println("maximum B error            = ", B_error)
println("FilterDDP iterations/status= ", data.k, "/", data.status)

# Actual FilterDDP objective trace for the three-method appendix plot.
trace = trace_case(solve_ddp, 3; ps_lo=0.0, pb_lo=PB_MIN_KKT, pb_hi=PB_MAX_KKT,
                   B_lo=B_MIN_KKT, B_hi=B_MAX_KKT,
                   c=cost, pL=demand, tol=1e-10)
figdir = normpath(joinpath(@__DIR__, "..", "..", "paper", "figures"))
open(joinpath(figdir, "user_ddp_t3_filter.csv"), "w") do io
    println(io, "iteration,objective")
    for row in trace
        row.k == 0 && continue  # k=0 is FilterDDP's unsolved initial guess
        println(io, "$(row.k),$(row.obj)")
    end
end
println("wrote FilterDDP objective trace with ", length(trace), " rows")

# Correct battery-energy trajectories for the appendix. Both arrays include
# B^0, so the final row is B^3 and can be checked directly against B_init.
open(joinpath(figdir, "user_ddp_t3_optimal_battery.csv"), "w") do io
    println(io, "t,centralized_kwh,filterddp_kwh,b_init_kwh,b_max_kwh")
    for i in eachindex(central.B)
        println(io, "$(i - 1),$(1000central.B[i]),$(1000filtered.B[i]),$(1000B_INIT_KKT),$(1000B_MAX_KKT)")
    end
end
println("wrote centralized and FilterDDP battery trajectories")

# Battery-power intervals in the same physical-action convention as Fig. 3:
# charging is positive/upward and discharging is negative/downward.
open(joinpath(figdir, "user_ddp_t3_optimal_battery_power.csv"), "w") do io
    println(io, "tc,charge_kw,discharge_kw,filter_action_kw")
    for t in eachindex(central.P_B)
        central_action = -1000central.P_B[t]
        filter_action = -1000filtered.pb[t]
        println(io, "$(t - 0.5),$(max(central_action, 0.0)),$(min(central_action, 0.0)),$filter_action")
    end
end
println("wrote centralized and FilterDDP battery-power actions")

gap < 1e-7 || error("FilterDDP and centralized Ipopt do not agree")
pb_error < 1e-6 || error("FilterDDP dispatch does not match centralized Ipopt")
