# Verify the four centralized MPOPF stationarity equations used in the appendix.
# Convention: lambda denotes equality multipliers and mu denotes inequality
# multipliers. B^t is the energy at the end of period t and B^0 = B_init.

using JuMP
using Ipopt
using Printf

include("tadmm_profiles.jl")

const C_B_KKT = 0.05
const DT_KKT = 1.0
const B_INIT_KKT = 2.0
const W_KKT = 5.0
const PB_MIN_KKT = -0.45
const PB_MAX_KKT = 0.45
const B_MIN_KKT = 1.20
const B_MAX_KKT = 1.95

"""Solve the centralized problem with explicit constraints and return its duals."""
function solve_explicit(T::Int; bounded::Bool = true, tol = 1e-10,
                        cost = tadmm_cost(T), pL = tadmm_pL(T))
    length(cost) == T == length(pL) || error("profile lengths must equal T")
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "constr_viol_tol", tol)

    @variable(model, P_subs[1:T])
    @variable(model, P_B[1:T])
    @variable(model, B[1:T+1]) # B[t+1] stores B^t; B[1] is B^0.
    @objective(model, Min,
        sum((cost[t] * P_subs[t] + C_B_KKT * P_B[t]^2) * DT_KKT for t in 1:T) +
        W_KKT * (B[T+1] - B_INIT_KKT)^2)

    @constraint(model, soc_init, B[1] == B_INIT_KKT)
    @constraint(model, balance[t in 1:T], P_subs[t] + P_B[t] - pL[t] == 0)
    @constraint(model, soc[t in 1:T], B[t+1] - B[t] + DT_KKT * P_B[t] == 0)
    @constraint(model, subs_min[t in 1:T], -P_subs[t] <= 0)
    if bounded
        @constraint(model, pb_min[t in 1:T], PB_MIN_KKT - P_B[t] <= 0)
        @constraint(model, pb_max[t in 1:T], P_B[t] - PB_MAX_KKT <= 0)
        @constraint(model, b_min[t in 1:T], B_MIN_KKT - B[t+1] <= 0)
        @constraint(model, b_max[t in 1:T], B[t+1] - B_MAX_KKT <= 0)
    end

    optimize!(model)
    termination_status(model) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL) ||
        error("centralized solve failed at T=$T: $(termination_status(model))")

    zeros_T = zeros(T)
    return (
        T=T, cost=cost, pL=pL, P_subs=value.(P_subs), P_B=value.(P_B),
        B=value.(B), objective=objective_value(model),
        lambda_bal=[-dual(balance[t]) for t in 1:T],
        lambda_soc=[-dual(soc[t]) for t in 1:T],
        mu_subs=[-dual(subs_min[t]) for t in 1:T],
        mu_Pmin=bounded ? [-dual(pb_min[t]) for t in 1:T] : zeros_T,
        mu_Pmax=bounded ? [-dual(pb_max[t]) for t in 1:T] : zeros_T,
        mu_Bmin=bounded ? [-dual(b_min[t]) for t in 1:T] : zeros_T,
        mu_Bmax=bounded ? [-dual(b_max[t]) for t in 1:T] : zeros_T,
    )
end

kkt1(r, t) = r.cost[t] * DT_KKT + r.lambda_bal[t] - r.mu_subs[t]
kkt2(r, t) = 2C_B_KKT * r.P_B[t] * DT_KKT + r.lambda_bal[t] +
             DT_KKT * r.lambda_soc[t] - r.mu_Pmin[t] + r.mu_Pmax[t]
kkt3(r, t) = r.lambda_soc[t] - r.lambda_soc[t+1] -
             r.mu_Bmin[t] + r.mu_Bmax[t]
kkt4(r) = 2W_KKT * (r.B[r.T+1] - B_INIT_KKT) + r.lambda_soc[r.T] -
          r.mu_Bmin[r.T] + r.mu_Bmax[r.T]

function report(r; label)
    T = r.T
    residuals = (
        maximum(abs(kkt1(r, t)) for t in 1:T),
        maximum(abs(kkt2(r, t)) for t in 1:T),
        T > 1 ? maximum(abs(kkt3(r, t)) for t in 1:T-1) : 0.0,
        abs(kkt4(r)),
    )
    println("\n", label, ", T = ", T)
    @printf("objective = %.12f; price range = %.6e\n",
            r.objective, maximum(r.cost) - minimum(r.cost))
    @printf("K1 dL/dP_Subs: %.3e\n", residuals[1])
    @printf("K2 dL/dP_B:    %.3e\n", residuals[2])
    @printf("K3 dL/dB^t:    %.3e\n", residuals[3])
    @printf("K4 dL/dB^T:    %.3e\n", residuals[4])
    T > 1 && println("successor lambda_soc terms: ", r.lambda_soc[2:end])
    worst = maximum(residuals)
    worst < 1e-6 || error("KKT verification failed at T=$T (worst=$worst)")
    return worst
end

function main()
    println("lambda = equality multipliers; mu = inequality multipliers")
    worst = 0.0
    for T in (3, 6)
        worst = max(worst, report(solve_explicit(T; bounded=true), label="all bounds"))
        worst = max(worst, report(solve_explicit(T; bounded=false), label="P_Subs floor only"))
    end
    @printf("\nAll four equations verified; worst residual %.3e\n", worst)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
