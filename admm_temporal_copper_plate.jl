############################################################
# Copper-Plate MPOPF in per-unit
# - ADMM with exactly T single-step blocks (no K)
# - Brute-force (monolithic) baseline
#
# Bases:
#   kV_B = 4.16/sqrt(3)  (phase-to-neutral; not used in power-only)
#   kVA_B = 1000  => P_BASE = 1000 kW, E_BASE = 1000 kWh (per 1h)
#
# Variables (optimization; all per-unit):
#   P_Subs[t], P_B[t] on P_BASE;  B[t] on E_BASE  (puh)
# Constraints:
#   B[t] = B[t-1] - P_B[t]*Δt
#   P_Subs[t] + P_B[t] = P_L[t]
#   soc_min*E_Rated_pu ≤ B[t] ≤ soc_max*E_Rated_pu
# Objective ($):  sum_t price[t] * (P_Subs_pu[t]*P_BASE) * Δt
#
# Requires: Julia ≥ 1.9, JuMP, Ipopt
############################################################

using JuMP
using Ipopt
using LinearAlgebra
using Printf

# ----------------------- Bases ---------------------------
const KV_B = 4.16 / sqrt(3)   # kV (unused here)
const KVA_B = 1000.0           # kVA
const P_BASE = KVA_B            # kW
const E_BASE = P_BASE * 1.0     # kWh per 1 hour

# ----------------------- Data ----------------------------

struct BatteryParamsPU
    E_Rated_pu::Float64          # puh on E_BASE
    soc_min::Float64             # pu
    soc_max::Float64             # pu
    P_B_R_pu::Float64            # pu on P_BASE
    Δt::Float64                  # hours
    B0_pu::Float64               # initial energy (puh)
    B_T_target_pu::Union{Nothing,Float64}  # terminal (puh), optional
end

struct InstancePU
    T::Int
    price::Vector{Float64}       # $/kWh
    P_L_pu::Vector{Float64}      # pu on P_BASE
    bat::BatteryParamsPU
end

Bmin(b::BatteryParamsPU) = b.soc_min * b.E_Rated_pu
Bmax(b::BatteryParamsPU) = b.soc_max * b.E_Rated_pu

"""
build_instance_pu(...)
Converts physical inputs (kW/kWh) to per-unit given the fixed bases.
"""
function build_instance_pu(T::Int, price::Vector{Float64}, P_L_kW::Vector{Float64},
    E_Rated_kWh::Float64, soc_min::Float64, soc_max::Float64,
    P_B_R_kW::Float64, Δt::Float64, B0_kWh::Float64;
    B_T_target_kWh::Union{Nothing,Float64}=nothing)
    @assert length(price) == T && length(P_L_kW) == T
    P_L_pu = P_L_kW ./ P_BASE
    E_Rated_pu = E_Rated_kWh / E_BASE
    P_B_R_pu = P_B_R_kW / P_BASE
    B0_pu = B0_kWh / E_BASE
    B_T_target_pu = isnothing(B_T_target_kWh) ? nothing : (B_T_target_kWh / E_BASE)
    bat = BatteryParamsPU(E_Rated_pu, soc_min, soc_max, P_B_R_pu, Δt, B0_pu, B_T_target_pu)
    return InstancePU(T, price, P_L_pu, bat)
end

# ---------------- Single-step block solve ----------------
# For a single time step t (i.e., Tloc = 1).
# Decision vars: P_B, P_Subs, B_t, B_L, B_R (scalars).
# Todo Check for SOC Trajectory constraints again (Name the variables as per agreed upon formulation)
function solve_single_step_block(inst::InstancePU, t::Int,
    Bhat_L::Float64, Bhat_R::Float64,
    uL::Float64, uR::Float64, ρ::Float64)
    b = inst.bat
    m = Model(Ipopt.Optimizer)
    set_silent(m)

    @variables(m, begin
        -b.P_B_R_pu <= P_B <= b.P_B_R_pu         # pu
        P_Subs                                     # pu (can be ±)
        Bmin(b) <= B_t <= Bmax(b)                  # puh at time t
        Bmin(b) <= B_L <= Bmax(b)                  # left boundary puh (B at t-1)
        Bmin(b) <= B_R <= Bmax(b)                  # right boundary puh (B at t)
    end)

    # NRPB (pu)
    @constraint(m, P_Subs + P_B == inst.P_L_pu[t])

    # SOC for single step: B_t = B_L - P_B*Δt, and B_R = B_t
    @constraint(m, B_t == B_L - P_B * b.Δt)
    @constraint(m, B_R == B_t)

    # Economic cost + ADMM penalties on boundaries
    energy_cost = inst.price[t] * (P_Subs * P_BASE) * b.Δt
    aug = (ρ / 2) * ((B_L - Bhat_L + uL)^2 + (B_R - Bhat_R + uR)^2)
    @objective(m, Min, energy_cost + aug)

    optimize!(m)

    return value(B_L), value(B_R), value(P_Subs), value(P_B), value(B_t), objective_value(m)
end

# ------------------- ADMM (T single-step blocks) ---------
"""
admm_solve_Tblocks(inst; ρ=5.0, max_iter=200, eps_pri=1e-3, eps_dual=1e-3)

- Creates exactly T single-step blocks (one per time index).
- Consensus vector Bhat has length T+1 (boundaries 0..T).
- Duals uL,uR have length T (one (L,R) pair per block).

Returns Dict(:P_B, :P_Subs, :B, :objective_history, :consensus_boundaries)
"""
function admm_solve_Tblocks(inst::InstancePU; ρ::Float64=5.0,
    max_iter::Int=200, eps_pri::Float64=1e-3, eps_dual::Float64=1e-3)
    T = inst.T
    b = inst.bat

    # Consensus boundaries Bhat[1..T+1] correspond to B(0), B(1), ..., B(T)
    Bhat = zeros(T + 1)
    Bhat[1] = b.B0_pu
    Bhat[end] = isnothing(b.B_T_target_pu) ? clamp(b.B0_pu, Bmin(b), Bmax(b)) : b.B_T_target_pu
    # initialize internal boundaries by interpolation
    for j in 2:T
        λ = (j - 1) / T
        Bhat[j] = (1 - λ) * Bhat[1] + λ * Bhat[end]
    end
    Bhat .= clamp.(Bhat, Bmin(b), Bmax(b))

    # Duals for each block’s (left,right) boundaries
    uL = zeros(T)
    uR = zeros(T)

    # Storage
    BL = zeros(T)
    BR = zeros(T)       # local boundary solutions
    PB = zeros(T)
    PS = zeros(T)
    B = zeros(T)
    obj_hist = Float64[]

    @printf "ADMM[T-blocks]: T=%d, ρ=%.3f\n" T ρ
    for it in 1:max_iter
        total = 0.0

        # 1) Solve each single-step block (can be parallelized if desired)
        for t in 1:T
            BL[t], BR[t], PS[t], PB[t], B[t], obj =
                solve_single_step_block(inst, t, Bhat[t], Bhat[t+1], uL[t], uR[t], ρ)
            total += obj
        end

        # 2) Consensus update on internal boundaries (pin the two ends)
        Bhat_old = copy(Bhat)
        Bhat[1] = b.B0_pu
        Bhat[end] = isnothing(b.B_T_target_pu) ? Bhat[end] : b.B_T_target_pu
        for j in 2:T  # internal interfaces between blocks (1..T-1), mapped to indices 2..T
            kL = j - 1   # left block index
            kR = j     # right block index
            # average of the two locals plus duals
            # Bhat[j] = 0.5 * ((BR[kL] + uR[kL]) + (BL[kR] + uL[kR]))
            Bhat[j] = (BR[kL] + uR[kL])
        end
        Bhat .= clamp.(Bhat, Bmin(b), Bmax(b))

        # 3) Dual updates
        for t in 1:T
            uL[t] += BL[t] - Bhat[t]
            uR[t] += BR[t] - Bhat[t+1]
        end

        # 4) Residuals and stopping
        r = vcat([BL[t] - Bhat[t] for t in 1:T], [BR[t] - Bhat[t+1] for t in 1:T])
        r_norm = norm(r)
        s_norm = ρ * norm(Bhat - Bhat_old)

        push!(obj_hist, total)
        @printf "it=%3d  obj=%10.4f  ‖r‖=%.2e  ‖s‖=%.2e\n" it total r_norm s_norm
        if r_norm ≤ eps_pri && s_norm ≤ eps_dual
            @printf "Converged at iteration %d\n" it
            break
        end
    end

    return Dict(
        :P_B => PB,
        :P_Subs => PS,
        :B => B,
        :objective_history => obj_hist,
        :consensus_boundaries => Bhat,  # length T+1 (B(0)..B(T))
    )
end

# ---------------- Brute-force (monolithic, pu) -----------

function brute_force_solve(inst::InstancePU)
    T = inst.T
    b = inst.bat
    m = Model(Ipopt.Optimizer)
    set_silent(m)

    @variables(m, begin
        -b.P_B_R_pu <= P_B[1:T] <= b.P_B_R_pu
        P_Subs[1:T]
        Bmin(b) <= B[1:T] <= Bmax(b)
    end)

    @constraint(m, [t = 1:T], P_Subs[t] + P_B[t] == inst.P_L_pu[t])
    @constraint(m, B[1] == b.B0_pu - P_B[1] * b.Δt) # Todo: Check for left SOC in ADMM t=1 problem
    @constraint(m, [t = 2:T], B[t] == B[t-1] - P_B[t] * b.Δt)
    if !isnothing(b.B_T_target_pu)
        @constraint(m, B[T] == b.B_T_target_pu)
    end

    @objective(m, Min, sum(inst.price[t] * (P_Subs[t] * P_BASE) * b.Δt for t in 1:T))
    optimize!(m)

    return Dict(:P_B => value.(P_B), :P_Subs => value.(P_Subs), :B => value.(B),
        :objective => objective_value(m))
end

# ----------------------- Example -------------------------

# if abspath(PROGRAM_FILE) == @__FILE__
    # Example inputs (physical)
    T = 24
    price_dollars_per_kWh = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2
    P_L_kW = 1000 .+ 250 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2

    E_Rated_kWh = 4000.0
    soc_min, soc_max = 0.30, 0.95
    P_B_R_kW = 800.0
    Δt_h = 1.0
    B0_kWh = 0.60 * E_Rated_kWh
    # B_T_target_kWh = 0.60 * E_Rated_kWh

    inst = build_instance_pu(T, price_dollars_per_kWh, P_L_kW,
        E_Rated_kWh, soc_min, soc_max,
        P_B_R_kW, Δt_h, B0_kWh;
        B_T_target_kWh=nothing)

    # Brute-force baseline
    mono = brute_force_solve(inst)
    @printf "\nBrute-force objective: %.4f\n" mono[:objective]

    # ADMM with exactly T single-step blocks
    sol = admm_solve_Tblocks(inst; ρ=5.0, max_iter=25, eps_pri=1e-3, eps_dual=1e-3)
    @printf "ADMM objective (last iterate): %.4f\n" last(sol[:objective_history])

    # Quick pu checks
    @printf "BF  B[min,max]=[%.4f, %.4f] pu  | PB[min,max]=[%.4f, %.4f] pu\n" minimum(mono[:B]) maximum(mono[:B]) minimum(mono[:P_B]) maximum(mono[:P_B])
    @printf "ADM B[min,max]=[%.4f, %.4f] pu  | PB[min,max]=[%.4f, %.4f] pu\n" minimum(sol[:B]) maximum(sol[:B]) minimum(sol[:P_B]) maximum(sol[:P_B])

    # Convert back to physical for sanity/plots if desired
    PB_kW = sol[:P_B] .* P_BASE
    PS_kW = sol[:P_Subs] .* P_BASE
    B_kWh = sol[:B] .* E_BASE
    @printf "(ADMM) PB[kW] in [%.1f, %.1f], B[kWh] in [%.1f, %.1f]\n" minimum(PB_kW) maximum(PB_kW) minimum(B_kWh) maximum(B_kWh)
# end
