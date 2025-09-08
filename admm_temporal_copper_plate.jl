############################################################
# ADMM Temporal Decomposition + Brute-Force MPOPF (Copper Plate)
# Notation:
#   B[t] (kWh), P_B[t] (kW), P_Subs[t] (kW)
#   B[t] = B[t-1] - P_B[t]*Δt
#   P_Subs[t] + P_B[t] = P_L[t]
# Bounds:
#   soc_min*E_Rated ≤ B[t] ≤ soc_max*E_Rated
#
# Requires: Julia ≥ 1.9, JuMP, Ipopt
############################################################

using JuMP
using Ipopt
using LinearAlgebra
using Printf

# ----------------------- Data ----------------------------

struct BatteryParams
    E_Rated::Float64      # kWh
    soc_min::Float64      # pu of capacity (e.g., 0.30)
    soc_max::Float64      # pu (e.g., 0.95)
    P_B_R::Float64        # symmetric kW power limit (|P_B| ≤ P_B_R)
    Δt::Float64           # hours per step
    B0::Float64           # initial energy (kWh)
    B_T_target::Union{Nothing,Float64}  # optional terminal energy (kWh)
end

struct Instance
    T::Int
    price::Vector{Float64}   # $/kWh (applied to P_Subs)
    P_L::Vector{Float64}     # net load (kW) at each t
    bat::BatteryParams
end

Bmin(b::BatteryParams) = b.soc_min * b.E_Rated
Bmax(b::BatteryParams) = b.soc_max * b.E_Rated

# ---------------- Block subproblem (ADMM) ----------------

"""
solve_block(inst, tidx, B̂L, B̂R, uL, uR, ρ) -> (B_L*, B_R*, P_Subs, P_B, B, obj)
"""
function solve_block(inst::Instance, tidx::Vector{Int},
    B̂L::Float64, B̂R::Float64,
    uL::Float64, uR::Float64, ρ::Float64)

    b = inst.bat
    Tloc = length(tidx)

    m = Model(Ipopt.Optimizer)
    set_silent(m)

    @variables(m, begin
        -b.P_B_R <= P_B[1:Tloc] <= b.P_B_R     # kW
        P_Subs[1:Tloc]                          # kW
        Bmin(b) <= B[1:Tloc] <= Bmax(b)         # kWh
        Bmin(b) <= B_L <= Bmax(b)               # kWh
        Bmin(b) <= B_R <= Bmax(b)               # kWh
    end)

    # NRPB
    @constraint(m, [τ = 1:Tloc], P_Subs[τ] + P_B[τ] == inst.P_L[tidx[τ]])

    # SOC dynamics: B[t] = B[t-1] - P_B[t]*Δt
    if Tloc == 1
        @constraint(m, B[1] == B_L - P_B[1] * b.Δt)
        @constraint(m, B_R == B[1])
    else
        @constraint(m, B[1] == B_L - P_B[1] * b.Δt)
        @constraint(m, [τ = 1:Tloc-1], B[τ+1] == B[τ] - P_B[τ+1] * b.Δt)
        @constraint(m, B_R == B[Tloc])
    end

    energy_cost = sum(inst.price[tidx[τ]] * P_Subs[τ] * b.Δt for τ = 1:Tloc)
    aug = (ρ / 2) * ((B_L - B̂L + uL)^2 + (B_R - B̂R + uR)^2)
    @objective(m, Min, energy_cost + aug)

    optimize!(m)

    return value(B_L), value(B_R), value.(P_Subs), value.(P_B), value.(B), objective_value(m)
end

# ------------------- ADMM driver -------------------------

"""
admm_solve(inst; K, ρ, max_iter, eps_pri, eps_dual)
Returns Dict(:P_B, :P_Subs, :B, :objective_history, :consensus, :blocks)
"""
function admm_solve(inst::Instance; K::Int=4, ρ::Float64=5.0,
    max_iter::Int=200, eps_pri::Float64=1e-3, eps_dual::Float64=1e-3)

    T = inst.T
    b = inst.bat

    # Partition
    sizes = fill(div(T, K), K)
    for i in 1:rem(T, K)
        sizes[i] += 1
    end
    blocks = Vector{Vector{Int}}()
    s = 1
    for k in 1:K
        e = s + sizes[k] - 1
        push!(blocks, collect(s:e))
        s = e + 1
    end

    # Consensus + duals
    B̂ = zeros(K + 1)
    uL = zeros(K)
    uR = zeros(K)

    B̂[1] = b.B0
    B̂[end] = isnothing(b.B_T_target) ? clamp(b.B0, Bmin(b), Bmax(b)) : b.B_T_target
    for j in 2:K
        λ = (j - 1) / K
        B̂[j] = (1 - λ) * B̂[1] + λ * B̂[end]
    end
    B̂ .= clamp.(B̂, Bmin(b), Bmax(b))

    # storage
    BL = zeros(K)
    BR = zeros(K)
    PBs = [zeros(length(blocks[k])) for k in 1:K]
    PSs = [zeros(length(blocks[k])) for k in 1:K]
    Bs = [zeros(length(blocks[k])) for k in 1:K]
    obj_hist = Float64[]

    @printf "ADMM(K=%d, ρ=%.3f, T=%d) starting...\n" K ρ T
    for it in 1:max_iter
        total = 0.0
        for k in 1:K
            tidx = blocks[k]
            BL[k], BR[k], PSs[k], PBs[k], Bs[k], obj =
                solve_block(inst, tidx, B̂[k], B̂[k+1], uL[k], uR[k], ρ)
            total += obj
        end

        B̂_old = copy(B̂)
        B̂[1] = b.B0
        B̂[end] = isnothing(b.B_T_target) ? B̂[end] : b.B_T_target
        for j in 2:K
            kL = j - 1
            kR = j
            B̂[j] = 0.5 * ((BR[kL] + uR[kL]) + (BL[kR] + uL[kR]))
        end
        B̂ .= clamp.(B̂, Bmin(b), Bmax(b))

        for k in 1:K
            uL[k] += BL[k] - B̂[k]
            uR[k] += BR[k] - B̂[k+1]
        end

        r = vcat([BL[k] - B̂[k] for k in 1:K], [BR[k] - B̂[k+1] for k in 1:K])
        r_norm = norm(r)
        s_norm = ρ * norm(B̂ - B̂_old)

        push!(obj_hist, total)
        @printf "it=%3d  obj=%10.4f  ‖r‖=%.2e  ‖s‖=%.2e\n" it total r_norm s_norm
        if r_norm ≤ eps_pri && s_norm ≤ eps_dual
            @printf "Converged at iteration %d\n" it
            break
        end
    end

    P_B = reduce(vcat, PBs)
    P_Subs = reduce(vcat, PSs)
    B_traj = reduce(vcat, Bs)

    return Dict(
        :P_B => P_B,
        :P_Subs => P_Subs,
        :B => B_traj,
        :objective_history => obj_hist,
        :consensus => B̂,
        :blocks => blocks
    )
end

# ---------------- Brute-force MPOPF (monolithic) ---------

"""
brute_force_solve(inst) -> Dict(:P_B, :P_Subs, :B, :objective)

Solves the *entire* horizon in one JuMP model:
  min  Σ_t price[t]*P_Subs[t]*Δt
  s.t. P_Subs[t] + P_B[t] = P_L[t]
       B[t] = B[t-1] - P_B[t]*Δt
       soc_min*E_Rated ≤ B[t] ≤ soc_max*E_Rated
       |P_B[t]| ≤ P_B_R
       B[0] = B0, (optional) B[T] = B_T_target
"""
function brute_force_solve(inst::Instance)
    T = inst.T
    b = inst.bat

    m = Model(Ipopt.Optimizer)
    set_silent(m)

    @variables(m, begin
        -b.P_B_R <= P_B[1:T] <= b.P_B_R
        P_Subs[1:T]
        Bmin(b) <= B[1:T] <= Bmax(b)
    end)

    # NRPB
    @constraint(m, [t = 1:T], P_Subs[t] + P_B[t] == inst.P_L[t])

    # SOC dynamics with fixed initial energy
    @constraint(m, B[1] == b.B0 - P_B[1] * b.Δt)
    @constraint(m, [t = 2:T], B[t] == B[t-1] - P_B[t] * b.Δt)

    # Optional terminal energy
    if !isnothing(b.B_T_target)
        @constraint(m, B[T] == b.B_T_target)
    end

    @objective(m, Min, sum(inst.price[t] * P_Subs[t] * b.Δt for t in 1:T))

    optimize!(m)

    return Dict(
        :P_B => value.(P_B),
        :P_Subs => value.(P_Subs),
        :B => value.(B),
        :objective => objective_value(m)
    )
end

# ----------------------- Example -------------------------

# if abspath(PROGRAM_FILE) == @__FILE__
    # Example 24h problem
    T = 24
    price = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2
    P_L = 1000 .+ 250 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2

    bat = BatteryParams(
        4000.0,    # E_Rated (kWh)
        0.30,      # soc_min
        0.95,      # soc_max
        800.0,     # P_B_R (kW)
        1.0,       # Δt (h)
        0.60 * 4000, # B0 (kWh)
        nothing    # B_T_target
    )
    inst = Instance(T, price, P_L, bat)

    # Brute-force baseline
    mono = brute_force_solve(inst)
    @printf "\nBrute-force objective: %.4f\n" mono[:objective]

    # ADMM (decomposed)
    sol = admm_solve(inst; K=4, ρ=5.0, max_iter=200, eps_pri=1e-4, eps_dual=1e-4)
    @printf "ADMM objective (last iterate): %.4f\n" last(sol[:objective_history])

    # Quick consistency checks
    @printf "BF  B[min,max]=[%.1f, %.1f] | PB[min,max]=[%.1f, %.1f]\n" minimum(mono[:B]) maximum(mono[:B]) minimum(mono[:P_B]) maximum(mono[:P_B])
    @printf "ADM B[min,max]=[%.1f, %.1f] | PB[min,max]=[%.1f, %.1f]\n" minimum(sol[:B]) maximum(sol[:B]) minimum(sol[:P_B]) maximum(sol[:P_B])
# end
