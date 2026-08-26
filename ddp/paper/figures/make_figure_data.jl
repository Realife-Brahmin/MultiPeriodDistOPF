# Emits figures/balance.csv for balance_fig.tex.
#
# The power balance eq:cp_balance, P_Subs^t + P_B^t - p_L^t = 0, can be read as
# two sides that must match once the battery is split by sign:
#
#   supply side   P_Subs^t + P_d^t        P_d = max(P_B, 0)   discharging
#   demand side   p_L^t    + P_c^t        P_c = max(-P_B, 0)  charging
#
# Those two sums are equal at every t by construction, which is exactly what the
# figure shows: two stacked bars per interval that reach the same height. The
# split into segments then reads directly as each asset's share of the kW mix.
#
# Data is NOT hardcoded here. It is read from the verified centralized reference
# (results/copper_plate/centralized_reference.csv, case 6e_T6), so the figure
# cannot drift from the solution the experiment actually certifies -- which is
# how the previous hand-copied schedule figure went stale.
#
# Run:  julia --startup-file=no ddp/paper/figures/make_balance.jl

const CASE  = "6e_T6"
const KVA_B = 1000.0            # 1 p.u. = 1000 kW, 1 p.u.h = 1000 kWh
const CSV   = joinpath(@__DIR__, "..", "..", "results", "copper_plate",
                       "centralized_reference.csv")

# Analytic benchmark profiles, needed because the reference CSV stores only
# P_Subs, P_B, and B.
benchmark_load(T) = 0.8 .+ 0.2 .* (sin.(range(0, 2pi, length = T) .- 0.8) .+ 1) ./ 2
benchmark_cost(T) = 0.08 .+ 0.12 .* (sin.(range(0, 2pi, length = T)) .+ 1) ./ 2
benchmark_pL(T)   = 2.0 .* benchmark_load(T)

isfile(CSV) || error("""
    $CSV not found. Generate it first:
      julia --startup-file=no --project=envs/ddp2026 \\
            ddp/examples/power_system/copper_plate_centralized.jl""")

psub = Dict{Int,Float64}(); pb = Dict{Int,Float64}(); B = Dict{Int,Float64}()
T = 0
for (i, line) in enumerate(eachline(CSV))
    i == 1 && continue
    f = split(strip(line), ',')
    (length(f) == 7 && f[1] == CASE) || continue
    global T = parse(Int, f[2])
    t = parse(Int, f[3])
    t == -1 && error("case $CASE is marked infeasible in the reference")
    f[4] == "NaN" || (psub[t] = parse(Float64, f[4]))
    f[5] == "NaN" || (pb[t]   = parse(Float64, f[5]))
    f[6] == "NaN" || (B[t]    = parse(Float64, f[6]))
end
T > 0 || error("case $CASE not found in $CSV")

pL = benchmark_pL(T)
c  = benchmark_cost(T)

# Verify the balance before drawing it -- the figure's whole claim is that the
# two stacks match, so assert it rather than trusting the plot to look right.
for t = 1:T
    r = psub[t] + pb[t] - pL[t]
    abs(r) < 1e-7 || error("balance residual $r at t=$t exceeds 1e-7")
end
println("balance verified at all $T intervals (max residual < 1e-7)")

out = joinpath(@__DIR__, "balance.csv")
open(out, "w") do io
    # tsup / tdem are the paired bar centres; the two stacks sit either side of t.
    println(io, "t,tsup,tdem,psub,pd,pl,pc,total")
    for t = 1:T
        pd = max(pb[t], 0.0)                 # discharging -> supply side
        pc = max(-pb[t], 0.0)                # charging    -> demand side
        tot = (psub[t] + pd) * KVA_B
        println(io, t, ",", t - 0.22, ",", t + 0.22, ",",
                psub[t] * KVA_B, ",", pd * KVA_B, ",",
                pL[t] * KVA_B, ",", pc * KVA_B, ",", tot)
    end
end
println("wrote ", out)

# ---------------------------------------------------------------------------
# schedule_interval.csv / schedule_soc.csv for schedule_fig.tex
# ---------------------------------------------------------------------------
# Everything in kW and kWh, never per unit: battery quantities are read as
# ratings, and p.u. hides how large the device actually is.
#
# Plotting convention: charging is drawn POSITIVE and
# discharging NEGATIVE, so a bar's direction is the device's action rather than
# the sign of P_B. Interval quantities sit at the interval centre t-1/2; the
# stored energy sits at the nodes 0..T, which is where B actually lives.

outi = joinpath(@__DIR__, "schedule_interval.csv")
open(outi, "w") do io
    println(io, "tc,psub,pl,cost,pc,pdneg")
    for t = 1:T
        pc    =  max(-pb[t], 0.0) * KVA_B          # charging, drawn positive
        pdneg = -max( pb[t], 0.0) * KVA_B          # discharging, drawn negative
        cost  = psub[t] * KVA_B * c[t]             # $ in this period (Δt = 1 h)
        println(io, t - 0.5, ",", psub[t] * KVA_B, ",", pL[t] * KVA_B, ",",
                cost, ",", pc, ",", pdneg)
    end
end
println("wrote ", outi)

outs = joinpath(@__DIR__, "schedule_soc.csv")
open(outs, "w") do io
    println(io, "node,B")
    for t = 1:(T+1)
        println(io, t - 1, ",", B[t] * KVA_B)      # node 0 is the fixed B_0
    end
end
println("wrote ", outs)

# Numbers for the caption.
tot_sub = sum(psub[t] for t = 1:T) * KVA_B
tot_dis = sum(max(pb[t], 0.0) for t = 1:T) * KVA_B
tot_chg = sum(max(-pb[t], 0.0) for t = 1:T) * KVA_B
tot_lod = sum(pL) * KVA_B
println()
println("caption figures (energy over the horizon, Δt = 1 h):")
println("  substation import   ", round(tot_sub, digits = 1), " kWh")
println("  battery discharge   ", round(tot_dis, digits = 1), " kWh")
println("  battery charge      ", round(tot_chg, digits = 1), " kWh")
println("  load served         ", round(tot_lod, digits = 1), " kWh")
println("  discharge share of supply at its peak interval: ",
        round(100 * maximum(max(pb[t], 0.0) / (psub[t] + max(pb[t], 0.0)) for t = 1:T),
              digits = 1), " %")
println("  B range ", round(KVA_B * minimum(values(B)), digits = 1), " – ",
        round(KVA_B * maximum(values(B)), digits = 1), " kWh")
