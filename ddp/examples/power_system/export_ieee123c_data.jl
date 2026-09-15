# Export the existing tADMM/OpenDSS parser output so the FilterDDP environment
# does not need to duplicate the OpenDSS dependency.
#
# Run from the repository root:
#   julia --startup-file=no --project=envs/tadmm \
#         ddp/examples/power_system/export_ieee123c_data.jl [system] [T]

using Printf
using Serialization

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO, "envs", "tadmm", "parse_opendss.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
dt = 24.0 / T
# Phase sampling. `range(0, 2pi, length=T)` includes BOTH endpoints, so samples
# 1 and T land on the same phase at every T -- one wasted sample and a skewed
# profile. Worse, at T=3 it samples sin at 0, pi, 2pi, all zero, so the price
# comes out CONSTANT and the instance carries no arbitrage signal whatsoever.
# Every T=3 result generated before 2026-09-15 is on such an instance, with the
# battery moving only 5-33% of rating and, at large10k, 5%.
#
# PROFILE_PERIODIC=1 uses proper periodic sampling, phase_k = 2*pi*(k-1)/T,
# which has no duplicated endpoint and is non-degenerate at every T (at T=3 it
# gives price spread ~118%). It is OPT-IN because it changes the instances that
# existing committed results -- including the centralized IPOPT timing sweep
# feeding the TPEC table -- were measured on.
periodic = get(ENV, "PROFILE_PERIODIC", "0") == "1"
phase = periodic ? [2pi * (k - 1) / T for k in 1:T] : collect(range(0, 2pi, length=T))
load_shape = T == 1 ? [0.97174] : 0.8 .+ 0.2 .* (sin.(phase .- 0.8) .+ 1) ./ 2
cost_shape = T == 1 ? [0.14] : 0.08 .+ 0.12 .* (sin.(phase) .+ 1) ./ 2
pv_shape = let
    pv = zeros(T)
    if T >= 4
        a = max(1, round(Int, 0.25T))
        b = min(T, round(Int, 0.75T))
        pv[a:b] .= sin.(range(0, pi, length=b-a+1))
    elseif T == 3
        pv[2] = 1.0
    elseif T == 2
        pv .= 0.5
    else
        pv .= 1.0
    end
    pv
end

network_kv_base = system == "ieee2522C_1ph" ? 7.2 :
    (system == "large10kC_1ph" ? 12.47 : 2.4018)

data = parse_system_from_dss(system, T;
    LoadShapeLoad=load_shape,
    LoadShapeCost=cost_shape,
    LoadShapePV=pv_shape,
    C_B=1e-6 * minimum(cost_shape),
    delta_t_h=dt,
    kV_B=network_kv_base)

outdir = joinpath(REPO, "ddp", "results", "network_filterddp")
mkpath(outdir)
# Degeneracy guard. A benchmark with no price spread cannot exercise a battery
# scheduler at all, and that went unnoticed for weeks. Report the spreads on
# every export and shout when one is flat, so a degenerate instance can never
# again be mistaken for a passing test.
let
    spread(v) = (lo = minimum(v); hi = maximum(v);
                 abs(lo) < 1e-12 ? (hi - lo) : (hi - lo) / abs(lo))
    ps, ls, vs = spread(cost_shape), spread(load_shape), spread(pv_shape)
    @printf("PROFILE T=%d periodic=%s  price spread=%.1f%%  load spread=%.1f%%  pv spread=%.1f%%
",
            T, periodic, 100ps, 100ls, 100vs)
    if ps < 0.05
        @warn """DEGENERATE INSTANCE: price spread is $(round(100ps, digits=3))%.
        With a flat price there is no arbitrage signal and the battery has almost no
        reason to move, so this instance cannot meaningfully test a scheduling
        algorithm. Re-export with PROFILE_PERIODIC=1, or use a larger T."""
    end
    if ls < 0.05
        @warn "DEGENERATE INSTANCE: load spread is $(round(100ls, digits=3))%."
    end
end

# Periodic instances get their OWN filename. Overwriting the default would make
# every earlier result silently incomparable -- the same class of mistake as the
# flat T=3 price itself, which is exactly what must not happen again.
suffix = periodic ? "_periodic" : ""
outfile = joinpath(outdir, "network_data_$(system)_T$(T)$(suffix).jls")
data[:profile_periodic] = periodic
data[:price_spread] = (maximum(cost_shape) - minimum(cost_shape)) / max(abs(minimum(cost_shape)), 1e-12)
serialize(outfile, data)

println("exported=$outfile")
println("T=$T buses=$(length(data[:Nset])) branches=$(length(data[:Lset])) " *
        "batteries=$(length(data[:Bset])) pv=$(length(data[:Dset])) loads=$(length(data[:NLset]))")
