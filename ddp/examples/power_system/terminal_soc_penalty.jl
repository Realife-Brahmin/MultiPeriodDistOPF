# Soft terminal state-of-charge penalty -- the single source of truth.
#
# Standing instruction (user, 2026-09-18): terminal SOC is ALWAYS a soft
# constraint, in every algorithm. The objective carries
#
#     + gamma * sum_j (B_j^T - B_j^0)^2        (B in p.u.h)
#
# with gamma FIXED PER SYSTEM and independent of the horizon T. Every solver
# (FilterDDP driver, centralized Ipopt, and tADMM / run_bf when they are next
# run) must read gamma from here, so no two codes can ever solve different
# problems. gamma is NOT C_B: C_B prices battery power in every period
# (C_B * P_B^2); gamma prices only where the energy ends up.
#
# Sizing rule: the marginal penalty 2*gamma*dB equals the mean market value of
# energy, cbar * kVA_B, at a deviation dB of 10% of the system's median battery
# energy rating B_R. So gamma = 5 * cbar * kVA_B / median(B_R):
#   - soft: small deviations cost less than the energy is worth, so the battery
#     may still end off B^0 when arbitrage pays for it;
#   - binding: a full 30%-SOC discharge from B^0 = 62.5% costs ~3.25x the energy
#     value at the margin, so it will not dump its store for free at T;
#   - not overpowering: even that worst case is a few percent of the objective;
#   - convex for any gamma > 0.
# Inputs measured 2026-09-18 from the periodic exports: cbar = 0.14 $/kWh at
# every T (periodic sampling makes the mean exact), kVA_B = 1000, and median
# B_R = 0.039998 / 0.021282 / 1.666800 p.u.h. Ratings do not change with T,
# hence neither does gamma.
const GAMMA_TERMINAL_BY_SYSTEM = Dict(
    "ieee123C_1ph"  => 1.75e4,
    "ieee2522C_1ph" => 3.29e4,
    "large10kC_1ph" => 4.20e2,
)

"""True when the soft terminal SOC penalty is on (TERMINAL_SOC_SOFT=1).
Opt-in only so that runs already in flight, which read the same driver file,
keep their original formulation bit for bit."""
terminal_soc_soft() = get(ENV, "TERMINAL_SOC_SOFT", "0") == "1"

function gamma_terminal(system::AbstractString)
    haskey(ENV, "GAMMA_TERMINAL_OVERRIDE") && return parse(Float64, ENV["GAMMA_TERMINAL_OVERRIDE"])
    haskey(GAMMA_TERMINAL_BY_SYSTEM, system) ||
        error("no terminal-SOC gamma fixed for system $system; add it to terminal_soc_penalty.jl")
    return GAMMA_TERMINAL_BY_SYSTEM[system]
end
