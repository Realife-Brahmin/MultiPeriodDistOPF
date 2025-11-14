#!/usr/bin/env julia
# test_parser_multi_poi.jl - Test the multi-POI OpenDSS parser

import Pkg
Pkg.activate(joinpath(@__DIR__, "envs", "multi_poi"))

using OpenDSSDirect
using Printf

# Include the parser
include(joinpath(@__DIR__, "envs", "multi_poi", "parse_opendss_multi_poi.jl"))

println("\n" * "="^80)
println("TESTING MULTI-POI OPENDSS PARSER")
println("="^80)

# System parameters
systemName = "ieee123_5poi_1ph"
T = 24
delta_t_h = 1.0
kVA_B = 1000.0
kV_B = 11.547  # 20kV LL -> 11.547kV LN

# Create time-varying profiles
LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2
LoadShapePV = zeros(T)

# Create different cost profiles for each substation (phase-shifted)
LoadShapeCost_dict = Dict(
    1 => 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .+ π/4) .+ 1) ./ 2,
    2 => 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .- π/4) .+ 1) ./ 2,
    3 => 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2,
    4 => 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .+ π/2) .+ 1) ./ 2,
    5 => 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .- π/2) .+ 1) ./ 2
)

C_B = 1e-6 * 0.08

# Parse system
try
    data = parse_system_from_dss_multi_poi(
        systemName,
        T;
        kVA_B=kVA_B,
        kV_B=kV_B,
        LoadShapeLoad=LoadShapeLoad,
        LoadShapePV=LoadShapePV,
        LoadShapeCost_dict=LoadShapeCost_dict,
        C_B=C_B,
        delta_t_h=delta_t_h
    )
    
    # Print summary
    println("\n" * "="^80)
    println("PARSING SUMMARY")
    println("="^80)
    
    println("\n--- Multi-POI Information ---")
    println("  Number of substations: $(data[:num_substations])")
    println("  Substation buses: $(data[:substation_buses])")
    println("  POI lines: $(data[:poi_lines])")
    println("  POI connections:")
    for (sub_bus, dist_bus) in data[:poi_connections]
        println("    $sub_bus → $dist_bus")
    end
    
    println("\n--- Network Topology ---")
    println("  Total distribution buses: $(data[:N])")
    println("  Total distribution branches: $(data[:m])")
    println("  Branches from POI connections: $(length(data[:L1set]))")
    println("  Other branches: $(length(data[:Lm1set]))")
    
    println("\n--- System Components ---")
    println("  Loads:")
    println("    Buses with loads: $(data[:N_L])")
    println("    Total base load: $(round(sum(values(data[:p_L_R_pu])) * kVA_B, digits=1)) kW")
    
    println("  PV Systems:")
    println("    Buses with PV: $(data[:n_D])")
    if data[:n_D] > 0
        println("    Total PV capacity: $(round(sum(values(data[:p_D_R_pu])) * kVA_B, digits=1)) kW")
    end
    
    println("  Battery Storage:")
    println("    Buses with batteries: $(data[:n_B])")
    if data[:n_B] > 0
        println("    Total battery capacity: $(round(sum(values(data[:B_R])), digits=1)) kWh")
        println("    Total battery power: $(round(sum(values(data[:P_B_R])), digits=1)) kW")
    end
    
    println("\n--- Impedances ---")
    println("  Base impedance: $(round(data[:Z_B], digits=4)) Ω")
    println("  Distribution lines (Lm1set): $(data[:m_dist])")
    println("  POI lines (L1set): $(data[:m_poi])")
    println("  Total lines (Lset): $(data[:m])")
    
    println("\n--- Time Series ---")
    println("  Time horizon: T = $T")
    println("  Time step: Δt = $(delta_t_h) hours")
    println("  Load shape range: [$(round(minimum(LoadShapeLoad), digits=3)), $(round(maximum(LoadShapeLoad), digits=3))]")
    
    println("\n--- Voltage Limits ---")
    println("  Distribution buses: $(data[:Vminpu][1]) to $(data[:Vmaxpu][1]) pu")
    sub_vmin = data[:Vminpu_sub]["1s"]
    sub_vmax = data[:Vmaxpu_sub]["1s"]
    println("  Substation buses: $sub_vmin to $sub_vmax pu")
    
    println("\n--- Network Sets for MPOPF ---")
    println("  Sset (substations): $(length(data[:Sset])) buses")
    println("  Nset (all buses): $(length(data[:Nset])) buses (substations + distribution)")
    println("  Nm1set (distribution only): $(length(data[:Nm1set])) buses")
    println("  Lset (all lines): $(length(data[:Lset])) lines")
    println("  L1set (POI lines): $(length(data[:L1set])) lines")
    println("  Lm1set (distribution lines): $(length(data[:Lm1set])) lines")
    
    println("\n" * "="^80)
    println("✓ PARSER TEST COMPLETE - ALL DATA SUCCESSFULLY LOADED")
    println("="^80)
    
    # Verify some key data structures exist
    required_keys = [:systemName, :T, :N, :Sset, :Nset, :Nm1set, :Lset, :L1set, :Lm1set,
                     :num_substations, :substation_buses, :poi_connections, 
                     :rdict_pu, :xdict_pu, :p_L_pu, :q_L_pu, :Vminpu, :Vmaxpu]
    
    missing_keys = filter(k -> !haskey(data, k), required_keys)
    if !isempty(missing_keys)
        println("\n⚠ WARNING: Missing required keys: $missing_keys")
    else
        println("\n✓ All required data structures present in data dictionary")
    end
    
    # Print sample of data for verification
    println("\n--- Sample Data Verification ---")
    println("  Sample branch impedance (first branch):")
    first_branch = data[:Lset][1]
    println("    Branch $(first_branch): r = $(round(data[:rdict_pu][first_branch], digits=6)) pu, x = $(round(data[:xdict_pu][first_branch], digits=6)) pu")
    
    println("\n  Sample POI impedance:")
    first_poi_bus = data[:substation_buses][1]
    if haskey(data[:poi_rdict], first_poi_bus)
        println("    POI $first_poi_bus: R = $(data[:poi_rdict][first_poi_bus]) Ω, X = $(data[:poi_xdict][first_poi_bus]) Ω")
    end
    
    println("\n  Sample load (bus with highest load):")
    if !isempty(data[:p_L_R_pu])
        max_load_bus = argmax(data[:p_L_R_pu])
        println("    Bus $max_load_bus: P = $(round(data[:p_L_R_pu][max_load_bus] * kVA_B, digits=1)) kW, Q = $(round(data[:q_L_R_pu][max_load_bus] * kVA_B, digits=1)) kVAr")
    end
    
catch e
    println("\n" * "="^80)
    println("❌ ERROR DURING PARSING")
    println("="^80)
    println("\nError message:")
    showerror(stdout, e)
    println("\n\nStacktrace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end
