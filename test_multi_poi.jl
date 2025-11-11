# multi_poi.jl
# Onetime script - uncomment if needed, comment once finished
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "envs", "multi_poi"))
# Pkg.add("Crayons")
# Pkg.add("JuMP")
# Pkg.add("Ipopt")
# Pkg.add("OpenDSSDirect")
# Pkg.instantiate()
# Pkg.precompile()
using JuMP
using Ipopt
using LinearAlgebra
using Crayons
using Printf

# Optional: Uncomment if running OpenDSS validation
# import OpenDSSDirect as dss
# ----------------------------------------------------------------------

data = Dict()
V_1_pu = 1.07
delta_1_deg = 0.0
V_2_pu = 1.07
delta_2_deg = 0.0
alpha_share = 0.75

r1_ohm = 0.025; r2_ohm=0.025; x1_ohm=0.005; x2_ohm=0.005;
P_L_kW = 500
Q_L_kW = P_L_kW*0.75
kVA_B = 1000
kV_B = 11.547 # 20kV line-to-line, 11.547kV line-to-neutral
C_1_dollar_per_kWh = 0.50
C_2_dollar_per_kWh = 0.50

data = Dict(
    :V_1_pu => V_1_pu,
    :delta_1_deg => delta_1_deg,
    :V_2_pu => V_2_pu,
    :delta_2_deg => delta_2_deg,
    :alpha_share => alpha_share,
    :r1_ohm => r1_ohm,
    :x1_ohm => x1_ohm,
    :r2_ohm => r2_ohm,
    :x2_ohm => x2_ohm,
    :P_L_kW => P_L_kW,
    :Q_L_kW => Q_L_kW,
    :kVA_B => kVA_B,
    :kV_B => kV_B,
    :C_1_dollar_per_kWh => C_1_dollar_per_kWh,
    :C_2_dollar_per_kWh => C_2_dollar_per_kWh
)

function process_data!(data)
    # Base values
    S_base = data[:kVA_B]
    V_base = data[:kV_B]

    # Per-unit conversions
    P_L_pu = data[:P_L_kW] / S_base
    Q_L_pu = data[:Q_L_kW] / S_base
    r1_pu = data[:r1_ohm] / ((V_base^2) / S_base)
    x1_pu = data[:x1_ohm] / ((V_base^2) / S_base)
    r2_pu = data[:r2_ohm] / ((V_base^2) / S_base)
    x2_pu = data[:x2_ohm] / ((V_base^2) / S_base)

    # Cost per unit energy (pu)
    C_1_dollar_pu = data[:C_1_dollar_per_kWh] / S_base
    C_2_dollar_pu = data[:C_2_dollar_per_kWh] / S_base

    # Angles in radians
    delta_1_rad = data[:delta_1_deg] * pi / 180
    delta_2_rad = data[:delta_2_deg] * pi / 180

    # Voltage limits
    Vminpu = 0.95
    Vmaxpu = 1.05

    # Update data dict with all required fields
    data[:P_L_pu] = P_L_pu
    data[:Q_L_pu] = Q_L_pu
    data[:r1_pu] = r1_pu
    data[:x1_pu] = x1_pu
    data[:r2_pu] = r2_pu
    data[:x2_pu] = x2_pu
    data[:C_1_dollar_pu] = C_1_dollar_pu
    data[:C_2_dollar_pu] = C_2_dollar_pu
    data[:delta_1_rad] = delta_1_rad
    data[:delta_2_rad] = delta_2_rad
    data[:Vminpu] = Vminpu
    data[:Vmaxpu] = Vmaxpu

end

process_data!(data)


"""
    compute_network_matrices(data)

Compute the incidence matrices B, B_T (tree), B_⊥ (chord) and angle vector β
for the 3-bus system with topology:
- Node 0 (substation): slack bus with V_0 = V_1_pu ∠ δ_1
- Node 1 (substation): PV bus with V_1 = V_2_pu ∠ δ_2  
- Node 2 (load): load bus
- Link 1: connects node 0 → node 2
- Link 2: connects node 1 → node 2

Returns a Dict with:
- :B - reduced incidence matrix (m × n)
- :B_T - tree part of B (n × n)
- :B_⊥ - chord part of B ((m-n) × n)
- :β - angle differences on all links (after OPF solve)
- :spanning_tree_links - indices of tree links
- :chord_links - indices of chord links
"""
function compute_network_matrices(data)
    # Network topology for 3-bus system:
    # Nodes: 0 (slack at grid1), 1 (PV at grid2), 2 (load)
    # We remove node 0 from the formulation → working with nodes 1, 2
    # Links: e=1 (0→2), e=2 (1→2)
    
    # After removing node 0, we have:
    # n = 2 nodes (node 1, node 2 in reduced system, originally nodes 1,2)
    # m = 2 links
    
    n_nodes_reduced = 2  # nodes 1 and 2 (after removing slack node 0)
    m_links = 2          # links 1 and 2
    
    # Full incidence matrix C (before removing node 0)
    # Rows: nodes 0, 1, 2
    # Cols: links 1, 2
    # C[i,e] = +1 if link e leaves node i
    #        = -1 if link e enters node i
    #        =  0 otherwise
    
    # Link 1: 0 → 2 (leaves node 0, enters node 2)
    # Link 2: 1 → 2 (leaves node 1, enters node 2)
    
    C_full = [
        1   0;   # node 0: link 1 leaves, link 2 not incident
        0   1;   # node 1: link 2 leaves, link 1 not incident  
       -1  -1    # node 2: both links enter
    ]
    
    # Reduced incidence matrix B: remove first row (node 0) and transpose
    B = C_full[2:end, :]'  # Remove row 1 (node 0), then transpose
    
    # Choose spanning tree T
    spanning_tree_links = [1, 2]  # Both links are in the tree
    chord_links = Int[]            # No chord links
    
    # B_T: tree part (should be n × n = 2 × 2, and invertible)
    B_T = B[spanning_tree_links, :]
    
    # Check if B_T is invertible
    det_B_T = det(B_T)
    
    # B_⊥: chord part ((m-n) × n matrix)
    B_chord = zeros(0, n_nodes_reduced)
    
    println("Network: 3 nodes, 2 links (radial), det(B_T) = $(round(det_B_T, digits=3))")
    
    return Dict(
        :B => B,
        :B_T => B_T,
        :B_chord => B_chord,
        :n_nodes => n_nodes_reduced,
        :m_links => m_links,
        :spanning_tree_links => spanning_tree_links,
        :chord_links => chord_links,
        :C_full => C_full
    )
end


"""
    compute_beta_angles(data, opf_result, network_matrices)

Compute the angle difference β_ij = ∠(v_i - z*_ij S_ij) for each link (i,j) ∈ E
after solving the OPF.

Arguments:
- data: system data dictionary
- opf_result: Dict with OPF solution (:P_1j, :Q_1j, :P_2j, :Q_2j, :v_j, etc.)
- network_matrices: Dict from compute_network_matrices()

Returns Dict with:
- :beta - vector of angle differences for all links (in radians)
- :beta_deg - vector of angle differences for all links (in degrees)
"""
function compute_beta_angles(data, opf_result, network_matrices)
    # Extract voltage magnitudes squared
    v_0 = data[:V_1_pu]^2  # Node 0 (grid1) - slack
    v_1 = data[:V_2_pu]^2  # Node 1 (grid2) - PV bus
    v_2 = opf_result[:v_j]  # Node 2 (load) - from OPF
    
    # Extract angles (in radians)
    δ_0 = data[:delta_1_rad]  # Node 0 angle
    δ_1 = data[:delta_2_rad]  # Node 1 angle
    
    # For node 2, we need to compute the angle from the OPF results
    # Using the branch flow equations and power flows
    # This is implicit in the BFM, but we can estimate it
    
    # Extract power flows
    P_12 = opf_result[:P_1j]  # Power from grid1 (node 0) to load (node 2) via link 1
    Q_12 = opf_result[:Q_1j]
    P_22 = opf_result[:P_2j]  # Power from grid2 (node 1) to load (node 2) via link 2
    Q_22 = opf_result[:Q_2j]
    
    # Extract line impedances
    r_1 = data[:r1_pu]
    x_1 = data[:x1_pu]
    z_1 = r_1 + im*x_1
    
    r_2 = data[:r2_pu]
    x_2 = data[:x2_pu]
    z_2 = r_2 + im*x_2
    
    # For a radial network in BFM, we can compute the angle at node 2
    # from the relationship: δ_j ≈ δ_i - (r*P + x*Q)/(V_i * V_j)
    # But since we have two sources, let's use the dominant one
    
    # Approximate angle at node 2 from path through link 1
    V_0_pu = sqrt(v_0)
    V_2_pu = sqrt(v_2)
    δ_2_from_link1 = δ_0 - (r_1 * P_12 + x_1 * Q_12) / (V_0_pu * V_2_pu)
    
    # Approximate angle at node 2 from path through link 2
    V_1_pu = sqrt(v_1)
    δ_2_from_link2 = δ_1 - (r_2 * P_22 + x_2 * Q_22) / (V_1_pu * V_2_pu)
    
    # Average them (or pick one if one dominates)
    δ_2 = (δ_2_from_link1 + δ_2_from_link2) / 2
    
    # Compute β for each link: β_ij = ∠(v_i - z*_ij S_ij)
    # where S_ij = P_ij + j*Q_ij is the complex power flow
    
    # Link 1: from node 0 to node 2
    S_12 = P_12 + im*Q_12  # Complex power on link 1
    V_0_complex = V_0_pu * exp(im * δ_0)
    V_0_minus_zS_1 = V_0_complex - conj(z_1) * S_12
    β_1 = angle(V_0_minus_zS_1)
    
    # Link 2: from node 1 to node 2  
    S_22 = P_22 + im*Q_22  # Complex power on link 2
    V_1_complex = V_1_pu * exp(im * δ_1)
    V_1_minus_zS_2 = V_1_complex - conj(z_2) * S_22
    β_2 = angle(V_1_minus_zS_2)
    
    # Also compute the angle differences directly
    Δδ_01 = δ_0 - δ_2  # Angle difference across link 1
    Δδ_12 = δ_1 - δ_2  # Angle difference across link 2
    
    β_vec = [β_1, β_2]
    β_deg_vec = rad2deg.(β_vec)
    
    println("β: [$(round(β_deg_vec[1], digits=3))°, $(round(β_deg_vec[2], digits=3))°], δ₂ = $(round(rad2deg(δ_2), digits=3))°")
    
    return Dict(
        :beta => β_vec,
        :beta_deg => β_deg_vec,
        :delta_2 => δ_2,
        :delta_2_deg => rad2deg(δ_2),
        :angle_diff_link1 => Δδ_01,
        :angle_diff_link2 => Δδ_12,
        :V_0_minus_zS_1 => V_0_minus_zS_1,
        :V_1_minus_zS_2 => V_1_minus_zS_2
    )
end


function solve_two_poi_opf(data)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    # Variables
    @variable(model, P_1j >= 0)
    @variable(model, Q_1j)
    @variable(model, l_1j >= 0)
    @variable(model, P_2j >= 0)
    @variable(model, Q_2j)
    @variable(model, l_2j >= 0)
    @variable(model, v_j >= data[:Vminpu]^2)
    @constraint(model, v_j <= data[:Vmaxpu]^2)

    # Extract parameters
    r1j = data[:r1_pu]
    x1j = data[:x1_pu]
    r2j = data[:r2_pu]
    x2j = data[:x2_pu]
    v_1 = data[:V_1_pu]^2
    v_2 = data[:V_2_pu]^2
    alpha = data[:alpha_share]
    P_L_j = data[:P_L_pu]
    Q_L_j = data[:Q_L_pu]
    C_1 = data[:C_1_dollar_pu]
    C_2 = data[:C_2_dollar_pu]

    # Constraints
    # 1. NRealPB
    @constraint(model, (P_1j - r1j * l_1j) + (P_2j - r2j * l_2j) == P_L_j)
    # 2. NReacPB
    @constraint(model, (Q_1j - x1j * l_1j) + (Q_2j - x2j * l_2j) == Q_L_j)
    # 3. KVL1
    @constraint(model, v_j == v_1 - 2 * (r1j * P_1j + x1j * Q_1j) + (r1j^2 + x1j^2) * l_1j)
    # 4. KVL2
    @constraint(model, v_j == v_2 - 2 * (r2j * P_2j + x2j * Q_2j) + (r2j^2 + x2j^2) * l_2j)
    # 5. BCPF1
    @constraint(model, P_1j^2 + Q_1j^2 >= v_1 * l_1j)
    # 6. BCPF2
    @constraint(model, P_2j^2 + Q_2j^2 >= v_2 * l_2j)
    # 7. Power sharing
    @constraint(model, P_2j == alpha * P_1j)
    # 8. Voltage limits (already in variable bounds)
    # 9. P_1j, P_2j >= 0 (already in variable bounds)

    # Objective
    @objective(model, Min, C_1 * P_1j + C_2 * P_2j)

    # Solve
    optimize!(model)

    # Collect results
    modelDict = Dict(
        :P_1j => value(P_1j),
        :Q_1j => value(Q_1j),
        :l_1j => value(l_1j),
        :P_2j => value(P_2j),
        :Q_2j => value(Q_2j),
        :l_2j => value(l_2j),
        :v_j => value(v_j)
    )
    return modelDict
end

# Example usage:
modelDict = solve_two_poi_opf(data)

println("\n" * "="^80)
println("NETWORK TOPOLOGY & OPF RESULTS")
println("="^80)

# Compute network matrices (topology-dependent, OPF-independent)
network_matrices = compute_network_matrices(data)

println()
println("Power dispatches: P₁=$(round(modelDict[:P_1j] * data[:kVA_B], digits=1)) kW, P₂=$(round(modelDict[:P_2j] * data[:kVA_B], digits=1)) kW")

# Compute beta angles (requires OPF results)
beta_results = compute_beta_angles(data, modelDict, network_matrices)

println("="^80)

# ==================================================================================
# SLACK BUS CONFIGURATION COMPARISON
# ==================================================================================

"""
    analyze_slack_configuration(slack_node, data)

Solve OPF with given slack bus and compute:
1. Objective function value
2. Phase-shifted angles for non-slack substations
3. Recommended median angle for each non-slack substation
4. Total angle deviation (grand RMS)

Arguments:
- slack_node: 0 or 1 (which substation is slack)
- data: system data dictionary

Returns Dict with:
- :objective - optimal cost
- :angles_deg - Dict of non-slack substation angles vs time
- :median_angles_deg - recommended angle for each non-slack substation
- :angle_ranges_deg - [min, max] for each non-slack substation
- :rmse_values_deg - RMSE for each non-slack substation
- :total_deviation_deg - grand RMS across all non-slack substations
"""
function analyze_slack_configuration(slack_node::Int, data::Dict)
    # For now, we have a single time period
    # In MPOPF, this will loop over multiple time periods
    
    # Setup model with specified slack bus
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    @variable(model, P_1j >= 0)
    @variable(model, Q_1j)
    @variable(model, l_1j >= 0)
    @variable(model, P_2j >= 0)
    @variable(model, Q_2j)
    @variable(model, l_2j >= 0)
    @variable(model, v_j >= data[:Vminpu]^2)
    @constraint(model, v_j <= data[:Vmaxpu]^2)
    
    # Extract parameters
    r1j = data[:r1_pu]
    x1j = data[:x1_pu]
    r2j = data[:r2_pu]
    x2j = data[:x2_pu]
    v_1 = data[:V_1_pu]^2
    v_2 = data[:V_2_pu]^2
    alpha = data[:alpha_share]
    P_L_j = data[:P_L_pu]
    Q_L_j = data[:Q_L_pu]
    C_1 = data[:C_1_dollar_pu]
    C_2 = data[:C_2_dollar_pu]
    
    # Constraints
    @constraint(model, (P_1j - r1j * l_1j) + (P_2j - r2j * l_2j) == P_L_j)
    @constraint(model, (Q_1j - x1j * l_1j) + (Q_2j - x2j * l_2j) == Q_L_j)
    @constraint(model, v_j == v_1 - 2 * (r1j * P_1j + x1j * Q_1j) + (r1j^2 + x1j^2) * l_1j)
    @constraint(model, v_j == v_2 - 2 * (r2j * P_2j + x2j * Q_2j) + (r2j^2 + x2j^2) * l_2j)
    @constraint(model, P_1j^2 + Q_1j^2 >= v_1 * l_1j)
    @constraint(model, P_2j^2 + Q_2j^2 >= v_2 * l_2j)
    @constraint(model, P_2j == alpha * P_1j)
    
    @objective(model, Min, C_1 * P_1j + C_2 * P_2j)
    
    optimize!(model)
    
    obj_value = objective_value(model)
    
    # Extract results
    P_1j_val = value(P_1j)
    Q_1j_val = value(Q_1j)
    P_2j_val = value(P_2j)
    Q_2j_val = value(Q_2j)
    v_j_val = value(v_j)
    
    # Compute angles using branch flow relationships
    # We have 2 substations: node 0 (grid1) and node 1 (grid2)
    # One is slack, the other is non-slack
    
    # Extract impedances
    r_1 = data[:r1_pu]
    x_1 = data[:x1_pu]
    r_2 = data[:r2_pu]
    x_2 = data[:x2_pu]
    
    V_0_pu = sqrt(v_1)
    V_1_pu = sqrt(v_2)
    V_2_pu = sqrt(v_j_val)
    
    # Get slack bus angle (always 0° for slack)
    δ_slack = 0.0
    
    # Compute non-slack substation angles (θ) from OPF solution
    # These are the phase-shifter angles computed POST-OPTIMIZATION
    # NOT from input delta values (those are fixed physical angles, irrelevant here)
    # We compute from power flows using voltage drop relationships
    
    if slack_node == 0
        # Substation 1 is slack at θ₁ = 0°
        # Compute θ₂ (angle at substation 2) from power flow solution
        # Using voltage drop: θ_to = θ_from - (r*P + x*Q) / (V_from * V_to)
        
        # Angle at load from sub1 path
        θ_load_from_sub1 = 0.0 - (r_1 * P_1j_val + x_1 * Q_1j_val) / (V_0_pu * V_2_pu)
        
        # Angle at sub2 from load (going backwards)
        θ_sub2 = θ_load_from_sub1 + (r_2 * P_2j_val + x_2 * Q_2j_val) / (V_1_pu * V_2_pu)
        
        angles_deg = Dict(2 => rad2deg(θ_sub2))
        
    else  # slack_node == 1
        # Substation 2 is slack at θ₂ = 0°
        # Compute θ₁ (angle at substation 1) from power flow solution
        
        # Angle at load from sub2 path
        θ_load_from_sub2 = 0.0 - (r_2 * P_2j_val + x_2 * Q_2j_val) / (V_1_pu * V_2_pu)
        
        # Angle at sub1 from load (going backwards)
        θ_sub1 = θ_load_from_sub2 + (r_1 * P_1j_val + x_1 * Q_1j_val) / (V_0_pu * V_2_pu)
        
        angles_deg = Dict(1 => rad2deg(θ_sub1))
    end
    
    # For single time period, median is just the value itself
    median_angles_deg = Dict()
    angle_ranges_deg = Dict()
    
    for (sub, angle) in angles_deg
        median_angles_deg[sub] = angle
        angle_ranges_deg[sub] = [angle, angle]  # Single time period: min = max = angle
    end
    
    # Compute RMSE for each non-slack substation
    # RMSE = sqrt(mean((angle - median)^2))
    # For single time period: RMSE = 0
    rmse_values_deg = Dict()
    for sub in keys(angles_deg)
        # In MPOPF: rmse = sqrt(mean((angles[:, t] .- median).^2))
        rmse_values_deg[sub] = 0.0  # Single time period
    end
    
    # Total angle deviation (grand RMS)
    # Grand RMS = sqrt(mean(all_rmse_values.^2))
    if length(rmse_values_deg) > 0
        total_deviation = sqrt(sum(v^2 for v in values(rmse_values_deg)) / length(rmse_values_deg))
    else
        total_deviation = 0.0
    end
    
    return Dict(
        :objective => obj_value,
        :angles_deg => angles_deg,
        :median_angles_deg => median_angles_deg,
        :angle_ranges_deg => angle_ranges_deg,
        :rmse_values_deg => rmse_values_deg,
        :total_deviation_deg => total_deviation,
        :P_1j => P_1j_val,
        :P_2j => P_2j_val,
        :Q_1j => Q_1j_val,
        :Q_2j => Q_2j_val
    )
end

println("\n" * "="^80)
println("SLACK BUS CONFIGURATION COMPARISON")
println("="^80)
println()

# Analyze each slack bus configuration
results_by_slack = Dict()

for slack_sub in [1, 2]
    println("Analyzing with Substation $slack_sub as slack bus...")
    results_by_slack[slack_sub] = analyze_slack_configuration(slack_sub - 1, data)
end

println("\n" * "="^80)
println("TABLE I: COMPARISON OF SLACK BUS CONFIGURATIONS")
println("="^80)

# Print header
header = @sprintf("%-15s | %-20s | %-20s", "Slack Bus", "Total Angle", "Operational")
header2 = @sprintf("%-15s | %-20s | %-20s", "(Substation)", "Deviation* (°)", "Cost (\$)")
separator = "-"^length(header)

println(separator)
println(header)
println(header2)
println(separator)

for slack_sub in sort(collect(keys(results_by_slack)))
    result = results_by_slack[slack_sub]
    
    # Convert to physical units
    cost_dollars = result[:objective] * data[:kVA_B]  # Convert from pu
    total_dev = result[:total_deviation_deg]
    
    # Format with colors
    if slack_sub == 1
        sub_str = Crayon(foreground = :green, bold = true)(@sprintf("%-15d", slack_sub))
    else
        sub_str = Crayon(foreground = :cyan, bold = true)(@sprintf("%-15d", slack_sub))
    end
    
    dev_str = @sprintf("%-20.3f", total_dev)
    cost_str = @sprintf("%-20.2f", cost_dollars)
    
    @printf("%s | %s | %s\n", sub_str, dev_str, cost_str)
end

println(separator)
println("*Sum of RMS angle deviations across all non-slack substations")
println()

# Print detailed angle coordination for each configuration
for slack_sub in sort(collect(keys(results_by_slack)))
    result = results_by_slack[slack_sub]
    
    println("\n" * "="^80)
    println("TABLE II: ANGLE COORDINATION DETAILS (SLACK: SUBSTATION $slack_sub)")
    println("="^80)
    
    header_detail = @sprintf("%-15s | %-20s | %-25s | %-15s", 
                            "Substation", "Median", "Angle Range", "Std Dev")
    header_detail2 = @sprintf("%-15s | %-20s | %-25s | %-15s",
                             "(Non-Slack)", "Angle (°)", "(°)", "(°)")
    sep_detail = "-"^length(header_detail)
    
    println(sep_detail)
    println(header_detail)
    println(header_detail2)
    println(sep_detail)
    
    # Get non-slack substations
    non_slack_subs = sort(collect(keys(result[:angles_deg])))
    
    for sub in non_slack_subs
        median_angle = result[:median_angles_deg][sub]
        angle_range = result[:angle_ranges_deg][sub]
        rmse = result[:rmse_values_deg][sub]
        
        sub_str = Crayon(foreground = :yellow)(@sprintf("%-15d", sub))
        median_str = @sprintf("θ%-2d,%-2d = %7.3f", sub, slack_sub, median_angle)
        range_str = @sprintf("[%7.3f, %7.3f]", angle_range[1], angle_range[2])
        rmse_str = @sprintf("%-15.3f", rmse)
        
        @printf("%s | %-20s | %-25s | %s\n", sub_str, median_str, range_str, rmse_str)
    end
    
    println(sep_detail)
    
    # Print power dispatch summary
    println("\nPower Dispatches:")
    println("  Substation 1: P = $(round(result[:P_1j] * data[:kVA_B], digits=1)) kW, Q = $(round(result[:Q_1j] * data[:kVA_B], digits=1)) kVAr")
    println("  Substation 2: P = $(round(result[:P_2j] * data[:kVA_B], digits=1)) kW, Q = $(round(result[:Q_2j] * data[:kVA_B], digits=1)) kVAr")
end

println("\n" * "="^80)
println("\nScript complete! Slack bus configuration analysis finished.")
println("For OpenDSS validation, see: envs/multi_poi/archive/opendss_validation.jl")
println("="^80)

# Uncomment below to run OpenDSS validation:
# include("envs/multi_poi/archive/opendss_validation.jl")
# validation_results = validate_opf_with_opendss(modelDict, data)