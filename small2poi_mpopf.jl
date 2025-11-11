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

# r1_ohm = 0.025; r2_ohm=0.025; x1_ohm=0.025; x2_ohm=0.025;
r1_ohm = 0.025; r2_ohm=0.075; x1_ohm=0.025; x2_ohm=0.025;
P_L_kW = 500
Q_L_kW = P_L_kW*0.75
kVA_B = 1000
kV_B = 11.547 # 20kV line-to-line, 11.547kV line-to-neutral
C_1_dollar_per_kWh = 1.5
C_2_dollar_per_kWh = 1.5

# Power limits for each substation (kW)
P_1_max_kW = 350.0
P_2_max_kW = 350.0

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
    :C_2_dollar_per_kWh => C_2_dollar_per_kWh,
    :P_1_max_kW => P_1_max_kW,
    :P_2_max_kW => P_2_max_kW
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
    
    # Power limits (pu)
    P_1_max_pu = data[:P_1_max_kW] / S_base
    P_2_max_pu = data[:P_2_max_kW] / S_base

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
    data[:P_1_max_pu] = P_1_max_pu
    data[:P_2_max_pu] = P_2_max_pu

end

process_data!(data)


"""
    compute_network_matrices(data, slack_node)

Compute the incidence matrices B, B_T (tree), B_⊥ (chord) for the 3-bus system
with the specified slack bus removed.

Topology:
- Node 0 (substation 1): can be slack
- Node 1 (substation 2): can be slack
- Node 2 (load): never slack
- Link 1: connects node 0 → node 2
- Link 2: connects node 1 → node 2

Arguments:
- data: system data dictionary
- slack_node: which node is slack (0 or 1)

Returns a Dict with:
- :B - reduced incidence matrix (m × n) after removing slack node
- :B_T - tree part of B (n × n)
- :B_⊥ - chord part of B ((m-n) × n)
- :spanning_tree_links - indices of tree links
- :chord_links - indices of chord links
"""
function compute_network_matrices(data, slack_node::Int)
    # Network topology for 3-bus system:
    # Nodes: 0 (substation 1), 1 (substation 2), 2 (load)
    # Links: e=1 (0→2), e=2 (1→2)
    
    # Full incidence matrix C (before removing slack)
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
    
    # Remove slack node row and transpose to get B
    if slack_node == 0
        # Remove node 0 (row 1), keep nodes 1, 2
        B = C_full[2:end, :]'
        n_nodes_reduced = 2
    elseif slack_node == 1
        # Remove node 1 (row 2), keep nodes 0, 2
        B = C_full[[1,3], :]'
        n_nodes_reduced = 2
    else
        error("Invalid slack node: $slack_node. Must be 0 or 1.")
    end
    
    m_links = 2
    
    # For a radial network, all links are in the spanning tree
    spanning_tree_links = [1, 2]  # Both links are in the tree
    chord_links = Int[]            # No chord links
    
    # B_T: tree part (should be n × n and invertible)
    B_T = B[spanning_tree_links, :]
    
    # Check if B_T is invertible
    det_B_T = det(B_T)
    
    # B_⊥: chord part ((m-n) × n matrix) - empty for radial network
    B_chord = zeros(0, n_nodes_reduced)
    
    println("  Network with slack=$slack_node: det(B_T) = $(round(det_B_T, digits=3))")
    
    return Dict(
        :B => B,
        :B_T => B_T,
        :B_chord => B_chord,
        :n_nodes => n_nodes_reduced,
        :m_links => m_links,
        :spanning_tree_links => spanning_tree_links,
        :chord_links => chord_links,
        :C_full => C_full,
        :slack_node => slack_node
    )
end


"""
    compute_beta_angles(data, opf_result, slack_node)

Compute the angle difference β_ij = ∠(v_i - z*_ij S_ij) for each link (i,j) ∈ E
after solving the OPF, as per equation (27) in the paper.

Arguments:
- data: system data dictionary
- opf_result: Dict with OPF solution (:P_1j, :Q_1j, :P_2j, :Q_2j, :v_j, etc.)
- slack_node: which node is slack (0 or 1)

Returns Dict with:
- :beta - vector of angle differences for all links (in radians)
- :beta_deg - vector of angle differences for all links (in degrees)
"""
function compute_beta_angles(data, opf_result, slack_node::Int)
    # Compute β_ij = ∠(v_i - z*_ij S_ij) for each link (i,j) ∈ E
    # as per equation (27) in the paper
    
    # Extract voltage magnitudes
    V_0_pu = data[:V_1_pu]  # Node 0 (substation 1)
    V_1_pu = data[:V_2_pu]  # Node 1 (substation 2)
    V_2_pu = sqrt(opf_result[:v_j])  # Node 2 (load) - from OPF
    
    # For slack node, angle is 0
    if slack_node == 0
        δ_0 = 0.0
    else
        δ_0 = data[:delta_1_rad]  # Non-slack: use input angle (for β computation only)
    end
    
    if slack_node == 1
        δ_1 = 0.0
    else
        δ_1 = data[:delta_2_rad]  # Non-slack: use input angle (for β computation only)
    end
    
    # Extract power flows
    P_02 = opf_result[:P_1j]  # Power on link 1: node 0 → node 2
    Q_02 = opf_result[:Q_1j]
    P_12 = opf_result[:P_2j]  # Power on link 2: node 1 → node 2
    Q_12 = opf_result[:Q_2j]
    
    # Extract line impedances
    r_1 = data[:r1_pu]
    x_1 = data[:x1_pu]
    z_1 = r_1 + im*x_1
    
    r_2 = data[:r2_pu]
    x_2 = data[:x2_pu]
    z_2 = r_2 + im*x_2
    
    # Compute β for each link: β_ij = ∠(v_i - z*_ij S_ij)
    # Link 1: from node 0 to node 2
    S_02 = P_02 + im*Q_02
    V_0_complex = V_0_pu * exp(im * δ_0)
    V_0_minus_zS = V_0_complex - conj(z_1) * S_02
    β_1 = angle(V_0_minus_zS)
    
    # Link 2: from node 1 to node 2  
    S_12 = P_12 + im*Q_12
    V_1_complex = V_1_pu * exp(im * δ_1)
    V_1_minus_zS = V_1_complex - conj(z_2) * S_12
    β_2 = angle(V_1_minus_zS)
    
    β_vec = [β_1, β_2]
    
    return Dict(
        :beta => β_vec,
        :beta_deg => rad2deg.(β_vec)
    )
end


function solve_two_poi_opf(data)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    # Variables
    @variable(model, 0 <= P_1j <= data[:P_1_max_pu])
    @variable(model, Q_1j)
    @variable(model, l_1j >= 0)
    @variable(model, 0 <= P_2j <= data[:P_2_max_pu])
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
    # alpha = data[:alpha_share]  # Not used when power limits are active
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
    # 7. Power limits (already in variable bounds)
    # 8. Voltage limits (already in variable bounds)
    # Note: Power sharing constraint (P_2j == alpha * P_1j) commented out
    #       to allow optimizer to choose based on costs and power limits

    # Objective
    @objective(model, Min, C_1 * P_1j + C_2 * P_2j)

    # Solve
    optimize!(model)
    
    # Check termination status
    status = termination_status(model)
    
    if status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL
        status_crayon = Crayon(foreground = :green, bold = true)
        println(status_crayon("✓ OPF Solved: $status"))
    else
        status_crayon = Crayon(foreground = :red, bold = true)
        println(status_crayon("✗ OPF Failed: $status"))
    end

    # Collect results
    modelDict = Dict(
        :P_1j => value(P_1j),
        :Q_1j => value(Q_1j),
        :l_1j => value(l_1j),
        :P_2j => value(P_2j),
        :Q_2j => value(Q_2j),
        :l_2j => value(l_2j),
        :v_j => value(v_j),
        :status => status
    )
    return modelDict
end

println("\n" * "="^80)
println("NETWORK TOPOLOGY")
println("="^80)
println("3-bus system: 2 substations (nodes 0, 1) + 1 load (node 2)")
println("Links: 0→2 (link 1), 1→2 (link 2)")
println("Network matrices B, B_T computed per slack configuration")
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
    
    @variable(model, 0 <= P_1j <= data[:P_1_max_pu])
    @variable(model, Q_1j)
    @variable(model, l_1j >= 0)
    @variable(model, 0 <= P_2j <= data[:P_2_max_pu])
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
    # Note: Power sharing constraint removed, using power limits instead
    
    @objective(model, Min, C_1 * P_1j + C_2 * P_2j)
    
    optimize!(model)
    
    # Check termination status
    status = termination_status(model)
    
    if status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL
        status_crayon = Crayon(foreground = :green, bold = true)
        println(status_crayon("  ✓ Slack=$slack_node: $status"))
    else
        status_crayon = Crayon(foreground = :red, bold = true)
        println(status_crayon("  ✗ Slack=$slack_node: $status"))
    end
    
    obj_value = objective_value(model)
    
    # Extract results
    P_1j_val = value(P_1j)
    Q_1j_val = value(Q_1j)
    P_2j_val = value(P_2j)
    Q_2j_val = value(Q_2j)
    v_j_val = value(v_j)
    
    # Compute angles using matrix formula from paper: θ_★ = P(B^(-1)β)
    # where β_ij = ∠(v_i - z*_ij S_ij) from equation (27)
    # and P projects to (-π, π]
    
    # Get network matrices for this slack configuration
    network_matrices = compute_network_matrices(data, slack_node)
    B = network_matrices[:B]
    B_T = network_matrices[:B_T]
    
    # Package OPF results
    opf_result = Dict(
        :P_1j => P_1j_val,
        :Q_1j => Q_1j_val,
        :P_2j => P_2j_val,
        :Q_2j => Q_2j_val,
        :v_j => v_j_val
    )
    
    # Compute β vector for all links
    beta_result = compute_beta_angles(data, opf_result, slack_node)
    β = beta_result[:beta]
    
    # Compute angles: θ_★ = P(B^(-1)β)
    # For radial network, B = B_T (all links are tree links)
    θ_star = B_T \ β  # Solve B_T * θ_★ = β
    
    # Project to (-π, π]
    θ_star = mod.(θ_star .+ π, 2π) .- π
    
    # Map back to substation indices
    # θ_star corresponds to non-slack nodes in the reduced system
    if slack_node == 0
        # Removed node 0, so θ_star corresponds to [node 1, node 2]
        # We care about node 1 (substation 2)
        angles_deg = Dict(2 => rad2deg(θ_star[1]))
    else  # slack_node == 1
        # Removed node 1, so θ_star corresponds to [node 0, node 2]
        # We care about node 0 (substation 1)
        angles_deg = Dict(1 => rad2deg(θ_star[1]))
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
    println("Solving OPF with Substation $slack_sub as slack bus...")
    results_by_slack[slack_sub] = analyze_slack_configuration(slack_sub - 1, data)
end

println()

# Print detailed angle coordination for each configuration (TABLE II)
for slack_sub in sort(collect(keys(results_by_slack)))
    result = results_by_slack[slack_sub]
    
    println("\n" * "="^80)
    title2_crayon = Crayon(foreground = :light_magenta, bold = true)
    println(title2_crayon("TABLE II: ANGLE COORDINATION DETAILS (SLACK: SUBSTATION $slack_sub)"))
    println("="^80)
    
    header2_crayon = Crayon(foreground = :white, bold = true)
    header_detail = @sprintf("%-15s | %-20s | %-25s | %-15s", 
                            "Substation", "Median", "Angle Range", "Std Dev")
    header_detail2 = @sprintf("%-15s | %-20s | %-25s | %-15s",
                             "(Non-Slack)", "Angle (°)", "(°)", "(°)")
    sep_detail = "-"^length(header_detail)
    
    println(sep_detail)
    println(header2_crayon(header_detail))
    println(header2_crayon(header_detail2))
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
    
    # Print computed substation angles θ₁ and θ₂
    println()
    if slack_sub == 1
        # Sub 1 is slack (θ₁ = 0°), compute θ₂
        θ₁ = 0.0
        θ₂ = result[:angles_deg][2]
        
        theta1_str = Crayon(foreground = :green, bold = true)("θ₁ (Sub 1, slack)")
        theta2_str = Crayon(foreground = :yellow, bold = true)("θ₂ (Sub 2)")
        
        println("Computed Substation Angles:")
        @printf("  %s = %7.3f°\n", theta1_str, θ₁)
        @printf("  %s = %7.3f°\n", theta2_str, θ₂)
        
        delta_str = Crayon(foreground = :magenta, bold = true)("Δθ (θ₂ - θ₁)")
        @printf("  %s = %7.3f°\n", delta_str, θ₂ - θ₁)
        
    else  # slack_sub == 2
        # Sub 2 is slack (θ₂ = 0°), compute θ₁
        θ₁ = result[:angles_deg][1]
        θ₂ = 0.0
        
        theta1_str = Crayon(foreground = :yellow, bold = true)("θ₁ (Sub 1)")
        theta2_str = Crayon(foreground = :cyan, bold = true)("θ₂ (Sub 2, slack)")
        
        println("Computed Substation Angles:")
        @printf("  %s = %7.3f°\n", theta1_str, θ₁)
        @printf("  %s = %7.3f°\n", theta2_str, θ₂)
        
        delta_str = Crayon(foreground = :magenta, bold = true)("Δθ (θ₂ - θ₁)")
        @printf("  %s = %7.3f°\n", delta_str, θ₂ - θ₁)
    end
    
    # Print power dispatch summary
    println()
    p1_crayon = Crayon(foreground = :light_green)
    p2_crayon = Crayon(foreground = :light_cyan)
    
    println("Power Dispatches:")
    p1_val = round(result[:P_1j] * data[:kVA_B], digits=1)
    q1_val = round(result[:Q_1j] * data[:kVA_B], digits=1)
    p2_val = round(result[:P_2j] * data[:kVA_B], digits=1)
    q2_val = round(result[:Q_2j] * data[:kVA_B], digits=1)
    
    println("  ", p1_crayon("Substation 1:"), " P = $(p1_val) kW, Q = $(q1_val) kVAr")
    println("  ", p2_crayon("Substation 2:"), " P = $(p2_val) kW, Q = $(q2_val) kVAr")
    
    # Print cost breakdown
    cost1 = p1_val * data[:C_1_dollar_per_kWh]
    cost2 = p2_val * data[:C_2_dollar_per_kWh]
    total_cost = cost1 + cost2
    
    cost_crayon = Crayon(foreground = :light_red, bold = true)
    println()
    println("Cost Breakdown:")
    println("  Sub 1 cost: \$$(round(cost1, digits=2)) (@\$$(data[:C_1_dollar_per_kWh])/kWh)")
    println("  Sub 2 cost: \$$(round(cost2, digits=2)) (@\$$(data[:C_2_dollar_per_kWh])/kWh)")
    println("  ", cost_crayon("Total:"), " \$$(round(total_cost, digits=2))")
end

# Now print summary comparison table (TABLE I)
println("\n" * "="^80)
title_crayon = Crayon(foreground = :light_blue, bold = true)
println(title_crayon("TABLE I: COMPARISON OF SLACK BUS CONFIGURATIONS"))
println("="^80)

# Print header
header_crayon = Crayon(foreground = :white, bold = true)
header = @sprintf("%-15s | %-20s | %-20s", "Slack Bus", "Total Angle", "Operational")
header2 = @sprintf("%-15s | %-20s | %-20s", "(Substation)", "Deviation* (°)", "Cost (\$)")
separator = "-"^length(header)

println(separator)
println(header_crayon(header))
println(header_crayon(header2))
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

println("\n" * "="^80)
success_crayon = Crayon(foreground = :light_green, bold = true)
info_crayon = Crayon(foreground = :light_blue)
println()
println(success_crayon("✓ Script complete!"), " Slack bus configuration analysis finished.")
println(info_crayon("ℹ For OpenDSS validation, see:"), " envs/multi_poi/archive/opendss_validation.jl")
println("="^80)

# Uncomment below to run OpenDSS validation:
# include("envs/multi_poi/archive/opendss_validation.jl")
# validation_results = validate_opf_with_opendss(modelDict, data)