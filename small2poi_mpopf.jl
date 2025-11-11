# small2poi_mpopf.jl - Multi-Period Optimal Power Flow for 2-POI system
# 
# Key feature: 'same_voltage_levels' parameter (line ~30)
#   - true (default): All substations have fixed voltages (1.05 pu)
#   - false: Non-slack substation voltage becomes a decision variable (optimized)
#            This allows different optimal voltage setpoints and breaks angle symmetry
#
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
using Statistics
using Plots  # For input curve plotting

# Optional: Uncomment if running OpenDSS validation
# import OpenDSSDirect as dss
# ----------------------------------------------------------------------

# ============================================================================
# SYSTEM PARAMETERS
# ============================================================================

data = Dict()

# User-defined parameter: Should substations have same voltage levels?
# - true: Non-slack substations have fixed voltage (V_1_pu, V_2_pu)
# - false: Non-slack substation voltages become decision variables
# same_voltage_levels = true
same_voltage_levels = false

# Plotting options
showPlots = false  # Set to true to display plots in GUI
savePlots = true   # Set to true to save plots to file

V_1_pu = 1.05
delta_1_deg = 0.0
V_2_pu = 1.05
delta_2_deg = 0.0

# Network impedances (ohms)
r1_ohm = 0.25; r2_ohm = 0.75
x1_ohm = 0.25; x2_ohm = 0.25

# Base load (peak values for scaling)
P_L_base_kW = 5000
Q_L_base_kW = P_L_base_kW * 0.75

# System base values
kVA_B = 1000
kV_B = 11.547 # 20kV line-to-line, 11.547kV line-to-neutral

# Power limits for each substation (kW)
P_1_max_kW = P_L_base_kW * 0.7
P_2_max_kW = P_L_base_kW * 0.7

# ============================================================================
# TIME HORIZON PARAMETERS
# ============================================================================

T = 24  # Number of time periods
delta_t_h = 24.0 / T  # Time step duration in hours

# Time-varying load profile (sinusoidal, similar to tadmm_copper_plate)
LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2

# Time-varying energy cost ($/kWh) - different for each substation
# Substation 1: Peak during early hours
LoadShapeCost_1 = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .+ π/4) .+ 1) ./ 2

# Substation 2: Peak during later hours (shifted by π/2)
LoadShapeCost_2 = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T) .- π/4) .+ 1) ./ 2

data = Dict(
    :same_voltage_levels => same_voltage_levels,
    :V_1_pu => V_1_pu,
    :delta_1_deg => delta_1_deg,
    :V_2_pu => V_2_pu,
    :delta_2_deg => delta_2_deg,
    :r1_ohm => r1_ohm,
    :x1_ohm => x1_ohm,
    :r2_ohm => r2_ohm,
    :x2_ohm => x2_ohm,
    :P_L_base_kW => P_L_base_kW,
    :Q_L_base_kW => Q_L_base_kW,
    :kVA_B => kVA_B,
    :kV_B => kV_B,
    :P_1_max_kW => P_1_max_kW,
    :P_2_max_kW => P_2_max_kW,
    :T => T,
    :delta_t_h => delta_t_h,
    :LoadShapeLoad => LoadShapeLoad,
    :LoadShapeCost_1 => LoadShapeCost_1,
    :LoadShapeCost_2 => LoadShapeCost_2
)

function process_data!(data)
    # Base values
    S_base_kVA = data[:kVA_B]
    S_base_MVA = S_base_kVA / 1000  # Convert kVA to MVA
    V_base = data[:kV_B]
    T = data[:T]

    # Per-unit impedances: Z_pu = Z_ohm * (MVA_base / kV_base^2)
    Z_base_ohm = (V_base^2) / S_base_MVA  # Ohms
    r1_pu = data[:r1_ohm] / Z_base_ohm
    x1_pu = data[:x1_ohm] / Z_base_ohm
    r2_pu = data[:r2_ohm] / Z_base_ohm
    x2_pu = data[:x2_ohm] / Z_base_ohm

    # Time-varying load in per-unit [T vector]
    P_L_base_pu = data[:P_L_base_kW] / S_base_kVA
    Q_L_base_pu = data[:Q_L_base_kW] / S_base_kVA
    
    # Create time-indexed load arrays
    P_L_pu = [P_L_base_pu * data[:LoadShapeLoad][t] for t in 1:T]
    Q_L_pu = [Q_L_base_pu * data[:LoadShapeLoad][t] for t in 1:T]

    # Angles in radians
    delta_1_rad = data[:delta_1_deg] * pi / 180
    delta_2_rad = data[:delta_2_deg] * pi / 180

    # Voltage limits
    Vminpu = 0.95
    Vmaxpu = 1.05
    
    # Power limits (pu)
    P_1_max_pu = data[:P_1_max_kW] / S_base_kVA
    P_2_max_pu = data[:P_2_max_kW] / S_base_kVA
    
    # Time set
    Tset = 1:T

    # Update data dict with all required fields
    data[:P_L_pu] = P_L_pu  # Now a time-indexed vector
    data[:Q_L_pu] = Q_L_pu  # Now a time-indexed vector
    data[:r1_pu] = r1_pu
    data[:x1_pu] = x1_pu
    data[:r2_pu] = r2_pu
    data[:x2_pu] = x2_pu
    data[:delta_1_rad] = delta_1_rad
    data[:delta_2_rad] = delta_2_rad
    data[:Vminpu] = Vminpu
    data[:Vmaxpu] = Vmaxpu
    data[:P_1_max_pu] = P_1_max_pu
    data[:P_2_max_pu] = P_2_max_pu
    data[:Tset] = Tset

end

process_data!(data)

println("\n" * "="^80)
println("MULTI-PERIOD OPTIMAL POWER FLOW (MPOPF)")
println("="^80)
println("Time periods: T = $(data[:T])")
println("Time step: Δt = $(data[:delta_t_h]) hours")
println("Base load: P = $(data[:P_L_base_kW]) kW, Q = $(data[:Q_L_base_kW]) kVAr")
println("Load variation: $(round(minimum(data[:LoadShapeLoad]), digits=3)) to $(round(maximum(data[:LoadShapeLoad]), digits=3))")
println("Cost variation (Sub 1): \$$(round(minimum(data[:LoadShapeCost_1]), digits=3)) to \$$(round(maximum(data[:LoadShapeCost_1]), digits=3)) per kWh")
println("Cost variation (Sub 2): \$$(round(minimum(data[:LoadShapeCost_2]), digits=3)) to \$$(round(maximum(data[:LoadShapeCost_2]), digits=3)) per kWh")
voltage_config_str = data[:same_voltage_levels] ? "FIXED (all substations at 1.05 pu)" : "VARIABLE (non-slack optimized)"
println("Voltage configuration: ", Crayon(foreground = :light_magenta, bold = true)(voltage_config_str))
println("="^80)


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
    compute_beta_angles(data, opf_result, slack_node, v_nonslack=nothing)

Compute the angle difference β_ij = ∠(v_i - z*_ij S_ij) for each link (i,j) ∈ E
after solving the OPF, as per equation (27) in the paper.

Arguments:
- data: system data dictionary
- opf_result: Dict with OPF solution (:P_1j, :Q_1j, :P_2j, :Q_2j, :v_j, etc.)
- slack_node: which node is slack (0 or 1)
- v_nonslack: voltage magnitude squared for non-slack substation (if variable)

Returns Dict with:
- :beta - vector of angle differences for all links (in radians)
- :beta_deg - vector of angle differences for all links (in degrees)
"""
function compute_beta_angles(data, opf_result, slack_node::Int, v_nonslack=nothing)
    # Compute β_ij = ∠(v_i - z*_ij S_ij) for each link (i,j) ∈ E
    # as per equation (27) in the paper
    
    # Extract voltage magnitudes
    same_voltage_levels = data[:same_voltage_levels]
    
    if same_voltage_levels
        # Fixed voltages
        V_0_pu = data[:V_1_pu]  # Node 0 (substation 1)
        V_1_pu = data[:V_2_pu]  # Node 1 (substation 2)
    else
        # Variable voltage for non-slack substation
        if slack_node == 0
            V_0_pu = data[:V_1_pu]  # Slack (fixed)
            V_1_pu = sqrt(v_nonslack)  # Non-slack (variable)
        else  # slack_node == 1
            V_0_pu = sqrt(v_nonslack)  # Non-slack (variable)
            V_1_pu = data[:V_2_pu]  # Slack (fixed)
        end
    end
    
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


"""
    solve_two_poi_mpopf(data, slack_node)

Solve multi-period optimal power flow for 2-POI system with specified slack bus.

Arguments:
- data: system data dictionary
- slack_node: which node is slack (0 or 1, corresponding to substations 1 and 2)

Returns Dict with time-indexed arrays for all decision variables.
"""
function solve_two_poi_mpopf(data, slack_node::Int)
    # Extract time set and parameters
    Tset = data[:Tset]
    T = data[:T]
    Δt = data[:delta_t_h]
    P_BASE = data[:kVA_B]
    same_voltage_levels = data[:same_voltage_levels]
    
    # Compute network topology matrices for this slack configuration
    println("  Computing network matrices for slack bus = node $slack_node...")
    network_matrices = compute_network_matrices(data, slack_node)
    
    # Create model
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "max_iter", 3000)
    set_optimizer_attribute(model, "tol", 1e-6)
    
    # ========== VARIABLES (time-indexed) ==========
    @variable(model, 0 <= P_1j[t in Tset] <= data[:P_1_max_pu])
    @variable(model, Q_1j[t in Tset])
    @variable(model, l_1j[t in Tset] >= 0)
    @variable(model, 0 <= P_2j[t in Tset] <= data[:P_2_max_pu])
    @variable(model, Q_2j[t in Tset])
    @variable(model, l_2j[t in Tset] >= 0)
    @variable(model, v_j[t in Tset])
    
    # Voltage variables for non-slack substations (if same_voltage_levels = false)
    if !same_voltage_levels
        if slack_node == 0
            # Node 0 is slack, node 1 voltage becomes decision variable
            @variable(model, data[:Vminpu]^2 <= v_1[t in Tset] <= data[:Vmaxpu]^2)
        else  # slack_node == 1
            # Node 1 is slack, node 0 voltage becomes decision variable
            @variable(model, data[:Vminpu]^2 <= v_0[t in Tset] <= data[:Vmaxpu]^2)
        end
    end
    
    # Extract parameters (time-independent)
    r1j = data[:r1_pu]
    x1j = data[:x1_pu]
    r2j = data[:r2_pu]
    x2j = data[:x2_pu]
    Vminpu = data[:Vminpu]
    Vmaxpu = data[:Vmaxpu]
    
    # Set substation voltages based on same_voltage_levels flag
    if same_voltage_levels
        # Fixed voltages (current behavior)
        v_1 = data[:V_1_pu]^2
        v_2 = data[:V_2_pu]^2
    else
        # Voltage depends on slack configuration - will be used in constraints below
        # (handled differently in constraint section)
    end

    # ========== CONSTRAINTS (for each time period) ==========
    for t in Tset
        # Load at time t
        P_L_t = data[:P_L_pu][t]
        Q_L_t = data[:Q_L_pu][t]
        
        # 1. Real power balance at load
        @constraint(model, (P_1j[t] - r1j * l_1j[t]) + (P_2j[t] - r2j * l_2j[t]) == P_L_t)
        
        # 2. Reactive power balance at load
        @constraint(model, (Q_1j[t] - x1j * l_1j[t]) + (Q_2j[t] - x2j * l_2j[t]) == Q_L_t)
        
        # 3. KVL constraint for line 1 (node 0 → node 2)
        # 4. KVL constraint for line 2 (node 1 → node 2)
        # 5. SOC relaxation for line 1
        # 6. SOC relaxation for line 2
        if same_voltage_levels
            # Fixed voltages (current behavior)
            v_1_val = data[:V_1_pu]^2
            v_2_val = data[:V_2_pu]^2
            
            @constraint(model, v_j[t] == v_1_val - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
            @constraint(model, v_j[t] == v_2_val - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
            @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_1_val * l_1j[t])
            @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_2_val * l_2j[t])
        else
            # Variable voltages for non-slack substation
            if slack_node == 0
                # Node 0 (Sub 1) is slack with fixed voltage
                # Node 1 (Sub 2) has variable voltage v_1[t]
                v_0_val = data[:V_1_pu]^2  # Slack voltage (fixed)
                
                @constraint(model, v_j[t] == v_0_val - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
                @constraint(model, v_j[t] == v_1[t] - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
                @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_0_val * l_1j[t])
                @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_1[t] * l_2j[t])
            else  # slack_node == 1
                # Node 1 (Sub 2) is slack with fixed voltage
                # Node 0 (Sub 1) has variable voltage v_0[t]
                v_1_val = data[:V_2_pu]^2  # Slack voltage (fixed)
                
                @constraint(model, v_j[t] == v_0[t] - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
                @constraint(model, v_j[t] == v_1_val - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
                @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_0[t] * l_1j[t])
                @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_1_val * l_2j[t])
            end
        end
        
        # 7. Voltage limits (load bus always constrained)
        @constraint(model, Vminpu^2 <= v_j[t] <= Vmaxpu^2)
    end

    # ========== OBJECTIVE ==========
    # Total energy cost over all time periods
    # Cost is time-varying and different for each substation
    LoadShapeCost_1 = data[:LoadShapeCost_1]
    LoadShapeCost_2 = data[:LoadShapeCost_2]
    
    @expression(model, energy_cost_total,
        sum(LoadShapeCost_1[t] * P_1j[t] * P_BASE * Δt + 
            LoadShapeCost_2[t] * P_2j[t] * P_BASE * Δt for t in Tset))
    
    @objective(model, Min, energy_cost_total)

    # ========== SOLVE ==========
    println("\nSolving MPOPF with T=$(T) time periods...")
    solve_start = time()
    optimize!(model)
    solve_time_val = time() - solve_start
    
    # Check termination status
    status = termination_status(model)
    
    if status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL
        status_crayon = Crayon(foreground = :green, bold = true)
        println(status_crayon("✓ MPOPF Solved: $status"))
        println("Solve time: $(round(solve_time_val, digits=2)) seconds")
    else
        status_crayon = Crayon(foreground = :red, bold = true)
        println(status_crayon("✗ MPOPF Failed: $status"))
    end

    # ========== COLLECT RESULTS ==========
    result = Dict(
        :status => status,
        :objective => has_values(model) ? objective_value(model) : NaN,
        :solve_time => solve_time_val,
        :slack_node => slack_node,
        :same_voltage_levels => same_voltage_levels,
        # Store time-indexed arrays
        :P_1j => [value(P_1j[t]) for t in Tset],
        :Q_1j => [value(Q_1j[t]) for t in Tset],
        :l_1j => [value(l_1j[t]) for t in Tset],
        :P_2j => [value(P_2j[t]) for t in Tset],
        :Q_2j => [value(Q_2j[t]) for t in Tset],
        :l_2j => [value(l_2j[t]) for t in Tset],
        :v_j => [value(v_j[t]) for t in Tset],
        # Store network matrices for angle computation
        :network_matrices => network_matrices
    )
    
    # Store non-slack substation voltages if variable
    if !same_voltage_levels
        if slack_node == 0
            # Node 1 voltage was variable
            result[:v_nonslack] = [value(v_1[t]) for t in Tset]
        else  # slack_node == 1
            # Node 0 voltage was variable
            result[:v_nonslack] = [value(v_0[t]) for t in Tset]
        end
    end
    
    return result
end

"""
    analyze_slack_configuration_mpopf(slack_node, data)

Solve MPOPF with specified slack bus and compute angle statistics across time.

Returns Dict with:
- :objective - total cost over all time periods
- :angles_deg - matrix of angles [substation × time]
- :median_angles_deg - median angle for each non-slack substation
- :angle_ranges_deg - [min, max] angle for each non-slack substation
- :rmse_values_deg - RMSE (temporal variation) for each non-slack substation
- :total_deviation_deg - grand RMS of all RMSE values
- Other power flow results
"""
function analyze_slack_configuration_mpopf(slack_node::Int, data)
    # Solve MPOPF
    result = solve_two_poi_mpopf(data, slack_node)
    
    if !(result[:status] == MOI.OPTIMAL || result[:status] == MOI.LOCALLY_SOLVED)
        @warn "MPOPF failed for slack_node=$slack_node"
        return result
    end
    
    # Extract time parameters
    T = data[:T]
    Tset = data[:Tset]
    
    # Compute angles for each time period
    # angles_matrix[substation_idx, time_idx]
    angles_matrix = Dict()
    
    for t in Tset
        # Package results for this time period
        opf_result_t = Dict(
            :P_1j => result[:P_1j][t],
            :Q_1j => result[:Q_1j][t],
            :P_2j => result[:P_2j][t],
            :Q_2j => result[:Q_2j][t],
            :v_j => result[:v_j][t]
        )
        
        # Get non-slack voltage if variable
        v_nonslack_t = haskey(result, :v_nonslack) ? result[:v_nonslack][t] : nothing
        
        # Compute β angles
        beta_result = compute_beta_angles(data, opf_result_t, slack_node, v_nonslack_t)
        β = beta_result[:beta]
        
        # Get network matrices
        B_T = result[:network_matrices][:B_T]
        
        # Solve for angles: θ_★ = B_T^(-1) * β
        θ_star = B_T \ β
        
        # Project to (-π, π]
        θ_star = mod.(θ_star .+ π, 2π) .- π
        
        # Map to substation indices
        if slack_node == 0
            # Removed node 0, so θ_star[1] corresponds to node 1 (substation 2)
            sub_idx = 2
            angle_deg = rad2deg(θ_star[1])
        else  # slack_node == 1
            # Removed node 1, so θ_star[1] corresponds to node 0 (substation 1)
            sub_idx = 1
            angle_deg = rad2deg(θ_star[1])
        end
        
        # Store angle
        if !haskey(angles_matrix, sub_idx)
            angles_matrix[sub_idx] = Float64[]
        end
        push!(angles_matrix[sub_idx], angle_deg)
    end
    
    # Compute statistics across time for each non-slack substation
    median_angles_deg = Dict()
    angle_ranges_deg = Dict()
    rmse_values_deg = Dict()
    
    for (sub, angles_vec) in angles_matrix
        # Median
        median_angles_deg[sub] = median(angles_vec)
        
        # Range [min, max]
        angle_ranges_deg[sub] = [minimum(angles_vec), maximum(angles_vec)]
        
        # RMSE (temporal variation around median)
        rmse_values_deg[sub] = sqrt(mean((angles_vec .- median_angles_deg[sub]).^2))
    end
    
    # Total angle deviation (grand RMS)
    if length(rmse_values_deg) > 0
        total_deviation = sqrt(sum(v^2 for v in values(rmse_values_deg)) / length(rmse_values_deg))
    else
        total_deviation = 0.0
    end
    
    # Add statistics to result
    result[:angles_matrix] = angles_matrix
    result[:median_angles_deg] = median_angles_deg
    result[:angle_ranges_deg] = angle_ranges_deg
    result[:rmse_values_deg] = rmse_values_deg
    result[:total_deviation_deg] = total_deviation
    
    return result
end


println("\n" * "="^80)
println("NETWORK TOPOLOGY")
println("="^80)
println("3-bus system: 2 substations (nodes 0, 1) + 1 load (node 2)")
println("Links: 0→2 (link 1), 1→2 (link 2)")
println("="^80)

# ==================================================================================
# SOLVE MPOPF FOR BOTH SLACK BUS CONFIGURATIONS
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :cyan, bold = true)("SOLVING MULTI-PERIOD OPF FOR BOTH SLACK CONFIGURATIONS"))
println("="^80)

# Analyze each slack bus configuration
results_by_slack = Dict()

for slack_sub in [1, 2]
    slack_node = slack_sub - 1  # Convert substation number to node index
    println("\n" * "-"^80)
    println(Crayon(foreground = :yellow, bold = true)("Solving MPOPF with Substation $slack_sub (node $slack_node) as slack bus..."))
    println("-"^80)
    results_by_slack[slack_sub] = analyze_slack_configuration_mpopf(slack_node, data)
end

# ==================================================================================
# PRINT COMPARISON TABLES (TABLE II FIRST, THEN TABLE I)
# ==================================================================================

# ==================================================================================
# DETAILED ANGLE COORDINATION (TABLE II) FOR EACH CONFIGURATION
# ==================================================================================

for slack_sub in sort(collect(keys(results_by_slack)))
    result = results_by_slack[slack_sub]
    
    println("\n" * "="^80)
    title2_crayon = Crayon(foreground = :light_magenta, bold = true)
    println(title2_crayon("TABLE II: ANGLE COORDINATION DETAILS (SLACK: SUBSTATION $slack_sub)"))
    println("="^80)
    
    header2_crayon = Crayon(foreground = :white, bold = true)
    header_detail = @sprintf("%-15s | %-20s | %-25s | %-15s", 
                            "Substation", "Median", "Angle Range", "Temporal RMSE")
    header_detail2 = @sprintf("%-15s | %-20s | %-25s | %-15s",
                             "(Non-Slack)", "Angle (°)", "(°)", "(°)")
    sep_detail = "-"^length(header_detail)
    
    println(sep_detail)
    println(header2_crayon(header_detail))
    println(header2_crayon(header_detail2))
    println(sep_detail)
    
    # Get non-slack substations
    non_slack_subs = sort(collect(keys(result[:angles_matrix])))
    
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
    
    # Print computed substation angles (showing median and range)
    println()
    if slack_sub == 1
        # Sub 1 is slack (θ₁ = 0°), compute θ₂
        θ₁_median = 0.0
        θ₂_median = result[:median_angles_deg][2]
        θ₂_range = result[:angle_ranges_deg][2]
        
        theta1_str = Crayon(foreground = :green, bold = true)("θ₁ (Sub 1, slack)")
        theta2_str = Crayon(foreground = :yellow, bold = true)("θ₂ (Sub 2)")
        
        println("Computed Substation Angles (across time):")
        @printf("  %s = %7.3f° (fixed)\n", theta1_str, θ₁_median)
        @printf("  %s = %7.3f° (median), range [%7.3f°, %7.3f°]\n", 
                theta2_str, θ₂_median, θ₂_range[1], θ₂_range[2])
        
        delta_str = Crayon(foreground = :magenta, bold = true)("Δθ (θ₂ - θ₁)")
        @printf("  %s = %7.3f° (median)\n", delta_str, θ₂_median - θ₁_median)
        
    else  # slack_sub == 2
        # Sub 2 is slack (θ₂ = 0°), compute θ₁
        θ₁_median = result[:median_angles_deg][1]
        θ₁_range = result[:angle_ranges_deg][1]
        θ₂_median = 0.0
        
        theta1_str = Crayon(foreground = :yellow, bold = true)("θ₁ (Sub 1)")
        theta2_str = Crayon(foreground = :cyan, bold = true)("θ₂ (Sub 2, slack)")
        
        println("Computed Substation Angles (across time):")
        @printf("  %s = %7.3f° (median), range [%7.3f°, %7.3f°]\n",
                theta1_str, θ₁_median, θ₁_range[1], θ₁_range[2])
        @printf("  %s = %7.3f° (fixed)\n", theta2_str, θ₂_median)
        
        delta_str = Crayon(foreground = :magenta, bold = true)("Δθ (θ₂ - θ₁)")
        @printf("  %s = %7.3f° (median)\n", delta_str, θ₂_median - θ₁_median)
    end
    
    # Print bus voltages (time-averaged)
    println()
    v_crayon = Crayon(foreground = :light_yellow)
    v_var_crayon = Crayon(foreground = :light_magenta, bold = true)
    println("Bus Voltages (time-averaged):")
    
    same_v_levels = data[:same_voltage_levels]  # Local copy to avoid soft scope warning
    slack_node_val = result[:slack_node]  # Get slack node from result
    
    if same_v_levels
        # Fixed voltages
        v_0_pu = data[:V_1_pu]  # Substation 1 (node 0) - fixed
        v_1_pu = data[:V_2_pu]  # Substation 2 (node 1) - fixed
        v_2_pu = mean(sqrt.(result[:v_j]))  # Load bus (node 2) - averaged over time
        
        println("  ", v_crayon("Node 0 (Sub 1):"), " V = $(round(v_0_pu, digits=4)) pu (fixed)")
        println("  ", v_crayon("Node 1 (Sub 2):"), " V = $(round(v_1_pu, digits=4)) pu (fixed)")
        println("  ", v_crayon("Node 2 (Load): "), " V = $(round(v_2_pu, digits=4)) pu (avg)")
    else
        # Variable voltage for non-slack
        v_2_pu = mean(sqrt.(result[:v_j]))  # Load bus (node 2) - averaged over time
        
        if slack_node_val == 0
            # Node 0 is slack (fixed), node 1 is variable
            v_0_pu = data[:V_1_pu]
            v_1_pu_avg = mean(sqrt.(result[:v_nonslack]))
            v_1_pu_min = minimum(sqrt.(result[:v_nonslack]))
            v_1_pu_max = maximum(sqrt.(result[:v_nonslack]))
            
            println("  ", v_crayon("Node 0 (Sub 1):"), " V = $(round(v_0_pu, digits=4)) pu (slack, fixed)")
            println("  ", v_var_crayon("Node 1 (Sub 2):"), " V = $(round(v_1_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_1_pu_min, digits=4)), $(round(v_1_pu_max, digits=4))] pu")
        else  # slack_node_val == 1
            # Node 1 is slack (fixed), node 0 is variable
            v_1_pu = data[:V_2_pu]
            v_0_pu_avg = mean(sqrt.(result[:v_nonslack]))
            v_0_pu_min = minimum(sqrt.(result[:v_nonslack]))
            v_0_pu_max = maximum(sqrt.(result[:v_nonslack]))
            
            println("  ", v_var_crayon("Node 0 (Sub 1):"), " V = $(round(v_0_pu_avg, digits=4)) pu (avg, optimized)")
            println("                    Range: [$(round(v_0_pu_min, digits=4)), $(round(v_0_pu_max, digits=4))] pu")
            println("  ", v_crayon("Node 1 (Sub 2):"), " V = $(round(v_1_pu, digits=4)) pu (slack, fixed)")
        end
        
        println("  ", v_crayon("Node 2 (Load): "), " V = $(round(v_2_pu, digits=4)) pu (avg)")
    end
    
    # Print power dispatch summary (time-averaged)
    println()
    p1_crayon = Crayon(foreground = :light_green)
    p2_crayon = Crayon(foreground = :light_cyan)
    
    P_BASE = data[:kVA_B]
    println("Power Dispatches (time-averaged):")
    p1_avg = mean(result[:P_1j]) * P_BASE
    q1_avg = mean(result[:Q_1j]) * P_BASE
    p2_avg = mean(result[:P_2j]) * P_BASE
    q2_avg = mean(result[:Q_2j]) * P_BASE
    
    println("  ", p1_crayon("Substation 1:"), " P = $(round(p1_avg, digits=1)) kW, Q = $(round(q1_avg, digits=1)) kVAr")
    println("  ", p2_crayon("Substation 2:"), " P = $(round(p2_avg, digits=1)) kW, Q = $(round(q2_avg, digits=1)) kVAr")
    
    cost_crayon = Crayon(foreground = :light_red, bold = true)
    println()
    println("  ", cost_crayon("Total Energy Cost:"), " \$$(round(result[:objective], digits=2))")
    println("  Solve time: $(round(result[:solve_time], digits=2)) seconds")
end

# ==================================================================================
# SUMMARY COMPARISON TABLE (TABLE I)
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :light_blue, bold = true)("TABLE I: COMPARISON OF SLACK BUS CONFIGURATIONS"))
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
    
    total_cost = result[:objective]
    total_dev = result[:total_deviation_deg]
    
    # Format with colors
    if slack_sub == 1
        sub_str = Crayon(foreground = :green, bold = true)(@sprintf("%-15d", slack_sub))
    else
        sub_str = Crayon(foreground = :cyan, bold = true)(@sprintf("%-15d", slack_sub))
    end
    
    dev_str = @sprintf("%-20.3f", total_dev)
    cost_str = @sprintf("%-20.2f", total_cost)
    
    @printf("%s | %s | %s\n", sub_str, dev_str, cost_str)
end

println(separator)
println("*RMS of temporal angle variations across all non-slack substations")
println("="^80)

println("\n" * "="^80)
success_crayon = Crayon(foreground = :light_green, bold = true)
info_crayon = Crayon(foreground = :light_blue)
println(success_crayon("✓ MPOPF SCRIPT COMPLETE!"), " Analyzed both slack bus configurations.")
println(info_crayon("ℹ Tables I & II show angle coordination and cost comparison."))
println("="^80)

# ==================================================================================
# PLOT INPUT CURVES
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :magenta, bold = true)("GENERATING INPUT CURVE PLOTS"))
println("="^80)

# Load plotter
include(joinpath(@__DIR__, "envs", "multi_poi", "Plotter.jl"))

# Generate plot in proper location
plots_dir = joinpath(@__DIR__, "envs", "multi_poi", "processedData", "plots")
mkpath(plots_dir)
plot_filename = joinpath(plots_dir, "mpopf_input_curves.png")
plot_input_curves(data, showPlots=showPlots, savePlots=savePlots, filename=plot_filename)
print_curve_statistics(data)

println("="^80)