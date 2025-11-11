# small2poi_mpopf.jl - Multi-Period Optimal Power Flow for 2-POI system
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
using Plots

# Optional: Uncomment if running OpenDSS validation
# import OpenDSSDirect as dss
# ----------------------------------------------------------------------

# ============================================================================
# SYSTEM PARAMETERS
# ============================================================================

data = Dict()
V_1_pu = 1.05
delta_1_deg = 0.0
V_2_pu = 1.05
delta_2_deg = 0.0

# Network impedances (ohms)
r1_ohm = 0.025; r2_ohm = 0.075
x1_ohm = 0.025; x2_ohm = 0.025

# Base load (peak values for scaling)
P_L_base_kW = 500
Q_L_base_kW = P_L_base_kW * 0.75

# System base values
kVA_B = 1000
kV_B = 11.547 # 20kV line-to-line, 11.547kV line-to-neutral

# Power limits for each substation (kW)
P_1_max_kW = 350.0
P_2_max_kW = 350.0

# ============================================================================
# TIME HORIZON PARAMETERS
# ============================================================================

T = 24  # Number of time periods
delta_t_h = 24.0 / T  # Time step duration in hours

# Time-varying load profile (sinusoidal, similar to tadmm_copper_plate)
LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2

# Time-varying energy cost ($/kWh)
LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2

data = Dict(
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
    :LoadShapeCost => LoadShapeCost
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
    Vminpu = 0.90
    Vmaxpu = 1.10
    
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
println("Cost variation: \$$(round(minimum(data[:LoadShapeCost]), digits=3)) to \$$(round(maximum(data[:LoadShapeCost]), digits=3)) per kWh")
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


"""
    solve_two_poi_mpopf(data)

Solve multi-period optimal power flow for 2-POI system.

Returns Dict with time-indexed arrays for all decision variables.
"""
function solve_two_poi_mpopf(data)
    # Extract time set and parameters
    Tset = data[:Tset]
    T = data[:T]
    Δt = data[:delta_t_h]
    P_BASE = data[:kVA_B]
    
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
    
    # Extract parameters (time-independent)
    r1j = data[:r1_pu]
    x1j = data[:x1_pu]
    r2j = data[:r2_pu]
    x2j = data[:x2_pu]
    v_1 = data[:V_1_pu]^2
    v_2 = data[:V_2_pu]^2
    Vminpu = data[:Vminpu]
    Vmaxpu = data[:Vmaxpu]

    # ========== CONSTRAINTS (for each time period) ==========
    for t in Tset
        # Load at time t
        P_L_t = data[:P_L_pu][t]
        Q_L_t = data[:Q_L_pu][t]
        
        # 1. Real power balance at load
        @constraint(model, (P_1j[t] - r1j * l_1j[t]) + (P_2j[t] - r2j * l_2j[t]) == P_L_t)
        
        # 2. Reactive power balance at load
        @constraint(model, (Q_1j[t] - x1j * l_1j[t]) + (Q_2j[t] - x2j * l_2j[t]) == Q_L_t)
        
        # 3. KVL constraint for line 1
        @constraint(model, v_j[t] == v_1 - 2 * (r1j * P_1j[t] + x1j * Q_1j[t]) + (r1j^2 + x1j^2) * l_1j[t])
        
        # 4. KVL constraint for line 2
        @constraint(model, v_j[t] == v_2 - 2 * (r2j * P_2j[t] + x2j * Q_2j[t]) + (r2j^2 + x2j^2) * l_2j[t])
        
        # 5. SOC relaxation for line 1
        @constraint(model, P_1j[t]^2 + Q_1j[t]^2 >= v_1 * l_1j[t])
        
        # 6. SOC relaxation for line 2
        @constraint(model, P_2j[t]^2 + Q_2j[t]^2 >= v_2 * l_2j[t])
        
        # 7. Voltage limits
        @constraint(model, Vminpu^2 <= v_j[t] <= Vmaxpu^2)
    end

    # ========== OBJECTIVE ==========
    # Total energy cost over all time periods
    # Cost is time-varying from LoadShapeCost
    LoadShapeCost = data[:LoadShapeCost]
    
    @expression(model, energy_cost_total,
        sum(LoadShapeCost[t] * (P_1j[t] + P_2j[t]) * P_BASE * Δt for t in Tset))
    
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
        # Store time-indexed arrays
        :P_1j => [value(P_1j[t]) for t in Tset],
        :Q_1j => [value(Q_1j[t]) for t in Tset],
        :l_1j => [value(l_1j[t]) for t in Tset],
        :P_2j => [value(P_2j[t]) for t in Tset],
        :Q_2j => [value(Q_2j[t]) for t in Tset],
        :l_2j => [value(l_2j[t]) for t in Tset],
        :v_j => [value(v_j[t]) for t in Tset]
    )
    
    return result
end

println("\n" * "="^80)
println("NETWORK TOPOLOGY")
println("="^80)
println("3-bus system: 2 substations (nodes 0, 1) + 1 load (node 2)")
println("Links: 0→2 (link 1), 1→2 (link 2)")
println("="^80)

# ==================================================================================
# SOLVE MPOPF
# ==================================================================================

println("\n" * "="^80)
println(Crayon(foreground = :cyan, bold = true)("SOLVING MULTI-PERIOD OPF"))
println("="^80)

sol_mpopf = solve_two_poi_mpopf(data)

# ==================================================================================
# REPORT RESULTS
# ==================================================================================

if sol_mpopf[:status] == MOI.OPTIMAL || sol_mpopf[:status] == MOI.LOCALLY_SOLVED
    println("\n" * "="^80)
    println(Crayon(foreground = :green, bold = true)("MPOPF SOLUTION SUMMARY"))
    println("="^80)
    
    # Convert to physical units
    P_BASE = data[:kVA_B]
    Δt = data[:delta_t_h]
    T = data[:T]
    
    # Extract arrays
    P_1j_kW = sol_mpopf[:P_1j] .* P_BASE
    P_2j_kW = sol_mpopf[:P_2j] .* P_BASE
    Q_1j_kVAr = sol_mpopf[:Q_1j] .* P_BASE
    Q_2j_kVAr = sol_mpopf[:Q_2j] .* P_BASE
    v_j_pu = sqrt.(sol_mpopf[:v_j])
    
    # Total substation power
    P_total_kW = P_1j_kW .+ P_2j_kW
    
    # Compute energy costs
    energy_costs_t = data[:LoadShapeCost] .* P_total_kW .* Δt
    total_cost = sum(energy_costs_t)
    
    println("\n--- OBJECTIVE VALUE ---")
    @printf "Total Cost: \$%.2f\n" total_cost
    @printf "Average Cost per hour: \$%.2f/h\n" (total_cost / (T * Δt))
    
    println("\n--- POWER SUMMARY ---")
    @printf "Substation 1 Power (kW): min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_1j_kW) maximum(P_1j_kW) mean(P_1j_kW)
    @printf "Substation 2 Power (kW): min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_2j_kW) maximum(P_2j_kW) mean(P_2j_kW)
    @printf "Total Power (kW): min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_total_kW) maximum(P_total_kW) mean(P_total_kW)
    
    println("\n--- VOLTAGE SUMMARY ---")
    @printf "Load Voltage (pu): min=%.4f, max=%.4f, avg=%.4f\n" minimum(v_j_pu) maximum(v_j_pu) mean(v_j_pu)
    
    println("\n--- ENERGY SUMMARY ---")
    total_energy_kWh = sum(P_total_kW) * Δt
    @printf "Total Energy Delivered: %.2f kWh\n" total_energy_kWh
    @printf "Sub 1 Energy: %.2f kWh (%.1f%%)\n" (sum(P_1j_kW) * Δt) (100 * sum(P_1j_kW) / sum(P_total_kW))
    @printf "Sub 2 Energy: %.2f kWh (%.1f%%)\n" (sum(P_2j_kW) * Δt) (100 * sum(P_2j_kW) / sum(P_total_kW))
    
    println("="^80)
    
    # ==================================================================================
    # PLOT RESULTS
    # ==================================================================================
    
    println("\n" * "="^80)
    println(Crayon(foreground = :magenta, bold = true)("GENERATING PLOTS"))
    println("="^80)
    
    time_hours = collect(1:T) .* Δt
    
    # Plot 1: Load profile and substation power dispatch
    p1 = plot(time_hours, P_total_kW, label="Total Load", linewidth=2, 
              xlabel="Time (hours)", ylabel="Power (kW)", title="Load and Power Dispatch",
              legend=:topright, grid=true)
    plot!(p1, time_hours, P_1j_kW, label="Substation 1", linewidth=2, linestyle=:dash)
    plot!(p1, time_hours, P_2j_kW, label="Substation 2", linewidth=2, linestyle=:dashdot)
    
    # Plot 2: Energy cost profile
    p2 = plot(time_hours, data[:LoadShapeCost], label="Energy Cost", linewidth=2,
              xlabel="Time (hours)", ylabel="Cost (\$/kWh)", title="Time-Varying Energy Cost",
              legend=:topright, grid=true, color=:red)
    
    # Plot 3: Load voltage profile
    p3 = plot(time_hours, v_j_pu, label="Load Bus Voltage", linewidth=2,
              xlabel="Time (hours)", ylabel="Voltage (pu)", title="Load Bus Voltage Profile",
              legend=:topright, grid=true, color=:green)
    hline!(p3, [data[:Vminpu]], label="V_min", linestyle=:dot, color=:red)
    hline!(p3, [data[:Vmaxpu]], label="V_max", linestyle=:dot, color=:red)
    
    # Plot 4: Reactive power
    p4 = plot(time_hours, Q_1j_kVAr, label="Substation 1", linewidth=2,
              xlabel="Time (hours)", ylabel="Reactive Power (kVAr)", title="Reactive Power Dispatch",
              legend=:topright, grid=true)
    plot!(p4, time_hours, Q_2j_kVAr, label="Substation 2", linewidth=2, linestyle=:dash)
    
    # Combined plot
    plot_combined = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 800))
    
    # Save plot
    mkpath("plots")
    savefig(plot_combined, "plots/mpopf_2poi_T$(T)_results.png")
    println(Crayon(foreground = :green)("✓ Plots saved to plots/mpopf_2poi_T$(T)_results.png"))
    println("="^80)
    
else
    println("\n" * "="^80)
    println(Crayon(foreground = :red, bold = true)("MPOPF FAILED"))
    println("="^80)
    println("Status: ", sol_mpopf[:status])
    println("="^80)
end

println("\n" * "="^80)
println(Crayon(foreground = :green, bold = true)("✓ MPOPF SCRIPT COMPLETE"))
println("="^80)