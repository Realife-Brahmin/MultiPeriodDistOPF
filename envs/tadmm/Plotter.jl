"""
Plotting utilities for tADMM Multi-Period OPF
Standalone plotting functions for battery actions and input curves
"""

using Plots
using Printf

"""
    plot_input_curves(data::Dict; showPlots::Bool=true, savePlots::Bool=false, filename::String="input_curves.png")

Plot LoadShapeLoad, LoadShapePV, and LoadShapeCost curves on the same figure.
Uses dark yellow for load, orange for PV, and dark green for cost.

# Arguments
- data: Dictionary containing :T, :LoadShapeLoad, :LoadShapePV, :LoadShapeCost, :kVA_B
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Plot object
"""
function plot_input_curves(data::Dict; showPlots::Bool=true, savePlots::Bool=false, filename::String="input_curves.png")
    T = data[:T]
    LoadShapeLoad = data[:LoadShapeLoad]
    LoadShapePV = data[:LoadShapePV]
    LoadShapeCost = data[:LoadShapeCost]
    
    # Prepare data for plotting
    time_steps = 1:T
    load_cost_cents = LoadShapeCost .* 100  # Convert from $/kWh to cents/kWh
    
    # Calculate y-axis limits for left axis (load and PV)
    left_min = -0.05
    left_max = 1.05
    right_min = floor(0.95 * minimum(load_cost_cents))
    right_max = ceil(1.05 * maximum(load_cost_cents))

    # Set theme and backend
    gr()
    theme(:mute)
    
    # Create main plot with LoadShape and PV (dark yellow and orange)
    p = plot(
        time_steps, LoadShapeLoad,
        dpi=600,
        label="Load Shape (λᵗ)",
        xlabel="Time Period (t)",
        ylabel="Normalized Profile [0-1]",
        legend=:topright,
        lw=4,
        color=:darkgoldenrod2,  # Dark yellow theme
        markershape=:square,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        ylims=(left_min, left_max),
        xticks=1:T,
        title="Input Curves: Load, PV, and Cost Profiles",
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        right_margin=10Plots.mm
    )

    # Add PV profile (orange theme)
    plot!(
        p, time_steps, LoadShapePV,
        label="PV Generation Profile",
        lw=4,
        color=:darkorange,  # Orange theme
        markershape=:circle,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0
    )

    # Add secondary y-axis for cost curve (dark green theme)
    ax2 = twinx()
    plot!(
        ax2, time_steps, load_cost_cents,
        label="Energy Cost (Cᵗ)",
        lw=4,
        color=:darkgreen,  # Dark green theme
        linestyle=:solid,
        markershape=:diamond,
        markersize=7,
        markerstrokecolor=:black,
        markerstrokewidth=2.0,
        ylabel="Cost [cents/kWh]",
        ylims=(right_min, right_max),
        legend=:bottomright,
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )

    # Show the plot if requested
    if showPlots
        display(p)
    end

    # Save the plot if requested
    if savePlots
        @printf "Saving input curves plot to: %s\n" filename
        savefig(p, filename)
    end
    
    return p
end


"""
    plot_battery_actions(solution::Dict, data::Dict, method_name::String; 
                        showPlots::Bool=true, savePlots::Bool=false, 
                        filename::String="battery_actions_lindistflow.png")

Plot battery charging/discharging power and SOC for LinDistFlow solution.
Creates separate charging/discharging bars and SOC plot.

# Arguments
- solution: Solution dictionary containing :P_B (per-unit), :B (per-unit SOC)
- data: Data dictionary containing battery parameters and bases
- method_name: Name of the solution method (e.g., "LinDistFlow", "LinDistFlow-Ipopt")
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Combined plot object with charging/discharging and SOC subplots
"""
function plot_battery_actions(solution::Dict, data::Dict, method_name::String; 
                               showPlots::Bool=true, savePlots::Bool=false, 
                               filename::String="battery_actions_lindistflow.png")
    
    T = data[:T]
    P_BASE = data[:kVA_B]
    E_BASE = P_BASE * 1.0  # 1 hour window
    time_steps = 1:T
    
    # Extract battery parameters
    Bset = data[:Bset]
    if isempty(Bset)
        @warn "No batteries in the system, skipping battery actions plot"
        return nothing
    end
    
    # For now, plot only the first battery (can be extended for multiple batteries)
    b = first(Bset)
    
    # Convert to physical units for plotting
    P_B_pu = nothing
    if haskey(solution, :P_B)
        P_B_var = solution[:P_B]
        # JuMP DenseAxisArray or Dict-like structure
        try
            P_B_pu = [P_B_var[b, t] for t in 1:T]
        catch e
            @warn "Failed to extract P_B values: $e"
            @warn "P_B type: $(typeof(P_B_var))"
            return nothing
        end
    else
        @warn "P_B not found in solution dictionary. Available keys: $(keys(solution))"
        return nothing
    end
    
    B_pu = nothing
    if haskey(solution, :B)
        B_var = solution[:B]
        # JuMP DenseAxisArray or Dict-like structure
        try
            B_pu = [B_var[b, t] for t in 1:T]
        catch e
            @warn "Failed to extract B values: $e"
            @warn "B type: $(typeof(B_var))"
            return nothing
        end
    else
        @warn "B not found in solution dictionary. Available keys: $(keys(solution))"
        return nothing
    end
    
    P_B_kW = P_B_pu .* P_BASE
    B_kWh = [data[:B0_pu][b] * E_BASE; B_pu .* E_BASE]  # Include initial SOC
    
    # Separate charging and discharging (P_B > 0 is discharging, P_B < 0 is charging)
    charging_power_kW = -min.(P_B_kW, 0.0)  # P_B < 0 made positive for green bars
    discharging_power_kW = max.(P_B_kW, 0.0)  # P_B > 0 kept positive for wine red bars
    
    # Calculate SOC percentages
    E_Rated_kWh = data[:B_R][b]
    soc_percent = B_kWh ./ E_Rated_kWh .* 100
    
    # Physical limits
    P_B_R_kW = data[:P_B_R][b]
    ylimit = (-P_B_R_kW * 1.1, P_B_R_kW * 1.1)
    
    # Set theme
    gr()
    theme(:mute)
    
    # Create charging/discharging bar plot with exact positioning
    power_x_positions = collect(0.5:1.0:(T-0.5))
    
    # Set common x-axis limits for both plots
    x_min, x_max = -1.0, T + 0.5
    
    charging_discharge_plot = bar(
        power_x_positions, charging_power_kW,
        dpi=600,
        label="Charging (P_c)",
        color=:green,
        bar_width=1.0,
        legend=:topleft,
        legendfontsize=8,
        legend_background_color=RGBA(1,1,1,0.7),
        xlabel="Time Interval (t)",
        ylabel="P_c/P_d [kW]",
        ylim=ylimit,
        xlim=(x_min, x_max),
        xticks=0:T,
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        title="Battery Actions - $(method_name) - Battery $(b)\nCharging and Discharging",
        titlefont=font(12, "Computer Modern"),
        guidefont=font(15, "Computer Modern"),
        tickfontfamily="Computer Modern",
        top_margin=5Plots.mm
    )
    
    bar!(charging_discharge_plot,
        power_x_positions, -discharging_power_kW,
        dpi=600,
        label="Discharging (P_d)", 
        color=:maroon,
        bar_width=1.0
    )
    
    hline!(charging_discharge_plot, [0], color=:black, lw=2, label=false)
    
    # Create SOC bar plot with exact positioning
    soc_x_positions = [0.0; collect(1.0:T)]
    
    # Plot B0 (initial SOC) with dark plum color to distinguish it as fixed
    soc_plot = bar(
        [soc_x_positions[1]], [soc_percent[1]],
        dpi=600,
        label="Initial SOC B₀ (Fixed)",
        color=:indigo,
        alpha=0.9,
        bar_width=1.0,
        linewidth=0,
        legend=:top,
        legendfontsize=8,
        legend_background_color=RGBA(1,1,1,0.7),
        xlabel="Time Interval (t)",
        ylabel="SOC [%]",
        ylim=(data[:soc_min][b] * 100 * 0.95, data[:soc_max][b] * 100 * 1.10),
        xlim=(x_min, x_max),
        xticks=0:T,
        yticks=5*(div(data[:soc_min][b] * 100 * 0.95, 5)-1):10:5*(div(data[:soc_max][b] * 100 * 1.05, 5)+1),
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.05,
        title="SOC",
        titlefont=font(12, "Computer Modern"),
        guidefont=font(15, "Computer Modern"),
        tickfontfamily="Computer Modern",
        top_margin=0Plots.mm
    )
    
    # Add B1, B2, ..., BT (optimized SOC values)
    if length(soc_percent) > 1
        bar!(soc_plot,
            soc_x_positions[2:end], soc_percent[2:end],
            dpi=600,
            label="SOC (B)",
            color=:purple,
            alpha=1.0,
            bar_width=1.0
        )
    end
    
    # Combine the two plots in vertical layout
    plot_combined = plot(charging_discharge_plot, soc_plot, layout=(2,1))
    
    # Show the plot if requested
    if showPlots
        display(plot_combined)
    end
    
    # Save the plot if requested
    if savePlots
        @printf "Saving %s battery actions plot to: %s\n" method_name filename
        savefig(plot_combined, filename)
    end
    
    return plot_combined
end


"""
    plot_battery_actions_comparison(solutions::Vector{Dict}, data::Dict, method_names::Vector{String}; 
                                    showPlots::Bool=true, savePlots::Bool=false, 
                                    filename::String="battery_actions_comparison.png")

Plot battery actions for multiple solution methods side-by-side for comparison.

# Arguments
- solutions: Vector of solution dictionaries
- data: Data dictionary containing battery parameters and bases
- method_names: Vector of method names corresponding to solutions
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Combined plot object with all methods side-by-side
"""
function plot_battery_actions_comparison(solutions::Vector{Dict}, data::Dict, method_names::Vector{String}; 
                                         showPlots::Bool=true, savePlots::Bool=false, 
                                         filename::String="battery_actions_comparison.png")
    
    n_methods = length(solutions)
    @assert length(method_names) == n_methods "Number of method names must match number of solutions"
    
    # Create individual plots for each method
    plots = []
    for (sol, name) in zip(solutions, method_names)
        p = plot_battery_actions(sol, data, name; showPlots=false, savePlots=false)
        if !isnothing(p)
            push!(plots, p)
        end
    end
    
    if isempty(plots)
        @warn "No valid plots generated for comparison"
        return nothing
    end
    
    # Combine plots horizontally
    combined_plot = plot(plots..., layout=(1, length(plots)), size=(600*length(plots), 600))
    
    # Show the plot if requested
    if showPlots
        display(combined_plot)
    end
    
    # Save the plot if requested
    if savePlots
        @printf "Saving battery actions comparison plot to: %s\n" filename
        savefig(combined_plot, filename)
    end
    
    return combined_plot
end
