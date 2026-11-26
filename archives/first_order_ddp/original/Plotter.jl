"""
Plotting utilities for DDP Copper Plate Problem
Simplified plotting functions (no network, single battery)
"""

using Plots
using Printf
using Statistics

"""
    plot_battery_actions(solution::Dict, data, method_name::String; kwargs...)

Plot battery charging/discharging power and SOC for copper plate solution.
Creates separate charging/discharging bars and SOC plot.

# Arguments
- solution: Solution dictionary containing :P_B, :B arrays
- data: CopperPlateData instance
- method_name: Name of the solution method (e.g., "Brute Force", "DDP")
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Combined plot object with charging/discharging and SOC subplots
"""
function plot_battery_actions(solution::Dict, data, method_name::String; 
                               showPlots::Bool=true, savePlots::Bool=false, 
                               filename::String="battery_actions.png")
    
    T = data.T
    P_BASE = data.P_BASE
    E_BASE = data.E_BASE
    bat = data.bat
    Δt = bat.Δt
    
    # Extract battery power and SOC
    P_B_pu = solution[:P_B]
    B_pu = solution[:B]
    
    # Convert to physical units
    P_B_kW = P_B_pu .* P_BASE
    B_kWh = B_pu .* E_BASE
    
    # Include initial SOC B0 for SOC plot
    B0_kWh = bat.B0_pu * E_BASE
    B_kWh_with_initial = [B0_kWh; B_kWh]
    
    # Separate charging and discharging (P_B > 0 is discharging, P_B < 0 is charging)
    charging_power_kW = -min.(P_B_kW, 0.0)  # P_B < 0 made positive for green bars
    discharging_power_kW = max.(P_B_kW, 0.0)  # P_B > 0 kept positive for wine red bars
    
    # Calculate SOC percentages
    E_Rated_kWh = bat.E_Rated_pu * E_BASE
    soc_percent = B_kWh_with_initial ./ E_Rated_kWh .* 100
    
    # Physical limits
    P_B_R_kW = bat.P_B_R_pu * P_BASE
    ylimit_power = (-P_B_R_kW * 1.1, P_B_R_kW * 1.1)
    
    # Set theme
    gr()
    theme(:mute)
    
    # Create title
    title_text = "Battery Actions - $(method_name)\nCharging and Discharging"
    
    # P_B bars: positioned at interval boundaries [0.5, 1.5, 2.5, ...] for intervals [0,1), [1,2), [2,3), ...
    power_x_positions = collect(0.5:1.0:(T-0.5))
    
    # Set common x-axis limits for both plots
    x_min, x_max = -0.5, T + 0.5
    
    # Create charging/discharging bar plot
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
        ylim=ylimit_power,
        xlim=(x_min, x_max),
        xticks=0:T,
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        title=title_text,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(15, "Computer Modern"),
        tickfontfamily="Computer Modern",
        top_margin=5Plots.mm
    )
    
    bar!(charging_discharge_plot,
        power_x_positions, -discharging_power_kW,
        label="Discharging (P_d)", 
        color=:maroon,
        bar_width=1.0
    )
    
    hline!(charging_discharge_plot, [0], color=:black, lw=2, label=false)
    
    # Create SOC bar plot
    # B bars: B₀ at x=0, B₁ at x=1, B₂ at x=2, etc.
    soc_x_positions = collect(0:T)
    
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
        ylim=(bat.soc_min * 100 * 0.95, bat.soc_max * 100 * 1.10),
        xlim=(x_min, x_max),
        xticks=0:T,
        yticks=5*(div(bat.soc_min * 100 * 0.95, 5)-1):10:5*(div(bat.soc_max * 100 * 1.05, 5)+1),
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
    
    # Add B1, B2, ..., BT (optimized SOC values) with normal styling
    if length(soc_percent) > 1
        bar!(soc_plot,
            soc_x_positions[2:end], soc_percent[2:end],
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
    plot_ddp_convergence(solution_ddp; showPlots=true, savePlots=false, filename="ddp_convergence.png")

Plot DDP convergence showing objective function and convergence error across iterations.

# Arguments
- solution_ddp: DDP solution dictionary with :convergence_history
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot
"""
function plot_ddp_convergence(solution_ddp; showPlots::Bool=true, savePlots::Bool=false, 
                              filename::String="ddp_convergence.png", bf_objective::Union{Nothing, Float64}=nothing)
    
    conv_hist = solution_ddp[:convergence_history]
    obj_history = conv_hist[:obj_history]
    convergence_error = conv_hist[:convergence_error]
    
    iterations = 1:length(obj_history)
    
    # Get DDP final objective
    ddp_final_obj = solution_ddp[:objective]
    
    # Set theme and colors
    gr()
    theme(:mute)
    
    # Create objective function plot with DDP trace
    ddp_label = @sprintf "DDP (final: \$%.2f)" ddp_final_obj
    obj_plot = plot(
        iterations, obj_history,
        dpi=600,
        label=ddp_label,
        xlabel="Iteration (k)",
        ylabel="Objective Function [\$]",
        legend=:topright,
        lw=3,
        color=:magenta,  # Magenta for DDP (blue reserved for tADMM)
        markershape=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        title="DDP Convergence - Objective Function",
        titlefont=font(11, "Computer Modern"),
        guidefont=font(11, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Add BF reference line if provided
    if !isnothing(bf_objective)
        bf_label = @sprintf "BF Optimal: \$%.2f" bf_objective
        plot!(obj_plot, 
              iterations, 
              fill(bf_objective, length(iterations)),
              label=bf_label,
              linestyle=:dash,
              lw=2.5,
              color=:orange)
    end
    
    # Create convergence error plot (log scale if range is large enough)
    use_log_scale = (maximum(convergence_error) / minimum(convergence_error[convergence_error .> 0]) > 10.0)
    
    error_plot = plot(
        iterations, convergence_error,
        dpi=600,
        label="Convergence Error (max(||ΔB||, ||Δμ||))",
        xlabel="Iteration (k)",
        ylabel=use_log_scale ? "max(||ΔB||, ||Δμ||) [log scale]" : "max(||ΔB||, ||Δμ||)",
        legend=:topright,
        lw=3,
        color=:darkgreen,
        markershape=:square,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=false,  # Disable minor grid to avoid NaN issues
        title="Convergence Error",
        titlefont=font(11, "Computer Modern"),
        guidefont=font(11, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Only use log scale if data range is suitable
    if use_log_scale
        plot!(error_plot, yscale=:log10)
    end
    
    # Combine plots in vertical layout
    combined_plot = plot(obj_plot, error_plot, layout=(2,1), size=(800, 600))
    
    # Show the plot if requested
    if showPlots
        display(combined_plot)
    end
    
    # Save the plot if requested
    if savePlots
        @printf "Saving DDP convergence plot to: %s\n" filename
        savefig(combined_plot, filename)
    end
    
    return combined_plot
end

"""
    plot_comparison(sol_bf, sol_ddp, data; showPlots=true, savePlots=false, filename="comparison.png")

Plot side-by-side comparison of Brute Force and DDP solutions.

# Arguments
- sol_bf: Brute force solution dictionary
- sol_ddp: DDP solution dictionary
- data: CopperPlateData instance
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot
"""
function plot_comparison(sol_bf, sol_ddp, data; showPlots::Bool=true, savePlots::Bool=false,
                        filename::String="comparison.png")
    
    T = data.T
    P_BASE = data.P_BASE
    
    # Convert to physical units
    P_Subs_bf_kW = sol_bf[:P_Subs] .* P_BASE
    P_Subs_ddp_kW = sol_ddp[:P_Subs] .* P_BASE
    P_L_kW = data.P_L_pu .* P_BASE
    
    time_steps = 1:T
    
    # Set theme
    gr()
    theme(:mute)
    
    # Create comparison plot
    p = plot(
        time_steps, P_Subs_bf_kW,
        dpi=600,
        label="Brute Force",
        xlabel="Time Period (t)",
        ylabel="Substation Power [kW]",
        legend=:topright,
        lw=3,
        color=:darkorange,
        markershape=:circle,
        markersize=5,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        title="Substation Power Comparison: Brute Force vs DDP",
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    plot!(p, time_steps, P_Subs_ddp_kW,
        label="DDP",
        lw=3,
        color=:magenta,
        markershape=:square,
        markersize=5,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )
    
    plot!(p, time_steps, P_L_kW,
        label="Load",
        lw=2,
        color=:darkgoldenrod2,
        linestyle=:dash,
        markershape=:diamond,
        markersize=4
    )
    
    # Show the plot if requested
    if showPlots
        display(p)
    end
    
    # Save the plot if requested
    if savePlots
        @printf "Saving comparison plot to: %s\n" filename
        savefig(p, filename)
    end
    
    return p
end
