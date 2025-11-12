# Plotter.jl - Plotting functions for Multi-POI MPOPF
# Plots input curves: LoadShapeCost_1, LoadShapeCost_2, and LoadShapeLoad

using Plots
using Printf
using LaTeXStrings

"""
    plot_input_curves(data; showPlots::Bool=true, savePlots::Bool=false, filename::String="mpopf_input_curves.png")

Plot LoadShapeLoad, LoadShapeCost_1, and LoadShapeCost_2 curves on the same figure.
Uses dark yellow for load (left axis) and sexy green shades for costs (right axis).

# Arguments
- data: Dictionary containing :T, :delta_t_h, :LoadShapeLoad, :LoadShapeCost_1, :LoadShapeCost_2
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Plot object
"""
function plot_input_curves(data; showPlots::Bool=true, savePlots::Bool=false, filename::String="mpopf_input_curves.png")
    T = data[:T]
    Δt = data[:delta_t_h]
    
    # Prepare data for plotting
    time_steps = 1:T
    LoadShapeLoad = data[:LoadShapeLoad]
    LoadShapeCost_1 = data[:LoadShapeCost_1]
    LoadShapeCost_2 = data[:LoadShapeCost_2]
    
    # Convert costs to cents/kWh for better readability
    cost_1_cents = LoadShapeCost_1 .* 100
    cost_2_cents = LoadShapeCost_2 .* 100
    
    # Calculate y-axis limits for left axis (load)
    left_min = -0.05
    left_max = 1.05
    right_min = floor(0.95 * min(minimum(cost_1_cents), minimum(cost_2_cents)))
    right_max = ceil(1.05 * max(maximum(cost_1_cents), maximum(cost_2_cents)))
    
    # Smart x-tick spacing (avoid overlap at endpoints)
    xtick_vals = if T <= 24
        1:T
    else
        step = max(1, div(T, 5))
        ticks = collect(step:step:T)
        # Only include T if it's far enough from the last tick
        if isempty(ticks) || (T - last(ticks)) > step * 0.3
            vcat(1, ticks, T)
        else
            vcat(1, ticks[1:end-1], T)  # Replace last tick with T
        end
    end
    xtick_vals = sort(unique(xtick_vals))
    
    # Smart time step label (hours or minutes)
    dt_label = if Δt >= 1.0
        "Δt = $(Δt) h"
    else
        dt_minutes = round(Int, Δt * 60)
        "Δt = $(dt_minutes) min"
    end

    # Set theme and backend
    gr()
    theme(:mute)
    
    # Create main plot with LoadShapeLoad (dark yellow theme) and TOP legend (like tADMM style)
    p = plot(
        time_steps, LoadShapeLoad,
        dpi=600,
        label=L"\lambda^t",
        xlabel="Time Period (t)",
        ylabel="Normalized Load Profile [0-1]",
        legend=:top,
        legendfontsize=9,
        lw=4,
        color=:darkgoldenrod2,
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
        xticks=xtick_vals,
        title="Input Curves: Load and Substation Costs ($(dt_label))",
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        right_margin=10Plots.mm,
        top_margin=4Plots.mm
    )

    # Add secondary y-axis for cost curves (sexy green shades)
    ax2 = twinx()
    
    # Substation 1 cost: Dark forest green (#228B22)
    plot!(
        ax2, time_steps, cost_1_cents,
        label=L"C^t_1",
        lw=4,
        color=RGB(34/255, 139/255, 34/255),
        linestyle=:solid,
        markershape=:circle,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0,
        ylabel="Cost [cents/kWh]",
        ylims=(right_min, right_max),
        legend=false,
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Substation 2 cost: Bright lime green (#32CD32)
    plot!(
        ax2, time_steps, cost_2_cents,
        label=L"C^t_2",
        lw=4,
        color=RGB(50/255, 205/255, 50/255),
        linestyle=:dash,
        markershape=:diamond,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0
    )

    # Consolidate legend entries to a single top legend on primary axis
    plot!(p, time_steps, fill(NaN, T), label=L"C^t_1", lw=4, color=RGB(34/255,139/255,34/255), linestyle=:solid, marker=:circle, markersize=6, markerstrokecolor=:black)
    plot!(p, time_steps, fill(NaN, T), label=L"C^t_2", lw=4, color=RGB(50/255,205/255,50/255), linestyle=:dash, marker=:diamond, markersize=6, markerstrokecolor=:black)
    # consolidated legend uses the primary axis legend with fontsize set above

    # Show the plot if requested
    if showPlots
        display(p)
    end

    # Save the plot if requested
    if savePlots
        savefig(p, filename)
    end
    
    return p
end


"""
    plot_angle_voltage_trajectories(result, data, slack_sub; showPlots::Bool=true, savePlots::Bool=false, filename::String="angle_voltage_trajectory.png")

Plot angle trajectory θ^t_i (left axis) and optimized non-slack substation voltage V_i (right axis) over time.
Reference line δ_i shows the median coordination angle for the non-slack substation.

# Arguments
- result: Dictionary containing optimization results with :angles_matrix, :median_angles_deg, :v_nonslack, :v_j
- data: Dictionary containing system data (:T, :delta_t_h, :V_1_pu, :V_2_pu, :same_voltage_levels)
- slack_sub: Which substation is slack (1 or 2)
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Plot object
"""
function plot_angle_voltage_trajectories(result, data, slack_sub; showPlots::Bool=true, savePlots::Bool=false, filename::String="angle_voltage_trajectory.png")
    T = data[:T]
    Δt = data[:delta_t_h]
    time_steps = 1:T
    
    # Extract angle data
    angles_matrix = result[:angles_matrix]
    median_angles_deg = result[:median_angles_deg]
    
    # Get the non-slack substation index
    non_slack_subs = sort(collect(keys(angles_matrix)))
    if isempty(non_slack_subs)
        @warn "No non-slack substations found in angles_matrix"
        return nothing
    end
    non_slack_sub = non_slack_subs[1]  # For 2-substation system, only one non-slack
    
    # Get angle trajectory and median
    θ_trajectory = angles_matrix[non_slack_sub]
    θ_median = median_angles_deg[non_slack_sub]
    
    # Compute angle deviation: φ = θ - θ_median
    φ_trajectory = θ_trajectory .- θ_median
    
    # Extract voltage data
    same_voltage_levels = data[:same_voltage_levels]
    slack_node_val = result[:slack_node]
    
    # Load bus voltage (always available)
    v_load_pu = sqrt.(result[:v_j])
    
    # Substation voltages
    if same_voltage_levels
        # Both substations fixed
        v_sub1 = fill(data[:V_1_pu], T)
        v_sub2 = fill(data[:V_2_pu], T)
    else
        # One is variable
        if slack_node_val == 0
            # Sub 1 is slack (fixed), Sub 2 is variable
            v_sub1 = fill(data[:V_1_pu], T)
            v_sub2 = sqrt.(result[:v_nonslack])
        else  # slack_node_val == 1
            # Sub 2 is slack (fixed), Sub 1 is variable
            v_sub1 = sqrt.(result[:v_nonslack])
            v_sub2 = fill(data[:V_2_pu], T)
        end
    end
    
    # Calculate y-axis limits for angle θ (left axis)
    θ_min = minimum(θ_trajectory)
    θ_max = maximum(θ_trajectory)
    θ_range = θ_max - θ_min
    θ_margin = max(0.5, 0.15 * θ_range)
    left_min = θ_min - θ_margin
    left_max = θ_max + θ_margin
    
    # Calculate y-axis limits for voltage (right axis) - only non-slack voltage
    v_nonslack_data = slack_node_val == 0 ? v_sub2 : v_sub1
    v_min = minimum(v_nonslack_data)
    v_max = maximum(v_nonslack_data)
    v_range = v_max - v_min
    v_margin = max(0.01, 0.1 * v_range)
    right_min = v_min - v_margin
    right_max = v_max + v_margin
    
    # Smart x-tick spacing
    xtick_vals = if T <= 24
        1:T
    else
        step = max(1, div(T, 5))
        ticks = collect(step:step:T)
        if isempty(ticks) || (T - last(ticks)) > step * 0.3
            vcat(1, ticks, T)
        else
            vcat(1, ticks[1:end-1], T)
        end
    end
    xtick_vals = sort(unique(xtick_vals))
    
    # Smart time step label
    dt_label = if Δt >= 1.0
        "Δt = $(Δt) h"
    else
        dt_minutes = round(Int, Δt * 60)
        "Δt = $(dt_minutes) min"
    end
    
    # Set theme and backend
    gr()
    theme(:mute)
    
    # Determine slack voltage and reference angle for title
    if slack_node_val == 0
        slack_v_str = "V_1 = $(data[:V_1_pu]) \\, \\mathrm{pu}"
        δ_slack_str = "\\delta_1 = 0.0^{\\circ}"
        title_str = "Substation Angle and Voltage Profile vs Time (\$$(slack_v_str), $(δ_slack_str)\$)"
    else
        slack_v_str = "V_2 = $(data[:V_2_pu]) \\, \\mathrm{pu}"
        δ_slack_str = "\\delta_2 = 0.0^{\\circ}"
        title_str = "Substation Angle and Voltage Profile vs Time (\$$(slack_v_str), $(δ_slack_str)\$)"
    end
    
    # Create main plot with angle θ (blue theme) - use LaTeX labels
    δ_text = @sprintf("%.2f", round(θ_median, digits=2))
    
    p = plot(
        time_steps, θ_trajectory,
        dpi=600,
        label=L"\theta^t_{%$non_slack_sub}",
        xlabel="Time Period (t)",
        ylabel="Angle [degrees]",
        legend=:top,
        legendfontsize=9,
        lw=4,
        color=:dodgerblue,
        markershape=:circle,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        ylims=(left_min, left_max),
        xticks=xtick_vals,
        title=title_str,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        right_margin=10Plots.mm,
        top_margin=4Plots.mm
    )
    
    # Add reference line at δ (the median coordination angle)
    plot!(
        p, time_steps, fill(θ_median, T),
        label=L"\delta_{%$non_slack_sub} = %$δ_text^{\circ}",
        lw=2.5,
        color=:gray40,
        linestyle=:dash,
        alpha=0.9
    )
    
    # Add secondary y-axis for voltages
    ax2 = twinx()
    
    # Non-slack substation voltage (the one being optimized)
    if slack_node_val == 0
        # Sub 1 is slack, plot Sub 2 (optimized)
        v_nonslack = v_sub2
        nonslack_label = L"V_{2}"
        v_color = :limegreen
    else
        # Sub 2 is slack, plot Sub 1 (optimized)
        v_nonslack = v_sub1
        nonslack_label = L"V_{1}"
        v_color = :forestgreen
    end
    
    plot!(
        ax2, time_steps, v_nonslack,
        label="$(nonslack_label)",
        lw=4,
        color=v_color,
        linestyle=:solid,
        markershape=:square,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0,
        ylabel="Voltage Magnitude [pu]",
        ylims=(right_min, right_max),
        legend=false,
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )

    # Consolidate legends: register voltage series label on primary axis with proper color/marker but no data
    plot!(p, time_steps, fill(NaN, T), label=nonslack_label, lw=4, color=v_color, linestyle=:solid,
        marker=:square, markersize=6, markerstrokecolor=:black)
    # Remove legend frame background for cleaner look
    # unified legend font size applied via legendfontsize above
    
    # Removed inline annotations now that legends are restored
    
    # Show the plot if requested
    if showPlots
        display(p)
    end
    
    # Save the plot if requested
    if savePlots
        savefig(p, filename)
    end
    
    return p
end


"""
    print_curve_statistics(data)

Print statistics about the input curves.
"""
function print_curve_statistics(data)
    println("\n" * "="^80)
    println("INPUT CURVE STATISTICS")
    println("="^80)
    
    # Load profile statistics
    λ_L = data[:LoadShapeLoad]
    println("\nLoad Profile (λ_L):")
    @printf("  Min:    %.4f pu\n", minimum(λ_L))
    @printf("  Max:    %.4f pu\n", maximum(λ_L))
    @printf("  Mean:   %.4f pu\n", sum(λ_L) / length(λ_L))
    @printf("  Range:  %.4f pu\n", maximum(λ_L) - minimum(λ_L))
    
    # Substation 1 cost statistics
    C₁ = data[:LoadShapeCost_1]
    println("\nSubstation 1 Cost (C₁):")
    @printf("  Min:    \$%.4f/kWh\n", minimum(C₁))
    @printf("  Max:    \$%.4f/kWh\n", maximum(C₁))
    @printf("  Mean:   \$%.4f/kWh\n", sum(C₁) / length(C₁))
    @printf("  Range:  \$%.4f/kWh\n", maximum(C₁) - minimum(C₁))
    
    # Substation 2 cost statistics
    C₂ = data[:LoadShapeCost_2]
    println("\nSubstation 2 Cost (C₂):")
    @printf("  Min:    \$%.4f/kWh\n", minimum(C₂))
    @printf("  Max:    \$%.4f/kWh\n", maximum(C₂))
    @printf("  Mean:   \$%.4f/kWh\n", sum(C₂) / length(C₂))
    @printf("  Range:  \$%.4f/kWh\n", maximum(C₂) - minimum(C₂))
    
    # Cost difference analysis
    println("\nCost Difference (C₁ - C₂):")
    cost_diff = C₁ .- C₂
    @printf("  Min:    \$%.4f/kWh\n", minimum(cost_diff))
    @printf("  Max:    \$%.4f/kWh\n", maximum(cost_diff))
    @printf("  Mean:   \$%.4f/kWh\n", sum(cost_diff) / length(cost_diff))
    
    # Find when each substation is cheaper
    t_sub1_cheaper = sum(C₁ .< C₂)
    t_sub2_cheaper = sum(C₂ .< C₁)
    t_equal = sum(C₁ .== C₂)
    
    println("\nCost Advantage:")
    @printf("  Sub 1 cheaper: %d periods (%.1f%%)\n", t_sub1_cheaper, 100*t_sub1_cheaper/length(C₁))
    @printf("  Sub 2 cheaper: %d periods (%.1f%%)\n", t_sub2_cheaper, 100*t_sub2_cheaper/length(C₁))
    @printf("  Equal cost:    %d periods (%.1f%%)\n", t_equal, 100*t_equal/length(C₁))
    
    println("="^80)
end

# Example usage (uncomment to test):
# include("../../small2poi_mpopf.jl")  # Load data dictionary
# plot_input_curves(data)
# plot_input_curves_detailed(data)
# print_curve_statistics(data)
