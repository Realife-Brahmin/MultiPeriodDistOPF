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
    theme(:dao)
    
    # Create main plot with LoadShapeLoad (dark yellow theme) and TOP legend (like tADMM style)
    p = plot(
        time_steps, LoadShapeLoad,
        dpi=600,
        label=L"Load Shape (\lambda^t)",
        xlabel="Time Period (t)",
        ylabel="Normalized Load Profile [0-1]",
        legend=:topright,
        legend_background_color=RGBA(1,1,1,0.75),
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
        top_margin=6Plots.mm
    )

    # Add secondary y-axis for cost curves (sexy green shades)
    ax2 = twinx()
    
    # Substation 1 cost: Dark forest green (#228B22)
    plot!(
        ax2, time_steps, cost_1_cents,
        label="Cost Subs 1 \$(C^t_1)\$",
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
        label=L"Cost Subs 2 (C^t_2)",
        lw=4,
        color=RGB(50/255, 205/255, 50/255),
        linestyle=:dash,
        markershape=:diamond,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0
    )

    # Consolidate legend entries to a single top legend on primary axis
    plot!(p, [NaN], [NaN], seriestype=:scatter, label=L"Cost Subs 1 (C^t_1)",
        marker=:circle, markercolor=RGB(34/255,139/255,34/255), markersize=7,
        markerstrokecolor=:black, markerstrokewidth=2.0, lw=0)
    plot!(p, [NaN], [NaN], seriestype=:scatter, label=L"Cost Subs 2 (C^t_2)",
        marker=:diamond, markercolor=RGB(50/255,205/255,50/255), markersize=7,
        markerstrokecolor=:black, markerstrokewidth=2.0, lw=0)
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
    """
    Publication-ready stacked subplot visualization for IEEE PESGM.
    Top panel: Phase angles (θ), Bottom panel: Voltage magnitudes (V).
    """
    T = data[:T]
    Δt = data[:delta_t_h]
    time_steps = 1:T
    
    # Extract angle data
    angles_matrix = result[:angles_matrix]
    median_angles_deg = result[:median_angles_deg]
    
    # Get the non-slack substation indices
    non_slack_subs = sort(collect(keys(angles_matrix)))
    if isempty(non_slack_subs)
        @warn "No non-slack substations found in angles_matrix"
        return nothing
    end
    
    # For multiple substations, we'll plot all non-slack angles
    # Get all angle trajectories and medians
    θ_trajectories = Dict()
    θ_medians = Dict()
    for sub in non_slack_subs
        θ_trajectories[sub] = angles_matrix[sub]
        θ_medians[sub] = median_angles_deg[sub]
    end
    
    # Extract voltage data
    same_voltage_levels = data[:same_voltage_levels]
    slack_node_val = result[:slack_node]
    
    # Load bus voltage (always available)
    v_load_pu = sqrt.(result[:v_j])
    
    # Determine number of substations (count V_*_pu keys in data)
    n_subs = 0
    for i in 1:10  # Check up to 10 substations
        if haskey(data, Symbol("V_$(i)_pu"))
            n_subs = i
        else
            break
        end
    end
    
    # Substation voltages
    v_subs = Dict{Int, Vector{Float64}}()
    if same_voltage_levels
        # All substations fixed
        for sub in 1:n_subs
            v_subs[sub] = fill(data[Symbol("V_$(sub)_pu")], T)
        end
    else
        # Non-slack substations are variable
        for sub in 1:n_subs
            if sub - 1 == slack_node_val  # Convert sub index to node index
                # This is the slack substation (fixed)
                v_subs[sub] = fill(data[Symbol("V_$(sub)_pu")], T)
            else
                # Non-slack substation (variable)
                if haskey(result, :v_nonslack) && haskey(result[:v_nonslack], sub)
                    v_subs[sub] = sqrt.(result[:v_nonslack][sub])
                else
                    @warn "Missing voltage data for substation $sub"
                    v_subs[sub] = fill(NaN, T)
                end
            end
        end
    end
    
    # Calculate y-axis limits for angles θ (left axis)
    # Use all non-slack angle trajectories
    all_angles = vcat([θ_trajectories[sub] for sub in non_slack_subs]...)
    θ_min = minimum(all_angles)
    θ_max = maximum(all_angles)
    θ_range = θ_max - θ_min
    θ_margin = max(0.5, 0.15 * θ_range)
    left_min = θ_min - θ_margin
    left_max = θ_max + θ_margin
    
    # Calculate y-axis limits for voltage (right axis) - only non-slack voltages
    all_v_nonslack = vcat([v_subs[sub] for sub in non_slack_subs]...)
    v_min = minimum(all_v_nonslack)
    v_max = maximum(all_v_nonslack)
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
    theme(:dao)
    
    # Determine slack information for title
    slack_sub_idx = slack_node_val + 1
    slack_v_str = "V_{$slack_sub_idx} = $(data[Symbol("V_$(slack_sub_idx)_pu")]) \\, \\mathrm{pu}"
    δ_slack_str = "\\delta_{$slack_sub_idx} = 0.0^{\\circ}"
    title_str = "Substation Voltage Profile vs Time (Slack: Subs $slack_sub_idx)"
    
    # Color scheme: one family per substation (from cove's palette)
    # Sub 1: deep blue (angle) + sky blue (voltage) - calm, baseline
    # Sub 2: amber/orange (angle) + light amber (voltage) - warm, active
    # Sub 3: teal green (angle) + mint green (voltage) - stable, cool
    angle_colors = [RGB(0x00/255, 0x72/255, 0xB2/255),    # #0072B2 deep blue
                    RGB(0xE6/255, 0x9F/255, 0x00/255),    # #E69F00 amber/orange
                    RGB(0x00/255, 0x9E/255, 0x73/255)]    # #009E73 teal green
    voltage_colors = [RGB(0x66/255, 0xB2/255, 0xFF/255),  # #66B2FF sky blue
                      RGB(0xFF/255, 0xD5/255, 0x80/255),  # #FFD580 light amber
                      RGB(0x7F/255, 0xD1/255, 0xAE/255)]  # #7FD1AE mint green
    
    # Smart x-tick spacing: fewer ticks for clarity (every 2-3 hours)
    xtick_sparse = if T <= 24
        collect(1:2:T)
    else
        step = max(2, div(T, 8))
        collect(1:step:T)
    end
    
    # ==== TOP PANEL: ANGLES ====
    p_angle = plot(
        dpi=400,
        ylabel="Angle [°]",
        legend=:topright,
        legend_columns=1,
        legendfontsize=8,
        framestyle=:box,
        ylims=(left_min, left_max),
        xticks=(xtick_sparse, string.(xtick_sparse)),
        xformatter=_->"",  # Hide x-labels on top panel
        grid=true,
        gridlinewidth=0.5,
        gridalpha=0.25,
        gridstyle=:solid,
        titlefont=font(11, "Computer Modern"),
        guidefont=font(10, "Computer Modern"),
        tickfont=font(9, "Computer Modern"),
        left_margin=3Plots.mm,
        right_margin=3Plots.mm,
        top_margin=2Plots.mm,
        bottom_margin=0Plots.mm,
        title=title_str
    )
    
    # Plot angle trajectories with δ reference lines
    for (idx, sub) in enumerate(non_slack_subs)
        θ_traj = θ_trajectories[sub]
        θ_med = θ_medians[sub]
        
        angle_col = angle_colors[mod1(idx, length(angle_colors))]
        
        # δ reference line (thicker, more visible)
        hline!(
            p_angle,
            [θ_med],
            linestyle=:dash,
            linewidth=2.0,
            linecolor=angle_col,
            alpha=0.8,
            label="δ$sub = $(round(θ_med, digits=2))°"
        )
        
        # Angle trajectory (solid line, lw=2.2, marker size=5, with black edge)
        plot!(
            p_angle, time_steps, θ_traj,
            color=angle_col,
            lw=2.2,
            markershape=:circle,
            markersize=5,
            markerstrokewidth=1.5,
            markerstrokecolor=:black,
            markeralpha=0.9,
            label="θ @ Sub $sub"
        )
    end
    
    # ==== BOTTOM PANEL: VOLTAGES ====
    p_voltage = plot(
        dpi=400,
        xlabel="Time Period (t)",
        ylabel="Voltage [pu]",
        legend=:topright,
        legend_columns=2,
        legendfontsize=8,
        framestyle=:box,
        ylims=(right_min, right_max),
        xticks=(xtick_sparse, string.(xtick_sparse)),
        grid=true,
        gridlinewidth=0.5,
        gridalpha=0.25,
        gridstyle=:solid,
        guidefont=font(10, "Computer Modern"),
        tickfont=font(9, "Computer Modern"),
        left_margin=3Plots.mm,
        right_margin=3Plots.mm,
        top_margin=0Plots.mm,
        bottom_margin=4Plots.mm
    )
    
    # Plot voltage trajectories (lw=2.0, dashed, marker size=4)
    for (idx, sub) in enumerate(non_slack_subs)
        v_traj = v_subs[sub]
        
        volt_col = voltage_colors[mod1(idx, length(voltage_colors))]
        
        # Voltage trajectory (dashed line, lw=2.0, marker size=4, with black edge)
        plot!(
            p_voltage, time_steps, v_traj,
            color=volt_col,
            lw=2.0,
            linestyle=:dash,
            markershape=:square,
            markersize=4,
            markerstrokewidth=1.5,
            markerstrokecolor=:black,
            markeralpha=0.9,
            label="V @ Sub $sub"
        )
    end
    
    # Stack the two panels vertically with IEEE 2-column dimensions
    # 7.1 inches × 4.6 inches at 400 DPI for full 2-column width
    p = plot(p_angle, p_voltage, 
             layout=(2,1), 
             size=(7.1*100, 4.6*100),  # Convert inches to pixels (100 px/inch for screen)
             plot_title="",
             link=:x)
    
    # Show the plot if requested
    if showPlots
        display(p)
    end
    
    # Save the plot if requested (already at 400 dpi from panel settings)
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
