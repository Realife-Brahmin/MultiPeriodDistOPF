# Plotter.jl - Plotting functions for Multi-POI MPOPF
# Plots input curves: LoadShapeCost_1, LoadShapeCost_2, and LoadShapeLoad

using Plots
using Printf

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
    
    # Create main plot with LoadShapeLoad (dark yellow theme)
    p = plot(
        time_steps, LoadShapeLoad,
        dpi=600,
        label="Load Shape (λᵗ)",
        xlabel="Time Period (t)",
        ylabel="Normalized Load Profile [0-1]",
        legend=:left,
        lw=4,
        color=:darkgoldenrod2,  # Dark yellow theme
        markershape=:square,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
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
        right_margin=10Plots.mm
    )

    # Add secondary y-axis for cost curves (sexy green shades)
    ax2 = twinx()
    
    # Substation 1 cost: Dark forest green (#228B22)
    plot!(
        ax2, time_steps, cost_1_cents,
        label="Cost Sub 1 (C₁ᵗ)",
        lw=4,
        color=RGB(34/255, 139/255, 34/255),  # Forest green
        linestyle=:solid,
        markershape=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        ylabel="Cost [cents/kWh]",
        ylims=(right_min, right_max),
        legend=:right,
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Substation 2 cost: Bright lime green (#32CD32)
    plot!(
        ax2, time_steps, cost_2_cents,
        label="Cost Sub 2 (C₂ᵗ)",
        lw=4,
        color=RGB(50/255, 205/255, 50/255),  # Lime green
        linestyle=:dash,
        markershape=:diamond,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )

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
