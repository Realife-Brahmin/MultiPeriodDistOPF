"""
Plotting utilities for tADMM Multi-Period OPF
Standalone plotting functions for battery actions and input curves
"""

using Plots
using Printf
using LaTeXStrings

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
                               filename::String="battery_actions_lindistflow.png",
                               battery_index::Union{Int,Nothing}=nothing,
                               plot_all_batteries::Bool=false)
    
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
    
    # Determine which battery to plot
    if plot_all_batteries
        # Plot all batteries - save separate files
        for (idx, b) in enumerate(Bset)
            # Modify filename to include battery index
            base, ext = splitext(filename)
            battery_filename = "$(base)_battery$(b)$(ext)"
            _plot_single_battery(solution, data, method_name, b, T, P_BASE, E_BASE, 
                               showPlots, savePlots, battery_filename)
        end
        return nothing
    else
        # Plot single battery (first or specified)
        b = battery_index === nothing ? first(Bset) : battery_index
        return _plot_single_battery(solution, data, method_name, b, T, P_BASE, E_BASE, 
                                   showPlots, savePlots, filename)
    end
end

function _plot_single_battery(solution::Dict, data::Dict, method_name::String, b::Int,
                              T::Int, P_BASE::Float64, E_BASE::Float64,
                              showPlots::Bool, savePlots::Bool, filename::String)
    
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
        savefig(combined_plot, filename)
    end
    
    return combined_plot
end

"""
    plot_substation_power_and_cost(solution, data, method_name; showPlots=true, savePlots=false, filename="substation_power_cost.png")

Plot substation power (black) and substation power cost (green) vs time on dual y-axes.
Power cost is the actual objective function contribution: P_Subs * Cost * dt

# Arguments
- solution: Solution dictionary containing :P_Subs and :status
- data: Dictionary containing :T, :LoadShapeCost, :kVA_B, :delta_t_h
- method_name: Name of the optimization method (for title)
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot
"""
function plot_substation_power_and_cost(solution, data, method_name; 
                                       showPlots::Bool=true, savePlots::Bool=false, 
                                       filename::String="substation_power_cost.png")
    
    T = data[:T]
    LoadShapeCost = data[:LoadShapeCost]
    P_BASE = data[:kVA_B]
    delta_t_h = data[:delta_t_h]
    
    # Extract substation power
    P_Subs_var = solution[:P_Subs]
    
    # Convert to actual values
    try
        P_Subs_kW = [P_Subs_var[t] * P_BASE for t in 1:T]
        
        # Calculate actual power cost: P_Subs (kW) * Cost ($/kWh) * dt (h) = $ per time step
        power_cost_dollars = [P_Subs_kW[t] * LoadShapeCost[t] * delta_t_h for t in 1:T]
        
        # Set theme and backend
        gr()
        theme(:mute)
        
        # Set up x-axis ticks (all time steps if T=24)
        xtick_vals = T == 24 ? (1:T) : :auto
        
        # Create plot with dual y-axes
        p = plot(
            1:T, P_Subs_kW,
            dpi=600,
            label="Substation Power",
            xlabel="Time Period (t)",
            ylabel="Substation Power (kW)",
            legend=:topright,
            color=:black,
            linewidth=3,
            fontfamily="Computer Modern",
            grid=true,
            gridstyle=:dot,
            gridalpha=0.5,
            minorgrid=true,
            minorgridstyle=:dot,
            minorgridalpha=0.2,
            title="$method_name: Substation Power & Cost",
            size=(900, 500),
            xticks=xtick_vals
        )
        
        # Add cost on secondary y-axis
        plot!(twinx(), 1:T, power_cost_dollars,
            label="Power Cost",
            ylabel="Power Cost (\$/period)",
            legend=:topleft,
            color=:darkgreen,
            linewidth=3,
            linestyle=:solid,
            xticks=xtick_vals
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
        
    catch e
        @warn "Could not create substation power & cost plot" exception=(e, catch_backtrace())
        return nothing
    end
end

"""
    plot_voltage_profile_one_bus(solution, data, method_name; bus=nothing, showPlots=true, savePlots=false, filename="voltage_one_bus.png")

Plot voltage magnitude vs time for a single bus (default: last bus in the system).

# Arguments
- solution: Solution dictionary containing :v (voltage squared)
- data: Dictionary containing :T, :Nset
- method_name: Name of the optimization method (for title)
- bus: Bus number to plot (default: last bus in Nset)
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot
"""
function plot_voltage_profile_one_bus(solution, data, method_name; 
                                      bus::Union{Int,Nothing}=nothing,
                                      showPlots::Bool=true, savePlots::Bool=false, 
                                      filename::String="voltage_one_bus.png")
    
    T = data[:T]
    Nset = data[:Nset]
    
    # Select bus (default: last bus)
    selected_bus = isnothing(bus) ? maximum(Nset) : bus
    
    if !(selected_bus in Nset)
        @warn "Bus $selected_bus not found in system. Available buses: $Nset"
        return nothing
    end
    
    # Extract voltage (squared)
    v_var = solution[:v]
    
    try
        # Convert v (squared) to voltage magnitude in p.u.
        V_pu = [sqrt(v_var[selected_bus, t]) for t in 1:T]
        
        # Set theme and backend
        gr()
        theme(:mute)
        
        # Set up x-axis ticks (all time steps if T=24)
        xtick_vals = T == 24 ? (1:T) : :auto
        
        # Create plot with LaTeX notation
        p = plot(
            1:T, V_pu,
            dpi=600,
            label=L"V^t_{%$selected_bus}",
            xlabel=L"t",
            ylabel=L"V^t_j \textrm{ [p.u.]}",
            legend=:topright,
            color=:blue,
            linewidth=3,
            fontfamily="Computer Modern",
            grid=true,
            gridstyle=:dot,
            gridalpha=0.5,
            minorgrid=true,
            minorgridstyle=:dot,
            minorgridalpha=0.2,
            title="$method_name: " * L"V^t_j" * " vs " * L"t" * " for Bus $selected_bus",
            size=(900, 500),
            ylims=(0.92, 1.08),
            xticks=xtick_vals
        )
        
        # Add voltage limit lines (typically 0.95 and 1.05)
        hline!([0.95, 1.05], color=:red, linestyle=:dash, linewidth=2, label="Limits (0.95/1.05 p.u.)")
        
        # Show the plot if requested
        if showPlots
            display(p)
        end
        
        # Save the plot if requested
        if savePlots
            savefig(p, filename)
        end
        
        return p
        
    catch e
        @warn "Could not create voltage profile plot for bus $selected_bus" exception=(e, catch_backtrace())
        return nothing
    end
end

"""
    plot_voltage_profile_all_buses(solution, data, method_name; time_step=nothing, showPlots=true, savePlots=false, filename="voltage_all_buses.png", create_gif=false, gif_filename="voltage_animation.gif")

Plot voltage magnitude for all buses at a specific time step (default: T÷2).
Optionally create an animated GIF showing voltage evolution over all time steps.

# Arguments
- solution: Solution dictionary containing :v (voltage squared)
- data: Dictionary containing :T, :Nset
- method_name: Name of the optimization method (for title)
- time_step: Time step to plot (default: T÷2)
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the static plot
- create_gif: Whether to create an animated GIF
- gif_filename: Filename for the animated GIF
"""
function plot_voltage_profile_all_buses(solution, data, method_name; 
                                       time_step::Union{Int,Nothing}=nothing,
                                       showPlots::Bool=true, savePlots::Bool=false, 
                                       filename::String="voltage_all_buses.png",
                                       create_gif::Bool=false,
                                       gif_filename::String="voltage_animation.gif")
    
    T = data[:T]
    Nset = sort(collect(data[:Nset]))  # Sort buses for consistent plotting
    
    # Select time step (default: middle of horizon)
    selected_t = isnothing(time_step) ? T ÷ 2 : time_step
    
    if !(1 <= selected_t <= T)
        @warn "Time step $selected_t out of range [1, $T]"
        return nothing
    end
    
    # Extract voltage (squared)
    v_var = solution[:v]
    
    try
        if create_gif
            # Create animation with LaTeX notation
            anim = @animate for t in 1:T
                V_pu = [sqrt(v_var[bus, t]) for bus in Nset]
                
                plot(
                    Nset, V_pu,
                    dpi=300,
                    xlabel=L"j",
                    ylabel=L"V^t_j \textrm{ [p.u.]}",
                    legend=false,
                    color=:blue,
                    linewidth=3,
                    marker=:circle,
                    markersize=5,
                    fontfamily="Computer Modern",
                    grid=true,
                    gridstyle=:dot,
                    gridalpha=0.5,
                    minorgrid=true,
                    minorgridstyle=:dot,
                    minorgridalpha=0.2,
                    title="$method_name: " * L"V^t_j" * " vs " * L"j" * " at t=$t/$T",
                    size=(900, 500),
                    ylims=(0.92, 1.08)
                )
                
                # Add voltage limit lines
                hline!([0.95, 1.05], color=:red, linestyle=:dash, linewidth=2, alpha=0.7)
            end
            
            # Save the GIF (suppress output messages)
            redirect_stdout(devnull) do
                redirect_stderr(devnull) do
                    gif(anim, gif_filename, fps=2)
                end
            end
        end
        
        # Create static plot for selected time step with LaTeX notation
        V_pu = [sqrt(v_var[bus, selected_t]) for bus in Nset]
        
        gr()
        theme(:mute)
        
        p = plot(
            Nset, V_pu,
            dpi=600,
            xlabel=L"j",
            ylabel=L"V^t_j \textrm{ [p.u.]}",
            legend=false,
            color=:blue,
            linewidth=3,
            marker=:circle,
            markersize=6,
            fontfamily="Computer Modern",
            grid=true,
            gridstyle=:dot,
            gridalpha=0.5,
            minorgrid=true,
            minorgridstyle=:dot,
            minorgridalpha=0.2,
            title="$method_name: " * L"V^t_j" * " vs " * L"j" * " at t=$selected_t/$T",
            size=(900, 500),
            ylims=(0.92, 1.08)
        )
        
        # Add voltage limit lines
        hline!([0.95, 1.05], color=:red, linestyle=:dash, linewidth=2, label="Limits")
        
        # Show the plot if requested
        if showPlots
            display(p)
        end
        
        # Save the plot if requested
        if savePlots
            savefig(p, filename)
        end
        
        return p
        
    catch e
        @warn "Could not create voltage profile plot for all buses" exception=(e, catch_backtrace())
        return nothing
    end
end

"""
    plot_pv_power(sol::Dict, data::Dict, method_label::String="Method"; 
                  pv_index::Int=1, showPlots::Bool=true, savePlots::Bool=false, 
                  filename::String="pv_power.png")

Plot PV real and reactive power output over time for a single PV unit.
- Top subplot: p_D (real power) in orange
- Bottom subplot: q_D (reactive power) in purple (bidirectional)

# Arguments
- sol: Solution dictionary containing :p_D (from data) and :q_D (from optimization)
- data: Data dictionary containing :Dset, :p_D_pu, :T, :kVA_B
- method_label: Label for the method used
- pv_index: Index of PV to plot (1 = first PV, -1 = last PV)
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Plot object
"""
function plot_pv_power(sol::Dict, data::Dict, method_label::String="Method"; 
                       pv_index::Int=1, showPlots::Bool=true, savePlots::Bool=false, 
                       filename::String="pv_power.png")
    try
        Dset = data[:Dset]
        
        if isempty(Dset)
            @warn "No PV units in system, cannot plot PV power"
            return nothing
        end
        
        # Select PV unit
        pv_bus = pv_index == -1 ? maximum(Dset) : minimum(Dset)
        if pv_index > 0 && length(Dset) >= pv_index
            pv_bus = sort(collect(Dset))[pv_index]
        end
        
        T = data[:T]
        P_BASE = data[:kVA_B]
        p_D_pu = data[:p_D_pu]
        q_D_vals = sol[:q_D]
        
        time_steps = 1:T
        
        # Extract p_D (from data) and q_D (from solution) for the selected PV
        p_D_kW = [p_D_pu[pv_bus, t] * P_BASE for t in 1:T]
        q_D_kvar = [q_D_vals[pv_bus, t] * P_BASE for t in 1:T]
        
        # Set theme and backend
        gr()
        theme(:mute)
        
        # Create subplot layout
        p1 = plot(
            time_steps, p_D_kW,
            dpi=600,
            label=L"p_D^t (Real Power)",
            xlabel="",
            ylabel="PV Real Power (kW)",
            title="PV Power Output at Bus $pv_bus - $method_label",
            linewidth=3,
            linecolor=:orange,
            legend=:topright,
            size=(900, 700),
            grid=:on,
            minorgrid=:on,
            gridlinewidth=1,
            gridalpha=0.3,
            minorgridalpha=0.15,
            framestyle=:box,
            ylims=(0, 1.1 * maximum(p_D_kW))
        )
        
        # Add zero line for q_D subplot
        y_max = maximum(abs.(q_D_kvar))
        y_lim = y_max > 0 ? 1.2 * y_max : 1.0
        
        p2 = plot(
            time_steps, q_D_kvar,
            dpi=600,
            label=L"q_D^t (Reactive Power)",
            xlabel="Time Period (t)",
            ylabel="PV Reactive Power (kvar)",
            linewidth=3,
            linecolor=:purple,
            legend=:topright,
            size=(900, 700),
            grid=:on,
            minorgrid=:on,
            gridlinewidth=1,
            gridalpha=0.3,
            minorgridalpha=0.15,
            framestyle=:box,
            ylims=(-y_lim, y_lim)
        )
        
        # Add zero reference line
        hline!(p2, [0], linestyle=:dash, linecolor=:gray, linewidth=1, label="", alpha=0.5)
        
        # Combine subplots
        p = plot(p1, p2, layout=(2, 1), size=(900, 700))
        
        # Save and/or show
        if savePlots
            savefig(p, filename)
        end
        
        if showPlots
            display(p)
        end
        
        return p
        
    catch e
        @warn "Could not create PV power plot" exception=(e, catch_backtrace())
        return nothing
    end
end

"""
    plot_pv_power_circle_gif(sol::Dict, data::Dict, method_label::String="Method";
                             pv_index::Int=1, showPlots::Bool=true, savePlots::Bool=false,
                             filename::String="pv_power_circle.gif")

Create an animated GIF showing PV power as a point moving in the P-Q plane.
- X-axis: p_D (real power) >= 0
- Y-axis: q_D (reactive power), bidirectional
- Circle shows the apparent power limit S_D_R
- Point traces path over time

# Arguments
- sol: Solution dictionary containing :q_D
- data: Data dictionary containing :Dset, :p_D_pu, :S_D_R, :T, :kVA_B
- method_label: Label for the method used
- pv_index: Index of PV to plot
- showPlots: Whether to display the final frame
- savePlots: Whether to save the GIF
- filename: Filename for the GIF

# Returns
- Animation object
"""
function plot_pv_power_circle_gif(sol::Dict, data::Dict, method_label::String="Method";
                                  pv_index::Int=1, showPlots::Bool=true, savePlots::Bool=false,
                                  filename::String="pv_power_circle.gif")
    try
        Dset = data[:Dset]
        
        if isempty(Dset)
            @warn "No PV units in system, cannot create PV power circle GIF"
            return nothing
        end
        
        # Select PV unit
        pv_bus = pv_index == -1 ? maximum(Dset) : minimum(Dset)
        if pv_index > 0 && length(Dset) >= pv_index
            pv_bus = sort(collect(Dset))[pv_index]
        end
        
        T = data[:T]
        P_BASE = data[:kVA_B]
        p_D_pu = data[:p_D_pu]
        S_D_R = data[:S_D_R]
        q_D_vals = sol[:q_D]
        
        # Extract PV power trajectory
        p_D_kW = [p_D_pu[pv_bus, t] * P_BASE for t in 1:T]
        q_D_kvar = [q_D_vals[pv_bus, t] * P_BASE for t in 1:T]
        
        # Apparent power rating
        S_rating_kVA = S_D_R[pv_bus] * P_BASE
        
        # Set theme
        gr()
        theme(:mute)
        
        # Create capability circle (semicircle for p_D >= 0)
        θ = range(-π/2, π/2, length=100)
        circle_p = S_rating_kVA .* cos.(θ)
        circle_q = S_rating_kVA .* sin.(θ)
        
        # Determine axis limits
        p_max = 1.1 * S_rating_kVA
        q_max = 1.1 * S_rating_kVA
        
        # Create animation
        anim = @animate for t in 1:T
            # Plot capability circle
            p = plot(
                circle_p, circle_q,
                dpi=300,
                label="PV Capability (S = $(round(S_rating_kVA, digits=1)) kVA)",
                xlabel=L"p_D (Real Power, kW)",
                ylabel=L"q_D (Reactive Power, kvar)",
                title="PV Power at Bus $pv_bus - Time t=$t/$T",
                linewidth=2,
                linecolor=:gray,
                linestyle=:dash,
                legend=:topright,
                size=(700, 700),
                grid=:on,
                minorgrid=:on,
                gridlinewidth=1,
                gridalpha=0.3,
                minorgridalpha=0.15,
                framestyle=:box,
                xlims=(0, p_max),
                ylims=(-q_max, q_max),
                aspect_ratio=:equal
            )
            
            # Add zero reference lines
            hline!(p, [0], linestyle=:solid, linecolor=:black, linewidth=1, label="", alpha=0.3)
            vline!(p, [0], linestyle=:solid, linecolor=:black, linewidth=1, label="", alpha=0.3)
            
            # Plot trajectory up to current time (trail)
            if t > 1
                plot!(p, p_D_kW[1:t], q_D_kvar[1:t],
                      linewidth=2, linecolor=:blue, alpha=0.5, label="Trajectory")
            end
            
            # Plot current point (larger marker)
            scatter!(p, [p_D_kW[t]], [q_D_kvar[t]],
                    markersize=8, markercolor=:red, markerstrokewidth=2,
                    markerstrokecolor=:darkred, label="Current (t=$t)",
                    marker=:circle)
            
            # Add text annotation with values
            annotate!(p, p_max * 0.7, -q_max * 0.9,
                     text(@sprintf("p_D = %.2f kW\nq_D = %.2f kvar\nS = %.2f kVA",
                                   p_D_kW[t], q_D_kvar[t],
                                   sqrt(p_D_kW[t]^2 + q_D_kvar[t]^2)),
                          :left, 10))
        end
        
        # Save GIF (suppress output messages)
        if savePlots
            redirect_stdout(devnull) do
                redirect_stderr(devnull) do
                    gif(anim, filename, fps=2)
                end
            end
        end
        
        if showPlots && T > 0
            # Show final frame
            display(plot(
                circle_p, circle_q,
                dpi=300,
                label="PV Capability",
                xlabel=L"p_D (Real Power, kW)",
                ylabel=L"q_D (Reactive Power, kvar)",
                title="PV Power at Bus $pv_bus - Full Trajectory",
                linewidth=2,
                linecolor=:gray,
                linestyle=:dash,
                size=(700, 700),
                aspect_ratio=:equal
            ))
            plot!(p_D_kW, q_D_kvar, linewidth=3, linecolor=:blue, label="Trajectory")
            scatter!([p_D_kW[1]], [q_D_kvar[1]], markersize=8, markercolor=:green, label="Start")
            scatter!([p_D_kW[end]], [q_D_kvar[end]], markersize=8, markercolor=:red, label="End")
        end
        
        return anim
        
    catch e
        @warn "Could not create PV power circle GIF" exception=(e, catch_backtrace())
        return nothing
    end
end

"""
    plot_tadmm_ldf_convergence(sol_tadmm, sol_bf, eps_pri, eps_dual; showPlots::Bool=true, savePlots::Bool=false, filename::String="tadmm_convergence.png")

Plot tADMM convergence with 3 horizontal subplots matching copper plate example style:
1. Objective function value
2. Primal residual (log scale)  
3. Dual residual (log scale)

Includes BF optimal line and tolerance thresholds.

# Arguments
- sol_tadmm: tADMM solution dictionary with :convergence_history
- sol_bf: Brute force solution dictionary with :objective
- eps_pri: Primal residual tolerance
- eps_dual: Dual residual tolerance
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot to file
- filename: Filename for saving the plot

# Returns
- Combined plot object with 3 subplots
"""
function plot_tadmm_ldf_convergence(sol_tadmm, sol_bf, eps_pri::Float64, eps_dual::Float64; 
                                    showPlots::Bool=true, savePlots::Bool=false, 
                                    filename::String="tadmm_convergence.png")
    
    # Extract convergence history
    hist = sol_tadmm[:convergence_history]
    obj_history = hist[:obj_history]
    r_norm_history = hist[:r_norm_history]
    s_norm_history = hist[:s_norm_history]
    
    # Set theme and colors (matching copper plate example)
    gr()
    theme(:mute)
    line_colour_obj = :dodgerblue       # Blue for objective
    line_colour_primal = :darkgreen     # Green for primal residual
    line_colour_dual = :darkorange2     # Orange for dual residual
    
    # Subplot 1: Objective function
    p1 = plot(
        obj_history,
        dpi=600,
        xlabel="Iteration",
        ylabel="Cost (\$)",
        title="Objective",
        lw=2,
        color=line_colour_obj,
        marker=:circle,
        markersize=2,
        markerstrokewidth=0,
        legend=false,
        grid=true,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        titlefont=font(11, "Computer Modern"),
        guidefont=font(10, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Add BF optimal line if available (check for successful optimization)
    if haskey(sol_bf, :objective)
        bf_obj = sol_bf[:objective]
        # Only add line if objective is finite (successful solve)
        if isfinite(bf_obj)
            hline!(p1, [bf_obj], 
                   color=:darkorange, lw=2, linestyle=:dash, alpha=0.7)
        end
    end
    
    # Subplot 2: Primal residual (log scale)
    p2 = plot(
        r_norm_history,
        dpi=600,
        xlabel="Iteration",
        ylabel="‖r‖",
        title="Primal Residual",
        lw=2,
        color=line_colour_primal,
        marker=:circle,
        markersize=2,
        markerstrokewidth=0,
        yscale=:log10,
        legend=false,
        grid=true,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        titlefont=font(11, "Computer Modern"),
        guidefont=font(10, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Add tolerance threshold line
    hline!(p2, [eps_pri], 
           color=:red, lw=2, linestyle=:dash, alpha=0.7)
    
    # Subplot 3: Dual residual (log scale)
    p3 = plot(
        s_norm_history,
        dpi=600,
        xlabel="Iteration",
        ylabel="‖s‖",
        title="Dual Residual",
        lw=2,
        color=line_colour_dual,
        marker=:circle,
        markersize=2,
        markerstrokewidth=0,
        yscale=:log10,
        legend=false,
        grid=true,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        titlefont=font(11, "Computer Modern"),
        guidefont=font(10, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Add tolerance threshold line
    hline!(p3, [eps_dual], 
           color=:red, lw=2, linestyle=:dash, alpha=0.7)
    
    # Combine into horizontal layout (1 row, 3 columns)
    p_combined = plot(p1, p2, p3, 
                     layout=(1, 3), 
                     size=(1200, 400),
                     plot_title="tADMM Convergence Summary",
                     plot_titlefont=font(14, "Computer Modern"))
    
    # Show the plot if requested
    if showPlots
        display(p_combined)
    end
    
    # Save the plot if requested
    if savePlots
        @printf "Saving tADMM convergence plot to: %s\n" filename
        savefig(p_combined, filename)
    end
    
    return p_combined
end

