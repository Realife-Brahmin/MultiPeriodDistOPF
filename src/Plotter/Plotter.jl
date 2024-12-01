module Plotter

using LaTeXStrings
using Measures
using Plots
using Parameters: @unpack
import JuMP: value  # Import JuMP's value function to extract values of decision variables
import Base.Filesystem: mkpath, isdir  # To create directories
include("../helperFunctions.jl")
using .helperFunctions: myprintln

export plot_battery_actions, plot_input_forecast_curves, plot_line_losses, plot_substation_power, plot_substation_power_cost

#good light themes: :bright, :dao, :gruvbox_light, :solarized_light, :vibrant, :wong, :wong2
common_theme = :mute
common_marker_face = :circle
common_marker_stroke_color = :black
common_marker_stroke_width = 2.0

function plot_battery_actions(modelDict;
    showPlots::Bool=false,
    savePlots::Bool=true,
    macroItrNum::Int=1,
    verbose::Bool=false)

    theme(common_theme)
    # Extract necessary parameters from the `data` dictionary
    @unpack data, modelVals = modelDict;
    @unpack Tset, Bset, kVA_B, B_R_pu, P_B_R, Bref_pu, systemName, numAreas, T, DER_percent, Batt_percent, alpha, alphaAppendix, gamma, gammaAppendix, soc_min, soc_max, gedAppendix, solver = data;
    Tset = sort(collect(Tset))

    P_c = modelVals[:P_c]
    P_d = modelVals[:P_d]
    B = modelVals[:B]

    # Set the base directory for saving plots
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)",
        "batteryActionPlots", "macroItr_$(macroItrNum)")
    if savePlots && !isdir(base_dir)
        myprintln(verbose, "Creating directory: $base_dir")
        mkpath(base_dir)
    end

    # Loop through all battery buses and create plots
    for j in Bset
        time_intervals = collect(Tset)

        # Collect charging, discharging power, and state of charge values
        charging_power_kW = [P_c[j, t] * kVA_B for t in Tset]
        discharging_power_kW = [P_d[j, t] * kVA_B for t in Tset]
        soc = [Bref_pu[j] / B_R_pu[j] * 100; [B[j, t] / B_R_pu[j] * 100 for t in Tset]]

        ylimit = (-P_B_R[j], P_B_R[j])

        gr()

        # Create a plot for charging and discharging
        charging_discharge_plot = bar(
            time_intervals, charging_power_kW,
            dpi=600,
            label="Charging",
            color=:green,
            legend=:bottomleft,
            xlabel="Time Interval " * L"t",
            ylabel=L"P_c/P_d \, [kW]",
            ylim=ylimit,
            xticks=1:T,
            gridstyle=:solid,
            gridlinewidth=1.0,
            gridalpha=0.2,
            minorgrid=true,
            minorgridstyle=:solid,
            minorgridalpha=0.15,
            title="Battery at Bus $(j)\nCharging and Discharging",
            titlefont=font(12, "Computer Modern"),
            guidefont=font(15, "Computer Modern"),
            tickfontfamily="Computer Modern",
            # left_margin=10mm,   # Adds space to the left
            # right_margin=10mm,  # Adds space to the right
            top_margin=5mm,    # Adds space at the top
            # bottom_margin=10mm  # Adds space at the bottom
        )

        bar!(charging_discharge_plot, 
            time_intervals, -discharging_power_kW, 
            dpi=600,
            label="Discharging", color=:darkred)

        hline!(charging_discharge_plot, [0], color=:black, lw=2, label=false)

        # Create a plot for SOC
        soc_plot = bar(
            0:T, soc,
            dpi=600,
            label=L"SOC",
            color=:purple,
            legend=:bottomleft,
            xlabel="Time Interval " * L"t",
            ylabel="SOC "*L"[\%]",
            ylim=(soc_min[j]*100*0.95, soc_max[j]*100*1.10),
            # ylim=(0, 100),
            xticks=0:T,
            yticks=5*(div(soc_min[j] * 100 * 0.95, 5)-1):10:5*(div(soc_max[j] * 100 * 1.05, 5)+1),
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
            # left_margin=10mm,   # Adds space to the left
            # right_margin=10mm,  # Adds space to the right
            top_margin=0mm,    # Adds space at the top
            # bottom_margin=10mm  # Adds space at the bottom
        )

        # Combine the two plots in a layout
        plot_combined = plot(charging_discharge_plot, soc_plot, layout=@layout([a; b]))

        # Show the plot if `showPlots` is true
        if showPlots
            display(plot_combined)
        end

        # Save the plot if `savePlots` is true
        if savePlots
            @unpack alphaAppendix, gammaAppendix = data;
            filename = joinpath(base_dir, "Battery_$(j)_alpha_$(alphaAppendix)_gamma_$(gammaAppendix)_$(solver).png")
            myprintln(verbose, "Saving plot to: $filename")
            savefig(plot_combined, filename)
        end

    end
end

# Function to plot Substation Power over time
function plot_substation_power(data;
    showPlots::Bool=false,
    savePlots::Bool=true,
    macroItrNum::Int=1,
    verbose::Bool=false)
    @unpack Tset, PSubs_vs_t_1toT_kW, T, simNatureString, gedString, objfunString, systemName, objfunPrefix, gedAppendix, solver = data

    yvalues = PSubs_vs_t_1toT_kW

    theme(common_theme)

    # Setup for saving plot
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_1")
    if savePlots && !isdir(base_dir)
        myprintln(verbose, "Creating directory: $base_dir")
        mkpath(base_dir)  # Create the directory and its parents if needed
    end
    @unpack objfunConciseDescription, alphaAppendix, gammaAppendix = data;
    filename = joinpath(base_dir, "Horizon_$(T)_$(solver)_SubstationRealPowers_vs_t_$(gedAppendix)_for_$(objfunConciseDescription)_alpha_$(alphaAppendix)_gamma_$(gammaAppendix).png")

    gr()

    outputPlot = plot(
        Tset, PSubs_vs_t_1toT_kW,
        dpi=600,
        label=L"(P^t_{Subs})",
        xlabel="Time Period " * L"(t)",
        ylabel=L"P_{Subs} \, [kW]",
        title="Substation Power " * L"(P^t_{Subs})" * " across the Horizon\n" *
        "using $(simNatureString) OPF\n" *
        "with $(gedString)\n" *
        "optimizing for $(objfunString)",
        legend=:topleft,
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.05,
        lw=4,
        marker=common_marker_face,
        markersize=4,
        markerstrokecolor=common_marker_stroke_color,
        markerstrokewidth=common_marker_stroke_width,
        xlims=(0, T + 1),
        xticks=1:T,  # Set xticks for every hour from 1 to T
        ylims=(minimum(yvalues) * 0.95, maximum(yvalues) * 1.05),
        titlefont=font(8, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        # left_margin=10mm,   # Adds space to the left
        # right_margin=10mm,  # Adds space to the right
        top_margin=5mm,    # Adds space at the top
        # bottom_margin=10mm  # Adds space at the bottom
    )

    # Show the plot if `showPlots` is true
    if showPlots
        display(outputPlot)
    end

    # Save the plot if `savePlots` is true
    if savePlots
        myprintln(verbose, "Saving plot to: $filename")
        savefig(outputPlot, filename)
    end
end

function plot_substation_power_cost(data;
    showPlots::Bool=false,
    savePlots::Bool=true,
    macroItrNum::Int=1,
    verbose::Bool=false)

    theme(common_theme)

    @unpack Tset, PSubsCost_vs_t_1toT_dollar, T, simNatureString, gedString, objfunString, systemName, gedAppendix, solver = data

    yvalues = PSubsCost_vs_t_1toT_dollar

    gr()

    outputPlot = plot(
        Tset, yvalues,
        dpi=600,
        label=L"(P^t_{SubsCost})",
        xlabel="Time Period " * L"(t)",
        ylabel="Substation Power Cost " * L"[$]",
        title="Substation Power Cost " * L"(P^t_{SubsCost})" * " across the Horizon\n" * "using $(simNatureString) OPF\n" * "with $(gedString)\n" * "optimizing for $(objfunString)",
        legend=:topleft,
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.05,
        lw=4,
        marker=common_marker_face,
        markersize=4,
        markerstrokecolor=common_marker_stroke_color,
        markerstrokewidth=common_marker_stroke_width,
        xlims=(0, T + 1),
        xticks=1:1:T,
        ylims=(minimum(yvalues) * 0.95, maximum(yvalues) * 1.05),
        titlefont=font(8, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        top_margin=5mm,    # Adds space at the top
    )

    if showPlots
        display(outputPlot)
    end

    # Saving plot if requested
    if savePlots
        base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_1")
        if !isdir(base_dir)
            myprintln(verbose, "Creating directory: $base_dir")
            mkpath(base_dir)
        end
        @unpack objfunConciseDescription, alphaAppendix, gammaAppendix = data;
        filename = joinpath(base_dir, "Horizon_$(T)_$(solver)_SubstationPowerCost_vs_t__$(gedAppendix)_for_$(objfunConciseDescription)_alpha_$(alphaAppendix)_gamma_$(gammaAppendix).png")
        myprintln(verbose, "Saving plot to: $filename")
        savefig(outputPlot, filename)
    end
end

function plot_line_losses(data;
    showPlots::Bool=false,
    savePlots::Bool=true,
    macroItrNum::Int=1,
    verbose::Bool=false)

    @unpack numAreas, Tset, PLoss_vs_t_1toT_kW, T, simNatureString, gedString, objfunString, systemName, gedAppendix, solver = data

    yvalues = PLoss_vs_t_1toT_kW

    theme(common_theme)

    gr()

    outputPlot = plot(
        Tset, yvalues,
        dpi=600,
        label=L"(P^t_{Loss})",
        xlabel="Time Period " * L"(t)",
        ylabel="Line Losses " * L"[kW]",
        title="Line Losses " * L"(P^t_{Loss})" * " across the Horizon\n" * "using $(simNatureString) OPF\n" * "with $(gedString)\n" * "optimizing for $(objfunString)",
        legend=:topleft,
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.05,
        lw=4,
        marker=common_marker_face,
        markersize=4,
        markerstrokecolor=common_marker_stroke_color,
        markerstrokewidth=common_marker_stroke_width,
        xlims=(0, T + 1),
        xticks=1:1:T,
        ylims=(minimum(yvalues) * 0.95, maximum(yvalues) * 1.05),
        titlefont=font(8, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        # left_margin=10mm,   # Adds space to the left
        # right_margin=10mm,  # Adds space to the right
        top_margin=5mm,    # Adds space at the top
        # bottom_margin=10mm  # Adds space at the bottom
    )

    if showPlots
        display(outputPlot)
    end

    # Saving plot if requested
    if savePlots
        base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)",)
        if !isdir(base_dir)
            myprintln(verbose, "Creating directory: $base_dir")
            mkpath(base_dir)
        end
        @unpack objfunConciseDescription, alphaAppendix, gammaAppendix = data;
        filename = joinpath(base_dir, "Horizon_$(T)_$(solver)_LineLosses_vs_t__$(gedAppendix)_for_$(objfunConciseDescription)_alpha_$(alphaAppendix)_gamma_$(gammaAppendix).png")
        myprintln(verbose, "Saving plot to: $filename")
        savefig(outputPlot, filename)
    end
end

function plot_input_forecast_curves(data; showPlots::Bool=false, savePlots::Bool=true, filename::String="input_forecast_curves.png",
filenameSuffix::String="nonspecific", verbose::Bool=false)

    @unpack LoadShapeLoad, LoadShapePV, LoadShapeCost, T = data

    # Prepare data for plotting
    time_steps = 1:T
    load_cost_cents = LoadShapeCost .* 100  # Convert from $/kWh to cents/kWh

    # Calculate y-axis limits
    # left_min = 0.95 * minimum([minimum(LoadShapeLoad), minimum(LoadShapePV)])
    left_min = -0.05
    left_max = 1.05 * maximum([maximum(LoadShapeLoad), maximum(LoadShapePV)])
    right_min = floor(0.95 * minimum(load_cost_cents))
    right_max = ceil(1.05 * maximum(load_cost_cents))

    gr()
    theme(common_theme)
    # Plot LoadShapeLoad and LoadShapePV on the primary (left) y-axis
    outputPlot = plot(
        time_steps, LoadShapeLoad,
        dpi=600,
        label="Loading Factor "*L"(\lambda^t)",
        xlabel="Time Period "*L"(t)",
        ylabel="Loading/Irradiance Factor [Dimensionless]",
        legend=:left,
        lw=3,
        color=:darkgoldenrod2,
        markershape=:square,
        markerstrokecolor=common_marker_stroke_color,
        markerstrokewidth=common_marker_stroke_width,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.05,
        ylims=(left_min, left_max),
        xticks=1:T,
        title="Forecast Curves for Load, PV, and Cost",
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        right_margin=1cm
    )

    plot!(time_steps, LoadShapePV, 
        label="Solar Irradiance "*L"(\lambda^t_{PV})", 
        lw=3, 
        color=:orangered,
        markershape=:utriangle,
        markersize=8,
        markerstrokecolor=common_marker_stroke_color,
        markerstrokewidth=common_marker_stroke_width
    )

    # Use twinx() for the secondary y-axis for LoadShapeCost
    ax2 = twinx()
    plot!(
        ax2, time_steps, load_cost_cents,
        label="Substation Power Cost "*L"(C^t)",
        lw=3,
        color=:darkgreen,
        linestyle=:solid,
        markershape=:diamond,
        markerstrokecolor=common_marker_stroke_color,
        markerstrokewidth=common_marker_stroke_width,
        markersize=7,
        ylabel="Cost [cents/kWh]",
        ylims=(right_min, right_max),
        legend=:right,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )

    # Show the plot if requested
    if showPlots
        display(outputPlot)
    end

    # Saving plot if requested
    if savePlots
        @unpack systemName, gedAppendix = data;
        base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)")
        if !isdir(base_dir)
            myprintln(verbose, "Creating directory: $base_dir")
            mkpath(base_dir)
        end
        filename = joinpath(base_dir, "Horizon_$(T)_InputForecastCurves")
        ext = ".png"
        filename = filename*"_"*filenameSuffix*ext
        myprintln(verbose, "Saving plot to: $filename")
        savefig(outputPlot, filename)
    end
end


end  # module Plotter
