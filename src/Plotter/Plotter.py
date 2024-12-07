# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator

# # Common configurations
# common_theme = 'muted'
# common_marker_face = 'o'
# common_marker_stroke_color = 'black'
# common_marker_stroke_width = 2.0

# def ensure_dir(directory, verbose=False):
#     if not os.path.exists(directory):
#         if verbose:
#             print(f"Creating directory: {directory}")
#         os.makedirs(directory)

# def plot_battery_actions(model_dict, show_plots=False, save_plots=True, macro_itr_num=1, verbose=False):
#     data, model_vals = model_dict['data'], model_dict['modelVals']
#     Tset = sorted(list(data['Tset']))
#     Bset = data['Bset']
#     kVA_B = data['kVA_B']
#     Bref_pu, B_R_pu = data['Bref_pu'], data['B_R_pu']
#     system_name, ged_appendix = data['systemName'], data['gedAppendix']
#     P_c, P_d, B = model_vals['P_c'], model_vals['P_d'], model_vals['B']

#     base_dir = os.path.join("processedData", system_name, ged_appendix, f"Horizon_{data['T']}", f"numAreas_{data['numAreas']}", "batteryActionPlots", f"macroItr_{macro_itr_num}")
#     if save_plots:
#         ensure_dir(base_dir, verbose)

#     for j in Bset:
#         time_intervals = list(Tset)
#         charging_power_kW = [P_c[j, t] * kVA_B for t in Tset]
#         discharging_power_kW = [-P_d[j, t] * kVA_B for t in Tset]
#         soc = [Bref_pu[j] / B_R_pu[j] * 100] + [B[j, t] / B_R_pu[j] * 100 for t in Tset]

#         y_limit = (-data['P_B_R'][j], data['P_B_R'][j])

#         # Charging and Discharging Plot
#         plt.figure()
#         plt.bar(time_intervals, charging_power_kW, label="Charging", color="green")
#         plt.bar(time_intervals, discharging_power_kW, label="Discharging", color="darkred")
#         plt.axhline(0, color="black", linewidth=2)
#         plt.xlabel("Time Interval t")
#         plt.ylabel("P_c/P_d [kW]")
#         plt.ylim(y_limit)
#         plt.xticks(range(1, data['T'] + 1))
#         plt.title(f"Battery at Bus {j}\nCharging and Discharging")
#         plt.legend(loc="lower left")

#         if save_plots:
#             filename = os.path.join(base_dir, f"Battery_{j}_ActionPlot.png")
#             if verbose:
#                 print(f"Saving plot to: {filename}")
#             plt.savefig(filename)
#         if show_plots:
#             plt.show()
#         plt.close()

#         # SOC Plot
#         plt.figure()
#         plt.bar([0] + time_intervals, soc, color="purple", label="SOC")
#         plt.xlabel("Time Interval t")
#         plt.ylabel("SOC [%]")
#         plt.xticks(range(0, data['T'] + 1))
#         plt.ylim(data['soc_min'][j] * 100 * 0.95, data['soc_max'][j] * 100 * 1.10)
#         plt.title("SOC")
#         plt.legend(loc="lower left")

#         if save_plots:
#             filename = os.path.join(base_dir, f"Battery_{j}_SOCPlot.png")
#             if verbose:
#                 print(f"Saving plot to: {filename}")
#             plt.savefig(filename)
#         if show_plots:
#             plt.show()
#         plt.close()

# def plot_substation_power(model_dict, show_plots=False, save_plots=True, macro_itr_num=1, verbose=False):
#     data = model_dict['data']
#     Tset, PSubs_vs_t_1toT_kW = data['Tset'], data['PSubs_vs_t_1toT_kW']

#     base_dir = os.path.join("processedData", data['systemName'], data['gedAppendix'], f"Horizon_{data['T']}", "numAreas_1")
#     if save_plots:
#         ensure_dir(base_dir, verbose)

#     plt.figure()
#     plt.plot(Tset, PSubs_vs_t_1toT_kW, label="P_t_Subs", marker=common_marker_face, color="blue")
#     plt.xlabel("Time Period t")
#     plt.ylabel("P_Subs [kW]")
#     plt.title("Substation Power across the Horizon")
#     plt.legend(loc="upper left")

#     if save_plots:
#         filename = os.path.join(base_dir, "SubstationPowerPlot.png")
#         if verbose:
#             print(f"Saving plot to: {filename}")
#         plt.savefig(filename)
#     if show_plots:
#         plt.show()
#     plt.close()

# def plot_substation_power_cost(model_dict, show_plots=False, save_plots=True, macro_itr_num=1, verbose=False):
#     data = model_dict['data']
#     Tset, PSubsCost_vs_t_1toT_dollar = data['Tset'], data['PSubsCost_vs_t_1toT_dollar']

#     base_dir = os.path.join("processedData", data['systemName'], data['gedAppendix'], f"Horizon_{data['T']}", "numAreas_1")
#     if save_plots:
#         ensure_dir(base_dir, verbose)

#     plt.figure()
#     plt.plot(Tset, PSubsCost_vs_t_1toT_dollar, label="P_t_SubsCost", marker=common_marker_face, color="green")
#     plt.xlabel("Time Period t")
#     plt.ylabel("Substation Power Cost [$]")
#     plt.title("Substation Power Cost across the Horizon")
#     plt.legend(loc="upper left")

#     if save_plots:
#         filename = os.path.join(base_dir, "SubstationPowerCostPlot.png")
#         if verbose:
#             print(f"Saving plot to: {filename}")
#         plt.savefig(filename)
#     if show_plots:
#         plt.show()
#     plt.close()

# def plot_line_losses(model_dict, show_plots=False, save_plots=True, macro_itr_num=1, verbose=False):
#     data = model_dict['data']
#     Tset, PLoss_vs_t_1toT_kW = data['Tset'], data['PLoss_vs_t_1toT_kW']

#     base_dir = os.path.join("processedData", data['systemName'], data['gedAppendix'], f"Horizon_{data['T']}", f"numAreas_{data['numAreas']}")
#     if save_plots:
#         ensure_dir(base_dir, verbose)

#     plt.figure()
#     plt.plot(Tset, PLoss_vs_t_1toT_kW, label="P_t_Loss", marker=common_marker_face, color="red")
#     plt.xlabel("Time Period t")
#     plt.ylabel("Line Losses [kW]")
#     plt.title("Line Losses across the Horizon")
#     plt.legend(loc="upper left")

#     if save_plots:
#         filename = os.path.join(base_dir, "LineLossesPlot.png")
#         if verbose:
#             print(f"Saving plot to: {filename}")
#         plt.savefig(filename)
#     if show_plots:
#         plt.show()
#     plt.close()

# import os
# import numpy as np
# import matplotlib.pyplot as plt

# # Function to plot input forecast curves
# def plot_input_forecast_curves(data, showPlots=False, savePlots=True, filename="input_forecast_curves.png",
#                                filenameSuffix="nonspecific", verbose=False):
#     LoadShapeLoad = data['LoadShapeLoad']
#     LoadShapePV = data['LoadShapePV']
#     LoadShapeCost = data['LoadShapeCost']
#     T = data['T']

#     # Prepare data for plotting
#     time_steps = np.arange(1, T + 1)
#     load_cost_cents = LoadShapeCost * 100  # Convert from $/kWh to cents/kWh

#     # Calculate y-axis limits
#     left_min = -0.05
#     left_max = 1.05 * max(np.max(LoadShapeLoad), np.max(LoadShapePV))
#     right_min = np.floor(0.95 * np.min(load_cost_cents))
#     right_max = np.ceil(1.05 * np.max(load_cost_cents))

#     # Plot LoadShapeLoad and LoadShapePV on the primary (left) y-axis
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     ax1.plot(time_steps, LoadShapeLoad, label="Loading Factor (λ^t)",
#              linewidth=3, color="goldenrod", marker="s", markersize=6, markeredgecolor="black", markeredgewidth=1.5)
#     ax1.plot(time_steps, LoadShapePV, label="Solar Irradiance (λ^t_{PV})",
#              linewidth=3, color="orangered", marker="^", markersize=8, markeredgecolor="black", markeredgewidth=1.5)

#     ax1.set_xlabel("Time Period (t)")
#     ax1.set_ylabel("Loading/Irradiance Factor [Dimensionless]")
#     ax1.set_ylim([left_min, left_max])
#     ax1.set_xticks(time_steps)
#     ax1.legend(loc="upper left")
#     ax1.grid(alpha=0.3)
#     ax1.set_title("Forecast Curves for Load, PV, and Cost")

#     # Use a secondary y-axis for LoadShapeCost
#     ax2 = ax1.twinx()
#     ax2.plot(time_steps, load_cost_cents, label="Substation Power Cost (C^t)",
#              linewidth=3, color="darkgreen", marker="D", markersize=7, markeredgecolor="black", markeredgewidth=1.5)
#     ax2.set_ylabel("Cost [cents/kWh]")
#     ax2.set_ylim([right_min, right_max])
#     ax2.legend(loc="upper right")

#     # Show the plot if requested
#     if showPlots:
#         plt.show()

#     # Saving plot if requested
#     if savePlots:
#         systemName = data.get('systemName', 'defaultSystem')
#         gedAppendix = data.get('gedAppendix', 'defaultAppendix')
#         base_dir = os.path.join("processedData", systemName, gedAppendix, f"Horizon_{T}")
#         if not os.path.isdir(base_dir):
#             if verbose:
#                 print(f"Creating directory: {base_dir}")
#             os.makedirs(base_dir)

#         full_filename = os.path.join(base_dir, f"Horizon_{T}_InputForecastCurves_{filenameSuffix}.png")
#         if verbose:
#             print(f"Saving plot to: {full_filename}")
#         plt.savefig(full_filename, dpi=600)

#     plt.close(fig)

# Plotter.py

import os
import matplotlib.pyplot as plt
from pathlib import Path
from src.helperFunctions import myprintln  # Assuming this exists and behaves like in Julia

# Constants as in Julia code
common_theme = "mute"  # Not directly translatable to matplotlib theme
common_marker_face = 'o'  # circle
common_marker_stroke_color = 'black'
common_marker_stroke_width = 2.0

def plot_battery_actions(modelDict, showPlots=False, savePlots=True, macroItrNum=1, verbose=False):
    data = modelDict['data']
    modelVals = modelDict['modelVals']

    # Unpack
    Tset = data['Tset']
    Bset = data['Bset']
    kVA_B = data['kVA_B']
    B_R_pu = data['B_R_pu']
    P_B_R = data['P_B_R']
    Bref_pu = data['Bref_pu']
    systemName = data['systemName']
    numAreas = data['numAreas']
    T = data['T']
    alphaAppendix = data['alphaAppendix']
    gammaAppendix = data['gammaAppendix']
    soc_min = data['soc_min']
    soc_max = data['soc_max']
    gedAppendix = data['gedAppendix']
    solver = data['solver']

    # Sort Tset
    Tset = sorted(list(Tset))

    P_c = modelVals['P_c']
    P_d = modelVals['P_d']
    B = modelVals['B']

    base_dir = os.path.join("processedData", systemName, gedAppendix, f"Horizon_{T}", f"numAreas_{numAreas}",
                            "batteryActionPlots", f"macroItr_{macroItrNum}")
    if savePlots and not os.path.isdir(base_dir):
        myprintln(verbose, f"Creating directory: {base_dir}")
        Path(base_dir).mkdir(parents=True, exist_ok=True)

    # For each battery
    for j in Bset:
        time_intervals = list(Tset)

        charging_power_kW = [P_c[(j, t)] * kVA_B for t in Tset]
        discharging_power_kW = [P_d[(j, t)] * kVA_B for t in Tset]
        soc = [Bref_pu[j]/B_R_pu[j]*100] + [B[(j, t)]/B_R_pu[j]*100 for t in Tset]

        ylimit = (-P_B_R[j], P_B_R[j])

        # Two subplots (charging/discharging and SOC)
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, dpi=600)
        fig.suptitle(f"Battery at Bus {j}")

        # Charging/Discharging plot
        ax1.bar(time_intervals, charging_power_kW, label="Charging", color='green')
        ax1.bar(time_intervals, [-x for x in discharging_power_kW], label="Discharging", color='darkred')
        ax1.axhline(0, color='black', lw=2)
        ax1.set_xlabel("Time Interval $t$")
        ax1.set_ylabel("$P_c/P_d$ [kW]")
        ax1.set_ylim(ylimit)
        ax1.set_xticks(range(1, T+1))
        ax1.grid(True, linestyle='solid', linewidth=1.0, alpha=0.2)
        ax1.legend(loc='lower left')
        ax1.set_title("Charging and Discharging", fontsize=12)

        # SOC plot
        ax2.bar(range(0, T+1), soc, label="$SOC$", color='purple')
        ax2.set_xlabel("Time Interval $t$")
        ax2.set_ylabel("SOC [%]")
        ax2.set_ylim(soc_min[j]*100*0.95, soc_max[j]*100*1.10)
        ax2.set_xticks(range(0, T+1))
        ax2.grid(True, linestyle='solid', linewidth=1.0, alpha=0.2)
        ax2.legend(loc='lower left')
        ax2.set_title("SOC", fontsize=12)

        fig.tight_layout()

        if showPlots:
            plt.show()

        if savePlots:
            filename = os.path.join(base_dir, f"Battery_{j}_alpha_{alphaAppendix}_gamma_{gammaAppendix}_{solver}.png")
            myprintln(verbose, f"Saving plot to: {filename}")
            plt.savefig(filename)

        plt.close(fig)

def plot_substation_power(modelDict, showPlots=False, savePlots=True, macroItrNum=1, verbose=False):
    data = modelDict['data']

    Tset = data['Tset']
    PSubs_vs_t_1toT_kW = data['PSubs_vs_t_1toT_kW']
    T = data['T']
    simNatureString = data['simNatureString']
    gedString = data['gedString']
    objfunString = data['objfunString']
    systemName = data['systemName']
    gedAppendix = data['gedAppendix']
    solver = data['solver']
    objfunConciseDescription = data['objfunConciseDescription']
    alphaAppendix = data['alphaAppendix']
    gammaAppendix = data['gammaAppendix']

    yvalues = PSubs_vs_t_1toT_kW

    base_dir = os.path.join("processedData", systemName, gedAppendix, f"Horizon_{T}", "numAreas_1")
    if savePlots and not os.path.isdir(base_dir):
        myprintln(verbose, f"Creating directory: {base_dir}")
        Path(base_dir).mkdir(parents=True, exist_ok=True)

    filename = os.path.join(base_dir, f"Horizon_{T}_{solver}_SubstationRealPowers_vs_t_{gedAppendix}_for_{objfunConciseDescription}_alpha_{alphaAppendix}_gamma_{gammaAppendix}.png")

    fig, ax = plt.subplots(dpi=600)
    ax.plot(Tset, yvalues, label="$P^t_{Subs}$", lw=4, marker=common_marker_face, ms=4,
            mfc='none', mec=common_marker_stroke_color, mew=common_marker_stroke_width)
    ax.set_xlabel("Time Period $(t)$")
    ax.set_ylabel("$P_{Subs} [kW]$")
    ax.set_title("Substation Power $(P^t_{Subs})$ across the Horizon\n"
                 f"using {simNatureString} OPF\n"
                 f"with {gedString}\n"
                 f"optimizing for {objfunString}", fontsize=8)
    ax.legend(loc='upper left')
    ax.grid(True, linestyle='solid', linewidth=1.0, alpha=0.2)
    ax.set_xlim(0, T+1)
    ax.set_xticks(range(1, T+1))
    ax.set_ylim(min(yvalues)*0.95, max(yvalues)*1.05)

    fig.tight_layout()

    if showPlots:
        plt.show()

    if savePlots:
        myprintln(verbose, f"Saving plot to: {filename}")
        plt.savefig(filename)

    plt.close(fig)

def plot_substation_power_cost(modelDict, showPlots=False, savePlots=True, macroItrNum=1, verbose=False):
    data = modelDict['data']

    Tset = data['Tset']
    PSubsCost_vs_t_1toT_dollar = data['PSubsCost_vs_t_1toT_dollar']
    T = data['T']
    simNatureString = data['simNatureString']
    gedString = data['gedString']
    objfunString = data['objfunString']
    systemName = data['systemName']
    gedAppendix = data['gedAppendix']
    solver = data['solver']
    objfunConciseDescription = data['objfunConciseDescription']
    alphaAppendix = data['alphaAppendix']
    gammaAppendix = data['gammaAppendix']

    yvalues = PSubsCost_vs_t_1toT_dollar

    base_dir = os.path.join("processedData", systemName, gedAppendix, f"Horizon_{T}", "numAreas_1")
    if savePlots and not os.path.isdir(base_dir):
        myprintln(verbose, f"Creating directory: {base_dir}")
        Path(base_dir).mkdir(parents=True, exist_ok=True)

    filename = os.path.join(base_dir, f"Horizon_{T}_{solver}_SubstationPowerCost_vs_t__{gedAppendix}_for_{objfunConciseDescription}_alpha_{alphaAppendix}_gamma_{gammaAppendix}.png")

    fig, ax = plt.subplots(dpi=600)
    ax.plot(Tset, yvalues, label="$(P^t_{SubsCost})$", lw=4, marker=common_marker_face, ms=4,
            mfc='none', mec=common_marker_stroke_color, mew=common_marker_stroke_width)
    ax.set_xlabel("Time Period $(t)$")
    ax.set_ylabel("Substation Power Cost [$]")
    ax.set_title("Substation Power Cost $(P^t_{SubsCost})$ across the Horizon\n"
                 f"using {simNatureString} OPF\n"
                 f"with {gedString}\n"
                 f"optimizing for {objfunString}", fontsize=8)
    ax.legend(loc='upper left')
    ax.grid(True, linestyle='solid', linewidth=1.0, alpha=0.2)
    ax.set_xlim(0, T+1)
    ax.set_xticks(range(1, T+1))
    ax.set_ylim(min(yvalues)*0.95, max(yvalues)*1.05)

    fig.tight_layout()

    if showPlots:
        plt.show()

    if savePlots:
        myprintln(verbose, f"Saving plot to: {filename}")
        plt.savefig(filename)

    plt.close(fig)

def plot_line_losses(modelDict, showPlots=False, savePlots=True, macroItrNum=1, verbose=False):
    data = modelDict['data']

    numAreas = data['numAreas']
    Tset = data['Tset']
    PLoss_vs_t_1toT_kW = data['PLoss_vs_t_1toT_kW']
    T = data['T']
    simNatureString = data['simNatureString']
    gedString = data['gedString']
    objfunString = data['objfunString']
    systemName = data['systemName']
    gedAppendix = data['gedAppendix']
    solver = data['solver']
    objfunConciseDescription = data['objfunConciseDescription']
    alphaAppendix = data['alphaAppendix']
    gammaAppendix = data['gammaAppendix']

    yvalues = PLoss_vs_t_1toT_kW

    base_dir = os.path.join("processedData", systemName, gedAppendix, f"Horizon_{T}", f"numAreas_{numAreas}")
    if savePlots and not os.path.isdir(base_dir):
        myprintln(verbose, f"Creating directory: {base_dir}")
        Path(base_dir).mkdir(parents=True, exist_ok=True)

    filename = os.path.join(base_dir, f"Horizon_{T}_{solver}_LineLosses_vs_t__{gedAppendix}_for_{objfunConciseDescription}_alpha_{alphaAppendix}_gamma_{gammaAppendix}.png")

    fig, ax = plt.subplots(dpi=600)
    ax.plot(Tset, yvalues, label="$(P^t_{Loss})$", lw=4, marker=common_marker_face, ms=4,
            mfc='none', mec=common_marker_stroke_color, mew=common_marker_stroke_width)
    ax.set_xlabel("Time Period $(t)$")
    ax.set_ylabel("Line Losses [kW]")
    ax.set_title("Line Losses $(P^t_{Loss})$ across the Horizon\n"
                 f"using {simNatureString} OPF\n"
                 f"with {gedString}\n"
                 f"optimizing for {objfunString}", fontsize=8)
    ax.legend(loc='upper left')
    ax.grid(True, linestyle='solid', linewidth=1.0, alpha=0.2)
    ax.set_xlim(0, T+1)
    ax.set_xticks(range(1, T+1))
    ax.set_ylim(min(yvalues)*0.95, max(yvalues)*1.05)

    fig.tight_layout()

    if showPlots:
        plt.show()

    if savePlots:
        myprintln(verbose, f"Saving plot to: {filename}")
        plt.savefig(filename)

    plt.close(fig)

def plot_input_forecast_curves(data, showPlots=False, savePlots=True, filename="input_forecast_curves.png",
                               filenameSuffix="nonspecific", verbose=False):

    LoadShapeLoad = data['LoadShapeLoad']
    LoadShapePV = data['LoadShapePV']
    LoadShapeCost = data['LoadShapeCost']
    T = data['T']

    systemName = data['systemName']
    gedAppendix = data['gedAppendix']

    time_steps = range(1, T+1)
    load_cost_cents = [c*100 for c in LoadShapeCost]

    left_min = -0.05
    left_max = 1.05 * max(max(LoadShapeLoad), max(LoadShapePV))
    right_min = int(min(load_cost_cents)*0.95)
    right_max = int(max(load_cost_cents)*1.05)

    fig, ax1 = plt.subplots(dpi=600)
    ax1.plot(time_steps, LoadShapeLoad, label="Loading Factor $(\\lambda^t)$",
             lw=3, marker='s', mfc='none', mec=common_marker_stroke_color, mew=common_marker_stroke_width)
    ax1.plot(time_steps, LoadShapePV, label="Solar Irradiance $(\\lambda^t_{PV})$",
             lw=3, marker='^', mfc='none', mec=common_marker_stroke_color, mew=common_marker_stroke_width)
    ax1.set_xlabel("Time Period $(t)$")
    ax1.set_ylabel("Loading/Irradiance Factor [Dimensionless]")
    ax1.set_ylim(left_min, left_max)
    ax1.set_xticks(range(1, T+1))
    ax1.grid(True, linestyle='solid', linewidth=1.0, alpha=0.3)
    ax1.legend(loc='lower left')
    ax1.set_title("Forecast Curves for Load, PV, and Cost", fontsize=12)

    ax2 = ax1.twinx()
    ax2.plot(time_steps, load_cost_cents, label="Substation Power Cost $(C^t)$",
             lw=3, marker='D', mfc='none', mec=common_marker_stroke_color, mew=common_marker_stroke_width)
    ax2.set_ylabel("Cost [cents/kWh]")
    ax2.set_ylim(right_min, right_max)
    ax2.legend(loc='lower right')

    fig.tight_layout()

    if showPlots:
        plt.show()

    if savePlots:
        base_dir = os.path.join("processedData", systemName, gedAppendix, f"Horizon_{T}")
        if not os.path.isdir(base_dir):
            myprintln(verbose, f"Creating directory: {base_dir}")
            Path(base_dir).mkdir(parents=True, exist_ok=True)
        filename = os.path.join(base_dir, f"Horizon_{T}_InputForecastCurves_{filenameSuffix}.png")
        myprintln(verbose, f"Saving plot to: {filename}")
        plt.savefig(filename)

    plt.close(fig)

