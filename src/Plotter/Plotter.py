import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pyomo.environ import value

# Common settings
COMMON_THEME = "mute"
COMMON_MARKER_FACE = "o"
COMMON_MARKER_STROKE_COLOR = "black"
COMMON_MARKER_STROKE_WIDTH = 2.0

def myprintln(verbose, message):
    if verbose:
        print(message)

def create_directory_if_not_exists(directory, verbose=False):
    if not os.path.isdir(directory):
        myprintln(verbose, f"Creating directory: {directory}")
        os.makedirs(directory, exist_ok=True)

def save_plot(figure, filename, verbose=False):
    myprintln(verbose, f"Saving plot to: {filename}")
    figure.savefig(filename, dpi=600, bbox_inches="tight")
    plt.close(figure)

def plot_battery_actions(model, data, show_plots=False, save_plots=True, macro_itr_num=1, verbose=False):
    Tset = sorted(data['Tset'])
    P_c = model.P_c
    P_d = model.P_d
    B = model.B

    base_dir = os.path.join(
        "processedData", data['system_name'], data['ged_appendix'],
        f"Horizon_{data['T']}", f"numAreas_{data['num_areas']}",
        "batteryActionPlots", f"macroItr_{macro_itr_num}"
    )
    if save_plots:
        create_directory_if_not_exists(base_dir, verbose)

    for j in data['Bset']:
        time_intervals = list(range(1, len(Tset) + 1))
        charging_power_kw = [value(P_c[j, t]) * data['kVA_B'] for t in Tset]
        discharging_power_kw = [-value(P_d[j, t]) * data['kVA_B'] for t in Tset]
        soc = [value(data['Bref_pu'][j]) / value(data['B_R_pu'][j]) * 100] + [
            value(B[j, t]) / value(data['B_R_pu'][j]) * 100 for t in Tset
        ]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Charging and Discharging Plot
        ax1.bar(time_intervals, charging_power_kw, color="green", label="Charging")
        ax1.bar(time_intervals, discharging_power_kw, color="darkred", label="Discharging")
        ax1.axhline(0, color="black", linewidth=2)
        ax1.set_xlabel("Time Interval (t)")
        ax1.set_ylabel("P_c/P_d [kW]")
        ax1.set_title(f"Battery at Bus {j}\nCharging and Discharging")
        ax1.legend(loc="upper left")
        ax1.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

        # SOC Plot
        ax2.bar(range(len(soc)), soc, color="purple", label="SOC")
        ax2.set_xlabel("Time Interval (t)")
        ax2.set_ylabel("SOC [%]")
        ax2.set_title("State of Charge")
        ax2.legend(loc="upper left")
        ax2.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

        if show_plots:
            plt.show()

        if save_plots:
            filename = os.path.join(base_dir, f"Battery_{j}_macroItr_{macro_itr_num}.png")
            save_plot(fig, filename, verbose)

def plot_substation_power(data, show_plots=False, save_plots=True, verbose=False):
    Tset = data['Tset']
    yvalues = data['PSubs_vs_t_1toT_kW']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(Tset, yvalues, marker=COMMON_MARKER_FACE, label="P^t_{Subs}")
    ax.set_xlabel("Time Period (t)")
    ax.set_ylabel("P_{Subs} [kW]")
    ax.set_title("Substation Power across the Horizon")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

    if show_plots:
        plt.show()

    if save_plots:
        base_dir = os.path.join("processedData", data['system_name'], data['ged_appendix'], f"Horizon_{data['T']}")
        create_directory_if_not_exists(base_dir, verbose)
        filename = os.path.join(base_dir, "SubstationPower_vs_t.png")
        save_plot(fig, filename, verbose)

def plot_substation_power_cost(data, show_plots=False, save_plots=True, verbose=False):
    Tset = data['Tset']
    yvalues = data['PSubsCost_vs_t_1toT_dollar']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(Tset, yvalues, marker=COMMON_MARKER_FACE, label="P^t_{SubsCost}")
    ax.set_xlabel("Time Period (t)")
    ax.set_ylabel("Substation Power Cost [$]")
    ax.set_title("Substation Power Cost across the Horizon")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

    if show_plots:
        plt.show()

    if save_plots:
        base_dir = os.path.join("processedData", data['system_name'], data['ged_appendix'], f"Horizon_{data['T']}")
        create_directory_if_not_exists(base_dir, verbose)
        filename = os.path.join(base_dir, "SubstationPowerCost_vs_t.png")
        save_plot(fig, filename, verbose)

def plot_line_losses(data, show_plots=False, save_plots=True, verbose=False):
    Tset = data['Tset']
    yvalues = data['PLoss_vs_t_1toT_kW']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(Tset, yvalues, marker=COMMON_MARKER_FACE, label="P^t_{Loss}")
    ax.set_xlabel("Time Period (t)")
    ax.set_ylabel("Line Losses [kW]")
    ax.set_title("Line Losses across the Horizon")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

    if show_plots:
        plt.show()

    if save_plots:
        base_dir = os.path.join("processedData", data['system_name'], data['ged_appendix'], f"Horizon_{data['T']}")
        create_directory_if_not_exists(base_dir, verbose)
        filename = os.path.join(base_dir, "LineLosses_vs_t.png")
        save_plot(fig, filename, verbose)

def plot_input_forecast_curves(data, show_plots=False, save_plots=True, filename_suffix="nonspecific", verbose=False):
    time_steps = list(range(1, data['T'] + 1))
    load_cost_cents = [cost * 100 for cost in data['LoadShapeCost']]

    left_min = -0.05
    left_max = 1.05 * max(max(data['LoadShapeLoad']), max(data['LoadShapePV']))
    right_min = 0.95 * min(load_cost_cents)
    right_max = 1.05 * max(load_cost_cents)

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot LoadShapeLoad and LoadShapePV on primary axis
    ax1.plot(time_steps, data['LoadShapeLoad'], color="darkgoldenrod", label="Loading Factor (lambda^t)", linewidth=3, marker="s")
    ax1.plot(time_steps, data['LoadShapePV'], color="orangered", label="Solar Irradiance (lambda^t_PV)", linewidth=3, marker="^")
    ax1.set_xlabel("Time Period (t)")
    ax1.set_ylabel("Loading/Irradiance Factor [Dimensionless]")
    ax1.set_ylim(left_min, left_max)
    ax1.legend(loc="upper left")
    ax1.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

    # Plot LoadShapeCost on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(time_steps, load_cost_cents, color="darkgreen", label="Substation Power Cost (C^t)", linewidth=3, marker="d")
    ax2.set_ylabel("Cost [cents/kWh]")
    ax2.set_ylim(right_min, right_max)
    ax2.legend(loc="upper right")

    plt.title("Forecast Curves for Load, PV, and Cost")

    if show_plots:
        plt.show()

    if save_plots:
        base_dir = os.path.join("processedData", data['system_name'], data['ged_appendix'], f"Horizon_{data['T']}")
        create_directory_if_not_exists(base_dir, verbose)
        filename = os.path.join(base_dir, f"Horizon_{data['T']}_InputForecastCurves_{filename_suffix}.png")
        save_plot(fig, filename, verbose)
