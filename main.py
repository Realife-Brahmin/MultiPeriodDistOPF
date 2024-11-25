
## new main
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))
import numpy as np
from src.MultiPeriodDistOPF import parse_all_data, optimize_MPOPF_1ph_NL

# from parameters import unpack

from src.exporter import (
    export_optimization_model,
    export_decision_variables,
    export_simulation_key_results_txt)
from src.computeOutputs import compute_output_values
from src.Plotter.Plotter import (plot_battery_actions,
    plot_input_forecast_curves,
    plot_substation_power,
    plot_substation_power_cost,
    plot_line_losses
)

# Input Configuration
system_name = "ads10_1ph"
# system_name = "ieee123_1ph"
T0 = 24
# factor = 1 / 8
factor = 1
T = int(T0 * factor)
num_areas = 1
temporal_decmp = False
# objfun0 = "powerflow"
# objfun0 = "lineLossMin"
objfun0 = "subsPowerCostMin"
# objfun2 = "none"
objfun2 = "scd"
input_forecast_description = "bilevelCosts"
alpha = 1e-3
PSubsMax_kW = float("inf")  # Inf means no limit
# solver = "Ipopt"
# solver = "EAGO"
solver = "Gurobi"
# solver = "Juniper"
# solver = "MadNLP"

# Parse all data
data = parse_all_data(
    system_name,
    T,
    num_areas=num_areas,
    alpha=alpha,
    objfun0=objfun0,
    objfun2=objfun2,
    temporal_decmp=temporal_decmp,
    PSubsMax_kW=PSubsMax_kW,
    input_forecast_description=input_forecast_description,
    solver=solver,
)

# Optimize
model = optimize_MPOPF_1ph_NL(data)

# Post-simulation computation, plotting, and logging
def post_simulation_tasks(model, data):
    verbose = False
    # verbose = True

    # Export optimization model
    export_optimization_model(model, data, verbose=verbose)

    # Compute output values
    data = compute_output_values(model, data, verbose=verbose)

    # Export decision variables
    export_decision_variables(model, data, verbose=verbose)

    # Export simulation key results
    export_simulation_key_results_txt(model, data, verbose=verbose)

    # Plot results
    save_plots = False
    save_plots = True

    plot_battery_actions(model, data, show_plots=False, save_plots=save_plots, verbose=verbose)
    plot_input_forecast_curves(
        data, filename_suffix=input_forecast_description, show_plots=False, verbose=verbose
    )
    plot_substation_power(data, save_plots=save_plots, verbose=verbose)
    plot_substation_power_cost(data, save_plots=save_plots, verbose=verbose)
    plot_line_losses(data, save_plots=save_plots, verbose=verbose)


# Run post-simulation tasks
post_simulation_tasks(model, data)