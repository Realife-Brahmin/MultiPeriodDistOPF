## Solve MPOPF
from src.MultiPeriodDistOPF import *
# Equivalent parameter setup as in Julia
systemName = "ieee123_1ph"
T0 = 24
factor = 1/8
T = int(T0 * factor)
numAreas = 1
temporal_decmp = False
objfun0 = "subsPowerCostMin"
objfun2 = "scd"
inputForecastDescription = "bilevelCosts"
tSOC_hard = False
PSubsMax_kW = float('inf')
solver = "Ipopt"

# Parse all data
data = parse_all_data(
    systemName,
    T,
    numAreas=numAreas,
    objfun0=objfun0,
    objfun2=objfun2,
    temporal_decmp=temporal_decmp,
    PSubsMax_kW=PSubsMax_kW,
    inputForecastDescription=inputForecastDescription,
    solver=solver,
    tSOC_hard=tSOC_hard
)

# Solve the model
modelDict = optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)

# In Julia, we do:
# @unpack model, data = modelDict
# In Python, just access them directly if needed:
model = modelDict['model']
data = modelDict['data']

# Post-simulation computation, plotting, logging
verbose = False

modelDict = compute_output_values(modelDict, verbose=verbose)
export_decision_variables(modelDict, verbose=verbose)
export_simulation_key_results_txt(modelDict, verbose=verbose)

savePlots = True

plot_battery_actions(modelDict, showPlots=False, savePlots=savePlots, verbose=verbose)

# plot_input_forecast_curves expects data separately
plot_input_forecast_curves(data, filenameSuffix=inputForecastDescription, showPlots=False, verbose=verbose)

plot_substation_power(modelDict, savePlots=savePlots, verbose=verbose)
plot_substation_power_cost(modelDict, savePlots=savePlots, verbose=verbose)
plot_line_losses(modelDict, savePlots=savePlots, verbose=verbose)

modelDict = validate_opf_against_opendss(modelDict, verbose=verbose)

export_validation_key_results(modelDict, verbose=verbose)
export_validation_decision_variables(modelDict, verbose=verbose)
