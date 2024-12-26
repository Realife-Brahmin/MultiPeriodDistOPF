# solveMPOPF.jl
using Revise
using MultiPeriodDistOPF
using Parameters

Revise.revise()

systemName = "ads10_1ph"
# systemName = "ieee123_1ph"
# T0 = 24
T0 = 6
# factor = 1/8
factor = 1
T = Int(T0*factor) 
numAreas = 1
temporal_decmp = false
temporal_decmp = true
savePlots = false
savePlots = true
# objfun0 = "powerflow"
# objfun0 = "lineLossMin"
objfun0 = "subsPowerCostMin"
objfun2 = "scd"
inputForecastDescription = "bilevelCosts"
tSOC_hard = false
# tSOC_hard = true
PSubsMax_kW = Inf # Inf means no limit
solver = "Ipopt"
# solver = "Gurobi"
# solver = "Juniper"

# Parse all data
data = parse_all_data(systemName, T, objfun0=objfun0, temporal_decmp=temporal_decmp, tSOC_hard=tSOC_hard)

if !temporal_decmp
    modelDict = optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)
    @unpack model, modelVals, data = modelDict
elseif temporal_decmp
    modelDict = optimize_MPOPF_1ph_NL_DDP(data) # modelDict is basically ddpModel
    @unpack modelVals, data = modelDict
else
    error("temporal_decmp must be either true or false")
end

# postsim computation, plotting, logging
begin

    verbose = false
    # verbose = true

    modelDict = compute_output_values(modelDict, verbose=verbose)
    # @unpack data = modelDict

    export_decision_variables(modelDict, verbose=verbose)

    # Todo: Maybe separately save the simulation times? It is annoying to have file content differences every single run (for same exact sim)
    export_simulation_key_results_txt(modelDict, verbose=verbose)

    plot_battery_actions(modelDict, showPlots=false, savePlots=savePlots, verbose=verbose)

    @unpack data = modelDict;
    plot_input_forecast_curves(data, filenameSuffix=inputForecastDescription, showPlots=false, verbose=verbose)
    plot_substation_power(modelDict, savePlots=savePlots, verbose=verbose)
    plot_substation_power_cost(modelDict, savePlots=savePlots, verbose=verbose)
    plot_line_losses(modelDict, savePlots=savePlots, verbose=verbose)

    modelDict = validate_opf_against_opendss(modelDict, verbose=verbose)

    export_validation_key_results(modelDict, verbose=verbose)

    export_validation_decision_variables(modelDict, verbose=verbose)
end

