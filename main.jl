# solveMPOPF.jl
using Revise
using MultiPeriodDistOPF
using Parameters

Revise.revise()

begin
    # systemName0 = "ads10_1ph"
    # systemName0 = "ieee123_1ph-A"
    # systemName0 = "ieee123_1ph-B"
    systemName0 = "ieee729_1ph"
    # systemName0 = "ieee730_1ph"
    # T0 = 3
    T0 = 24
    # T0 = 11
    factor = 1
    # factor = 1/2
    temporal_decmp = false
    # temporal_decmp = true
end;

begin
    factor = 1
    T = Int(T0*factor) 
    numAreas = 1
    linearizedModel = false
    # linearizedModel = true
    maxiter_ddp = 8
    savePlots = false
    savePlots = true
    # objfun0 = "lineLossMin"
    objfun0 = "subsPowerCostMin"
    objfun2 = "scd"
    inputForecastDescription = "bilevelCosts"
    relax_terminal_soc_constraint = false
    relax_terminal_soc_constraint = true
    tSOC_hard = false
    # tSOC_hard = true
    solver = "Ipopt"
    # solver = "Gurobi"
    # solver = "Juniper"
    if systemName0 == "ieee123_1ph-A"
        systemName = "ieee123_1ph"
        DER_Percent_ud = 20
        DER_Rating_factor_ud = 1/3
        Batt_Percent_ud = 30
        Batt_Rating_factor_ud = 1/3
    elseif systemName0 == "ieee123_1ph-B"
        systemName = "ieee123_1ph"
        DER_Percent_ud = 40
        DER_Rating_factor_ud = 1
        Batt_Percent_ud = 50
        Batt_Rating_factor_ud = 1
    elseif systemName0 == "ieee730_1ph" 
        systemName = "ieee730_1ph"
        DER_Percent_ud = 20
        DER_Rating_factor_ud = 1
        Batt_Percent_ud = 30
        Batt_Rating_factor_ud = 1
    elseif systemName0 == "ieee729_1ph"
        systemName = "ieee729_1ph"
        DER_Percent_ud = 40
        DER_Rating_factor_ud = 1
        Batt_Percent_ud = 40
        Batt_Rating_factor_ud = 1
    elseif systemName0 == "ads10_1ph"
        systemName = "ads10_1ph"
        DER_Percent_ud = 25
        DER_Rating_factor_ud = 1/3
        Batt_Percent_ud = 25
        Batt_Rating_factor_ud = 1/3
    else
        error("systemName must be either ieee123_1ph, ieee730_1ph, or ads10_1ph")
    end
    gedDict_ud = Dict(:DER_Percent_ud=>DER_Percent_ud, :DER_Rating_factor_ud=>DER_Rating_factor_ud, :Batt_Percent_ud=>Batt_Percent_ud, :Batt_Rating_factor_ud=>Batt_Rating_factor_ud)
end;

# Parse all data
data = parse_all_data(systemName, T, temporal_decmp=temporal_decmp, linearizedModel=linearizedModel, relax_terminal_soc_constraint=relax_terminal_soc_constraint, gedDict_ud=gedDict_ud)

@unpack kVA_B_dict, MVA_B_dict, kV_B_dict, rdict, xdict, rdict_pu, xdict_pu, Z_B_dict, Lset, Nset = data;

if !temporal_decmp
    if !linearizedModel 
        modelDict = optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)
    elseif linearizedModel
        modelDict = optimize_MPOPF_1ph_L(data)
    else
        error("linearizedModel must be either true or false")
    end
    
    @unpack model, modelVals, data = modelDict
    # print_mu(modelDict)

elseif temporal_decmp
    modelDict = optimize_MPOPF_1ph_NL_DDP(data, maxiter=maxiter_ddp)
    @unpack modelVals, data = modelDict
    # print_mu(modelDict)
else
    error("temporal_decmp must be either true or false")
end

if temporal_decmp
    plot_substation_power_cost_allT_vs_k(modelDict, savePlots=savePlots)
end

# postsim computation, plotting, logging
begin

    verbose = false
    # verbose = true

    modelDict = compute_output_values(modelDict, verbose=verbose)

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

