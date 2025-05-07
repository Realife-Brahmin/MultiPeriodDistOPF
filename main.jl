# main.jl
include("./src/setupMultiPeriodDistOPF.jl") 

begin
    # systemName0 = "ads3_1ph"
    # systemName0 = "ads10_1ph"
    systemName0 = "ieee123_1ph-A"
    # systemName0 = "ieee123_1ph-B"
    # systemName0 = "ieee729_1ph"
    # systemName0 = "ieee730_1ph"
    # T0 = 8
    # T0 = 1
    # T0 = 96
    # T0 = 3
    # T0 = 6
    T0 = 12
    # T0 = 24
    # T0 = 48
    factor = 1
    # factor = 1/2
    # linearizedModel = false
    linearizedModel = true
    # temporal_decmp = false
    temporal_decmp = true
    # gamma_fpi = 0.5    
    gamma_fpi = 1.0
    warmStart_mu = "none"
    # warmStart_mu = "nonlinear"
    # warmStart_mu = "linear"
    # savePlots = false
    savePlots = true
end;

begin
    alpha_fpi = 1.00
    # alpha_fpi = 0.001
    T = Int(T0*factor) 
    numAreas = 1
    maxiter_ddp = 50
    threshold_conv_iters = 3
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
    elseif systemName0 == "ads3_1ph"
        systemName = "ads3_1ph"
        DER_Percent_ud = 50
        DER_Rating_factor_ud = 1/3
        Batt_Percent_ud = 50
        Batt_Rating_factor_ud = 1/3
    else
        error("systemName must be either ieee123_1ph, ieee730_1ph, or ads10_1ph")
    end
    gedDict_ud = Dict(:systemName0=>systemName0, :DER_Percent_ud=>DER_Percent_ud, :DER_Rating_factor_ud=>DER_Rating_factor_ud, :Batt_Percent_ud=>Batt_Percent_ud, :Batt_Rating_factor_ud=>Batt_Rating_factor_ud)
end;

begin
    # Parse all data
    data = Parser.parse_all_data(systemName, T, temporal_decmp=temporal_decmp, linearizedModel=linearizedModel, relax_terminal_soc_constraint=relax_terminal_soc_constraint, gedDict_ud=gedDict_ud, alpha_fpi=alpha_fpi, gamma_fpi=gamma_fpi, warmStart_mu=warmStart_mu,
    threshold_conv_iters=threshold_conv_iters)

    @unpack kVA_B_dict, MVA_B_dict, kV_B_dict, rdict, xdict, rdict_pu, xdict_pu, Z_B_dict, Lset, Nset, NLset = data;
end;

#region Optimization
if !temporal_decmp
    if !linearizedModel
        modelDict = Playbook.optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)
    elseif linearizedModel
        modelDict = Playbook.optimize_MPOPF_1ph_L(data)
    else
        error("linearizedModel must be either true or false")
    end
    
    @unpack model, modelVals, data = modelDict
    dualVariablesStateDict = Playbook.get_dual_variables_state_fullMPOPF(modelDict, verbose=true)
elseif temporal_decmp
    if !linearizedModel
        if warmStart_mu == "nonlinear"     
            modelDictBF = Playbook.optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)
            muDict = Playbook.get_soc_dual_variables_fullMPOPF(modelDictBF)
        elseif warmStart_mu == "linear"
            modelDictBF = Playbook.optimize_MPOPF_1ph_L(data)
            muDict = Playbook.get_soc_dual_variables_fullMPOPF(modelDictBF)
        else
            muDict = nothing
        end
        modelDict = DDP.optimize_MPOPF_1ph_NL_DDP(data, maxiter=maxiter_ddp, muDict=muDict)
    elseif linearizedModel
        muDict = nothing
        modelDict = DDPLinear.optimize_MPOPF_1ph_L_DDP(data, maxiter=maxiter_ddp, muDict=muDict)
    else
        error("linearizedModel must be either true or false")
    end

    @unpack modelVals, data = modelDict
    dualVariablesStateDict = Playbook.get_dual_variables_state_fullMPOPF(modelDict, verbose=false)
else
    error("temporal_decmp must be either true or false")
end

if temporal_decmp
    Plotter.plot_substation_power_cost_allT_vs_k(modelDict, savePlots=true)
end
#endregion

# postsim computation, plotting, logging
begin
    # savePlots = true
    verbose = false
    # verbose = true

    modelDict = CO.compute_output_values(modelDict, verbose=verbose)

    Exporter.export_decision_variables(modelDict, verbose=verbose)

    # Todo: Maybe separately save the simulation times? It is annoying to have file content differences every single run (for same exact sim)
    Exporter.export_simulation_key_results_txt(modelDict, verbose=verbose)

    Plotter.plot_battery_actions(modelDict, showPlots=false, savePlots=savePlots, verbose=verbose, maxBatts=2)

    @unpack data = modelDict;
    Plotter.plot_input_forecast_curves(data, filenameSuffix=inputForecastDescription, showPlots=false, verbose=verbose)
    Plotter.plot_substation_power(modelDict, savePlots=savePlots, verbose=verbose)
    Plotter.plot_substation_power_cost(modelDict, savePlots=savePlots, verbose=verbose)
    Plotter.plot_line_losses(modelDict, savePlots=savePlots, verbose=verbose)

    modelDict = validate_opf_against_opendss(modelDict, verbose=verbose)
    @unpack valdVals = modelDict;
    Exporter.export_validation_key_results(modelDict, verbose=verbose)

    Exporter.export_validation_decision_variables(modelDict, verbose=verbose)
end

