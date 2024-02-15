function ac_power_flow_opf_int(nBus,
                               nLines,
                               nNcurt_gen,
                               rdata_buses,
                               nSc,
                               nTP,
                               prob_scs,
                               time_step,
                               stoch_model,
                               Pg_curt_value,
                               p_fl_inc_value,
                               p_fl_dec_value,
                               p_strg_value,
                               Qg_curt_value,
                               nw_buses,
                               nw_lines,
                               rdata_loads,
                               bus_data_lsheet,
                               bus_data_Ssheet,
                               node_data,
                               nd_fl,
                               nd_curt_gen,
                               nd_ncurt_gen,
                               p_load,
                               q_load,
                               idx_Gs_lsheet,
                               idx_Bs_lsheet,
                               yii_sh,
                               i_curt_gens,
                               sbase,
                               flex_oltc,
                               vbase,
                               Ibase,
                               pg_min,
                               pg_max,
                               qg_min,
                               qg_max,
                               yij_line,
                               dLines,
                               i_ncurt_gens,
                               error_msg,
                               pgen_tol,
                               qgen_tol,
                               oltc_ratio_pf,
                               vol_cstr_tol)

    (acpf,e,f,Pg,Qg) = opf_model_initialization_pf_int(nBus,nNcurt_gen,rdata_buses)

    slack_flow = opf_model_objective_pf_int(acpf,nSc,nTP,nNcurt_gen,prob_scs,time_step,stoch_model,Pg,Qg)
    @objective(acpf,Min,slack_flow)

    # ACPF is the name of julia file or julia file, Min is the sense, slack_flow is the quantity to minimize.
    # basically, he's not minimizing a cost value in $, but minimizing power-import from the slack (grid)

    opf_model_power_balance_cons_pf_int(acpf,
                                        Pg_curt_value,
                                        p_fl_inc_value,
                                        p_fl_dec_value,
                                        p_strg_value,
                                        Qg_curt_value,
                                        Pg,
                                        Qg,
                                        e,
                                        f,
                                        nw_buses,
                                        rdata_buses,
                                        rdata_loads,
                                        bus_data_lsheet,
                                        bus_data_Ssheet,
                                        node_data,
                                        nd_fl,
                                        nd_curt_gen,
                                        nd_ncurt_gen,
                                        p_load,
                                        q_load,
                                        idx_Gs_lsheet,
                                        idx_Bs_lsheet,
                                        pg_max,
                                        yii_sh,
                                        i_curt_gens,
                                        sbase,
                                        nSc,
                                        nTP,
                                        nBus,
                                        flex_oltc,
                                        oltc_ratio_pf,
                                        vol_cstr_tol)

    network_constraints_no_function_new_pf_int(acpf,
                                               node_data,
                                               nw_buses,
                                               nw_lines,
                                               nBus,
                                               nLines,
                                               nNcurt_gen,
                                               rdata_buses,
                                               vbase,
                                               Ibase,
                                               Pg,
                                               Qg,
                                               e,
                                               f,
                                               pg_min,
                                               pg_max,
                                               qg_min,
                                               qg_max,
                                               yij_line,
                                               dLines,
                                               i_ncurt_gens,
                                               nSc,
                                               nTP)

    term_status_pf = opf_model_solve_pf_int(acpf,error_msg)

    (vm_ac_pf,va_ac_pf,p_slack_pf,q_slack_pf,v_rect_check) = opf_model_solution_pf_int(term_status_pf,e,f,Pg,Qg,nNcurt_gen,pgen_tol,qgen_tol,nTP,nSc,nBus,sbase)

    vol_nodes_mag_pf   = zeros(nSc_pf,nTP_pf,size(nw_buses,1))                      # Initialization of vol magnitue and vol angle terms. Everytime SC-PF is run, vol magnitudes should be set to zero
    vol_nodes_theta_pf = zeros(nSc_pf,nTP_pf,size(nw_buses,1))
    vol_rect_pf        = zeros(ComplexF64,(nSc_pf,nTP_pf,size(nw_buses,1)))

    for s in 1:nSc
        for t in 1:nTP
            for i in 1:size(nw_buses_pf,1)
                vol_nodes_mag_pf[s,t,i]   = vm_ac_pf[i,t,s]
                # vol_nodes_theta_pf[:,t,i] .= rad2deg(v_angle_pf[i])
                vol_nodes_theta_pf[s,t,i] = (va_ac_pf[i,t,s])                         # Angles are in Radian Because input of Cosine and Sine functions are in radian (i.e., Julia consideres them as radian)
                vol_rect_pf[s,t,i] = vm_ac_pf[i,t,s]*(cos(va_ac_pf[i,t,s]) + sin(va_ac_pf[i,t,s])im)
            end
        end
    end

ObjValue = 0 #value.(slack_flow)
return vol_nodes_mag_pf,vol_nodes_theta_pf,vol_rect_pf,ObjValue
end
