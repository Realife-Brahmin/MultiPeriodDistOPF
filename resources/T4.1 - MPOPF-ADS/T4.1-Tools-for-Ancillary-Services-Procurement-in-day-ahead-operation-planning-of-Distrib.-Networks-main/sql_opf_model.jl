function sql_opf_model(vol_nodes_mag_pf,vol_nodes_theta_pf,nSc,nTP,nBus,nNcurt_gen,nCurt_gen,nTrsf,nStr_active,nFl,nd_fl,flex_apc,flex_oltc,flex_adpf,flex_str,flex_fl,str_bin,fl_bin,rdata_buses,nLines,trsf_fict_nodes,tap_ratio_range,
    oltc_bin,prob_scs,time_step,tdata,bus_data_Ssheet,bus_data_lsheet,cost_a_str,cost_b_str,cost_c_str,cost_load_inc,cost_load_dec,nw_buses,rdata_loads,node_data,nd_curt_gen,nd_ncurt_gen,p_load,q_load,idx_Gs_lsheet,idx_Bs_lsheet,pg_max,
    yii_sh,i_curt_gens,sbase,dg_pf,iStr_active,yij_line,dLines,i_ncurt_gens,nw_lines,bus_data_gsheet,pg_min,qg_min,qg_max,pgen_tol,qgen_tol,p_tol,q_tol,rcvr_branch,rcvr_inj,rcvr_tso_dso,Ibase,idx_slack_bus_pf,solver_ipopt,solver_cbc,
    solver_bonmin,sql_itr,sql_itr_max,s_inj_error_rel,s_error_tol,angle_lb,angle_ub,cycleList,idx_slack_from,idx_slack_to,num_slack_cnctd_nodes,slack_cnctd_nodes,slack_nd,stoch_model,tratio_init,rdata_storage,load_theta,nd_Str_active,
    max_err_per_bus,norm_s_err_max_power,p_err_max_power,q_err_max_power,avg_err_per_bus,norm_s_err_avg,p_err_avg,q_err_avg,rel_max_error,rel_avg_error,sql_obj,sql_time,lin_itr,v_initial_pf,min_vol_limit,max_vol_limit,max_crnt_limit,oltc_ratio,
    vm_ac,va_ac,vm_plr,va_plr,vm_min,vm_max,gen_Ncurt_P,gen_Ncurt_Q,gen_Curt_P_sc,gen_Curt_P_gen,total_p_curt,fl_inc_sc,fl_dec_sc,fl_inc_nd,fl_dec_nd,p_ch_nd,p_dis_nd,soc_nd,p_strg_nd,
    trsf_ratio,dgpf,gen_Curt_Q,br_crnt_pu,br_crnt_si,br_pwr_ij_pu,br_pwr_ij_SI,br_crnt_ij_pu_ex,br_crnt_ij_SI_ex,br_pwr_ij_pu_ex,br_pwr_ij_SI_ex,app_err,p_inj,q_inj,p_imp_td,q_imp_td,
    str_bv,fl_bv,oltc_bv,oltc_aux_bv,ratio_sql,p_inj_nl,q_inj_nl,avg_gen_curt,avg_str_chr,avg_str_dch,avg_str_soc,avg_fl_od,avg_fl_ud,avg_tratio,line_alpha,line_eta,line_beta,dg_pf_value,v_mag_int_prv_itr,v_ang_int_prv_itr,
    term_status,solver_cplex,int_tol,opt_tol,num_thread,rdata_trsf,nTrsf_s1,nw_trsf,n_tap)

    v_mag_int   = zeros(Float64,(nSc,nTP,nBus))
    v_ang_int   = zeros(Float64,(nSc,nTP,nBus))

    allowableGap = 0
    ratioGap     = 0
    sql_itr      = 1
    feasb_status = 2      # feasb_status = 1 is OK; feasb_status = 2 means the OPF model is not solved, keep looking
    infsb_model  = 0
    lin_itr_old  = deepcopy(lin_itr)
    term_status  = 0

    while (s_inj_error_rel > s_error_tol) && (sql_itr <= sql_itr_max)
        println(repeat("\n",3))
        println("SQL_iteration: $sql_itr")
        if sql_itr==1
            v_mag_int  = deepcopy(vol_nodes_mag_pf)
            v_ang_int  = deepcopy(vol_nodes_theta_pf)
            oltc_ratio = deepcopy(tratio_init)                                   # For first iteration, oltc_ratio is set to the off-nominal tap ratio mentioned in the 'ratio' column
            if flex_oltc==1
                for i in 1:size(trsf_fict_nodes,1)
                    trsf_from_node = trsf_fict_nodes[i,1]
                    idx_trsf_from = (findall(x->x==trsf_from_node,rdata_buses[:,1]))
                    idx_trsf_from = idx_trsf_from[1,1]
                    # vol_from = vol_nodes_mag_pf[:,:,trsf_fict_nodes[i,1]]
                    vol_from = vol_nodes_mag_pf[:,:,idx_trsf_from]
                    v_mag_int = cat(v_mag_int,vol_from,dims=3)
                end
            end
        end
        if (s_inj_error_rel > s_error_tol) && (sql_itr <= sql_itr_max) # same condition of the WHILE loop above
            feasb_status = 2 # i guess feasb_status=2 means we have not YET found the solution, keep going. this doesn't mean FAILURE
        end

        while (feasb_status != 1) && (infsb_model < 1) # second condition added on June 27 by baraa. feasb_status is updated below as feasb_status=term_status
            println(repeat("\n",3))
            (acopf,v_mag,v_ang,Pg,Qg,Pg_curt,aux_tratio,tratio,Qg_curt,p_ch,p_dis,soc,p_strg,bin_ch,p_fl_inc,p_fl_dec,bin_fl,obj_exp,con_exp,aux_tap,bin_tap,v_ang_node,tratio_aux) = opf_model_initialization(v_mag_int,v_ang_int,nBus,nNcurt_gen,nCurt_gen,nTrsf,nStr_active,nFl,nSc,nTP,flex_oltc,flex_adpf,flex_str,flex_fl,str_bin,fl_bin,rdata_buses,nLines,trsf_fict_nodes,tap_ratio_range,oltc_bin,solver_ipopt,solver_cbc,solver_bonmin,idx_slack_from,idx_slack_to,num_slack_cnctd_nodes,solver_cplex,int_tol,opt_tol,num_thread,nTrsf_s1)
##----------------- Optimal Redispatch of Flexible Resources -----------------##
            cost_expr = opf_model_objective(acopf,flex_oltc,flex_str,flex_fl,flex_apc,nTrsf_s1,nStr_active,nCurt_gen,nFl,nd_fl,nSc,nTP,prob_scs,time_step,tdata,bus_data_Ssheet,bus_data_lsheet,cost_a_str,cost_b_str,cost_c_str,cost_load_inc,cost_load_dec,aux_tratio,p_ch,p_dis,p_fl_inc,p_fl_dec,Pg_curt,sbase,stoch_model)
            @objective(acopf,Min,cost_expr)
## ----------------------- Setting of Constraints -------------------------#####
            opf_model_power_balance_cons(acopf,Pg_curt,Pg,Qg,p_strg,p_fl_inc,p_fl_dec,tratio,v_mag,v_ang,nw_buses,rdata_buses,rdata_loads,bus_data_lsheet,bus_data_Ssheet,node_data,nd_fl,nd_curt_gen,nd_ncurt_gen,p_load,q_load,idx_Gs_lsheet,idx_Bs_lsheet,pg_max,yii_sh,i_curt_gens,sbase,dg_pf,nSc,nTP,nBus,flex_oltc,flex_adpf,v_mag_int,v_ang_int,con_exp,trsf_fict_nodes,v_initial_pf,Qg_curt,min_vol_limit,max_vol_limit,lin_itr_old,p_dis,p_ch,rdata_trsf,nTrsf_s1,nw_trsf)
            network_constraints_no_function_new_milp(acopf,node_data,nw_buses,nBus,nNcurt_gen,nCurt_gen,nFl,nStr_active,bus_data_gsheet,Pg,Qg,Pg_curt,p_strg,p_ch,p_dis,soc,p_fl_inc,p_fl_dec,v_mag,v_ang,v_mag_int,v_ang_int,pg_min,pg_max,qg_min,qg_max,p_load,nw_storage,iStr_active,yij_line,dLines,rdata_buses,tdata,tratio,aux_tratio,aux_tap,bin_tap,bin_fl,bin_ch,flex_apc,flex_oltc,flex_adpf,flex_str,flex_fl,fl_bin,str_bin,i_ncurt_gens,i_curt_gens,trsf_fict_nodes,tap_ratio_range,oltc_bin,angle_lb,angle_ub,cycleList,v_ang_node,slack_cnctd_nodes,slack_nd,tratio_aux,oltc_ratio,Qg_curt,max_crnt_limit,rdata_trsf,nTrsf_s1,nw_trsf)
            # print(acopf)
##--------------------- Solving the Optimization Model---------------------#####
            (term_status,s_obj,s_time) = opf_model_solve(acopf,error_msg,sql_itr)
            global term_status = term_status
            if term_status == 3
                println("The problem is infeasible")
                if lin_itr_old == 1
                    println("Problem is IF in Linear Iteration 1. No bounds tightening is DONE. Check the MODEL again")
                elseif lin_itr_old > 1
                    lin_itr_old = lin_itr_old-1
                    v_mag_int    = deepcopy(v_mag_int_prv_itr[:,:,:,lin_itr_old])
                    v_ang_int    = deepcopy(v_ang_int_prv_itr[:,:,:,lin_itr_old])
                end
                infsb_model  = infsb_model+1
                # println("Term status in sql_opf_modle now after opf_model_solve is $term_status, infsb_model = $infsb_model.") # inserted by Baraa

            elseif term_status == 1
                println("The OPF model is Feasible and is solved SUCEESSFULLY!")
                (vm_ac,
                 va_ac,
                 vm_plr,
                 va_plr,
                 vm_min,
                 vm_max,
                 gen_Ncurt_P,
                 gen_Ncurt_Q,
                 gen_Curt_P_sc,
                 gen_Curt_P_gen,
                 total_p_curt,
                 fl_inc_sc,
                 fl_dec_sc,
                 fl_inc_nd,
                 fl_dec_nd,
                 p_ch_nd,
                 p_dis_nd,
                 soc_nd,
                 p_strg_nd,
                 trsf_ratio,
                 gen_Curt_Q,
                 br_crnt_pu,
                 br_crnt_si,
                 br_pwr_ij_pu,
                 br_pwr_ij_SI,
                 br_crnt_ij_pu_ex,
                 br_crnt_ij_SI_ex,
                 br_pwr_ij_pu_ex,
                 br_pwr_ij_SI_ex,
                 app_err,
                 p_inj,
                 q_inj,
                 p_imp_td,
                 q_imp_td,
                 str_bv,
                 fl_bv,
                 oltc_bv,
                 oltc_aux_bv,
                 ratio_sql,
                 p_inj_nl,
                 q_inj_nl,
                 avg_gen_curt,
                 avg_str_chr,
                 avg_str_dch,
                 avg_str_soc,
                 avg_fl_od,
                 avg_fl_ud,
                 avg_tratio,
                 oltc_ratio_tmp)=
                opf_model_solution(term_status,
                                   v_mag,
                                   v_ang,
                                   Pg,
                                   Qg,
                                   Pg_curt,
                                   p_ch,
                                   p_dis,
                                   soc,
                                   p_strg,
                                   p_fl_inc,
                                   p_fl_dec,
                                   tratio,
                                   bin_ch,
                                   bin_fl,
                                   bin_tap,
                                   aux_tap,
                                   flex_str,
                                   flex_fl,
                                   flex_oltc,
                                   flex_adpf,
                                   nNcurt_gen,
                                   nCurt_gen,
                                   nFl,
                                   nTrsf,
                                   nLines,
                                   nw_lines,
                                   node_data,
                                   rdata_buses,
                                   yij_line,
                                   pg_max,
                                   pgen_tol,
                                   qgen_tol,
                                   p_tol,
                                   q_tol,
                                   prob_scs,
                                   nTP,
                                   nSc,
                                   rcvr_branch,
                                   rcvr_inj,
                                   rcvr_tso_dso,
                                   i_curt_gens,
                                   sbase,
                                   v_mag_int,
                                   v_ang_int,
                                   Ibase,
                                   dLines,
                                   idx_slack_bus_pf,
                                   trsf_fict_nodes,
                                   vm_ac,
                                   va_ac,
                                   vm_plr,
                                   va_plr,
                                   vm_min,
                                   vm_max,
                                   gen_Ncurt_P,
                                   gen_Ncurt_Q,
                                   gen_Curt_P_sc,
                                   gen_Curt_P_gen,
                                   total_p_curt,
                                   fl_inc_sc,
                                   fl_dec_sc,
                                   fl_inc_nd,
                                   fl_dec_nd,
                                   p_ch_nd,
                                   p_dis_nd,
                                   soc_nd,
                                   p_strg_nd,
                                   trsf_ratio,
                                   dgpf,
                                   gen_Curt_Q,
                                   br_crnt_pu,
                                   br_crnt_si,
                                   br_pwr_ij_pu,
                                   br_pwr_ij_SI,
                                   br_crnt_ij_pu_ex,
                                   br_crnt_ij_SI_ex,
                                   br_pwr_ij_pu_ex,
                                   br_pwr_ij_SI_ex,
                                   app_err,
                                   p_inj,
                                   q_inj,
                                   p_imp_td,
                                   q_imp_td,
                                   str_bv,
                                   fl_bv,
                                   oltc_bv,
                                   oltc_aux_bv,
                                   ratio_sql,
                                   p_inj_nl,
                                   q_inj_nl,
                                   avg_gen_curt,
                                   avg_str_chr,
                                   avg_str_dch,
                                   avg_str_soc,
                                   avg_fl_od,
                                   avg_fl_ud,
                                   avg_tratio,
                                   oltc_ratio,
                                   dg_pf_value,
                                   dg_pf,
                                   Qg_curt,
                                   rdata_trsf,
                                   nTrsf_s1,
                                   nw_trsf,
                                   n_tap)

                v_mag_int = deepcopy(vm_plr)
                v_ang_int = deepcopy(va_plr)

                if flex_oltc == 1
                    # oltc_ratio = deepcopy(ratio_sql)
                    oltc_ratio = deepcopy(oltc_ratio_tmp)                        # OLTC_ratio is updated with the New tap ratio values (For OLTC with status = 1, tap ratio values are determined
                end                                                              # by SLA, whereas for OLTC with status = 0, tap ratio values are equal to the off-nominal ratio values)

                (max_err_per_bus,norm_s_err_max_power,p_err_max_power,q_err_max_power) = s_inj_power_max_power(p_inj,q_inj,p_inj_nl,q_inj_nl,max_err_per_bus,norm_s_err_max_power,p_err_max_power,q_err_max_power,sql_itr,lin_itr)
                (avg_err_per_bus,norm_s_err_avg,p_err_avg,q_err_avg) = s_inj_power_per_bus(p_inj,q_inj,p_inj_nl,q_inj_nl,avg_err_per_bus,norm_s_err_avg,p_err_avg,q_err_avg,sql_itr,lin_itr)

                rel_max_error[sql_itr,lin_itr] = maximum(norm_s_err_max_power[:,sql_itr,lin_itr])*sbase # % pu changed to % MW
                rel_avg_error[sql_itr,lin_itr] = maximum(norm_s_err_avg[:,sql_itr,lin_itr])*sbase       # % pu changed to % MW

                s_inj_error_rel = maximum(norm_s_err_max_power[:,sql_itr,lin_itr])./100      # Converting from % to normal value (value is in pu)

                sql_obj[sql_itr,lin_itr]  = s_obj
                sql_time[sql_itr,lin_itr] = s_time

                feasb_status = term_status
                empty!(acopf)
                GC.gc()
            end
            # println("We reached the END of the inner WHILE (in SQL_OPF_MODEL), line #250, and feasb_status=$feasb_status, term_status=$term_status")
        end
    sql_itr = sql_itr+1
    global term_status = term_status
    # println("We exited the inner WHILE (in SQL_OPF_MODEL), and we're in the outer WHILE now, line #253, and feasb_status=$feasb_status, term_status=$term_status")
    end
global term_status = term_status
# println("We exited all WHILE loops (in SQL_OPF_MODEL), line #255, and feasb_status=$feasb_status, term_status=$term_status")
return  vm_ac,
        va_ac,
        vm_plr,
        va_plr,
        vm_min,
        vm_max,
        gen_Ncurt_P,
        gen_Ncurt_Q,
        gen_Curt_P_sc,
        gen_Curt_P_gen,
        total_p_curt,
        fl_inc_sc,
        fl_dec_sc,
        fl_inc_nd,
        fl_dec_nd,
        p_ch_nd,
        p_dis_nd,
        soc_nd,
        p_strg_nd,
        trsf_ratio,
        gen_Curt_Q,
        br_crnt_pu,
        br_crnt_si,
        br_pwr_ij_pu,
        br_pwr_ij_SI,
        br_crnt_ij_pu_ex,
        br_crnt_ij_SI_ex,
        br_pwr_ij_pu_ex,
        br_pwr_ij_SI_ex,
        app_err,
        p_inj,
        q_inj,
        p_inj_nl,
        q_inj_nl,
        p_imp_td,
        q_imp_td,
        str_bv,
        fl_bv,
        oltc_bv,
        oltc_aux_bv,
        ratio_sql,
        v_mag_int,
        v_ang_int,
        max_err_per_bus,
        norm_s_err_max_power,
        p_err_max_power,
        q_err_max_power,
        avg_err_per_bus,
        norm_s_err_avg,
        p_err_avg,
        q_err_avg,
        rel_max_error,
        rel_avg_error,
        sql_itr,
        sql_obj,
        sql_time,
        avg_gen_curt,
        avg_str_chr,
        avg_str_dch,
        avg_str_soc,
        avg_fl_od,
        avg_fl_ud,
        avg_tratio,
        oltc_ratio,
        infsb_model,
        term_status
end
