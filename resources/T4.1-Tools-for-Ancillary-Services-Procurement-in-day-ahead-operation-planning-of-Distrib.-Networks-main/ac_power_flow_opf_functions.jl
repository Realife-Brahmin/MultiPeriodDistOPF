function  opf_model_initialization_pf(nBus,nNcurt_gen,rdata_buses)
#-------------------------- Model Initialization -------------------------------
#---------------------------Decalration of OPF Variables------------------------
acpf = Model()
set_optimizer(acpf,Ipopt.Optimizer)                                             # Here the solver is being set
set_optimizer_attributes(acpf,"tol"=>1e-5,"max_iter"=>6000)  # "max_iter"=>100

#-------------------- AmplNLWriter (Bonmin) Solver -----------------------------
# optimizer  = AmplNLWriter.Optimizer
# acopf = Model(() -> AmplNLWriter.Optimizer("C:/bonmin-win64/bonmin.exe",["print_level=1","bonmin.file_solution=yes","file_print_level=5"]))
# acopf = Model(() -> AmplNLWriter.Optimizer("C:/bonmin-win64/bonmin.exe"))
################################################################################
# @timeit md_init "md_init_1"
@variables(acpf, begin
    (e[i=1:nSc, j=1:nTP,k=1:nBus],start = 1.00)                                  # This procedure can be used to set the start values of all variables
    (f[i=1:nSc, j=1:nTP,k=1:nBus],start = 0.00)
#---------- nNcurt_gen and nCurt_gen now forms the complete nGens --------------
    Pg[i=1:nSc, j=1:nTP,k=1:nNcurt_gen]                                          # Pg = non-curtailable Gen
    Qg[i=1:nSc, j=1:nTP,k=1:nNcurt_gen]
    end
    )
## ------------------ Setting values of voltage ----------------------------- ##
####--------------------- Flat Start Voltage -------------------------------####
idx_slack = findall(x->x==3,rdata_buses[:,2])
idx_slack = idx_slack[1,1]
fix.(f[:,:,idx_slack],0.0)
##---------------------------------------------------------------------------###
return acpf,e,f,Pg,Qg
end

function opf_model_objective_pf(acpf,nSc,nTP,nNcurt_gen,prob_scs,time_step,stoch_model,Pg,Qg)

######################### Slack Bus Power Flow #################################
##------------- Minimizing the active and reactive power flow ------------------
################################################################################
app_power_slack = []
for i in 1:nNcurt_gen
    for s in 1:nSc
        for t in 1:nTP
            if stoch_model == 0
                flow_slack = time_step*sbase*(Pg[s,t,i]+Qg[s,t,i])
            elseif stoch_model == 1
                flow_slack = prob_scs[s,i]*time_step*sbase*(Pg[s,t,i]+Qg[s,t,i])

            end
        push!(app_power_slack,flow_slack)
        end
    end
end
a = sum(app_power_slack)
slack_flow = a
return slack_flow
end


################################################################################
##---------------------Script for Power Balance Constraint-----------------------
################################################################################
function opf_model_power_balance_cons_pf(acpf,Pg_curt_value,p_fl_inc_value,p_fl_dec_value,p_strg_value,Qg_curt_value,Pg,Qg,e,f,nw_buses,rdata_buses,rdata_loads,bus_data_lsheet,bus_data_Ssheet,node_data,nd_fl,nd_curt_gen,nd_ncurt_gen,p_load,q_load,idx_Gs_lsheet,idx_Bs_lsheet,pg_max,yii_sh,i_curt_gens,sbase,nSc,nTP,nBus,flex_oltc,oltc_ratio_pf,vol_cstr_tol)
    for i in 1:nBus
        nd      = node_data[i,1].node_num                                            # Node num is independant of time period and scenario
        nd_num  = unique(nd)
        nd_num  = nd_num[1,1]
        idx_nd_nw_buses = findall(x->x==nd_num,rdata_buses[:,1])
        idx_nd_nw_buses = idx_nd_nw_buses[1,1]
        idx_bus_lsheet  = findall(x->x==nd_num,bus_data_lsheet)                      # Index of buses in load sheet
        idx_bus_Stsheet = findall(x->x==nd_num,bus_data_Ssheet)                      # Index of buses in Storage sheet
        idx_fl          = findall(x->x==nd_num,nd_fl)
        idx_curt_gen    = findall(x->x==nd_num,nd_curt_gen)
        idx_ncurt_gen   = findall(x->x==nd_num,nd_ncurt_gen)
## --------------------------------------------------------------------------###
# Accessing any shunt element connected at load/other buses (shunt capacitance etc.)
        Pii_sh  = rdata_loads[idx_bus_lsheet,idx_Gs_lsheet]/sbase                    # Active power injected or absorbed by the shunt element
        Qii_sh  = rdata_loads[idx_bus_lsheet,idx_Bs_lsheet]/sbase                    # Reactive power injected (capacitve) or absorbed (inductive) by the shunt elements
        gii_sh  = check_shunt_acpf(Pii_sh,idx_nd_nw_buses,nw_buses)
        bii_sh  = check_shunt_acpf(Qii_sh,idx_nd_nw_buses,nw_buses)
        push!(yii_sh,gii_sh+(bii_sh)im)

        for s in 1:nSc
            for t in 1:nTP
##------------- If gens are both curtailable and  non-curtailable --------------
                (pg,qg)   = nodal_components_values_acpf(Pg_curt_value,Pg,Qg,s,t,idx_ncurt_gen,idx_curt_gen,Qg_curt_value)
##-------------------- If gens are only non-curtailable ------------------------
            # pg = nodal_components_values(Pg,s,t,idx_bus_gsheet)
            # qg = nodal_components_values(Qg,s,t,idx_bus_gsheet)
##--------- Saving flexible load data, fixed load and storage data -----------##
                p_flx_inc = nodal_components_values_acpf(p_fl_inc_value,s,t,idx_fl)
                p_flx_dec = nodal_components_values_acpf(p_fl_dec_value,s,t,idx_fl)
                pd        = nodal_components_values_acpf(p_load,s,t,idx_bus_lsheet)
                qd        = nodal_components_values_acpf(q_load,s,t,idx_bus_lsheet)
                pstr      = nodal_components_values_acpf(p_strg_value,s,t,idx_bus_Stsheet)
                v_sq_out  = @NLexpression(acpf,(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2))
##----------------Nodal active and reactive power injections--------------------
                pinj_expr = []
                qinj_expr = []
                for j in 1:size(nd,1)                                                # length of each node vector in 'node_data' variable
                    gij_line     = node_data[i,1].node_gij_sr[j,1]
                    bij_line     = node_data[i,1].node_bij_sr[j,1]
                    gij_line_sh  = node_data[i,1].node_gij_sh[j,1]
                    bij_line_sh  = node_data[i,1].node_bij_sh[j,1]
                    cnctd_nd     = node_data[i,1].node_cnode[j,1]
                    idx_cnctd_nd = findall(x->x==cnctd_nd,rdata_buses[:,1])
                    idx_cnctd_nd = idx_cnctd_nd[1,1]
                    v_r  =  @NLexpression(acpf,e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
                    v_im =  @NLexpression(acpf,f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
                    v_sq =  @NLexpression(acpf,(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2))
                    if  (node_data[i,1].node_tap_ratio[j,1]!=0.0)  &&  (node_data[i,1].node_tap_ratio_min[j,1]!=0.0 && node_data[i,1].node_tap_ratio_max[j,1]!=0.0)                          # Transformer branch
                        if  nd_num==node_data[i,1].node_from[j,1]
                            if flex_oltc == 1
                                tap = oltc_ratio_pf[s,t,node_data[i,1].node_idx_trsf[j,1]]
                                # tap =  tratio[s,t,node_data[i,1].node_idx_trsf[j,1]]       # Node is connected to the sending end of the transformer (variable reference for the transformer)
                            else
                                tap =  node_data[i,1].node_tap_ratio[j,1]                # The hard value is used to validate the code against Carleton model!
                            end
                            # (pij,qij) = power_flow_line_trsf_minlp(acopf,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
                            (pij,qij) = power_flow_line_trsf_acpf(acpf,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,v_r,v_im,v_sq)
                        else
                            tap1 = 1.0
                            if flex_oltc == 1
                                tap2 = oltc_ratio_pf[s,t,node_data[i,1].node_idx_trsf[j,1]]
                                # tap2 =  tratio[s,t,node_data[i,1].node_idx_trsf[j,1]]   # Node is connected to the receiving end of the transformer (variable reference for the transformer)
                            else
                                tap2 =  node_data[i,1].node_tap_ratio[j,1]                # The hard value is used to validate the code against Carleton model!
                            end
                            # (pij,qij) = power_flow_line_trsf_minlp(acopf,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
                            (pij,qij) = power_flow_line_trsf_acpf(acpf,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,v_r,v_im,v_sq)
                        end
                    elseif (node_data[i,1].node_tap_ratio[j,1]!=0.0)  &&  (node_data[i,1].node_tap_ratio_min[j,1]==0.0 && node_data[i,1].node_tap_ratio_max[j,1]==0.0)
                            println("WARNING! The transformer tap ratio is non-zero but both minimum and maximum tap ratios are set to zero")
                            println("This transformer branch is treated as a simple line while setting the power injection constraint!")
                            tap = 1.0
                            (pij,qij) = power_flow_line_acpf(acpf,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
                    else                                                                      # Transmission/Distribution line
                            tap = 1.0
                            # (pij,qij) = power_flow_line_minlp(acopf,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
                            (pij,qij) = power_flow_line_acpf(acpf,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
                    end
                    push!(pinj_expr,pij)
                    push!(qinj_expr,qij)
                end
                p_inj = @NLexpression(acpf,sum(pinj_expr[j,1] for j in 1:size(pinj_expr,1)))
                q_inj = @NLexpression(acpf,sum(qinj_expr[j,1] for j in 1:size(qinj_expr,1)))

## ------ Power balance constraints when DG curtailement, Storage and ----------
##------------- load flexibility variables are considered-----------------------

                if isempty(idx_ncurt_gen) && ~isempty(idx_curt_gen)
################################################################################
##----------Constant power factor of DGs is considered!-------------------------
################################################################################
                    if flex_adpf == 0                                                # The power factor can be 1 (No Reactive power) or other than 1 (Reactive power is injected)
                            p_avl    = pg_max[s,t,i_curt_gens[idx_curt_gen][1]]
                            @NLconstraint(acpf,(sum(p_avl[j,1] for j in 1:size(p_avl,1))-sum(pg[j,1] for j in 1:size(pg,1)))+sum(pstr[j,1] for j in 1:size(pstr,1))-(sum(pd[j,1] for j in 1:size(pd,1))+sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(gii_sh[j,1] for j in 1:size(gii_sh,1))*((e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2)))==p_inj)
                            @NLconstraint(acpf,tan(acos(dg_pf))*(sum(p_avl[j,1] for j in 1:size(p_avl,1))-sum(pg[j,1] for j in 1:size(pg,1)))-(sum(qd[j,1] for j in 1:size(qd,1))+load_theta*sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-load_theta*sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(bii_sh[j,1] for j in 1:size(bii_sh,1))*((e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2)))==q_inj) # Q_inj of DG DOES NOT depend upon P_curt
###############################################################################
##-------------- Adaptive power factor of DGs is considered --------------------
################################################################################
                    elseif flex_adpf == 1
                        p_avl = pg_max[s,t,i_curt_gens[idx_curt_gen][1]]
                        idx   = idx_curt_gen[1]
                        @NLconstraint(acpf,sum(p_avl[j,1] for j in 1:size(p_avl,1))-sum(pg[j,1] for j in 1:size(pg,1))+sum(pstr[j,1] for j in 1:size(pstr,1))-(sum(pd[j,1] for j in 1:size(pd,1))+sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(gii_sh[j,1] for j in 1:size(gii_sh,1))*((e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2)))==p_inj)
                        @NLconstraint(acpf,sum(qg[j,1] for j in 1:size(qg,1))-(sum(qd[j,1] for j in 1:size(qd,1))+load_theta*sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-load_theta*sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(bii_sh[j,1] for j in 1:size(bii_sh,1))*((e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2)))==q_inj)

                    end
##------------------------------------------------------------------------------
                else
                    @NLconstraint(acpf,sum(pg[j,1] for j in 1:size(pg,1))+sum(pstr[j,1] for j in 1:size(pstr,1))-(sum(pd[j,1] for j in 1:size(pd,1))+sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(gii_sh[j,1] for j in 1:size(gii_sh,1))*((e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2)))==p_inj)
                    @NLconstraint(acpf,sum(qg[j,1] for j in 1:size(qg,1))-(sum(qd[j,1] for j in 1:size(qd,1))+load_theta*sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-load_theta*sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(bii_sh[j,1] for j in 1:size(bii_sh,1))*((e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2)))==q_inj)
                end
################################################################################
############## Setting voltage constraints on network buses ####################
################################################################################
                if nw_buses[i].bus_type==3                                       # Fixing the voltage at that bus whose type is 3 (Assuming this to be the Point of Common Coupling!)
                    # @constraint(acpf,[i,s,t],((nw_buses[idx_nd_nw_buses].bus_vmin)^2)-vol_cstr_tol<=(e[s,t,idx_nd_nw_buses]^2)<=((nw_buses[idx_nd_nw_buses].bus_vmax)^2)+vol_cstr_tol)
                    # @constraint(acopf,[i,s,t],(f[s,t,idx_nd_nw_buses]==0.0))
                    # @constraint(acpf,[i,s,t],(e[s,t,idx_nd_nw_buses]^2)==1.00^2)
                    if nTP == 24
                        @constraint(acpf,[i,s,t],(e[s,t,idx_nd_nw_buses]^2)==v_initial_pf[1,t]^2)
                    else
                        @constraint(acpf,[i,s,t],(e[s,t,idx_nd_nw_buses]^2)==v_initial_pf^2)
                    end
                else                                                             # Constraint on the remaining buses
                    @constraint(acpf,((nw_buses[idx_nd_nw_buses].bus_vmin)^2)-vol_cstr_tol<=(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2)<=((nw_buses[idx_nd_nw_buses].bus_vmax)^2)+vol_cstr_tol)
                end
################################################################################
################################################################################
            end
        end
    end
return acpf
end

################################################################################
##------------------------- Script for Network Constraint-----------------------
################################################################################
function network_constraints_no_function_new_pf(acpf,node_data,nw_buses,nw_lines,nBus,nLines,nNcurt_gen,rdata_buses,vbase,Ibase,Pg,Qg,e,f,pg_min,pg_max,qg_min,qg_max,yij_line,dLines,i_ncurt_gens,nSc,nTP)
################################################################################
####################### Non-Curtailable Gen Constraints ########################
################################################################################
for i in 1:nNcurt_gen
    for s in 1:nSc
        for t in 1:nTP
        # @constraint(acpf,pg_min[s,t,i_ncurt_gens[i][1]]<=Pg[s,t,i]<=pg_max[s,t,i_ncurt_gens[i][1]])
        # @constraint(acpf,qg_min[s,t,i_ncurt_gens[i][1]]<=Qg[s,t,i]<=qg_max[s,t,i_ncurt_gens[i][1]])

        @constraint(acpf,pg_min[s,t,i_ncurt_gens[i][1]]-p_cstr_fctr<=Pg[s,t,i]<=pg_max[s,t,i_ncurt_gens[i][1]]+p_cstr_fctr)
        @constraint(acpf,qg_min[s,t,i_ncurt_gens[i][1]]-q_cstr_fctr<=Qg[s,t,i]<=qg_max[s,t,i_ncurt_gens[i][1]]+q_cstr_fctr)
        end
    end
end
################################################################################
####################### Logitudinal Current Constraint #########################
################################################################################
    # for i in 1:nLines
    #     gij_line  = real(yij_line[i,1])
    #     bij_line  = imag(yij_line[i,1])
    #     tap       = nw_lines[i].line_tap_ratio
    #     idx_tap   = Int64(dLines[i,end])
    #     f_bus     = nw_lines[i].line_from
    #     t_bus     = nw_lines[i].line_to
    #     idx_t_bus = findall(x->x==t_bus,rdata_buses[:,1])
    #     idx_f_bus = findall(x->x==f_bus,rdata_buses[:,1])
    #     idx_f_bus = idx_f_bus[1,1]
    #     idx_t_bus = idx_t_bus[1,1]
    #
    #     Imax_from = (nw_lines[i].line_Smax_A*1000)/vbase[idx_f_bus,1]            # Imax_from is in Amperes
    #     Imax_to   = (nw_lines[i].line_Smax_A*1000)/vbase[idx_t_bus,1]            # Imax_to is in Amperes
    #     Imax_lng  = (nw_lines[i].line_Smax_A*1000)/vbase[idx_f_bus,1]            # Imax_lng is in Amperes
    #
    #     Imax_fpu  = Imax_from/Ibase[idx_f_bus,1]
    #     Imax_tpu  = Imax_to/Ibase[idx_t_bus,1]
    #     Imax_lpu  = Imax_lng/Ibase[idx_f_bus,1]
    #
    #     for s in 1:nSc
    #         for t in 1:nTP
    #             if tap!=0.0
    # # This is to compare the longitiudinal current constraint with the Imax constraint
    #                 if flex_oltc == 1
    #                     # t_tap = tratio[s,t,idx_tap]                                      # Variable ref for the transformer
    #                     t_tap = oltc_ratio[s,t,idx_tap]
    #                 elseif flex_oltc == 0
    #                     t_tap = tap                                                    # Hard value is set to validate the model against PowerModels results
    #                 end
    #             @NLconstraint(acpf,((gij_line^2+bij_line^2)*(t_tap^2*(e[s,t,idx_f_bus]^2+f[s,t,idx_f_bus]^2) + (e[s,t,idx_t_bus]^2+f[s,t,idx_t_bus]^2) + t_tap*(-2(e[s,t,idx_f_bus]*e[s,t,idx_t_bus]+f[s,t,idx_f_bus]*f[s,t,idx_t_bus]))))<=Imax_tpu^2)
    #
    #             else
    #             @NLconstraint(acpf, ((gij_line^2+bij_line^2)*((e[s,t,idx_f_bus]^2+f[s,t,idx_f_bus]^2) + (e[s,t,idx_t_bus]^2+f[s,t,idx_t_bus]^2) + (-2(e[s,t,idx_f_bus]*e[s,t,idx_t_bus]+f[s,t,idx_f_bus]*f[s,t,idx_t_bus]))))<=Imax_lpu^2)
    #             end
    #         end
    #     end
    # end
return acpf
end

################################################################################
############################ Solve OPF Model ###################################
################################################################################
function opf_model_solve_pf(acpf,error_msg)
    if isempty(error_msg)    # There is no error message
        status = optimize!(acpf)
    else
        println("The model is not formulated Successfully!")
    end

    term_status = JuMP.termination_status(acpf)                                 # Status of the resulting problem
    if isequal(termination_status(acpf),MOI.OPTIMAL) || isequal(termination_status(acpf),MOI.LOCALLY_SOLVED)
        println("Objective value: ", JuMP.objective_value(acpf))
        println("Solver Time: ", JuMP.solve_time(acpf))

        JuMP.primal_status(acpf)
        JuMP.dual_status(acpf)
        file  = string("results_",filename)
        term_status = 1
    elseif termination_status(acpf) == MOI.TIME_LIMIT && has_values(acpf)
        suboptimal_objective = objective_value(acpf)
        term_status = 2
    else
        println("The model has not been solved Successfully!")
        term_status = 3
    end
    return term_status
end

################################################################################
##################### Recovering the OPF Solution ##############################
################################################################################
function opf_model_solution_pf(term_status,e,f,Pg,Qg,nNcurt_gen,pgen_tol,qgen_tol,nTP,nSc,nBus,sbase)
    # println("Entering recovery_voltages_pf, term_status     is [$term_status]") # Inserted by baraa
    if term_status == 1 || term_status == 2
        vr            = JuMP.value.(e[:,:,:])
        vim           = JuMP.value.(f[:,:,:])
        Pgen_nCurt_pu = JuMP.value.(Pg[:,:,:])
        Qgen_nCurt_pu = JuMP.value.(Qg[:,:,:])

################################################################################
######################### Recovering Voltages ##################################
################################################################################
        (vm_ac,va_ac,v_rect,vm_min,vm_max) = recovery_voltages_pf(nTP,nSc,nBus,vr,vim)
        ################################################################################
        ########### Recovering Generator Injection Values (Non-Curtailable) ############
        ################################################################################
        (gen_Ncurt_P,gen_Ncurt_Q) = recovery_ncurt_gen_pf(nTP,nSc,nNcurt_gen,Pgen_nCurt_pu,Qgen_nCurt_pu,pgen_tol,qgen_tol)
        gen_Ncurt_P = gen_Ncurt_P*sbase
        gen_Ncurt_Q = gen_Ncurt_Q*sbase
    else
        # println("Not sure what to do here. need to fix this")
    end
return vm_ac,va_ac,gen_Ncurt_P,gen_Ncurt_Q,v_rect
end
