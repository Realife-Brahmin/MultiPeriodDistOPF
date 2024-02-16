## ---------------------- Main AC Power Flow Function ----------------------####
function ac_power_flow_model(slack_bus_type,vol_ctrl_bus_type,load_bus_type,rdata_buses,rdata_loads,rdata_gens,rdata_storage,nw_pPrf_data_load,
    nw_qPrf_data_load,nw_pPrf_data_gen_max,nw_qPrf_data_gen_max,scenario_data_p_min,scenario_data_p_max,scenario_data_q_min,scenario_data_q_max,nw_buses_pf,
    nw_lines_pf,nw_loads_pf,nw_gens_pf,nw_gcost_pf,nw_sbase_pf,v_initial_pf,v_mag_pf,v_angle_pf,max_mismatch_pf,epsilon,iteration,itr_max,ordata_buses_pf,
    p_curt_pf,p_dis_pf,p_ch_pf,p_od_pf,p_ud_pf,load_theta,nd_curt_gen,nd_fl,nd_Str_active,dg_ac_pf,q_dg_opf,flex_adpf,trsf_ratio_opf_lin,flex_oltc,term_status,vm_ac,va_ac)

    vol_nodes_mag_pf   = zeros(nSc_pf,nTP_pf,size(nw_buses,1))                      # Initialization of vol magnitue and vol angle terms. Everytime SC-PF is run, vol magnitudes should be set to zero
    vol_nodes_theta_pf = zeros(nSc_pf,nTP_pf,size(nw_buses,1))


    vol_rect_pf      = zeros(ComplexF64,(nSc_pf,nTP_pf,size(nw_buses,1)))
    idx_slack_bus_pf = findall(x->x==slack_bus_type,rdata_buses[:,2])
    idx_PV_buses_pf  = findall(x->x==vol_ctrl_bus_type,rdata_buses[:,2])
    idx_PQ_buses_pf  = findall(x->x==load_bus_type,rdata_buses[:,2])

    idx_oPV_buses_pf = findall(x->x==vol_ctrl_bus_type,rdata_buses[:,2])

    num_slack_bus_pf = length(idx_slack_bus_pf)
    num_PQ_buses_pf  = length(idx_PQ_buses_pf)
    num_PV_buses_pf  = length(idx_PV_buses_pf)

    old_num_PV_buses_pf = num_PV_buses_pf
    new_num_PV_buses_pf = num_PV_buses_pf

    idx_cmb_pf = vcat(idx_slack_bus_pf,idx_PV_buses_pf)
    ac_pf_con  = []

    for i in 1:length(idx_cmb_pf)
        nd_num = nw_gens_pf[i].gen_bus_num
        sp_vol = nw_gens_pf[i].gen_V_set
        idx = findall(x->x==nd_num,rdata_buses[:,1])
        v_mag_pf[idx[1,1]] = sp_vol
    end

    # if term_status == 1       # OPF model has converged to a solution
    #     v_mag_pf    = vm_ac[1:nBus,:,:]
    #     v_angle_pf  = va_ac[1:nBus,:,:]
    # else
    #     v_mag_pf    = v_mag_pf
    #     v_angle_pf  = v_angle_pf
    # end
### ------------------- Average Scenario Generation ------------------------####
    if scenario == 0 && nTP_pf == 1                                              # Deterministic Single-Period AC-Power Flow
        active_gen   = rdata_gens[:,9]./sbase                                    # P_max column in 'Gens' sheet
        # reactive_gen = rdata_gens[:,4]./sbase                                  # Q_max column in 'Gens' sheet
        reactive_gen = tan(acos(dg_ac_pf)).*active_gen                            # Q_max column in 'Gens' sheet

    elseif scenario == 0 && nTP_pf > 1                                           # Deterministic Multi-Period AC-Power Flow
        active_gen   = nw_pPrf_data_gen_max[:,2:end]./sbase                      # P_max values in 'P_Profiles_Gen_Max'
        # reactive_gen = nw_qPrf_data_gen_max[:,2:end]./sbase                    # Q_max values in 'Q_Profiles_Gen_Max'
        reactive_gen = tan(acos(dg_ac_pf)).*active_gen                           # Q_max values in 'Q_Profiles_Gen_Max'
    elseif scenario == 1 && nTP_pf>1
        if nSc_pf == 1                                                              # Average Sochastic Multi-Period AC-Power Flow
            active_gen   = avg_scenario_gen(prob_scs,scenario_data_p_max,nTP_pf,nGens)    # nGens = total number of generators
            # reactive_gen = avg_scenario_gen(prob_scs,scenario_data_q_max,nTP_pf,nGens)    # nGens = total number of generators
            reactive_gen = tan(acos(dg_ac_pf)).*active_gen
        elseif nSc_pf>1                                                          # Full Stocahstic Multi-Period Problem
            active_gen   = scenario_data_p_max                                   # Full Stochastic Multi-Period AC-Power Flow
            reactive_gen = scenario_data_q_max
        end
    end

### ------------------------------------------------------------------------####
    J11 = zeros(size(nw_buses_pf,1)-num_slack_bus_pf,size(nw_buses_pf,1)-num_slack_bus_pf)
    J12 = zeros(size(nw_buses_pf,1)-num_slack_bus_pf,size(nw_buses_pf,1)-num_slack_bus_pf)
    J21 = zeros(size(nw_buses_pf,1)-num_slack_bus_pf,size(nw_buses_pf,1)-num_slack_bus_pf)
    J22 = zeros(size(nw_buses_pf,1)-num_slack_bus_pf,size(nw_buses_pf,1)-num_slack_bus_pf)

## ---------------- Formation of Admittance matrix -----------------------------
    y_bus = Y_bus(nw_lines_pf,nw_buses_pf,trsf_ratio_opf_lin,1,1,flex_oltc) # For non OLTC option, s & t are replaced by 1
##-------------------------- AC-Power Flow Iteration ------------------------ ##
v_mag_pf_old    = deepcopy(v_mag_pf)
v_angle_pf_old  = deepcopy(v_angle_pf)
for s in 1:nSc_pf
    for t in 1:nTP_pf
        if t>=1
            max_mismatch_pf = 1
            iteration       = 0
            v_mag_pf    = deepcopy(v_mag_pf_old)
            v_angle_pf  = deepcopy(v_angle_pf_old)
        end
        if flex_oltc==1
            y_bus = Y_bus(nw_lines_pf,nw_buses_pf,trsf_ratio_opf_lin,s,t,flex_oltc)
        end

        println("")
        if num_slack_bus_pf>1
            println("ERROR! There are more than 1 swing buses in the data. Please correct the data provided!")
            println("WARNING! The load flow program is terminated without running!")
        else
            p_sch_pf   = zeros(size(nw_buses_pf,1)-num_slack_bus_pf)                              # Net schedule active power except slack bus
            q_sch_pf   = zeros(size(nw_buses_pf,1)-num_slack_bus_pf)                              # Net schedule reactive power except slack bus
            p_node_pf  = zeros(size(nw_buses_pf,1)-num_slack_bus_pf)                              # Injected active power vector (dimension is one less than the number of nodes; excludes slack bus)                                    # Injected active power
            q_node_pf  = zeros(size(nw_buses_pf,1)-num_slack_bus_pf)                              # Injected reactive power vector (dimension is one less than the number of nodes; excludes slack bus)
####------------ Schedule (specified) Power at PQ/PV buses ------------------ ##
#### ------- Determining P_Gen-P_load at PQ buses as well as PV buses ----------
            idx_sch = 0
            (p_sch_pf,q_sch_pf) = schd_power(nw_buses_pf,nw_gens_pf,nw_loads_pf,rdata_buses,rdata_gens,rdata_loads,nw_pPrf_data_load,nw_qPrf_data_load,active_gen,reactive_gen,p_sch_pf,q_sch_pf,idx_sch,nw_sbase_pf,t,s,p_curt_pf,p_dis_pf,p_ch_pf,p_od_pf,p_ud_pf,rdata_storage,load_theta,nd_curt_gen,nd_fl,nd_Str_active,dg_ac_pf,q_dg_opf,flex_adpf)
            idx_sch = 0
##------- Conversion from Cartesian coordiantes to Polar coordinates ---------##
            idx   = 0
            idx_i = 0
            idx_j = 0

            org_idx_PV_buses_pf = zeros(size(idx_oPV_buses_pf))
            mdf_idx_PV_buses_pf = zeros(size(idx_oPV_buses_pf))
            org_idx_PV_buses_pf = idx_oPV_buses_pf
            println("Scenario: S = $s")
            println("Time Period: T = $t")
            while max_mismatch_pf > epsilon && iteration<=itr_max
                # println("Iteration: $iteration")
#### ---------Information about the indices of PV and PQ buses ------------ ###
##-----Indices of PV buses to be deleted from Jacobian and mismatch vector -----
                (idx_dPV_buses,idx_dPQ_buses) = pv_pq_data(idx_PV_buses_pf,idx_slack_bus_pf,idx_PQ_buses_pf) # Modified indices of PV and PQ buses which will be used while deleting specific rows and columns in Jacobian Matrix
##------------- Calculation of vector y-f(x) (mismatch vector) --------------###
                delta_mismatch_inj = mismatch_vector(nw_buses_pf,rdata_buses,v_angle_pf,v_mag_pf,y_bus,p_node_pf,q_node_pf,idx_dPV_buses,p_sch_pf,q_sch_pf,s,t)
##-------------------- Formation of Jacobian Matrix -------------------------###
                J = jacobian(nw_buses_pf,rdata_buses,v_angle_pf,v_mag_pf,y_bus,J11,J12,J21,J22,idx_dPV_buses,s,t)
##------------------ Solving for delta x = inv(J)*delta y -------------------###
                delta_x_pf = J\delta_mismatch_inj
##--------------------- Retrieving voltage magnitude and angles -------------###
                (v_angle_sv,v_mag_sv) = solution_voltage(slack_bus_type,rdata_buses,nw_buses_pf,idx_PQ_buses_pf,idx_PV_buses_pf,idx_dPQ_buses,idx_dPV_buses,v_angle_pf,v_mag_pf,delta_x_pf,idx_slack_bus_pf,s,t)
                v_mag_pf   = v_mag_sv                                                  # sv = solution voltage (Here, voltage is updated)
                v_angle_pf = v_angle_sv
##------------------------- PV-PQ Switching ---------------------------------###
                (nw_buses_pf,rdata_buses,nw_gens_pf) = pv_pq_switching(y_bus,nw_buses_pf,nw_loads_pf,nw_gens_pf,v_mag_pf,v_angle_pf,nw_sbase_pf,rdata_buses,rdata_loads,rdata_gens,idx_oPV_buses_pf,load_bus_type,vol_ctrl_bus_type,mdf_idx_PV_buses_pf,s,t)
                check_pv_buses_pf = length(findall(x->x==vol_ctrl_bus_type,rdata_buses[:,2]))
                if check_pv_buses_pf!=old_num_PV_buses_pf                            # Switching of PV bus to PQ has taken place
                    println("Is there any PQ-PV switching?")
                    (p_sch_pf,q_sch_pf) = schd_power(nw_buses_pf,nw_gens_pf,nw_loads_pf,rdata_buses,rdata_gens,rdata_loads,nw_pPrf_data_load,nw_qPrf_data_load,active_gen,reactive_gen,p_sch_pf,q_sch_pf,idx_sch,nw_sbase_pf,t,s,p_curt_pf,p_dis_pf,p_ch_pf,p_od_pf,p_ud_pf,rdata_storage,load_theta,nd_curt_gen,nd_fl,nd_Str_active,dg_ac_pf,q_dg_opf,flex_adpf)
                end
                old_num_PV_buses_pf = check_pv_buses_pf
####--------------------------------------------------------------------#######
                max_mismatch_pf   = maximum(abs.(delta_mismatch_inj))
                max_mismatch_pf_P = maximum(abs.(delta_mismatch_inj[1:(num_PV_buses_pf+num_PQ_buses_pf)]))
                max_mismatch_pf_Q = maximum(abs.(delta_mismatch_inj[(num_PV_buses_pf+num_PQ_buses_pf)+1:end]))
                iteration = iteration+1
                # println("Maximum mismatch P: $max_mismatch_pf_P")
                # println("Maximum mismatch Q: $max_mismatch_pf_Q")
            end # End while loop
        end     # End If condition
####--------- Priniting the OUTCOME of AC Power Flow Solution --------------####
        oPV_bus_type = ordata_buses_pf
        nPV_bus_type = rdata_buses[:,2]
        switch_PV    = oPV_bus_type-nPV_bus_type
        idx_sPV      = findall(x->x!=0,switch_PV)
        # println("")
        for i in 1:length(idx_sPV)
            nd_num = Int64(rdata_buses[idx_sPV[i,1],1])
            println("WARNING! Generator G$nd_num bus has not switched back as PV!")
        end
## Print message whether the power flow algorithm is solved correctly or not ###
            # println("")
            if max_mismatch_pf<epsilon && iteration < itr_max
                println("The power flow has converged successfully")
                # println("The maximum mismatch is $max_mismatch_pf")
                # println("The converged solution is obtained in $iteration iterations" )
                push!(ac_pf_con,"converged")
            elseif max_mismatch_pf>epsilon && iteration >= itr_max
                println("The power flow has NOT converged in the required number of iteration")
                # println("The maximum mismatch is $max_mismatch_pf")
                # println("The iteration has hit the limit:: $iteration")
                push!(ac_pf_con,"Not converged")
            end
            println("")
            for i in 1:size(nw_buses_pf,1)
                vol_nodes_mag_pf[s,t,i]   = v_mag_pf[i,t,s]
                # vol_nodes_theta_pf[:,t,i] .= rad2deg(v_angle_pf[i])
                vol_nodes_theta_pf[s,t,i] = (v_angle_pf[i,t,s])                         # Angles are in Radian Because input of Cosine and Sine functions are in radian (i.e., Julia consideres them as radian)
                vol_rect_pf[s,t,i] = v_mag_pf[i,t,s]*(cos(v_angle_pf[i,t,s]) + sin(v_angle_pf[i,t,s])im)
            end
    end  # End Time loop
end      # End Scneario Loop
################################################################################
##--------------- Retrieve the value of line loadings ------------------------##
################################################################################


################################################################################
##---------- Retrieve the value of injections and voltages -------------------##
################################################################################
p_node_pf_calc_pf  = zeros(Float64,(nSc_pf,nTP_pf,size(nw_buses_pf,1)))
q_node_pf_calc_pf  = zeros(Float64,(nSc_pf,nTP_pf,size(nw_buses_pf,1)))
p_gen_node_pf      = zeros(Float64,(nSc_pf,nTP_pf,size(nw_gens_pf,1)))
q_gen_node_pf      = zeros(Float64,(nSc_pf,nTP_pf,size(nw_gens_pf,1)))
for s in 1:nSc_pf
    for t in 1:nTP_pf
        v_mag_aux   = vol_nodes_mag_pf[s,t,:]
        v_angle_aux = vol_nodes_theta_pf[s,t,:]
        for i in 1:size(nw_buses_pf,1)
            nd_num        = nw_buses_pf[i].bus_num
            idx_nd_bsheet = findall(x->x==nd_num,rdata_buses[:,1])
            idx_nd_bsheet = idx_nd_bsheet[1,1]
            idx_nd_lsheet = findall(x->x==nd_num,rdata_loads[:,1])
            idx_nd_gsheet = findall(x->x==nd_num,rdata_gens[:,1])
            theta_diff    = vol_nodes_theta_pf[s,t,idx_nd_bsheet,1].-v_angle_aux
            p_inj = vol_nodes_mag_pf[s,t,idx_nd_bsheet,1]*sum(v_mag_aux.*((real(y_bus[idx_nd_bsheet,:])).*cos.(theta_diff)+(imag(y_bus[idx_nd_bsheet,:])).*sin.(theta_diff)),dims=1)
            q_inj = vol_nodes_mag_pf[s,t,idx_nd_bsheet,1]*sum(v_mag_aux.*((real(y_bus[idx_nd_bsheet,:])).*sin.(theta_diff)-(imag(y_bus[idx_nd_bsheet,:])).*cos.(theta_diff)),dims=1)
            p_node_pf_calc_pf[s,t,idx_nd_bsheet,1] = p_inj[1,1]
            q_node_pf_calc_pf[s,t,idx_nd_bsheet,1] = q_inj[1,1]
        end

        for i in 1:size(nw_gens_pf,1)
            nd_num   = nw_gens_pf[i].gen_bus_num
            idx_nd_bsheet = findall(x->x==nd_num,rdata_buses[:,1])
            idx_nd_gsheet = findall(x->x==nd_num,rdata_gens[:,1])
            idx_nd_lsheet = findall(x->x==nd_num,rdata_loads[:,1])
            bus_type = rdata_buses[idx_nd_bsheet[1,1],2]
            if !isempty(idx_nd_lsheet)
                idx_nd_lsheet = idx_nd_lsheet[1,1]
            end
            # if bus_type != slack_bus_type
                if !isempty(idx_nd_lsheet)
                    p_gen_node_pf[s,t,i] = p_node_pf_calc_pf[s,t,idx_nd_bsheet[1,1]]-nw_loads_pf[idx_nd_lsheet[1,1]].load_G/nw_sbase_pf[1].sbase+nw_loads_pf[idx_nd_lsheet[1,1]].load_P/nw_sbase_pf[1].sbase
                    q_gen_node_pf[s,t,i] = q_node_pf_calc_pf[s,t,idx_nd_bsheet[1,1]]-nw_loads_pf[idx_nd_lsheet[1,1]].load_B/nw_sbase_pf[1].sbase+nw_loads_pf[idx_nd_lsheet[1,1]].load_Q/nw_sbase_pf[1].sbase
                else
                    p_gen_node_pf[s,t,i] = p_node_pf_calc_pf[s,t,idx_nd_bsheet[1,1]]
                    q_gen_node_pf[s,t,i] = q_node_pf_calc_pf[s,t,idx_nd_bsheet[1,1]]
                end
            # end
        end
    end
end
p_gen_node_pf = p_gen_node_pf.*nw_sbase_pf[1].sbase
q_gen_node_pf = q_gen_node_pf.*nw_sbase_pf[1].sbase


    return vol_nodes_mag_pf,vol_nodes_theta_pf,vol_rect_pf,idx_slack_bus_pf,p_gen_node_pf,q_gen_node_pf,ac_pf_con
end
