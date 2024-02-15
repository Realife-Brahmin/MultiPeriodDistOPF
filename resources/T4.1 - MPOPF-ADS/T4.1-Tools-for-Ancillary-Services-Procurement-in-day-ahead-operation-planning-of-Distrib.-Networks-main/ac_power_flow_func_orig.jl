##-----------------  Schedule Power Calculation Function --------------------###
function schd_power(nw_buses,nw_gens,nw_loads,rdata_buses,rdata_gens,rdata_loads,nw_pPrf_data_load,nw_qPrf_data_load,active_gen,reactive_gen,p_sch,q_sch,idx_sch,
    nw_sbase,t,s,p_curt_pf,p_dis_pf,p_ch_pf,p_od_pf,p_ud_pf,rdata_storage,load_theta,nd_curt_gen,nd_fl,nd_Str_active,dg_ac_pf,q_dg_opf,flex_adpf)
    idx_sch = 0
    for i in 1:size(nw_buses,1)
        nd_num   = nw_buses[i].bus_num
        bus_type = nw_buses[i].bus_type

        if bus_type == 3
            idx_sch = idx_sch+1
        elseif bus_type !=3
##----------- Indices of Nodes in Buses, Loads and Generator Sheets ---------###
            if nTP_pf == 1
                idx_nd_bsheet = findall(x->x==nd_num,rdata_buses[:,1])               # Index of node in bus sheet
                idx_nd_lsheet = findall(x->x==nd_num,rdata_loads[:,1])               # Index of node in load sheet
                idx_nd_gsheet = findall(x->x==nd_num,rdata_gens[:,1])                # Index of node in generator sheet
            elseif nTP_pf>1
                idx_nd_bsheet = findall(x->x==nd_num,rdata_buses[:,1])               # Index of node in bus sheet
                idx_nd_lsheet = findall(x->x==nd_num,rdata_pProfile_load[:,1])       # Index of node in load sheet
                idx_nd_gsheet = findall(x->x==nd_num,nw_pPrf_data_gen_max[:,1])      # Index of node in generator sheet
            end

            idx_nd_Fl     = findall(x->x==nd_num,nd_fl)                             # Index of node in the flexible Load data (nd_fl)
            idx_nd_Ssheet = findall(x->x==nd_num,nd_Str_active)                     # Index of node in Storage sheet
            idx_nd_nCurt  = findall(x->x==nd_num,nd_curt_gen)                       # Index of node in Curtailable generator data!

            if !isempty(idx_nd_bsheet)
                idx_nd_bsheet = idx_nd_bsheet[1,1]-idx_sch
            end
            if !isempty(idx_nd_lsheet)
                idx_nd_lsheet = idx_nd_lsheet[1,1]
            end
            if !isempty(idx_nd_gsheet)
                idx_nd_gsheet = idx_nd_gsheet[1,1]
            end
            if !isempty(idx_nd_Ssheet)
                idx_nd_Ssheet = idx_nd_Ssheet[1,1]
            end
            if !isempty(idx_nd_Fl)
                idx_nd_Fl = idx_nd_Fl[1,1]
            end
            if !isempty(idx_nd_nCurt)
                idx_nd_nCurt = idx_nd_nCurt[1,1]
            end

##--- Selection of Load Profiles based upon Single-Period OR Multi-Period ----##
            if nTP_pf>1
                if !isempty(idx_nd_lsheet) && !isempty(idx_nd_Fl)
                    active_load_profile   = (nw_pPrf_data_load[idx_nd_lsheet,t+1]./nw_sbase[1].sbase)+(p_od_pf[s,t,idx_nd_Fl]-p_ud_pf[s,t,idx_nd_Fl])   # Data is in MW
                    reactive_load_profile = (nw_qPrf_data_load[idx_nd_lsheet,t+1]./nw_sbase[1].sbase)+load_theta*(p_od_pf[s,t,idx_nd_Fl]-p_ud_pf[s,t,idx_nd_Fl])   # Data is in MVAr
                elseif !isempty(idx_nd_lsheet) && isempty(idx_nd_Fl)
                    active_load_profile   = (nw_pPrf_data_load[idx_nd_lsheet,t+1]./nw_sbase[1].sbase)   # Data is in MW
                    reactive_load_profile = (nw_qPrf_data_load[idx_nd_lsheet,t+1]./nw_sbase[1].sbase)  # Data is in MVAr
                end

                if scenario == 0 || (scenario == 1 && nSc_pf == 1)
                    if !isempty(idx_nd_gsheet) && !isempty(idx_nd_nCurt)
                        active_gen_profile    = (active_gen[idx_nd_gsheet,t]-p_curt_pf[s,t,idx_nd_nCurt])          # Data is in pu
                        reactive_gen_profile  = q_dg_opf[s,t,idx_nd_nCurt]

                    elseif !isempty(idx_nd_gsheet) && isempty(idx_nd_nCurt)
                        active_gen_profile    = (active_gen[idx_nd_gsheet,t])          # Data is in pu
                        reactive_gen_profile  = (reactive_gen[idx_nd_gsheet,t])        # Data is in pu
                    end
                elseif scenario == 1 && nSc_pf > 1
                    if !isempty(idx_nd_gsheet) && !isempty(idx_nd_nCurt)
                        active_gen_profile    = active_gen[s,t,idx_nd_gsheet]-p_curt_pf[s,t,idx_nd_nCurt]         # Data is in pu
                        reactive_gen_profile  = q_dg_opf[s,t,idx_nd_nCurt]
                        # reactive_gen_profile  = reactive_gen[s,t,idx_nd_nCurt]
                    elseif !isempty(idx_nd_gsheet) && isempty(idx_nd_nCurt)
                        active_gen_profile    = active_gen[s,t,idx_nd_gsheet]       # Data is in pu
                        reactive_gen_profile  = reactive_gen[s,t,idx_nd_gsheet]        # Data is in pu
                    end
                end
            elseif nTP_pf == 1
                if !isempty(idx_nd_lsheet) && !isempty(idx_nd_Fl)
                    active_load_profile   = ((nw_loads[idx_nd_lsheet].load_P)./nw_sbase[1].sbase)+(p_od_pf[s,t,idx_nd_Fl]-p_ud_pf[s,t,idx_nd_Fl]) # Data is in MW
                    reactive_load_profile = ((nw_loads[idx_nd_lsheet].load_Q)./nw_sbase[1].sbase)+load_theta*(p_od_pf[s,t,idx_nd_Fl]-p_ud_pf[s,t,idx_nd_Fl]) # Data is in MVAr
                elseif !isempty(idx_nd_lsheet) && isempty(idx_nd_Fl)
                    active_load_profile   = ((nw_loads[idx_nd_lsheet].load_P)./nw_sbase[1].sbase)   # Data is in MW
                    reactive_load_profile = ((nw_loads[idx_nd_lsheet].load_Q)./nw_sbase[1].sbase)  # Data is in MVAr
                end

                if !isempty(idx_nd_gsheet) && !isempty(idx_nd_nCurt)
                    active_gen_profile    = active_gen[idx_nd_gsheet,1]-p_curt_pf[s,t,idx_nd_nCurt] # Data is in pu
                    reactive_gen_profile  = q_dg_opf[s,t,idx_nd_nCurt]

                    # println("Active Gen: Gen_num: $nd_num, Idx_gSheet: $idx_nd_gsheet, Idx_curt: $idx_nd_nCurt, Active_profile: $active_gen_profile")
                else
                    active_gen_profile    = active_gen[idx_nd_gsheet,1]          # Data is in pu
                    reactive_gen_profile  = reactive_gen[idx_nd_gsheet,1]        # Data is in pu
                end
            end

            if !isempty(idx_nd_Ssheet)
                actv_pwr_str =  (p_dis_pf[s,t,idx_nd_Ssheet]-p_ch_pf[s,t,idx_nd_Ssheet])
            else
                actv_pwr_str = 0
            end

            if !isempty(idx_nd_lsheet)
                active_load_profile = active_load_profile
                reactive_load_profile = reactive_load_profile
                active_shunt_pwr = nw_loads[idx_nd_lsheet].load_G./nw_sbase[1].sbase
                react_shunt_pwr  = nw_loads[idx_nd_lsheet].load_B./nw_sbase[1].sbase
            else
                active_load_profile = 0
                reactive_load_profile = 0
                active_shunt_pwr = 0
                react_shunt_pwr = 0
            end

            if !isempty(idx_nd_gsheet)
                active_gen_profile = active_gen_profile
                reactive_gen_profile = reactive_gen_profile
            else
                active_gen_profile = 0
                reactive_gen_profile = 0
            end


##------------------- Active and Reactive Power Balance Constraints----------###
            # if !isempty(idx_nd_lsheet) && !isempty(idx_nd_gsheet)                # Node contains both load and generation (P_gen-P_load)
            #     p_sch[idx_nd_bsheet,1] = (active_gen_profile+actv_pwr_str+(nw_loads[idx_nd_lsheet].load_G./nw_sbase[1].sbase))-active_load_profile
            #     q_sch[idx_nd_bsheet,1] = (reactive_gen_profile+(nw_loads[idx_nd_lsheet].load_B./nw_sbase[1].sbase))-reactive_load_profile
            #
            # elseif !isempty(idx_nd_lsheet) && isempty(idx_nd_gsheet)             # Node contains only load (0-P_load)
            #     p_sch[idx_nd_bsheet,1] = (actv_pwr_str+(nw_loads[idx_nd_lsheet].load_G./nw_sbase[1].sbase))-active_load_profile
            #     q_sch[idx_nd_bsheet,1] = (nw_loads[idx_nd_lsheet].load_B./nw_sbase[1].sbase)-reactive_load_profile
            #
            # elseif isempty(idx_nd_lsheet) && !isempty(idx_nd_gsheet)             # Node contains only generation (P_gen-0)
            #     p_sch[idx_nd_bsheet,1] = active_gen_profile+actv_pwr_str
            #     q_sch[idx_nd_bsheet,1] = reactive_gen_profile
            # else                                                                 # Node does not contains any load and generation
            #     p_sch[idx_nd_bsheet,1]  = 0.0
            #     q_sch[idx_nd_bsheet,1]  = 0.0
            # end

            p_sch[idx_nd_bsheet,1] = (active_gen_profile+actv_pwr_str+active_shunt_pwr)-active_load_profile
            q_sch[idx_nd_bsheet,1] = (reactive_gen_profile+react_shunt_pwr)-reactive_load_profile

        end
    end
    return p_sch,q_sch
end
###---------------------- Average Scenario Generation -----------------------###
function avg_scenario_gen(prob_scs,data,nTP_pf,nGens)
    data_aux = reshape(sum(prob_scs.*data[:,1:nTP_pf,2:end],dims=1),(nTP_pf,nGens-1))
    data_aux = data_aux'
    data_aux = vcat(data[1,nTP_pf,1],data_aux)
    return  data_aux
end
## -------------------- Formation of Y-bus matrix ---------------------------###
function Y_bus(nw_lines,nw_buses,trsf_ratio_opf_lin,s,t,flex_oltc)
    Ypr_l    = zeros(ComplexF64,size(nw_lines,1),size(nw_lines,1))
    Ypr_r    = zeros(ComplexF64,size(nw_lines,1),size(nw_lines,1))
    Ypr_lng  = zeros(ComplexF64,size(nw_lines,1),size(nw_lines,1))

    incd_l   = zeros(size(nw_lines,1),size(nw_buses,1))
    incd_r   = zeros(size(nw_lines,1),size(nw_buses,1))
    incd_lng = zeros(size(nw_lines,1),size(nw_buses,1))

    for i in 1:size(nw_lines,1)
        tap_ratio = nw_lines[i].line_tap_ratio
        if tap_ratio ==0
            sdg_node = nw_lines[i].line_from
            rcv_node = nw_lines[i].line_to
            zlong    = nw_lines[i].line_r + (nw_lines[i].line_x)im
            y_trsv   = nw_lines[i].line_g_shunt + (nw_lines[i].line_b_shunt)im

            idx_sdg_bsheet = findall(x->x==sdg_node,rdata_buses[:,1])            # Sending node position in rdata_buses sheet (Index of sending node in buses sheet)
            idx_rcv_bsheet = findall(x->x==rcv_node,rdata_buses[:,1])            # Receiving node position in rdata_buses sheet (Index of receiving node in buses sheet)

            Ypr_lng[i,i]  = inv(zlong)
            Ypr_l[i,i]    = y_trsv/2
            Ypr_r[i,i]    = y_trsv/2

            incd_l[i,idx_sdg_bsheet[1,1]]   = 1
            incd_r[i,idx_rcv_bsheet[1,1]]   = 1
            incd_lng[i,idx_sdg_bsheet[1,1]] = 1
            incd_lng[i,idx_rcv_bsheet[1,1]] = -1
        end
    end

    y_long = transpose(incd_lng)*Ypr_lng*incd_lng                                # This command is taking a lot of time to execute!
    y_l    = transpose(incd_l)*Ypr_l*incd_l
    y_r    = transpose(incd_r)*Ypr_r*incd_r

    y_bus = y_l.+y_long.+y_r

    idx_tap_line_sheet = findall(x->x=="ratio",rheader_lines)
    idx_trsf = findall(x->x!=0,rdata_lines[:,idx_tap_line_sheet[1,1]])

    for i in 1:size(idx_trsf,1)
        sdg_node = nw_lines[idx_trsf[i]].line_from
        rcv_node = nw_lines[idx_trsf[i]].line_to
        zlong    = nw_lines[idx_trsf[i]].line_r + (nw_lines[idx_trsf[i]].line_x)im
        y_pr     = inv(zlong)
        y_trsv   = nw_lines[idx_trsf[i]].line_g_shunt + (nw_lines[idx_trsf[i]].line_b_shunt)im
        idx_sdg_bsheet = findall(x->x==sdg_node,rdata_buses[:,1])
        idx_rcv_bsheet = findall(x->x==rcv_node,rdata_buses[:,1])
        if flex_oltc == 0
            tap_ratio = nw_lines[idx_trsf[i]].line_tap_ratio
        else
            tap_ratio = trsf_ratio_opf_lin[s,t,i]
        end
        y_bus[idx_sdg_bsheet[1,1],idx_sdg_bsheet[1,1]] = y_bus[idx_sdg_bsheet[1,1],idx_sdg_bsheet[1,1]]+tap_ratio^2*(y_pr+y_trsv/2)
        y_bus[idx_sdg_bsheet[1,1],idx_rcv_bsheet[1,1]] = y_bus[idx_sdg_bsheet[1,1],idx_rcv_bsheet[1,1]]+tap_ratio*(-y_pr)
        y_bus[idx_rcv_bsheet[1,1],idx_sdg_bsheet[1,1]] = y_bus[idx_rcv_bsheet[1,1],idx_sdg_bsheet[1,1]]+tap_ratio*(-y_pr)
        y_bus[idx_rcv_bsheet[1,1],idx_rcv_bsheet[1,1]] = y_bus[idx_rcv_bsheet[1,1],idx_rcv_bsheet[1,1]]+(y_pr+y_trsv/2)
    end
return y_bus
end
## -------------------- Formation of Jacobian Matrix -----------------------###
function jacobian(nw_buses,rdata_buses,v_angle,v_mag,y_bus,J11,J12,J21,J22,idx_dPV_buses,s,t)
    idx_i = 0
    idx_j = 0

    for i in 1:size(nw_buses,1)                                                 # Outer loop controls the indexing over rows of each sub-matrix in Jacobian matrix
        nd_num_i   = nw_buses[i].bus_num
        bus_type_i = nw_buses[i].bus_type

        if bus_type_i==3
            idx_i = idx_i+1
        elseif bus_type_i!=3
            idx_nd_bsheet_i = findall(x->x==nd_num_i,rdata_buses[:,1])
            idx_nd_bsheet_i = idx_nd_bsheet_i[1,1]-idx_i
            for j in 1:size(nw_buses,1)                                          # Inner loop controls the indexing over columns of each sub-matrix in Jacobian matrix
                nd_num_j   = nw_buses[j].bus_num
                bus_type_j = nw_buses[j].bus_type

                if bus_type_j==3
                    idx_j = idx_j+1

                elseif bus_type_j!=3
                    if nd_num_i == nd_num_j                                      # Diagonal elements of Jacobian sub-matrices
                        idx_nd_bsheet_j = findall(x->x==nd_num_j,rdata_buses[:,1])
                        idx_nd_bsheet_j = idx_nd_bsheet_j[1,1]-idx_j
                        theta_diff_jj   = v_angle[idx_nd_bsheet_j+idx_j,t,s].-v_angle[:,t,s]
                        p_inj = v_mag[idx_nd_bsheet_j+idx_j,t,s]*sum(v_mag[:,t,s].*((real(y_bus[idx_nd_bsheet_j+idx_j,:])).*cos.(theta_diff_jj)+(imag(y_bus[idx_nd_bsheet_j+idx_j,:])).*sin.(theta_diff_jj)),dims=1)
                        q_inj = v_mag[idx_nd_bsheet_j+idx_j,t,s]*sum(v_mag[:,t,s].*((real(y_bus[idx_nd_bsheet_j+idx_j,:])).*sin.(theta_diff_jj)-(imag(y_bus[idx_nd_bsheet_j+idx_j,:])).*cos.(theta_diff_jj)),dims=1)
                        J11[idx_nd_bsheet_j,idx_nd_bsheet_j] = -q_inj[1,1]-v_mag[idx_nd_bsheet_j+idx_j,t,s]^2*imag(y_bus[idx_nd_bsheet_j+idx_j,idx_nd_bsheet_j+idx_j])
                        J21[idx_nd_bsheet_j,idx_nd_bsheet_j] =  p_inj[1,1]-v_mag[idx_nd_bsheet_j+idx_j,t,s]^2*real(y_bus[idx_nd_bsheet_j+idx_j,idx_nd_bsheet_j+idx_j])
                        J12[idx_nd_bsheet_j,idx_nd_bsheet_j] =  p_inj[1,1]+v_mag[idx_nd_bsheet_j+idx_j,t,s]^2*real(y_bus[idx_nd_bsheet_j+idx_j,idx_nd_bsheet_j+idx_j])
                        J22[idx_nd_bsheet_j,idx_nd_bsheet_j] =  q_inj[1,1]-v_mag[idx_nd_bsheet_j+idx_j,t,s]^2*imag(y_bus[idx_nd_bsheet_j+idx_j,idx_nd_bsheet_j+idx_j])
                    else                                                         # Off-diagonal elements of Jacobian sub-matrices
                        idx_nd_bsheet_j = findall(x->x==nd_num_j,rdata_buses[:,1])
                        idx_nd_bsheet_j = idx_nd_bsheet_j[1,1]-idx_j
                        theta_diff_ij   = v_angle[(idx_nd_bsheet_i+idx_i),t,s]-v_angle[(idx_nd_bsheet_j+idx_j),t,s]
                        partial_p_inj   = +1*(v_mag[idx_nd_bsheet_i+idx_i,t,s]*v_mag[idx_nd_bsheet_j+idx_j,t,s])*(real(y_bus[idx_nd_bsheet_i+idx_i,idx_nd_bsheet_j+idx_j])*sin(theta_diff_ij)-imag(y_bus[idx_nd_bsheet_i+idx_i,idx_nd_bsheet_j+idx_j])*cos(theta_diff_ij))
                        partial_q_inj   = -1*(v_mag[idx_nd_bsheet_i+idx_i,t,s]*v_mag[idx_nd_bsheet_j+idx_j,t,s])*(real(y_bus[idx_nd_bsheet_i+idx_i,idx_nd_bsheet_j+idx_j])*cos(theta_diff_ij)+imag(y_bus[idx_nd_bsheet_i+idx_i,idx_nd_bsheet_j+idx_j])*sin(theta_diff_ij))
                        J11[idx_nd_bsheet_i,idx_nd_bsheet_j] = partial_p_inj
                        J21[idx_nd_bsheet_i,idx_nd_bsheet_j] = partial_q_inj
                        J12[idx_nd_bsheet_i,idx_nd_bsheet_j] = -partial_q_inj
                        J22[idx_nd_bsheet_i,idx_nd_bsheet_j] = partial_p_inj
                    end
                end
            end                 # End inner for-loop
            idx_j = 0
        end
    end                         # End outer for-loop
    # idx_i = 0
    if isempty(idx_dPV_buses)   # No PV bus in a network
        J = [J11 J12;J21 J22]
    else                        # Deleting rows and columns which corresponds to PV buses in J12, J21 and J22 submatrices.
        J12_tmp = J12[:,setdiff(1:end,idx_dPV_buses)]
        J21_tmp = J21[setdiff(1:end,idx_dPV_buses),:]
        J22_tmp = J22[setdiff(1:end,idx_dPV_buses),setdiff(1:end,idx_dPV_buses)]
        J = [J11 J12_tmp;J21_tmp J22_tmp]
    end
    return J
end

## -------------------------- Formation of Mismatch Vector ------------------###
function mismatch_vector(nw_buses,rdata_buses,v_angle,v_mag,y_bus,p_node,q_node,idx_dPV_buses,p_sch_pf,q_sch_pf,s,t)
    idx = 0
    for i in 1:size(nw_buses,1)                                                  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses
        bus_type = nw_buses[i].bus_type
        nd_num   = nw_buses[i].bus_num

        if bus_type ==3
            idx = idx+1
        elseif bus_type!=3                                                       # Injection for PQ and PV buses
            idx_nd_bsheet = findall(x->x==nd_num,rdata_buses[:,1])
            idx_nd_bsheet = idx_nd_bsheet[1,1]-idx
            theta_diff    = v_angle[idx_nd_bsheet+idx,t,s].-v_angle[:,t,s]
            p_inj = v_mag[idx_nd_bsheet+idx,t,s]*sum(v_mag[:,t,s].*((real(y_bus[idx_nd_bsheet+idx,:])).*cos.(theta_diff)+(imag(y_bus[idx_nd_bsheet+idx,:])).*sin.(theta_diff)),dims=1)
            q_inj = v_mag[idx_nd_bsheet+idx,t,s]*sum(v_mag[:,t,s].*((real(y_bus[idx_nd_bsheet+idx,:])).*sin.(theta_diff)-(imag(y_bus[idx_nd_bsheet+idx,:])).*cos.(theta_diff)),dims=1)
            p_node[idx_nd_bsheet,1] = p_inj[1,1]
            q_node[idx_nd_bsheet,1] = q_inj[1,1]
        end
    end
    delta_p = p_sch_pf.-p_node    # Mismatch active power
    delta_q = q_sch_pf.-q_node    # Mismatch reactive power

    if isempty(idx_dPV_buses)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p,delta_q)
    else                                                                         # Existence of PV buses in a system
        delta_q_tmp = deleteat!(delta_q,idx_dPV_buses)                           # Remove the corresponding reactive power mismatch entry from the delta_Y vector
        delta_mismatch_inj = vcat(delta_p,delta_q_tmp)
    end
return delta_mismatch_inj
end

## -------------------- Handling PV and PQ buses data -----------------------###
##-----Indices of PV buses to be deleted from Jacobian and mismatch vector----##
function pv_pq_data(idx_PV_buses,idx_slack_bus,idx_PQ_buses)
            idx_dPV_buses = idx_PV_buses.-idx_slack_bus                          # Indices of deleted PV buses

            if all(idx_dPV_buses .>0)                                            # All PV buses are placed below the slack bus in the buses sheet
                idx_dPV_buses = idx_PV_buses.-1
            elseif all(idx_dPV_buses .<0)                                        # All PV buses are placed above the slack bus in the buses sheet
                idx_dPV_buses = idx_PV_buses
            else                                                                 # Few PV buses are below and few are above the slack bus in the buses sheet
                for i in 1:size(idx_dPV_buses,1)
                    if idx_dPV_buses[i]<0
                        idx_dPV_buses[i] = idx_PV_buses[i]
                    else
                        idx_dPV_buses[i] =  idx_PV_buses[i]-1
                    end
                end
            end
##---Indices of modified PQ buses obtained with respect to the slack bus ----##
            idx_dPQ_buses = idx_PQ_buses.-idx_slack_bus                          # Indices of modified PQ buses

            if all(idx_dPQ_buses .>0)                                            # All PQ buses are placed below the slack bus in the buses sheet
                idx_dPQ_buses = idx_PQ_buses.-1
            elseif all(idx_dPQ_buses .<0)                                        # All PQ buses are placed below the slack bus in the buses sheet
                idx_dPQ_buses = idx_PQ_buses
            else                                                                 # Few PQ buses are below and few are above the slack bus in the buses sheet
                for i in 1:size(idx_dPQ_buses,1)
                    if idx_dPQ_buses[i]<0
                        idx_dPQ_buses[i] = idx_PQ_buses[i]
                    else
                        idx_dPQ_buses[i] =  idx_PQ_buses[i]-1
                    end
                end
            end
return idx_dPV_buses, idx_dPQ_buses
end

## ------ Code to retrieve the voltage magnitude and angle data ------------ ###
##----------- based upon the slack bus and the presence of PV buses -------- ###
function solution_voltage(slack_bus_type,rdata_buses,nw_buses,idx_PQ_buses,idx_PV_buses,idx_dPQ_buses,idx_dPV_buses,v_angle,v_mag,delta_x,idx_slack_bus,s,t)
    if isempty(idx_PV_buses)
        idx_slack_bus = findall(x->x==slack_bus_type,rdata_buses[:,2])           # All buses are PQ buses and there is no PV bus

            if !isempty(idx_slack_bus)
                idx_slack_bus = idx_slack_bus[1,1]
            end

        if isempty(idx_slack_bus)
            v_angle[1:end,t,s] = v_angle[1:end,t,s]+delta_x[1:Int64(size(delta_x,1)/2),1]
            v_mag[1:end,t,s]   = v_mag[1:end,t,s].*(1.0 .+delta_x[Int64((size(delta_x,1)/2)+1):end,1])
        elseif idx_slack_bus == 1                                                # First bus in the Buses sheet is slack bus
            v_angle[idx_slack_bus+1:end,t,s] = v_angle[idx_slack_bus+1:end,t,s]+delta_x[1:Int64(size(delta_x,1)/2),1]
            v_mag[idx_slack_bus+1:end,t,s]   = v_mag[idx_slack_bus+1:end,t,s].*(1.0 .+delta_x[Int64((size(delta_x,1)/2)+1):end,1])
        elseif idx_slack_bus == size(nw_buses,1)                                 # Last bus in the Buses sheet is slack bus
            v_angle[1:idx_slack_bus-1,t,s] = v_angle[1:idx_slack_bus-1,t,s]+delta_x[1:Int64(size(delta_x,1)/2),1]
            v_mag[1:idx_slack_bus-1,t,s]   = v_mag[1:idx_slack_bus-1,t,s].*(1.0 .+delta_x[Int64((size(delta_x,1)/2)+1):end,1])
        elseif idx_slack_bus>1 && idx_slack_bus < size(nw_buses,1)               # Few buses are placed above and few are placed below the slack bus in the buses sheet
            before_length = length(1:idx_slack_bus-1)
            after_length  = length(idx_slack_bus+1:size(nw_buses,1))
            v_angle[1:idx_slack_bus-1,t,s]       = v_angle[1:idx_slack_bus-1,t,s]+delta_x[1:before_length,1]
            v_angle[idx_slack_bus+1:end,t,s]     = v_angle[idx_slack_bus+1:end,t,s]+delta_x[before_length+1:Int64(size(delta_x,1)/2),1]
            v_mag[1:idx_slack_bus-1,t,s]         = v_mag[1:idx_slack_bus-1,t,s].*(1.0 .+delta_x[Int64((size(delta_x,1)/2)+1):Int64(size(delta_x,1)/2)+before_length,1])
            v_mag[idx_slack_bus+1:end,t,s]       = v_mag[idx_slack_bus+1:end,t,s].*(1.0 .+delta_x[Int64(size(delta_x,1)/2)+before_length+1:end,1])
        end
    else                                                        # Network contains both PQ and PV buses. Voltage magnitude at PV buses is kept constant.
        delta_angle = delta_x[1:size(nw_buses,1)-1]
        delta_mag   = delta_x[size(nw_buses,1):end]
        v_angle[idx_PQ_buses,t,s] = v_angle[idx_PQ_buses,t,s]+delta_angle[idx_dPQ_buses]
        v_angle[idx_PV_buses,t,s] = v_angle[idx_PV_buses,t,s]+delta_angle[idx_dPV_buses]
        v_mag[idx_PQ_buses,t,s]   = v_mag[idx_PQ_buses,t,s].*(1.0.+delta_mag)
    end
return v_angle,v_mag
end

## -------------- Determining reactive generation at PV buses -------------####
# and switching PV bus to PQ bus in case Q-Gen hits the minimum or maximum limit#
function pv_pq_switching(y_bus,nw_buses_pf,nw_loads_pf,nw_gens_pf,v_mag_pf,v_angle_pf,nw_sbase_pf,rdata_buses,rdata_loads,rdata_gens,idx_oPV_buses_pf,load_bus_type,vol_ctrl_bus_type,mdf_idx_PV_buses_pf,s,t)

    q_gen_pv_pf = zeros(size(idx_oPV_buses_pf,1))

    for i in 1:size(idx_oPV_buses_pf,1)                                          # Reactive power injection at PV buses
        nd_num        = nw_buses_pf[idx_oPV_buses_pf[i]].bus_num
        idx_nd_bsheet = findall(x->x==nd_num,rdata_buses[:,1])
        idx_nd_lsheet = findall(x->x==nd_num,rdata_loads[:,1])
        idx_nd_gsheet = findall(x->x==nd_num,rdata_gens[:,1])
        q_gen_max     = nw_gens_pf[idx_nd_gsheet[1,1]].gen_Qg_max/nw_sbase_pf[1].sbase
        q_gen_min     = nw_gens_pf[idx_nd_gsheet[1,1]].gen_Qg_min/nw_sbase_pf[1].sbase
        theta_diff    = v_angle_pf[idx_oPV_buses_pf[i],t,s].-v_angle_pf[:,t,s]
        q_inj         = v_mag_pf[idx_oPV_buses_pf[i],t,s]*sum(v_mag_pf[:,t,s].*((real(y_bus[idx_oPV_buses_pf[i],:])).*sin.(theta_diff)-(imag(y_bus[idx_oPV_buses_pf[i],:])).*cos.(theta_diff)),dims=1)
            if !isempty(idx_nd_lsheet)                                           # Load connected to the PV bus
                q_gen_pv_pf[i,1] = q_inj[1,1] + (nw_loads_pf[idx_nd_lsheet[1,1]].load_Q-nw_loads_pf[idx_nd_lsheet[1,1]].load_B)/nw_sbase_pf[1].sbase
            else                                                                 # No load connected to the PV bus and hence there is no mentioning of PV bus number in the load sheet
                q_gen_pv_pf[i,1] = q_inj[1,1]
            end
            if q_gen_pv_pf[i,1] > q_gen_max                                      # Generated power becomes greator than the upper reactive power generation limit
                pwr = q_gen_pv_pf[i,1]
                nw_buses_pf[idx_nd_bsheet[1,1]].bus_type = load_bus_type         # Switch to PQ bus
                rdata_buses[idx_nd_bsheet[1,1],2] = load_bus_type                # Switch to PQ bus
                nw_gens_pf[idx_nd_gsheet[1,1]].gen_Qg_avl = q_gen_max*nw_sbase_pf[1].sbase   # Setting Q_gen to the maximum limit
                mdf_idx_PV_buses_pf[i] = load_bus_type
                println("Generator G$nd_num: Q = $pwr >= Qmax = $q_gen_max => switched as PQ!")
            elseif q_gen_pv_pf[i,1] < q_gen_min                                 # Generated power becomes lower than the lower reactive power generation limit
                pwr = q_gen_pv_pf[i,1]
                nw_buses_pf[idx_nd_bsheet[1,1]].bus_type = load_bus_type        # Switch to PQ bus
                rdata_buses[idx_nd_bsheet[1,1],2] = load_bus_type            # Switch to PQ bus
                nw_gens_pf[idx_nd_gsheet[1,1]].gen_Qg_avl = q_gen_min*nw_sbase_pf[1].sbase    # Setting Q_gen to the minimum limit
                mdf_idx_PV_buses_pf[i] = load_bus_type
                println("Generator G$nd_num: Q = $pwr <= Qmin = $q_gen_min => switched as PQ!")
                q_gen_min <= q_gen_pv_pf[i,1] <=q_gen_max && mdf_idx_PV_buses_pf[i] == 1   # Qgen is within limits but the bus has switched from PV to PQ in the previous iteration
                pwr = q_gen_pv_pf[i,1]
                mdf_idx_PV_buses_pf[i] = 0
                nw_buses_pf[idx_nd_bsheet[1,1]].bus_type = vol_ctrl_bus_type            # Switch the PQ bus back to PV bus in case there is no violation of reactive power generation limit
                rdata_buses[idx_nd_bsheet[1,1],2] = vol_ctrl_bus_type
                println("Generator G$nd_num: Qmin = $q_gen_min <= Q = $pwr <= Qmax = $q_gen_max => switched back to PV!")
            end
    end
return nw_buses_pf,rdata_buses,nw_gens_pf
end
