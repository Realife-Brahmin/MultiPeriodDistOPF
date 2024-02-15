###-------- This file contains sub-functions of OPF model module ------------###
################################################################################
###-------- Full expression of active/reactive power flow at ----- #############
###------- Sending end of transformer active/reactive power flows -----------###
################################################################################
function power_flow_linear_theta_sdg_trsf_milp(acopf,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,trsf_fict_nodes,flex_oltc,status_oltc)

    line    = [nd_num cnctd_nd]
    id_ft   = line.-ft_line
    id_tf   = line.-tf_line
    line_ft = findall(all(id_ft.==idx_zero,dims=2))
    line_tf = findall(all(id_tf.==idx_zero,dims=2))

    if !isempty(line_ft) line_ft = line_ft[1][1] end
    if !isempty(line_tf) line_tf = line_tf[1][1] end

    if flex_oltc == 1 && status_oltc == 1                                        # OLTC transformer with status set to 1 (Need to find the index of fictitious node)
        if !isempty(line_ft) && isempty(line_tf)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,1])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,3])
        elseif !isempty(line_tf) && isempty(line_ft)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,3])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,1])
        end
        idx_fict_nd  = intersect(from_nd_trsf,to_nd_trsf)
        idx_fict_nd  = idx_fict_nd[1,1]
        # fict_nd      = trsf_fict_nodes[idx_fict_nd,2]
        fict_nd      = trsf_fict_nodes[idx_fict_nd,4]
        v_fbus_var   = v_mag[s,t,fict_nd]
        v_fbus       = v_mag_int[s,t,fict_nd]
    elseif (flex_oltc == 0 || flex_oltc == 1) && status_oltc == 0
        v_fbus_var = v_mag[s,t,idx_nd_nw_buses]
        v_fbus     = v_mag_int[s,t,idx_nd_nw_buses]
    end

    # v_fbus     = v_mag_int[s,t,idx_nd_nw_buses]
    v_tbus     = v_mag_int[s,t,idx_cnctd_nd]
    v_tbus_var = v_mag[s,t,idx_cnctd_nd]

    if !isempty(line_ft)&&isempty(line_tf)
        theta_diff = @expression(acopf,+1*v_ang[s,t,line_ft])
    elseif isempty(line_ft)&& !isempty(line_tf)
        theta_diff = @expression(acopf,-1*v_ang[s,t,line_tf])
    end

    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]  # v_ang_int has to be in radians!

    c_int = cos(theta_diff_int)    # cosine argument must be in radians!
    s_int = sin(theta_diff_int)    # sine argument must be in radians!

    alpha_i_P = tap*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)
    alpha_j_P = tap*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)

    alpha_i_Q = tap*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)
    alpha_j_Q = tap*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)

    beta_ij_P = tap*v_fbus*v_tbus*(s_int*gij_line-c_int*bij_line)
    beta_ij_Q = tap*v_fbus*v_tbus*(-s_int*bij_line-c_int*gij_line)

    gamma_P = tap*gij_line*((-v_fbus*v_tbus*theta_diff_int*s_int)-(c_int*(v_fbus-v_tbus)^2)/2)+tap*bij_line*((v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)
    gamma_Q = tap*bij_line*((v_fbus*v_tbus*theta_diff_int*s_int)+(c_int*(v_fbus-v_tbus)^2)/2)+tap*gij_line*((+v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)

    # offset_p = v_fbus*v_tbus*theta_diff_int*(c_int*bij_line-s_int*gij_line)+(((v_fbus-v_tbus)^2)/2)*(-c_int*gij_line-s_int*bij_line)
    # offset_q = v_fbus*v_tbus*theta_diff_int*(s_int*bij_line+c_int*gij_line)+(((v_fbus-v_tbus)^2)/2)*(c_int*bij_line-s_int*gij_line)

    pij = @expression(acopf,+tap^2*(gij_line_sh+gij_line)*v_fbus_var+(alpha_i_P*v_fbus_var)+(alpha_j_P*v_tbus_var)+(beta_ij_P*theta_diff)+gamma_P)
    qij = @expression(acopf,-tap^2*(bij_line_sh+bij_line)*v_fbus_var+(alpha_i_Q*v_fbus_var)+(alpha_j_Q*v_tbus_var)+(beta_ij_Q*theta_diff)+gamma_Q)

return pij, qij
end
################################################################################
###-------- Full expression of active/reactive power flow at ----- #############
### ------ Receiving end of transformer active/reactive power flows ---------###
################################################################################
function power_flow_linear_theta_rcv_trsf_milp(acopf,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,trsf_fict_nodes,flex_oltc,status_oltc)

    line    = [nd_num cnctd_nd]
    id_ft   = line.-ft_line
    id_tf   = line.-tf_line
    line_ft = findall(all(id_ft.==idx_zero,dims=2))
    line_tf = findall(all(id_tf.==idx_zero,dims=2))

    if !isempty(line_ft) line_ft = line_ft[1][1] end
    if !isempty(line_tf) line_tf = line_tf[1][1] end

    if flex_oltc == 1 && status_oltc == 1
        if !isempty(line_ft) && isempty(line_tf)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,1])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,3])
        elseif !isempty(line_tf) && isempty(line_ft)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,3])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,1])
        end
        idx_fict_nd  = intersect(from_nd_trsf,to_nd_trsf)
        idx_fict_nd  = idx_fict_nd[1,1]
        # fict_nd      = trsf_fict_nodes[idx_fict_nd,2]
        fict_nd      = trsf_fict_nodes[idx_fict_nd,4]
        v_tbus_var   = v_mag[s,t,fict_nd]
        v_tbus       = v_mag_int[s,t,fict_nd]
    elseif (flex_oltc == 0 || flex_oltc == 1) && status_oltc == 0
        v_tbus_var   = v_mag[s,t,idx_cnctd_nd]
        v_tbus     = v_mag_int[s,t,idx_cnctd_nd]
    end

    v_fbus     = v_mag_int[s,t,idx_nd_nw_buses]
    v_fbus_var = v_mag[s,t,idx_nd_nw_buses]
    # v_tbus_var = v_mag[s,t,idx_cnctd_nd]

    if !isempty(line_ft)&&isempty(line_tf)
        theta_diff = @expression(acopf,+1*v_ang[s,t,line_ft])
    elseif isempty(line_ft)&& !isempty(line_tf)
        theta_diff = @expression(acopf,-1*v_ang[s,t,line_tf])
    end

    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]

    c_int = cos(theta_diff_int)
    s_int = sin(theta_diff_int)

    alpha_i_P = tap2*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)
    alpha_j_P = tap2*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)

    alpha_i_Q = tap2*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)
    alpha_j_Q = tap2*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)

    beta_ij_P = tap2*v_fbus*v_tbus*(s_int*gij_line-c_int*bij_line)
    beta_ij_Q = tap2*v_fbus*v_tbus*(-s_int*bij_line-c_int*gij_line)

    gamma_P = tap2*gij_line*((-v_fbus*v_tbus*theta_diff_int*s_int)-(c_int*(v_fbus-v_tbus)^2)/2)+tap2*bij_line*((v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)
    gamma_Q = tap2*bij_line*((v_fbus*v_tbus*theta_diff_int*s_int)+(c_int*(v_fbus-v_tbus)^2)/2)+tap2*gij_line*((+v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)

    # offset_p = v_fbus*v_tbus*theta_diff_int*(c_int*bij_line-s_int*gij_line)+(((v_fbus-v_tbus)^2)/2)*(-c_int*gij_line-s_int*bij_line)
    # offset_q = v_fbus*v_tbus*theta_diff_int*(s_int*bij_line+c_int*gij_line)+(((v_fbus-v_tbus)^2)/2)*(c_int*bij_line-s_int*gij_line)

    pij = @expression(acopf,+tap1^2*(gij_line)*v_fbus_var+(alpha_i_P*v_fbus_var)+(alpha_j_P*v_tbus_var)+(beta_ij_P*theta_diff)+gamma_P)
    qij = @expression(acopf,-tap1^2*(bij_line)*v_fbus_var+(alpha_i_Q*v_fbus_var)+(alpha_j_Q*v_tbus_var)+(beta_ij_Q*theta_diff)+gamma_Q)

return pij, qij
end

################################################################################
###----- Full expression of active/reactive power flow for simple branch ----- #
################################################################################
function power_flow_linear_branch_theta_milp(acopf,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero)

    line    = [nd_num cnctd_nd]
    id_ft   = line.-ft_line
    id_tf   = line.-tf_line
    line_ft = findall(all(id_ft.==idx_zero,dims=2))
    line_tf = findall(all(id_tf.==idx_zero,dims=2))

    if !isempty(line_ft) line_ft = line_ft[1][1] end
    if !isempty(line_tf) line_tf = line_tf[1][1] end

    v_fbus     = v_mag_int[s,t,idx_nd_nw_buses]
    v_tbus     = v_mag_int[s,t,idx_cnctd_nd]

    v_fbus_var = v_mag[s,t,idx_nd_nw_buses]
    v_tbus_var = v_mag[s,t,idx_cnctd_nd]

    if !isempty(line_ft)&&isempty(line_tf)
        theta_diff = @expression(acopf,+1*v_ang[s,t,line_ft])
    elseif isempty(line_ft)&& !isempty(line_tf)
        theta_diff = @expression(acopf,-1*v_ang[s,t,line_tf])
    end

    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]

    c_int = cos(theta_diff_int)
    s_int = sin(theta_diff_int)

    alpha_i_P = ((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)
    alpha_j_P = ((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)

    alpha_i_Q = ((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)
    alpha_j_Q = ((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)

    beta_ij_P = v_fbus*v_tbus*(s_int*gij_line-c_int*bij_line)
    beta_ij_Q = v_fbus*v_tbus*(-s_int*bij_line-c_int*gij_line)

    gamma_P = gij_line*((-v_fbus*v_tbus*theta_diff_int*s_int)-(c_int*(v_fbus-v_tbus)^2)/2)+bij_line*((v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)
    gamma_Q = bij_line*((v_fbus*v_tbus*theta_diff_int*s_int)+(c_int*(v_fbus-v_tbus)^2)/2)+gij_line*((+v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)

    # offset_p = v_fbus*v_tbus*theta_diff_int*(c_int*bij_line-s_int*gij_line)+(((v_fbus-v_tbus)^2)/2)*(-c_int*gij_line-s_int*bij_line)
    # offset_q = v_fbus*v_tbus*theta_diff_int*(s_int*bij_line+c_int*gij_line)+(((v_fbus-v_tbus)^2)/2)*(c_int*bij_line-s_int*gij_line)

    pij = @expression(acopf,+(gij_line_sh/2+gij_line)*v_fbus_var+(alpha_i_P*v_fbus_var)+(alpha_j_P*v_tbus_var)+(beta_ij_P*theta_diff)+gamma_P)
    qij = @expression(acopf,-(bij_line_sh/2+bij_line)*v_fbus_var+(alpha_i_Q*v_fbus_var)+(alpha_j_Q*v_tbus_var)+(beta_ij_Q*theta_diff)+gamma_Q)
return pij, qij
end
################################################################################
################## Linear Longitudnal Current Expression #######################
################################################################################
function longitudinal_current_theta_milp(acopf,tap,v_mag,v_mag_int,v_ang,v_ang_int,gij_line,bij_line,s,t,idx_nd_nw_buses,idx_cnctd_nd,idx_line)
    v_fbus = v_mag_int[s,t,idx_nd_nw_buses]
    v_tbus = v_mag_int[s,t,idx_cnctd_nd]

    v_fbus_var = v_mag[s,t,idx_nd_nw_buses]
    v_tbus_var = v_mag[s,t,idx_cnctd_nd]

    theta_diff     = @expression(acopf,v_ang[s,t,idx_line])
    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]

    c_int = cos(theta_diff_int)
    s_int = sin(theta_diff_int)

    c_cnst = (tap*c_int)/(v_fbus+v_tbus)

    alpha_i_I  = (tap^2+c_cnst*(v_fbus-3*v_tbus))
    alpha_j_I  = (1+c_cnst*(v_tbus-3*v_fbus))
    beta_ij_I  = 2*tap*s_int*v_fbus*v_tbus
    gamma_ij_I = tap*(-2*s_int*v_fbus*v_tbus*theta_diff_int)+tap*(-c_int*(v_fbus-v_tbus)^2)

    Iij = @expression(acopf,(gij_line^2+bij_line^2)*(alpha_i_I*v_fbus_var+alpha_j_I*v_tbus_var+beta_ij_I*theta_diff+gamma_ij_I))
return Iij
end
################################################################################
################## Linear Longitudnal Current Recovery #########################
################################################################################
function longitudinal_current_theta_recovery_milp(tap,v_mag,v_mag_int,v_ang,v_ang_int,gij_line,bij_line,s,t,idx_nd_nw_buses,idx_cnctd_nd,idx_line)
    v_fbus = v_mag_int[s,t,idx_nd_nw_buses]
    v_tbus = v_mag_int[s,t,idx_cnctd_nd]

    v_fbus_var = v_mag[s,t,idx_nd_nw_buses].^2
    v_tbus_var = v_mag[s,t,idx_cnctd_nd].^2

    theta_diff     = v_ang[s,t,idx_line]
    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]

    c_int = cos(theta_diff_int)
    s_int = sin(theta_diff_int)

    c_cnst = (tap*c_int)/(v_fbus+v_tbus)
    alpha_i_I = (tap^2+c_cnst*(v_fbus-3*v_tbus))
    alpha_j_I = (1+c_cnst*(v_tbus-3*v_fbus))
    beta_ij_I = 2*tap*s_int*v_fbus*v_tbus
    gamma_i_I = tap*(-2*s_int*v_fbus*v_tbus*theta_diff_int)+tap*(-c_int*(v_fbus-v_tbus)^2)

    Iij = (gij_line^2+bij_line^2)*(alpha_i_I*v_fbus_var+alpha_j_I*v_tbus_var+beta_ij_I*theta_diff+gamma_i_I)
return Iij
end
###############################################################################
#################### Angle Retrieving Function ################################
###############################################################################
function angle_recovery(nw_buses,nw_lines,node_data,v_ang,nSc,nTP,nBus,nTrsf,idx_slack_bus_pf,rdata_buses)
    node_nxt_itr  = []
    node_process  = []       # Node which has been processed
    line_process  = []       # Line which has been processed
    if flex_oltc==1
        # dim_bus = nBus+nTrsf
        dim_bus = nBus
    else
        dim_bus = nBus
    end

    v_ang_rcvr    = zeros(nSc,nTP,dim_bus)
    # nBus    = size(nw_buses,1)
    v_acp   = zeros(Float64,dim_bus,length(1:nTP),(length(1:nSc)))
    va_acp  = zeros(Float64,dim_bus,length(1:nTP),(length(1:nSc)))
    node_cn_node = zeros(dim_bus,dim_bus)
    for i in 1:dim_bus
        idx_c_nodes_tmp = []
        # idx_c_node = node_data[i].node_cnode
        # node_cn_node[idx_c_node,i] .= 1
        c_nodes    = node_data[i].node_cnode
        for j in 1:size(c_nodes,1)
            tmp_data = findall(x->x==c_nodes[j,1],rdata_buses[:,1])
            tmp_data = tmp_data[1,1]
            push!(idx_c_nodes_tmp,tmp_data)
        end
        node_cn_node[idx_c_nodes_tmp,i] .= 1
    end
    va = v_ang
    idx_slack_bus = idx_slack_bus_pf[1,1]
    ang_slck_bus  = 0.0

    for i in 1:dim_bus
        if i==1                                                                      # Slack Node
            slack_node = idx_slack_bus                                               # Nodes which are connected to Slack bus
            (node_nxt_itr,cn_nodes) = node_next_iteration(node_cn_node,node_process,node_nxt_itr,idx_slack_bus)
            (v_ang_rcvr,line_process,node_process) = node_angle(slack_node,node_data,line_process,node_process,nw_lines,cn_nodes,v_ang_rcvr,va,i,rdata_buses)
        else                                                                         # Nodes other than Slack Node
            node_crnt_itr_tmp = deepcopy(node_nxt_itr)
            node_nxt_itr      = []
            global cn_nodes
############# Nodes to be processed in next iterartion #########################
            for l in 1:length(node_crnt_itr_tmp)
                idx_node = node_crnt_itr_tmp[l,1]
                idx_node = idx_node[1,1]
                (node_nxt_itr,cn_nodes) = node_next_iteration(node_cn_node,node_process,node_nxt_itr,idx_node)
            end
            (v_ang_rcvr,line_process,node_process) = node_angle(node_crnt_itr_tmp,node_data,line_process,node_process,nw_lines,cn_nodes,v_ang_rcvr,va,i,rdata_buses)
        end
    end
return v_ang_rcvr
end

function node_angle(crnt_node,node_data,line_process,node_process,nw_lines,node_nxt_itr,v_ang_rcvr,va,i,rdata_buses)

    for j in 1:size(crnt_node,1)
        idx_node  = crnt_node                                                    # This is the index of the node
        idx_node  = idx_node[j,1]
        line      = node_data[idx_node].node_iline

        for k in 1:length(line)
            line_tmp = findall(x->x==line[k,1],line_process[:,1])
            if isempty(line_tmp)
                from_node = nw_lines[line[k,1]].line_from
                to_node   = nw_lines[line[k,1]].line_to
                cn_node   = node_data[idx_node].node_cnode[k]

                idx_from_node = findall(x->x==from_node,rdata_buses[:,1])
                idx_to_node   = findall(x->x==to_node,rdata_buses[:,1])

                if cn_node == from_node
                    node_angle =  v_ang_rcvr[:,:,idx_node].+va[:,:,line[k,1]]
                    v_ang_rcvr[:,:,idx_from_node] = node_angle
                elseif cn_node == to_node
                    node_angle =  v_ang_rcvr[:,:,idx_node].-va[:,:,line[k,1]]
                    v_ang_rcvr[:,:,idx_to_node]   = node_angle
                end
                push!(line_process,line[k,1])
            end
        end
        push!(node_process,idx_node)
    end
    return v_ang_rcvr,line_process,node_process
end

function node_next_iteration(node_cn_node,node_process,node_nxt_itr,idx_node)
    cn_nodes = findall(x->x!=0,node_cn_node[:,idx_node])
    for j in 1:length(cn_nodes)
        if !isempty(node_process)
        node_tmp = findall(x->x==cn_nodes[j,1],node_process[:,1])
            if isempty(node_tmp)
                push!(node_nxt_itr,cn_nodes[j,1])
            end
        else
            push!(node_nxt_itr,cn_nodes[j,1])
        end
    end
    return node_nxt_itr,cn_nodes
end
################################################################################
############ Recovery: Power Flow Sending end of transformer branch ############
################################################################################
function power_flow_linear_theta_sdg_trsf_recovery_milp(tap,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,flex_oltc,trsf_fict_nodes,tap_nl,status_oltc)

    line    = [nd_num cnctd_nd]
    id_ft   = line.-ft_line
    id_tf   = line.-tf_line
    line_ft = findall(all(id_ft.==idx_zero,dims=2))
    line_tf = findall(all(id_tf.==idx_zero,dims=2))

    if !isempty(line_ft) line_ft = line_ft[1][1] end
    if !isempty(line_tf) line_tf = line_tf[1][1] end


    if flex_oltc == 1 && status_oltc == 1

        if !isempty(line_ft) && isempty(line_tf)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,1])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,3])
        elseif !isempty(line_tf) && isempty(line_ft)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,3])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,1])
        end
        idx_fict_nd  = intersect(from_nd_trsf,to_nd_trsf)
        idx_fict_nd  = idx_fict_nd[1,1]
        # fict_nd      = trsf_fict_nodes[idx_fict_nd,2]
        fict_nd      = trsf_fict_nodes[idx_fict_nd,4]                            # Index of Fictitious node

        v_fbus_var   = v_mag[s,t,fict_nd]
        v_fbus       = v_mag_int[s,t,fict_nd]

        # v_fbus_var   = v_mag[s,t,idx_nd_nw_buses]
        # v_fbus       = v_mag_int[s,t,idx_nd_nw_buses]
    elseif (flex_oltc == 0 || flex_oltc == 1) && status_oltc == 0

        v_fbus_var   = v_mag[s,t,idx_nd_nw_buses]
        v_fbus       = v_mag_int[s,t,idx_nd_nw_buses]
    end

    v_tbus     = v_mag_int[s,t,idx_cnctd_nd]
    v_tbus_var = v_mag[s,t,idx_cnctd_nd]

    if !isempty(line_ft)&&isempty(line_tf)
        theta_diff = +1*v_ang[s,t,line_ft]
    elseif isempty(line_ft)&& !isempty(line_tf)
        theta_diff = -1*v_ang[s,t,line_tf]
    end

    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]

    c_int = cos(theta_diff_int)
    s_int = sin(theta_diff_int)
    node_i_pconst = tap*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)
    node_j_pconst = tap*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)

    node_i_qconst = tap*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)
    node_j_qconst = tap*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)

    theta_ij_pconst = tap*v_fbus*v_tbus*(s_int*gij_line-c_int*bij_line)
    theta_ij_qconst = tap*v_fbus*v_tbus*(-s_int*bij_line-c_int*gij_line)

    offset_p = tap*gij_line*((-v_fbus*v_tbus*theta_diff_int*s_int)-(c_int*(v_fbus-v_tbus)^2)/2)+tap*bij_line*((v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)
    offset_q = tap*bij_line*((v_fbus*v_tbus*theta_diff_int*s_int)+(c_int*(v_fbus-v_tbus)^2)/2)+tap*gij_line*((+v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)

    pij = +tap^2*(gij_line_sh+gij_line)*v_fbus_var+(node_i_pconst*v_fbus_var)+(node_j_pconst*v_tbus_var)+(theta_ij_pconst*theta_diff)+offset_p
    qij = -tap^2*(bij_line_sh+bij_line)*v_fbus_var+(node_i_qconst*v_fbus_var)+(node_j_qconst*v_tbus_var)+(theta_ij_qconst*theta_diff)+offset_q

    pij_nl = (+tap_nl^2*(gij_line_sh+gij_line)*v_fbus_var)-(tap_nl*gij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*cos(theta_diff))-(tap_nl*bij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*sin(theta_diff))
    qij_nl = (-tap_nl^2*(bij_line_sh+bij_line)*v_fbus_var)-(tap_nl*gij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*sin(theta_diff))+(tap_nl*bij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*cos(theta_diff))

return pij,qij,pij_nl,qij_nl
end
################################################################################
############ Recovery: Power Flow Receiving end of transformer branch ##########
################################################################################
function power_flow_linear_theta_rcv_trsf_recovery_milp(tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,flex_oltc,trsf_fict_nodes,tap_nl,status_oltc)

    line    = [nd_num cnctd_nd]
    id_ft   = line.-ft_line
    id_tf   = line.-tf_line
    line_ft = findall(all(id_ft.==idx_zero,dims=2))
    line_tf = findall(all(id_tf.==idx_zero,dims=2))

    if !isempty(line_ft) line_ft = line_ft[1][1] end
    if !isempty(line_tf) line_tf = line_tf[1][1] end


    if flex_oltc == 1 && status_oltc == 1
        if !isempty(line_ft) && isempty(line_tf)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,1])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,3])
        elseif !isempty(line_tf) && isempty(line_ft)
            from_nd_trsf = findall(x->x==nd_num,trsf_fict_nodes[:,3])
            to_nd_trsf   = findall(x->x==cnctd_nd,trsf_fict_nodes[:,1])
        end
        idx_fict_nd  = intersect(from_nd_trsf,to_nd_trsf)
        idx_fict_nd  = idx_fict_nd[1,1]
        # fict_nd      = trsf_fict_nodes[idx_fict_nd,2]
        fict_nd      = trsf_fict_nodes[idx_fict_nd,4]
        v_tbus_var   = v_mag[s,t,fict_nd]
        v_tbus       = v_mag_int[s,t,fict_nd]

        # v_tbus_var   = v_mag[s,t,idx_cnctd_nd]
        # v_tbus       = v_mag_int[s,t,idx_cnctd_nd]
    elseif (flex_oltc == 0 || flex_oltc == 1) && status_oltc == 0
        v_tbus_var   = v_mag[s,t,idx_cnctd_nd]
        v_tbus       = v_mag_int[s,t,idx_cnctd_nd]
    end

    v_fbus     = v_mag_int[s,t,idx_nd_nw_buses]
    v_fbus_var = v_mag[s,t,idx_nd_nw_buses]

    if !isempty(line_ft)&&isempty(line_tf)
        theta_diff = +1*v_ang[s,t,line_ft]
    elseif isempty(line_ft)&& !isempty(line_tf)
        theta_diff = -1*v_ang[s,t,line_tf]
    end

    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]

    c_int = cos(theta_diff_int)
    s_int = sin(theta_diff_int)
    node_i_pconst = tap2*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)
    node_j_pconst = tap2*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)

    node_i_qconst = tap2*((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)
    node_j_qconst = tap2*((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)

    theta_ij_pconst = tap2*v_fbus*v_tbus*(s_int*gij_line-c_int*bij_line)
    theta_ij_qconst = tap2*v_fbus*v_tbus*(-s_int*bij_line-c_int*gij_line)

    offset_p = tap2*gij_line*((-v_fbus*v_tbus*theta_diff_int*s_int)-(c_int*(v_fbus-v_tbus)^2)/2)+tap2*bij_line*((v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)
    offset_q = tap2*bij_line*((v_fbus*v_tbus*theta_diff_int*s_int)+(c_int*(v_fbus-v_tbus)^2)/2)+tap2*gij_line*((+v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)

    pij = +tap1^2*(gij_line)*v_fbus_var+(node_i_pconst*v_fbus_var)+(node_j_pconst*v_tbus_var)+(theta_ij_pconst*theta_diff)+offset_p
    qij = -tap1^2*(bij_line)*v_fbus_var+(node_i_qconst*v_fbus_var)+(node_j_qconst*v_tbus_var)+(theta_ij_qconst*theta_diff)+offset_q

    pij_nl = (+(gij_line)*v_fbus_var)-(tap_nl*gij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*cos(theta_diff))-(tap_nl*bij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*sin(theta_diff))
    qij_nl = (-(bij_line)*v_fbus_var)-(tap_nl*gij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*sin(theta_diff))+(tap_nl*bij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*cos(theta_diff))

return pij,qij,pij_nl,qij_nl
end
################################################################################
############ ------ Recovery: Power Flow of simple Branch --------- ############
################################################################################
function power_flow_linear_branch_theta_recovery_milp(gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero)
    v_fbus     = v_mag_int[s,t,idx_nd_nw_buses]
    v_tbus     = v_mag_int[s,t,idx_cnctd_nd]

    # v_fbus_var = v_mag[s,t,idx_nd_nw_buses].^2
    # v_tbus_var = v_mag[s,t,idx_cnctd_nd].^2

    v_fbus_var = v_mag[s,t,idx_nd_nw_buses]
    v_tbus_var = v_mag[s,t,idx_cnctd_nd]

    line    = [nd_num cnctd_nd]
    id_ft   = line.-ft_line
    id_tf   = line.-tf_line
    line_ft = findall(all(id_ft.==idx_zero,dims=2))
    line_tf = findall(all(id_tf.==idx_zero,dims=2))

    if !isempty(line_ft) line_ft = line_ft[1][1] end
    if !isempty(line_tf) line_tf = line_tf[1][1] end

    if !isempty(line_ft)&&isempty(line_tf)
        theta_diff = +1*v_ang[s,t,line_ft]
    elseif isempty(line_ft)&& !isempty(line_tf)
        theta_diff = -1*v_ang[s,t,line_tf]
    end

    theta_diff_int = v_ang_int[s,t,idx_nd_nw_buses]-v_ang_int[s,t,idx_cnctd_nd]

    c_int = cos(theta_diff_int)
    s_int = sin(theta_diff_int)
    # s_int = 0
    node_i_pconst = ((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)
    node_j_pconst = ((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(-c_int*gij_line-s_int*bij_line)

    node_i_qconst = ((3*v_tbus-v_fbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)
    node_j_qconst = ((3*v_fbus-v_tbus)/(2*v_fbus+2*v_tbus))*(c_int*bij_line-s_int*gij_line)

    theta_ij_pconst = v_fbus*v_tbus*(s_int*gij_line-c_int*bij_line)
    theta_ij_qconst = v_fbus*v_tbus*(-s_int*bij_line-c_int*gij_line)

    offset_p = gij_line*((-v_fbus*v_tbus*theta_diff_int*s_int)-(c_int*(v_fbus-v_tbus)^2)/2)+bij_line*((v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)
    offset_q = bij_line*((v_fbus*v_tbus*theta_diff_int*s_int)+(c_int*(v_fbus-v_tbus)^2)/2)+gij_line*((+v_fbus*v_tbus*theta_diff_int*c_int)-(s_int*(v_fbus-v_tbus)^2)/2)

    pij = (gij_line_sh/2+gij_line)*v_fbus_var+(node_i_pconst*v_fbus_var)+(node_j_pconst*v_tbus_var)+(theta_ij_pconst*theta_diff)+offset_p
    qij = -(bij_line_sh/2+bij_line)*v_fbus_var+(node_i_qconst*v_fbus_var)+(node_j_qconst*v_tbus_var)+(theta_ij_qconst*theta_diff)+offset_q

    pij_nl = (+(gij_line_sh/2+gij_line)*v_fbus_var)-(gij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*cos(theta_diff))-(bij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*sin(theta_diff))
    qij_nl = (-(bij_line_sh/2+bij_line)*v_fbus_var)-(gij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*sin(theta_diff))+(bij_line*sqrt(v_fbus_var)*sqrt(v_tbus_var)*cos(theta_diff))

return pij,qij,pij_nl,qij_nl
end
################################################################################
################################################################################
function opf_results_variables(nTP,nSc,nBus,nLines,nNcurt_gen,nCurt_gen,nFl,nStr_active,nTrsf,nTrsf_s1,n_tap)

    if flex_oltc==0
        vm_ac          = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
        va_ac          = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
        vm_plr         = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
        va_plr         = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
    elseif flex_oltc==1
        # vm_ac          = zeros(Float64,(nBus+nTrsf),length(1:nTP),(length(1:nSc)))
        # va_ac          = zeros(Float64,(nBus+nTrsf),length(1:nTP),(length(1:nSc)))
        # vm_plr         = zeros(Float64,(nBus+nTrsf),length(1:nTP),(length(1:nSc)))
        # va_plr         = zeros(Float64,(nBus+nTrsf),length(1:nTP),(length(1:nSc)))

        vm_ac          = zeros(Float64,(nBus+nTrsf_s1),length(1:nTP),(length(1:nSc)))
        va_ac          = zeros(Float64,(nBus+nTrsf_s1),length(1:nTP),(length(1:nSc)))
        vm_plr         = zeros(Float64,(nBus+nTrsf_s1),length(1:nTP),(length(1:nSc)))
        va_plr         = zeros(Float64,(nBus+nTrsf_s1),length(1:nTP),(length(1:nSc)))
    end
    vm_min         = zeros(Float64,nSc,nTP)
    vm_max         = zeros(Float64,nSc,nTP)
    gen_Ncurt_P    = zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))
    gen_Ncurt_Q    = zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))
    gen_Curt_P_sc  = zeros(Float64,(nCurt_gen,length(1:nTP),length(1:nSc)))
    gen_Curt_P_gen = zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    total_p_curt   = zeros(nSc,1)
    fl_inc_sc      = zeros(Float64,(nFl,length(1:nTP),length(1:nSc)))
    fl_dec_sc      = zeros(Float64,(nFl,length(1:nTP),length(1:nSc)))
    fl_inc_nd      = zeros(Float64,(length(1:nSc),length(1:nTP),nFl))
    fl_dec_nd      = zeros(Float64,(length(1:nSc),length(1:nTP),nFl))
    p_ch_nd        = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    p_dis_nd       = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    p_dis_nd       = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    soc_nd         = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    p_strg_nd      = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    trsf_ratio     = zeros(Float64,(length(1:nTP),length(1:nSc),nTrsf_s1))
    dgpf           = zeros(Float64,(length(1:nTP),length(1:nSc),nCurt_gen))
    gen_Curt_Q     = zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    dg_pf_value    = zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    p_inj          = zeros((length(1:nSc),length(1:nTP),nBus))                             # Active power injection at each node
    q_inj          = zeros((length(1:nSc),length(1:nTP),nBus))
    p_inj_nl       = zeros((length(1:nSc),length(1:nTP),nBus))                             # Active power injection at each node
    q_inj_nl       = zeros((length(1:nSc),length(1:nTP),nBus))
    p_imp_td       = zeros(nTP,1)
    q_imp_td       = zeros(nTP,1)
    str_bv         = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    fl_bv          = zeros(Float64,(length(1:nSc),length(1:nTP),nFl))

    # oltc_bv        = zeros(Float64,(length(1:nSc),length(1:nTP),nTrsf_s1))
    # oltc_aux_bv    = zeros(Float64,(length(1:nSc),length(1:nTP),nTrsf_s1))
    # ratio_sql      = zeros(Float64,(length(1:nSc),length(1:nTP),nTrsf_s1))

    oltc_bv        = zeros(Float64,(length(1:nSc),length(1:nTP),n_tap,nTrsf_s1))
    oltc_aux_bv    = zeros(Float64,(length(1:nSc),length(1:nTP),n_tap,nTrsf_s1))
    ratio_sql      = zeros(Float64,(length(1:nSc),length(1:nTP),n_tap,nTrsf_s1))

    avg_gen_curt   = zeros(nCurt_gen,nTP)
    avg_str_chr    = zeros(nStr_active,nTP)
    avg_str_dch    = zeros(nStr_active,nTP)
    avg_str_soc    = zeros(nStr_active,nTP)
    avg_fl_od      = zeros(nFl,nTP)
    avg_fl_ud      = zeros(nFl,nTP)
    avg_tratio     = zeros(nTrsf,nTP)

    br_crnt_pu     = zeros((length(1:nSc),length(1:nTP),nLines))
    br_crnt_si     = zeros((length(1:nSc),length(1:nTP),nLines))
    br_pwr_ij_pu   = zeros((length(1:nSc),length(1:nTP),nLines))
    br_pwr_ij_SI   = zeros((length(1:nSc),length(1:nTP),nLines))
    br_crnt_ij_pu_ex  = zeros((length(1:nSc),length(1:nTP),nLines))
    br_crnt_ij_SI_ex  = zeros((length(1:nSc),length(1:nTP),nLines))
    br_pwr_ij_pu_ex   = zeros((length(1:nSc),length(1:nTP),nLines))
    br_pwr_ij_SI_ex   = zeros((length(1:nSc),length(1:nTP),nLines))

    line_alpha     = zeros((length(1:nSc),length(1:nTP),nLines))
    line_eta       = zeros((length(1:nSc),length(1:nTP),nLines))
    app_err        = zeros((length(1:nSc),length(1:nTP),nLines))
    line_beta      = zeros((length(1:nSc),length(1:nTP),nLines))
    return vm_ac,va_ac,vm_plr,va_plr,vm_min,vm_max,gen_Ncurt_P,gen_Ncurt_Q,gen_Curt_P_sc,gen_Curt_P_gen,
    total_p_curt,fl_inc_sc,fl_dec_sc,fl_inc_nd,fl_dec_nd,p_ch_nd,p_dis_nd,soc_nd,p_strg_nd,trsf_ratio,dgpf,
    gen_Curt_Q,br_crnt_pu,br_crnt_si,br_pwr_ij_pu,br_pwr_ij_SI,br_crnt_ij_pu_ex,br_crnt_ij_SI_ex,br_pwr_ij_pu_ex,br_pwr_ij_SI_ex,app_err,
    p_inj,q_inj,p_imp_td,q_imp_td,str_bv,fl_bv,oltc_bv,oltc_aux_bv,ratio_sql,
    p_inj_nl,q_inj_nl,avg_gen_curt,avg_str_chr,avg_str_dch,avg_str_soc,avg_fl_od,avg_fl_ud,avg_tratio,line_alpha,line_eta,line_beta,dg_pf_value
end
################################################################################
############## Calculation of Apparent Power Injection Error ###################
################################################################################
function s_inj_power_max_power(p_inj,q_inj,p_inj_nl,q_inj_nl,max_err_per_bus,norm_max_err_per_bus,p_err_max_power,q_err_max_power,sql_itr,lin_itr)
    p_inj_avg    = reshape(sum(sum(p_inj,dims=1),dims=2),(size(p_inj,3),1))
    q_inj_avg    = reshape(sum(sum(q_inj,dims=1),dims=2),(size(q_inj,3),1))
    p_inj_nl_avg = reshape(sum(sum(p_inj_nl,dims=1),dims=2),(size(p_inj_nl,3),1))
    q_inj_nl_avg = reshape(sum(sum(q_inj_nl,dims=1),dims=2),(size(q_inj_nl,3),1))

    p_inj_avg = p_inj_avg/(nTP*nSc)       # power in pu
    q_inj_avg = q_inj_avg/(nTP*nSc)       # power in pu
    p_inj_nl_avg = p_inj_nl_avg/(nTP*nSc) # power in pu
    q_inj_nl_avg = q_inj_nl_avg/(nTP*nSc) # power in pu

    s_inj_avg    = sqrt.((p_inj_avg).^2+(q_inj_avg).^2)
    s_inj_nl_avg = sqrt.((p_inj_nl_avg).^2+(q_inj_nl_avg).^2)

    p_err_max_power[:,sql_itr,lin_itr] = (abs.(p_inj_nl_avg.-p_inj_avg)./maximum(abs.(p_inj_nl_avg)))*100 # % in pu
    q_err_max_power[:,sql_itr,lin_itr] = (abs.(q_inj_nl_avg.-q_inj_avg)./maximum(abs.(q_inj_nl_avg)))*100 # % in pu

    max_err_per_bus[:,sql_itr,lin_itr]       = abs.(s_inj_avg.-s_inj_nl_avg)*100 # % in pu
    norm_max_err_per_bus[:,sql_itr,lin_itr]  = abs.(s_inj_avg.-s_inj_nl_avg)./maximum(s_inj_nl_avg)*100 # % in pu
    return max_err_per_bus,norm_max_err_per_bus,p_err_max_power,q_err_max_power
end

function s_inj_power_per_bus(p_inj,q_inj,p_inj_nl,q_inj_nl,avg_err_per_bus,norm_avg_err_per_bus,p_err_avg,q_err_avg,sql_itr,lin_itr)
    # p_inj_avg    = sum(p_inj)
    # q_inj_avg    = sum(q_inj)
    # p_inj_nl_avg = sum(p_inj_nl)
    # q_inj_nl_avg = sum(q_inj_nl)

    p_inj_avg    = reshape(sum(sum(p_inj,dims=1),dims=2),(size(p_inj,3),1))
    q_inj_avg    = reshape(sum(sum(q_inj,dims=1),dims=2),(size(q_inj,3),1))
    p_inj_nl_avg = reshape(sum(sum(p_inj_nl,dims=1),dims=2),(size(p_inj_nl,3),1))
    q_inj_nl_avg = reshape(sum(sum(q_inj_nl,dims=1),dims=2),(size(q_inj_nl,3),1))

    p_inj_avg = p_inj_avg/(nTP*nSc)          # power in pu
    q_inj_avg = q_inj_avg/(nTP*nSc)          # power in pu
    p_inj_nl_avg = p_inj_nl_avg/(nTP*nSc)    # power in pu
    q_inj_nl_avg = q_inj_nl_avg/(nTP*nSc)    # power in pu

    s_inj_avg    = sqrt.((p_inj_avg).^2+(q_inj_avg).^2)
    s_inj_nl_avg = sqrt.((p_inj_nl_avg).^2+(q_inj_nl_avg).^2)

    p_err_avg[:,sql_itr,lin_itr] = (abs.(p_inj_nl_avg.-p_inj_avg)/size(p_inj,3))*100  # % in pu
    q_err_avg[:,sql_itr,lin_itr] = (abs.(q_inj_nl_avg.-q_inj_avg)/size(p_inj,3))*100  # % in pu

    avg_err_per_bus[:,sql_itr,lin_itr]       = abs.(s_inj_avg-s_inj_nl_avg)*100       # % in pu
    norm_avg_err_per_bus[:,sql_itr,lin_itr]  = abs.(s_inj_avg-s_inj_nl_avg)/size(p_inj,3)*100 # % in pu


    return avg_err_per_bus,norm_avg_err_per_bus,p_err_avg,q_err_avg
end

function br_crnt_aux(Iij,Iij_ex,error,i)
    err = abs(error)
    if abs(error)<=1e-3 && Iij<0 # App current is negative
        # println("branch $i app current is LESS than 0.001 but is NEGATIVE")
        Iij = abs(Iij_ex)
    elseif abs(error)<=1e-3 && Iij >=0 # App current is Positive
        # println("branch $i app current is LESS than 0.001 and is POSITIVE")
        # do nothing
    elseif abs(error)>=1e-3 && Iij >=0 #
        # println("branch $i app current is GREATOR than 0.001 and is POSITIVE" )
    elseif abs(error)>=1e-3 && Iij <=0
        # println("WARNING: branch $i app current is GREATOR than 0.001 and is NEGATIVE")
        Iij = abs(Iij_ex)
    end
        return br_crnt  = sqrt(Iij[1,1])
end

function br_crnt_aux(Iij_ex,i)
    if Iij_ex<0
        # println("Exact branch current for line $i is NEGATIVE and has the value $Iij_ex")
        br_crnt  = sqrt(abs(Iij_ex))
    else
        br_crnt = sqrt(Iij_ex)
    end
    return br_crnt
end

# Constraint violation check after MILP OPF
function constraint_violation(nBus,nw_buses,nLines,nw_lines,nTP,nSc,vm_ac_pf,br_crnt,br_crnt_rat,lin_itr,min_vol_limit,max_vol_limit,max_crnt_limit,vol_viol_tol,crnt_viol_tol,vol_cnst_fctr,crnt_cnst_fctr,lin_itr_max)
##------------------------------------------------------------------------------
##################### Voltage Violation Check ##################################
##------------------------------------------------------------------------------
vol_nodes     = vm_ac_pf
vol_viol      = []
crnt_viol     = []
max_vol_viol  = []
max_crnt_viol = []
viol_nodes    = []
vol_viol_cn   = []
viol_num      = 0
nVol_viol     = 0
avg_vol_viol  = 0
above_avg_vol_viol = 0
below_avg_vol_viol = 0

for i in 1:nBus
    nd = nw_buses[i].bus_num
    vol_min = nw_buses[i].bus_vmin
    vol_max = nw_buses[i].bus_vmax
    for t in 1:nTP
        for s in 1:nSc
            vol = vm_ac_pf[s,t,i]
            if vol>vol_max || vol<vol_min
                if vol>vol_max
                    viol = (vol-vol_max)*100
                    if viol >=(vol_viol_tol) && lin_itr < lin_itr_max
                        # max_vol_limit[s,t,i,lin_itr+1] = max_vol_limit[s,t,i,lin_itr]*(1-((vol_max-max_vol_limit[s,t,i,lin_itr])+vol_cnst_fctr))
                        max_vol_limit[s,t,i,lin_itr+1] = vol_max-((vol_max-max_vol_limit[s,t,i,lin_itr])+vol_cnst_fctr)
                        # max_vol_limit[s,t,i,lin_itr+1] = max_vol_limit[s,t,i,lin_itr+1]*(1-vol_cnst_fctr*lin_itr)
                    elseif  viol <(vol_viol_tol) && lin_itr < lin_itr_max
                        max_vol_limit[s,t,i,lin_itr+1] = max_vol_limit[s,t,i,lin_itr]
                    end
                elseif vol<vol_min
                    viol = (vol_min-vol)*100
                    if viol >=(vol_viol_tol) && lin_itr < lin_itr_max
                        # min_vol_limit[s,t,i,lin_itr+1] = min_vol_limit[s,t,i,lin_itr+1]*(1+vol_cnst_fctr*lin_itr)
                        # min_vol_limit[s,t,i,lin_itr+1] = min_vol_limit[s,t,i,lin_itr]*(1-((min_vol_limit[s,t,i,lin_itr]-vol_min)+vol_cnst_fctr))
                        min_vol_limit[s,t,i,lin_itr+1] = vol_min+((min_vol_limit[s,t,i,lin_itr]-vol_min)+vol_cnst_fctr)
                    elseif viol <(vol_viol_tol) && lin_itr < lin_itr_max
                        min_vol_limit[s,t,i,lin_itr+1] = min_vol_limit[s,t,i,lin_itr]
                    end
                end
                viol_num = viol_num+1
                push!(max_vol_viol,viol)
                push!(vol_viol,("Itr:$lin_itr","viol_num:$viol_num",s,t,nd,vol,viol))
                push!(viol_nodes,nd)
                if viol >= 1
                    push!(vol_viol_cn,("Itr:$lin_itr","viol_num:$viol_num",s,t,nd,vol,viol))
                end
            else
                if lin_itr < lin_itr_max
                    max_vol_limit[s,t,i,lin_itr+1] = max_vol_limit[s,t,i,lin_itr]
                    min_vol_limit[s,t,i,lin_itr+1] = min_vol_limit[s,t,i,lin_itr]
                end
            end
        end
    end
end

nVol_viol    = viol_num
if !isempty(max_vol_viol)
    avg_vol_viol = sum(max_vol_viol)/size(max_vol_viol,1)
    above_avg_vol_viol = size(findall(x->x>avg_vol_viol,max_vol_viol),1)
    below_avg_vol_viol = size(findall(x->x<avg_vol_viol,max_vol_viol),1)
end

##------------------------------------------------------------------------------
####################### Current Violation Check ################################
##------------------------------------------------------------------------------
crnt_branch = br_crnt
viol_num = 0
nCrnt_viol     = 0
avg_crnt_viol  = 0
above_avg_crnt_viol = 0
below_avg_crnt_viol = 0

for i in 1:nLines
    line_from = nw_lines[i].line_from
    line_to   = nw_lines[i].line_to
    i_rat     = br_crnt_rat[i]
    for t in 1:nTP
        for s in 1:nSc
            i_line = br_crnt[s,t,i]
            if i_line>i_rat
                viol_num = viol_num+1
                prct_loading = ((i_line-i_rat)/i_rat)*100
                if prct_loading >= crnt_viol_tol
                    max_crnt_limit[s,t,i] = max_crnt_limit[s,t,i]*(1-crnt_cnst_fctr*(crnt_viol_tol / 100))
                else
                    max_crnt_limit[s,t,i] = max_crnt_limit[s,t,i]
                end

                push!(max_crnt_viol,prct_loading)
                push!(crnt_viol,("Itr:$lin_itr","viol_num:$viol_num",s,t,line_from,line_to,i_line,i_rat,prct_loading))
            else
                # do nothing
            end
        end
    end
end
nCrnt_viol    = viol_num
if !isempty(max_crnt_viol)
    avg_crnt_viol = sum(max_crnt_viol)/size(max_crnt_viol,1)
    above_avg_crnt_viol = size(findall(x->x>avg_crnt_viol,max_crnt_viol),1)
    below_avg_crnt_viol = size(findall(x->x<avg_crnt_viol,max_crnt_viol),1)
end


return vol_viol,crnt_viol,max_vol_viol,max_crnt_viol,viol_nodes,nVol_viol,avg_vol_viol,above_avg_vol_viol,below_avg_vol_viol,
nCrnt_viol,avg_crnt_viol,above_avg_crnt_viol,below_avg_crnt_viol,min_vol_limit,max_vol_limit,max_crnt_limit,vol_viol_cn
end

# Constraint violation check in Initial AC-PF
function constraint_violation_pf(nBus,nw_buses,nLines,nw_lines,nTP,nSc,vm_ac_pf,br_crnt,br_crnt_rat,lin_itr)
##------------------------------------------------------------------------------
##################### Voltage Violation Check ##################################
##------------------------------------------------------------------------------
vol_nodes     = vm_ac_pf
vol_viol      = []
vol_viol_cn   = []
crnt_viol     = []
max_vol_viol  = []
max_crnt_viol = []
viol_nodes    = []
viol_num      = 0
nVol_viol     = 0
avg_vol_viol  = 0
above_avg_vol_viol = 0
below_avg_vol_viol = 0

for i in 1:nBus
    nd = nw_buses[i].bus_num
    vol_min = nw_buses[i].bus_vmin
    vol_max = nw_buses[i].bus_vmax
    for t in 1:nTP
        for s in 1:nSc
            vol = vm_ac_pf[s,t,i]
            if vol>vol_max || vol<vol_min
                if vol>vol_max
                    viol = (vol-vol_max)#*100 # the *100 was deactivated by baraa
                elseif vol<vol_min
                    viol = (vol_min-vol)#*100 # the *100 was deactivated by baraa
                end
                viol_num = viol_num+1
                push!(max_vol_viol,viol)
                push!(vol_viol,("Itr:$lin_itr","viol_num:$viol_num",s,t,nd,vol,viol))
                push!(viol_nodes,nd)
                if viol >= 1
                    push!(vol_viol_cn,("Itr:$lin_itr","viol_num:$viol_num",s,t,nd,vol,viol))
                end
            else
                # do nothing
            end
        end
    end
end

nVol_viol    = viol_num
if !isempty(max_vol_viol)
    avg_vol_viol = sum(max_vol_viol)/size(max_vol_viol,1)
    above_avg_vol_viol = size(findall(x->x>avg_vol_viol,max_vol_viol),1)
    below_avg_vol_viol = size(findall(x->x<avg_vol_viol,max_vol_viol),1)
end

##------------------------------------------------------------------------------
####################### Current Violation Check ################################
##------------------------------------------------------------------------------
crnt_branch = br_crnt
viol_num = 0
nCrnt_viol     = 0
avg_crnt_viol  = 0
above_avg_crnt_viol = 0
below_avg_crnt_viol = 0

for i in 1:nLines
    line_from = nw_lines[i].line_from
    line_to   = nw_lines[i].line_to
    i_rat     = br_crnt_rat[i]
    for t in 1:nTP
        for s in 1:nSc
            i_line = br_crnt[s,t,i]
            if  i_line>i_rat
                viol_num = viol_num+1
                prct_loading = ((i_line-i_rat)/i_rat)*100
                push!(max_crnt_viol,prct_loading)
                push!(crnt_viol,("Itr:$lin_itr","viol_num:$viol_num",s,t,line_from,line_to,i_line,i_rat,prct_loading))
            else
                # do nothing
            end
        end
    end
end
nCrnt_viol    = viol_num
if !isempty(max_crnt_viol)
    avg_crnt_viol = sum(max_crnt_viol)/size(max_crnt_viol,1)
    above_avg_crnt_viol = size(findall(x->x>avg_crnt_viol,max_crnt_viol),1)
    below_avg_crnt_viol = size(findall(x->x<avg_crnt_viol,max_crnt_viol),1)
end


return vol_viol,crnt_viol,max_vol_viol,max_crnt_viol,viol_nodes,nVol_viol,avg_vol_viol,above_avg_vol_viol,below_avg_vol_viol,
nCrnt_viol,avg_crnt_viol,above_avg_crnt_viol,below_avg_crnt_viol,vol_viol_cn
end
