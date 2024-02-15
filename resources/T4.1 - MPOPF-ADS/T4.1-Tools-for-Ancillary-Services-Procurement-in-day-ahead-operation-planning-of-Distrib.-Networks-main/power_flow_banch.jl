function power_flow_branch(acopf,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero)

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
        # theta_diff = 1*v_ang[s,t,line_ft]
    elseif isempty(line_ft)&& !isempty(line_tf)
        theta_diff = @expression(acopf,-1*v_ang[s,t,line_tf])
        # theta_diff = -1*v_ang[s,t,line_tf]
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

    pij = @expression(acopf,(gij_line_sh/2+gij_line)*v_fbus_var+(alpha_i_P*v_fbus_var)+(alpha_j_P*v_tbus_var)+(beta_ij_P*theta_diff)+gamma_P)
    qij = @expression(acopf,(bij_line_sh/2+bij_line)*v_fbus_var+(alpha_i_Q*v_fbus_var)+(alpha_j_Q*v_tbus_var)+(beta_ij_Q*theta_diff)+gamma_Q)

    # pij = @expression(acopf,(gij_line_sh/2+gij_line)*v_fbus_var+(alpha_i_P*v_fbus_var)+(alpha_j_P*v_tbus_var))
    # qij = @expression(acopf,(bij_line_sh/2+bij_line)*v_fbus_var+(alpha_i_Q*v_fbus_var)+(alpha_j_Q*v_tbus_var))

return pij, qij
end
