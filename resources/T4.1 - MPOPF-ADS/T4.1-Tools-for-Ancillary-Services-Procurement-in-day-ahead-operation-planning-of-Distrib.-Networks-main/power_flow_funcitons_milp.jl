############################# Shunt Status ####################################
function check_shunt_milp(shunt_power,nd_num,nw_buses)

    if size(shunt_power,1)==1
        shunt_power = shunt_power[1,1]
        shunt_elem  = shunt_power/(nw_buses[nd_num].bus_vmag)^2                       # The shunt active power injected/abosorbed is based on the V_Mag (equal to1 )value provided in the MATPOWER.
    else
        shunt_elem  = shunt_power[:,1]/(nw_buses[nd_num].bus_vmag)^2
    end

    if isempty(shunt_elem) shunt_elem = 0.0 end
    return shunt_elem
end
######################### Nodal Components ####################################
function nodal_components_values_milp(var,s,t,idx)
    output = var[s,t,idx]
    if isempty(output) output = 0.0 end
    return output
end

function nodal_components_values_milp(var1,var2,var3,s,t,idx1,idx2,var4)
    if isempty(idx1) && ~isempty(idx2)
        out1 = var1[s,t,idx2]
        out2 = var4[s,t,idx2]
    elseif  ~isempty(idx1) && isempty(idx2)
        out1 = var2[s,t,idx1]     #Pg = Pg_nCurt (Non-curtailable DG)
        out2 = var3[s,t,idx1]
    elseif  isempty(idx1) && isempty(idx2)
        out1 = 0.0
        out2 = 0.0
    end
return out1, out2
end
############################## Power Flow Line #################################
# function power_flow_line_minlp(model,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
#     pij = @NLexpression(model,+tap^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
#     qij = @NLexpression(model,-tap^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
# return pij, qij
# end

function power_flow_line_minlp(model,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
    pij = @NLexpression(model,+tap1^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap2*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
    qij = @NLexpression(model,-tap1^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap2*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
return pij, qij
end

function power_flow_line_minlp(model,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,v_r,v_im,v_sq)
    pij = @NLexpression(model,+tap^2*(gij_line_sh/2+gij_line)*v_sq - tap*gij_line*v_r - tap*bij_line*v_im)
    qij = @NLexpression(model,-tap^2*(bij_line_sh/2+bij_line)*v_sq + tap*bij_line*v_r - tap*gij_line*v_im)
return pij, qij
end


######################## Power Flow Transformer ##############################
# function power_flow_line_trsf_minlp(model,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
#     pij = @NLexpression(model,+tap^2*(gij_line_sh+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
#     qij = @NLexpression(model,-tap^2*(bij_line_sh+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
# return pij, qij
# end
#
# function power_flow_line_trsf_minlp(model,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
#     pij = @NLexpression(model,+tap1^2*(gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap2*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
#     qij = @NLexpression(model,-tap1^2*(bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap2*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
# return pij, qij
# end

function power_flow_line_trsf_minlp(model,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,v_r,v_im,v_sq)
    pij = @NLexpression(model,+tap^2*(gij_line_sh+gij_line)*v_sq - tap*gij_line*v_r - tap*bij_line*v_im)
    qij = @NLexpression(model,-tap^2*(bij_line_sh+bij_line)*v_sq + tap*bij_line*v_r - tap*gij_line*v_im)
return pij, qij
end

function power_flow_line_trsf_minlp(model,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,v_r,v_im,v_sq)
    pij = @NLexpression(model,+tap1^2*(gij_line)*v_sq - tap2*gij_line*v_r - tap2*bij_line*v_im)
    qij = @NLexpression(model,-tap1^2*(bij_line)*v_sq + tap2*bij_line*v_r - tap2*gij_line*v_im)
return pij, qij
end
######################## Power Flow Line recovery ##############################
function power_flow_line_recovery_minlp(tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
    pij = +tap^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
    qij = -tap^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
return pij, qij
end


function power_flow_line_recovery_minlp(tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
    pij = +tap1^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap2*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
    qij = -tap1^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap2*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
return pij, qij
end
