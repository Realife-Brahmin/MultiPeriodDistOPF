############################# Shunt Status ####################################
function check_shunt_acpf(shunt_power,nd_num,nw_buses)

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
function nodal_components_values_acpf(var,s,t,idx)
    output = var[s,t,idx]
    if isempty(output) output = 0.0 end
    return output
end

function nodal_components_values_acpf(var1,var2,var3,s,t,idx1,idx2,var4)
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

function power_flow_line_acpf(model,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
    pij = @NLexpression(model,+tap1^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap2*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
    qij = @NLexpression(model,-tap1^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap2*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
return pij, qij
end

function power_flow_line_acpf(model,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
    pij = @NLexpression(model,+tap^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
    qij = @NLexpression(model,-tap^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]))
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

function power_flow_line_trsf_acpf(model,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,v_r,v_im,v_sq)
    pij = @NLexpression(model,+tap^2*(gij_line_sh+gij_line)*v_sq - tap*gij_line*v_r - tap*bij_line*v_im)
    qij = @NLexpression(model,-tap^2*(bij_line_sh+bij_line)*v_sq + tap*bij_line*v_r - tap*gij_line*v_im)
return pij, qij
end

function power_flow_line_trsf_acpf(model,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,v_r,v_im,v_sq)
    pij = @NLexpression(model,+tap1^2*(gij_line)*v_sq - tap2*gij_line*v_r - tap2*bij_line*v_im)
    qij = @NLexpression(model,-tap1^2*(bij_line)*v_sq + tap2*bij_line*v_r - tap2*gij_line*v_im)
return pij, qij
end
######################## Power Flow Line recovery ##############################
function power_flow_line_recovery_acpf(tap,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
    pij = +tap^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
    qij = -tap^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
return pij, qij
end


function power_flow_line_recovery_acpf(tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,e,f,s,t,idx_nd_nw_buses,idx_cnctd_nd)
    pij = +tap1^2*(gij_line_sh/2+gij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) - tap2*gij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*bij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
    qij = -tap1^2*(bij_line_sh/2+bij_line)*(e[s,t,idx_nd_nw_buses]^2+f[s,t,idx_nd_nw_buses]^2) + tap2*bij_line*(e[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]+f[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd]) - tap2*gij_line*(f[s,t,idx_nd_nw_buses]*e[s,t,idx_cnctd_nd]-e[s,t,idx_nd_nw_buses]*f[s,t,idx_cnctd_nd])
return pij, qij
end

################################################################################
########### Recovering Generator Injection Values (Non-Curtailable) ############
################################################################################
function recovery_ncurt_gen_pf(nTP,nSc,nNcurt_gen,Pgen_nCurt_pu,Qgen_nCurt_pu,pgen_tol,qgen_tol)
    gen_Ncurt_P = zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))
    gen_Ncurt_Q = zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))
    for i in 1:size(gen_Ncurt_P,3)
        for s in 1:nSc
            for t in 1:nTP
                if abs.(Pgen_nCurt_pu[s,t,i]) < pgen_tol
                    gen_Ncurt_P[s,t,i] = 0.0
                else
                    gen_Ncurt_P[s,t,i] = Pgen_nCurt_pu[s,t,i]
                end
                if abs.(Qgen_nCurt_pu[s,t,i]) < qgen_tol
                    gen_Ncurt_Q[s,t,i] = 0.0
                else
                    gen_Ncurt_Q[s,t,i] = Qgen_nCurt_pu[s,t,i]
                end
            end
        end
    end
    return gen_Ncurt_P,gen_Ncurt_Q
end
################################################################################
######################### Recovering Voltages ##################################
################################################################################
function recovery_voltages_pf(nTP,nSc,nBus,vr,vim)
    v_rect = vr+(vim)im
    v_mag  = abs.(v_rect)
    v_tht  = angle.(v_rect)   # angle(z) determines the phase angle in radians of a complex number z
    vm_ac  = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
    va_ac  = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
    min_vm = zeros(Float64,nSc,nTP)
    max_vm = zeros(Float64,nSc,nTP)
    for s in 1:nSc
        for i in 1:nBus
            vm_ac[i,:,s] = v_mag[s,:,i]
            va_ac[i,:,s] = v_tht[s,:,i]
        end
        for t in 1:nTP
            min_vm[s,t] = minimum(vm_ac[:,t,s])
            max_vm[s,t] = maximum(vm_ac[:,t,s])
        end
    end
    return vm_ac,va_ac,v_rect,min_vm,max_vm
end
