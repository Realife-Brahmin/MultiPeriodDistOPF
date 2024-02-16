#################### Recovering Branch Current Values AC-OPF ###################
function recovery_branch_current(nSc,nTP,nLines,nw_lines,rdata_buses,vm,v_mag_int,va,v_ang_int,yij_line,dLines,Ibase,br_crnt_ij_pu,br_crnt_ij_SI,oltc_ratio,br_pwr_ij_pu,br_pwr_ij_SI,br_crnt_ij_pu_ex,br_crnt_ij_SI_ex,br_pwr_ij_pu_ex,br_pwr_ij_SI_ex,app_err)

    # br_crnt_ij_pu  = zeros((length(1:nSc),length(1:nTP),nLines))
    # br_crnt_ij_SI  = zeros((length(1:nSc),length(1:nTP),nLines))
    vm = vm.^2
    for i in 1:nLines
        gij_line  = real(yij_line[i,1])
        bij_line  = imag(yij_line[i,1])
        tap       = nw_lines[i].line_tap_ratio
        idx_tap   = dLines[i,end]
        f_bus     = nw_lines[i].line_from
        t_bus     = nw_lines[i].line_to
        idx_f_bus = findall(x->x==f_bus,rdata_buses[:,1])
        idx_t_bus = findall(x->x==t_bus,rdata_buses[:,1])
        idx_f_bus = idx_f_bus[1,1]
        idx_t_bus = idx_t_bus[1,1]
        idx_line  = i
        for s in 1:nSc
            for t in 1:nTP
                theta_diff_int = v_ang_int[s,t,idx_f_bus]-v_ang_int[s,t,idx_t_bus]
                theta_ij       = va[s,t,idx_f_bus]-va[s,t,idx_t_bus]
                c_int = cos(theta_diff_int)
                vfrom = vm[s,t,idx_f_bus]   # Vm is the square of the voltage
                vto   = vm[s,t,idx_t_bus]

                vfrom_si = (sqrt.(vfrom)).*vbase[idx_f_bus]*1000                   # SI voltage is in Volts
                vto_si   = (sqrt.(vto)).*vbase[idx_t_bus]*1000                     # SI voltage is in Volts

                vfrom_int = v_mag_int[s,t,idx_f_bus]   # Vm is the square of the voltage
                vto_int   = v_mag_int[s,t,idx_t_bus]

                coeff_vf = (((1+c_int)*vfrom_int)+((1-3*c_int)*vto_int))/(vfrom_int+vto_int)
                coeff_vt = (((1+c_int)*vto_int)+((1-3*c_int)*vfrom_int))/(vfrom_int+vto_int)
                cons_ij  = c_int*(vfrom_int-vto_int)^2

                if tap!=0.0
# This is to compare the longitiudinal current constraint with the Imax constraint
                    # t_tap = oltc_ratio[s,t,i]
                    t_tap = tap
####--------------------- Linear constraint ------------------------------######
                # Constraint using alpha-beta

                    # Iij    = (gij_line^2+bij_line^2)*(((t_tap^2*vfrom)+vto)-((2*t_tap)*(c_int*(line_alpha[s,t,i]-line_beta[s,t,i])))) # Iij is the square of the longitudinal current
                    # Iij_ex = (gij_line^2+bij_line^2)*(((t_tap^2*vfrom)+vto)-((2*t_tap)*(sqrt(vfrom)*sqrt(vto)*cos(theta_ij))))        # Iij is the square of the longitudinal current

                    # Cosntraint # 01
                    # Iij    = (gij_line^2+bij_line^2)*(((t_tap^2+c_int*t_tap)*vfrom)+((1+c_int*t_tap)*vto)-(c_int*t_tap*((vfrom_int+vto_int).^2)))
                    # Iij_ex = (gij_line^2+bij_line^2)*(((t_tap^2*vfrom)+vto)-((2*t_tap)*(sqrt(vfrom)*sqrt(vto)*cos(theta_ij))))      # Iij is the square of the longitudinal current

                    # Cosntraint # 02
                    # Iij    = (gij_line^2+bij_line^2)*(((t_tap^2-c_int*t_tap)*vfrom)+((1-c_int*t_tap)*vto)+(c_int*t_tap*((vfrom_int-vto_int).^2)))
                    # Iij_ex = (gij_line^2+bij_line^2)*(((t_tap^2*vfrom)+vto)-((2*t_tap)*(sqrt(vfrom)*sqrt(vto)*cos(theta_ij))))      # Iij is the square of the longitudinal current

                    # Constraint by putting the value of (V_i-V_j)^2
                    Iij = (gij_line^2+bij_line^2)*((coeff_vf*vfrom)+(coeff_vt*vto)-cons_ij)
                    Iij_ex = (gij_line^2+bij_line^2)*(((t_tap^2*vfrom)+vto)-((2*t_tap)*(sqrt(vfrom)*sqrt(vto)*cos(theta_ij))))        # Iij is the square of the longitudinal current

                    # app_err[s,t,i] = 2*(gij_line^2+bij_line^2)*(((sqrt(vfrom)*sqrt(vto)*cos(theta_ij))-(c_int*(line_alpha[s,t,i]-line_beta[s,t,i]))))
                    # app_err[s,t,i] = 2*(gij_line^2+bij_line^2)*((sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij))-(c_int*((vfrom+vto)/2)-(((vfrom_int-vto_int).^2)/2)))

                    # br_crnt    = br_crnt_aux(Iij,Iij_ex,app_err[s,t,i],i)
                    br_crnt_ex = br_crnt_aux(Iij_ex,i)

                    br_crnt = (Iij)
                    # br_crnt_ex = sqrt.(Iij_ex)
                    br_crnt_ij_pu[s,t,i] = br_crnt
                    br_crnt_ij_SI[s,t,i] = br_crnt*Ibase[idx_t_bus,1]
                    br_pwr_ij_pu[s,t,i]  = (sqrt(vto))*br_crnt
                    br_pwr_ij_SI[s,t,i]  = (vto_si*(br_crnt*Ibase[idx_t_bus,1]))/1e6

                    br_crnt_ij_pu_ex[s,t,i] = br_crnt_ex
                    br_crnt_ij_SI_ex[s,t,i] = br_crnt_ex*Ibase[idx_t_bus,1]
                    br_pwr_ij_pu_ex[s,t,i]  = (sqrt(vto))*br_crnt_ex
                    br_pwr_ij_SI_ex[s,t,i]  = (vto_si*(br_crnt_ex*Ibase[idx_t_bus,1]))/1e6

                else
                    # Iij    = (gij_line^2+bij_line^2)*((vfrom+vto)-(2*c_int*(line_alpha[s,t,i]-line_beta[s,t,i])))
                    # Iij_ex = (gij_line^2+bij_line^2)*((vfrom+vto)-(2*sqrt(vfrom)*sqrt(vto)*cos(theta_ij)))      # Iij is the square of the longitudinal current

                    # Constraint # 01
                    # Iij    = (gij_line^2+bij_line^2)*(((1+c_int)*(vfrom+vto))-(c_int*((vfrom_int+vto_int).^2)))
                    # Iij_ex = (gij_line^2+bij_line^2)*((vfrom+vto)-(2*sqrt(vfrom)*sqrt(vto)*cos(theta_ij)))      # Iij is the square of the longitudinal current

                    # Constraint # 02
                    # Iij    = (gij_line^2+bij_line^2)*(((1-c_int)*(vfrom+vto))+(c_int*((vfrom_int-vto_int).^2)))
                    # Iij_ex = (gij_line^2+bij_line^2)*((vfrom+vto)-(2*sqrt(vfrom)*sqrt(vto)*cos(theta_ij)))      # Iij is the square of the longitudinal current

                    # Iij    = (gij_line^2+bij_line^2)*(((1-c_int)*(vfrom+vto)))

                    Iij    = (gij_line^2+bij_line^2)*((coeff_vf*vfrom)+(coeff_vt*vto)-cons_ij)
                    Iij_ex = (gij_line^2+bij_line^2)*((vfrom+vto)-(2*sqrt(vfrom)*sqrt(vto)*cos(theta_ij)))      # Iij is the square of the longitudinal current

                    # app_err[s,t,i] = ((sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij))-(c_int*(line_alpha[s,t,i]-line_beta[s,t,i])))
                    # app_err[s,t,i] = 2*(gij_line^2+bij_line^2)*(((sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij))-(c_int*(line_alpha[s,t,i]-line_beta[s,t,i]))))
                    # app_err[s,t,i] = 2*(gij_line^2+bij_line^2)*((sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij))-(c_int*((vfrom+vto)/2)-(((vfrom_int-vto_int).^2)/2)))

                    br_crnt_ex = br_crnt_aux(Iij_ex,i)
                    # br_crnt    = br_crnt_aux(Iij,Iij_ex,app_err[s,t,i],i)
                    if Iij<0
                        Iij = Iij
                    else
                        Iij =sqrt(Iij)
                    end
                    br_crnt = (Iij)
                    # br_crnt_ex = sqrt.(Iij_ex)
                    br_crnt_ij_pu[s,t,i] = br_crnt
                    br_crnt_ij_SI[s,t,i] = br_crnt*Ibase[idx_f_bus,1]
                    br_pwr_ij_pu[s,t,i]  = (sqrt(vfrom))*br_crnt
                    br_pwr_ij_SI[s,t,i]  = (vfrom_si*(br_crnt*Ibase[idx_f_bus,1]))/1e6

                    br_crnt_ij_pu_ex[s,t,i] = br_crnt_ex
                    br_crnt_ij_SI_ex[s,t,i] = br_crnt_ex*Ibase[idx_f_bus,1]
                    br_pwr_ij_pu_ex[s,t,i]  = (sqrt(vfrom))*br_crnt_ex
                    br_pwr_ij_SI_ex[s,t,i]  = (vfrom_si*(br_crnt_ex*Ibase[idx_f_bus,1]))/1e6
                end
            end
        end
    end
    return br_crnt_ij_pu,br_crnt_ij_SI,br_pwr_ij_pu,br_pwr_ij_SI,br_crnt_ij_pu_ex,br_crnt_ij_SI_ex,br_pwr_ij_pu_ex,br_pwr_ij_SI_ex,app_err
end

#################### Recovering Branch Current EXACT Values AC-OPF ###################
function recovery_branch_current_exact(nSc,nTP,nLines,nBus,nw_lines,rdata_buses,vm_ac,v_mag_int,va_ac,v_ang_int,yij_line,dLines,Ibase,oltc_ratio)
       # br_crnt_ij_pu  = zeros((length(1:nSc),length(1:nTP),nLines))
    # br_crnt_ij_SI  = zeros((length(1:nSc),length(1:nTP),nLines))

    va = reshape(va_ac,(nSc,nTP,nBus))
    vm = reshape(vm_ac,(nSc,nTP,nBus))
    br_crnt_ij_pu_ex  = zeros((length(1:nSc),length(1:nTP),nLines))
    br_crnt_ij_SI_ex  = zeros((length(1:nSc),length(1:nTP),nLines))
    br_pwr_ij_pu_ex   = zeros((length(1:nSc),length(1:nTP),nLines))
    br_pwr_ij_SI_ex   = zeros((length(1:nSc),length(1:nTP),nLines))
    line_alpha_ex     = zeros((length(1:nSc),length(1:nTP),nLines))
    line_eta_ex       = zeros((length(1:nSc),length(1:nTP),nLines))
    exact_value_ex    = zeros((length(1:nSc),length(1:nTP),nLines))
    line_beta_ex      = zeros((length(1:nSc),length(1:nTP),nLines))
    vm = vm.^2
    for i in 1:nLines
        gij_line  = real(yij_line[i,1])
        bij_line  = imag(yij_line[i,1])
        tap       = nw_lines[i].line_tap_ratio
        idx_tap   = dLines[i,end]
        f_bus     = nw_lines[i].line_from
        t_bus     = nw_lines[i].line_to
        idx_f_bus = findall(x->x==f_bus,rdata_buses[:,1])
        idx_t_bus = findall(x->x==t_bus,rdata_buses[:,1])
        idx_f_bus = idx_f_bus[1,1]
        idx_t_bus = idx_t_bus[1,1]
        idx_line  = i
        for s in 1:nSc
            for t in 1:nTP
                theta_diff_int = v_ang_int[s,t,idx_f_bus]-v_ang_int[s,t,idx_t_bus]
                theta_ij       = va[s,t,idx_f_bus]-va[s,t,idx_t_bus]
                c_int = cos(theta_diff_int)
                vfrom = vm[s,t,idx_f_bus]
                vto   = vm[s,t,idx_t_bus]

                vfrom_si = (sqrt.(vfrom)).*vbase[idx_f_bus]*1000                   # SI voltage is in Volts
                vto_si   = (sqrt.(vto)).*vbase[idx_t_bus]*1000                     # SI voltage is in Volts

                if tap!=0.0
# This is to compare the longitiudinal current constraint with the Imax constraint
                    # index_trsf = Int64(dLines[i,end])
                    # t_tap = oltc_ratio[s,t,index_trsf]
                    t_tap = oltc_ratio[s,t,i]
####--------------------- Linear constraint ------------------------------######
                    # Iij = (gij_line^2+bij_line^2)*(t_tap^2*vfrom+vto-2*t_tap*c_int*(line_alpha[s,t,i]-line_beta[s,t,i])) # Iij is the square of the longitudinal current
                    # Iij = (gij_line^2+bij_line^2)*((t_tap^2*vfrom)+vto-(2*t_tap*1*(line_alpha[s,t,i]))) # Iij is the square of the longitudinal current
                    # Iij = (gij_line^2+bij_line^2)*(t_tap^2*vfrom+vto-2*t_tap*sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij)) # Iij is the square of the longitudinal current
                    Iij = (t_tap^2*vfrom+vto-2*t_tap*sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij)) # Iij is the square of the longitudinal current

                    br_crnt = br_crnt_aux(Iij)
                    br_crnt_ij_pu_ex[s,t,i] = br_crnt
                    br_crnt_ij_SI_ex[s,t,i] = br_crnt*Ibase[idx_t_bus,1]
                    br_pwr_ij_pu_ex[s,t,i]  = (sqrt.(vto))*br_crnt
                    br_pwr_ij_SI_ex[s,t,i]  = (vto_si*(br_crnt*Ibase[idx_t_bus,1]))/1e6

                    # exact_value[s,t,i]     = vfrom*vto*cos(theta_ij)
                    exact_value_ex[s,t,i]     = sqrt.(vfrom)*sqrt.(vto)

                else
                    # Iij = (gij_line^2+bij_line^2)*(vfrom+vto-2*c_int*(line_alpha[s,t,i]-line_beta[s,t,i]))
                    # Iij = (gij_line^2+bij_line^2)*(vfrom+vto-2*1*(line_alpha[s,t,i]))
                    # Iij = (gij_line^2+bij_line^2)*(vfrom+vto-2*sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij)) # Iij is the square of the longitudinal current
                    Iij = (vfrom+vto-2*sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij)) # Iij is the square of the longitudinal current

                    br_crnt = br_crnt_aux(Iij)
                    br_crnt_ij_pu_ex[s,t,i] = br_crnt
                    br_crnt_ij_SI_ex[s,t,i] = br_crnt*Ibase[idx_f_bus,1]
                    br_pwr_ij_pu_ex[s,t,i]  = (sqrt.(vfrom))*br_crnt
                    br_pwr_ij_SI_ex[s,t,i]  = (vfrom_si*(br_crnt*Ibase[idx_f_bus,1]))/1e6


                    exact_value_ex[s,t,i]     = sqrt.(vfrom)*sqrt.(vto)*cos(theta_ij)
                    # exact_value[s,t,i]     = vfrom*vto
                end
            end
        end
    end
    return br_crnt_ij_pu_ex,br_crnt_ij_SI_ex,br_pwr_ij_pu_ex,br_pwr_ij_SI_ex,exact_value_ex
end
############## Recovering Branch Current Values _AC-Power Flow #################
function recovery_branch_current_pf(nSc,nTP,nLines,nw_lines,v_rect,yij_line,Ibase,sbase,ybase,vbase,oltc_ratio,dLines)
    br_crnt_ij_pu   = zeros(Float64,(length(1:nSc),length(1:nTP),nLines))
    br_crnt_ij_SI   = zeros(Float64,(length(1:nSc),length(1:nTP),nLines))

    # br_power_ij_pu  = zeros(Float64,(length(1:nSc),length(1:nTP),nLines))
    # br_power_ij_SI  = zeros(Float64,(length(1:nSc),length(1:nTP),nLines))

    br_power_ij_pu  = zeros(ComplexF64,(length(1:nSc),length(1:nTP),nLines))
    br_power_ij_SI  = zeros(ComplexF64,(length(1:nSc),length(1:nTP),nLines))

    for i in 1:nLines
        f_bus     = nw_lines[i].line_from
        t_bus     = nw_lines[i].line_to
        tap       = nw_lines[i].line_tap_ratio

        idx_f_bus = findall(x->x==f_bus,rdata_buses[:,1])
        idx_t_bus = findall(x->x==t_bus,rdata_buses[:,1])
        idx_f_bus = idx_f_bus[1,1]
        idx_t_bus = idx_t_bus[1,1]

        for s in 1:nSc
            for t in 1:nTP
                vfrom = v_rect[s,t,idx_f_bus]
                vto   = v_rect[s,t,idx_t_bus]

                vfrom_si = vfrom.*vbase[idx_f_bus]*1000                          # SI voltage is in Volts
                vto_si   = vto.*vbase[idx_t_bus]*1000                            # SI voltage is in Volts


                if tap!=0.0                                                      # Transformer Branch
                    index_trsf = Int64(dLines[i,end])
                    tap_trsf = oltc_ratio[s,t,index_trsf]

                    iline_sdg = (tap_trsf^2*yij_line[i,1]*vfrom)-(tap_trsf*yij_line[i,1]*vto)  # Sending End complex Current
                    iline_pu   = iline_sdg
                    iline_SI   = iline_sdg*Ibase[idx_f_bus,1]                         # I line is in A
                    # Sij_sdg_pu  = abs.(vfrom*iline_pu)
                    # Sij_sdg_si  = abs.(vfrom_si*iline_SI)/1e6                      # VA/1e6 = MVA

                    Sij_sdg_pu  = (vfrom*conj(iline_pu))
                    Sij_sdg_si  = (vfrom_si*conj(iline_SI))/1e6                        # VA/1e6 = MVA

                    # iline_rcv = (-tap*yij_line[i,1]*vfrom)+(yij_line[i,1]*vto)     # Receiving End Complex Curent
                    # iline_pu   = iline_rcv
                    # iline_SI   = iline_rcv*Ibase[idx_t_bus,1]                         # I line is in A

                    # Sij_sdg_pu  = abs.(vto*iline_pu)
                    # Sij_sdg_si  = abs.(vto_si*iline_SI)/1e6                      # VA/1e6 = MVA
                else
                    iline_pu    = yij_line[i,1]*(vfrom-vto)                      # This is a complex line current
                    iline_SI    = iline_pu*Ibase[idx_f_bus,1]
                    # Sij_sdg_pu  = abs.(vfrom*iline_pu)
                    # Sij_sdg_si  = abs.(vfrom_si*iline_SI)/1e6                    # Transforming the VA into MVA

                    Sij_sdg_pu  = (vfrom*conj(iline_pu))
                    Sij_sdg_si  = (vfrom_si*conj(iline_SI))/1e6                           # Transforming the VA into MVA
                end
                br_crnt_ij_pu[s,t,i]   = abs.(iline_pu)                           # iline for transformer is taken from the sending end
                br_crnt_ij_SI[s,t,i]   = abs.(iline_SI)                           # From branch is used for the calculating the current
                br_power_ij_pu[s,t,i]  = Sij_sdg_pu
                br_power_ij_SI[s,t,i]  = Sij_sdg_si
            end
        end
    end

    return br_crnt_ij_pu,br_crnt_ij_SI,br_power_ij_pu,br_power_ij_SI
end
################################################################################
################ Recovering Adaptive Power Factor Values #######################
################################################################################
function recovery_adpf(nTP,nSc,nCurt_gen,dg_pf_value,pg_max,i_curt_gens,sbase,dgpf,gen_Curt_Q)
    # dgpf       = zeros(Float64,(length(1:nTP),length(1:nSc),nCurt_gen))
    # gen_Curt_Q = zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    for i in 1:size(dg_pf_value,3)
        for s in 1:nSc
            for t in 1:nTP
                pcnst = pg_max[s,t,i_curt_gens[i][1]]
                dgpf[t,s,i] = atan(dg_pf_value[s,t,i])
                gen_Curt_Q[s,t,i]  = pcnst*dgpf[t,s,i]*sbase
            end
        end
    end
    return dgpf,gen_Curt_Q
end
################################################################################
##################### Recovering DGs Reactive Power ############################
################################################################################
function recovery_react_dg(nTP,nSc,nCurt_gen,Qgen_Curt_pu,sbase)
    gen_Curt_Q = zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    for i in 1:nCurt_gen
        for s in 1:nSc
            for t in 1:nTP
                gen_Curt_Q[s,t,i] = Qgen_Curt_pu[s,t,i]*sbase
            end
        end
    end
    return gen_Curt_Q
end
################################################################################
################ Recovering Reactive Power Injection from DGs ##################
################################################################################
function recovery_reactive_dg(nTP,nSc,nCurt_gen,dg_pf,pg_max,i_curt_gens,sbase,Pgen_Curt_pu)
    dgpf       = dg_pf
    gen_Curt_Q = zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    for i in 1:nCurt_gen
        for s in 1:nSc
            for t in 1:nTP
                p_dg = pg_max[s,t,i_curt_gens[i][1]]-Pgen_Curt_pu[s,t,i]
                gen_Curt_Q[s,t,i]  = p_dg*tan(acos(dgpf))*sbase
            end
        end
    end
    return gen_Curt_Q
end
################################################################################
######################## Recovering OLTC Values ################################
################################################################################
function recovery_oltc(nTP,nSc,nTrsf,trsf_tap,trsf_ratio)
    # trsf_ratio = zeros(Float64,(length(1:nTP),length(1:nSc),nTrsf))
    for i in 1:size(trsf_tap,3)
        for s in 1:nSc
            for t in 1:nTP
                trsf_ratio[t,s,i] = trsf_tap[s,t,i]
            end
        end
    end
    return trsf_ratio
end
################################################################################
##################### Recovering Flexible Load Values ##########################
################################################################################
function recovery_flexible(nTP,nSc,nFl,pfl_inc,pfl_dec,fl_inc_sc,fl_dec_sc,fl_inc_nd,fl_dec_nd)
    # fl_inc_sc = zeros(Float64,(nFl,length(1:nTP),length(1:nSc)))
    # fl_dec_sc = zeros(Float64,(nFl,length(1:nTP),length(1:nSc)))
    #
    # fl_inc_nd = zeros(Float64,(length(1:nSc),length(1:nTP),nFl))
    # fl_dec_nd = zeros(Float64,(length(1:nSc),length(1:nTP),nFl))
##------------------ Extracting information for each scenario --------------####
        for s in 1:nSc
            for i in 1:nFl
                fl_inc_sc[i,:,s] = pfl_inc[s,:,i]
                fl_dec_sc[i,:,s] = pfl_dec[s,:,i]
            end
        end
##------------- Extracting information for each Flexible Load --------------####
        for i in 1:nFl
            for s in 1:nSc
                fl_inc_nd[s,:,i] = pfl_inc[s,:,i]
                fl_dec_nd[s,:,i] = pfl_dec[s,:,i]
            end
        end
    return fl_inc_sc,fl_dec_sc,fl_inc_nd,fl_dec_nd
end
################################################################################
######################### Recovering Storage Values ############################
################################################################################
function recovery_storage(nTP,nSc,nStr_active,p_ch,p_dis,soc,p_strg,p_ch_nd,p_dis_nd,soc_nd,p_strg_nd)
    # p_ch_nd   = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    # p_dis_nd  = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    # soc_nd    = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    # p_strg_nd = zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
##------------- Extracting information for each Flexible Load --------------####
        for i in 1:nStr_active
            for s in 1:nSc
                p_ch_nd[s,:,i]   = p_ch[s,:,i]
                p_dis_nd[s,:,i]  = p_dis[s,:,i]
                soc_nd[s,:,i]    = soc[s,:,i]
                p_strg_nd[s,:,i] = p_strg[s,:,i]
            end
        end
    return p_ch_nd,p_dis_nd,soc_nd,p_strg_nd
end
################################################################################
########### Recovering Generator Injection Values (Curtailable) ################
################################################################################
function recovery_curt_gen(nTP,nSc,nCurt_gen,Pgen_Curt_pu,pgen_tol,qgen_tol,gen_Curt_P_sc,gen_Curt_P_gen,total_p_curt,avg_gen_curt)
    # gen_Curt_P_sc  = zeros(Float64,(nCurt_gen,length(1:nTP),length(1:nSc)))
    # gen_Curt_P_gen = zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    for s in 1:nSc
        for i in 1:nCurt_gen
            gen_Curt_P_sc[i,:,s] = Pgen_Curt_pu[s,:,i]                           # Power curtailed for all DGs in each scenario!
        end
    end
    for i in 1:nCurt_gen
            for s in 1:nSc
                for t in 1:nTP
                    if abs(Pgen_Curt_pu[s,t,i]) < pgen_tol
                        gen_Curt_P_gen[s,t,i] = 0.0                              # Power curtailed of each DG in all scenarios!
                    else
                        gen_Curt_P_gen[s,t,i] = Pgen_Curt_pu[s,t,i]
                    end
                end
            end
    end
    avg_gen_curt_aux = reshape(sum(prob_scs.*gen_Curt_P_gen,dims=1),(nTP,nCurt_gen))
    for i in 1:nCurt_gen
        for t in 1:nTP
            avg_gen_curt[i,t] = avg_gen_curt_aux[t,i]
        end
    end
    # avg_gen_curt     = permutedims(avg_gen_curt_aux)
    total_p_curt = sum(sum(Pgen_Curt_pu,dims=2),dims=3)
    return gen_Curt_P_sc,gen_Curt_P_gen,total_p_curt,avg_gen_curt
end
################################################################################
########### Recovering Generator Injection Values (Non-Curtailable) ############
################################################################################
function recovery_ncurt_gen(nTP,nSc,nNcurt_gen,Pgen_nCurt_pu,Qgen_nCurt_pu,pgen_tol,qgen_tol,gen_Ncurt_P,gen_Ncurt_Q)
    # gen_Ncurt_P = zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))
    # gen_Ncurt_Q = zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))
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
function recovery_voltages(nTP,nSc,nBus,vr,vim,vm_ac,va_ac,vm_min,vm_max,flex_oltc)
    v_rect  = vr+(vim)im
    v_mag  = abs.(v_rect)
    v_tht  = angle.(v_rect)     # angle(z) determines the phase angle in radians of a complex number z
    # vm_ac  = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
    # va_ac  = zeros(Float64,nBus,length(1:nTP),(length(1:nSc)))
    # min_vm = zeros(Float64,nSc,nTP)
    # max_vm = zeros(Float64,nSc,nTP)
    if flex_oltc==1
        dim_bus = nBus+nTrsf_s1
    else
        dim_bus = nBus
    end
    for s in 1:nSc
        for i in 1:dim_bus
            vm_ac[i,:,s] = v_mag[s,:,i]
            va_ac[i,:,s] = v_tht[s,:,i]
        end
        for t in 1:nTP
            vm_min[s,t] = minimum(vm_ac[:,t,s])
            vm_max[s,t] = maximum(vm_ac[:,t,s])
        end
    end
    return vm_ac,va_ac,v_rect,vm_min,vm_max
end
################################################################################
################# Recovering Power Injeciton Values ############################
################################################################################
function recovery_injection(rdata_lines,rdata_buses,nSc,nTP,nBus,node_data,flex_oltc,tratio,vm,va,v_mag_int,v_ang_int,trsf_fict_nodes,ratio_sql,p_inj,q_inj,p_inj_nl,q_inj_nl)
    vm = vm.^2
    ft_line  = [rdata_lines[:,1] rdata_lines[:,2]]
    tf_line  = [rdata_lines[:,2] rdata_lines[:,1]]
    idx_zero = [0.0 0.0]
    # p_inj    = zeros((length(1:nSc),length(1:nTP),nBus))                             # Active power injection at each node
    # q_inj    = zeros((length(1:nSc),length(1:nTP),nBus))                             # Reactive power injection at each node
    iline = []
    p_tol = 1e-4
    q_tol = 1e-4
    for i in 1:nBus
        nd      = node_data[i,1].node_num                                            # Node num is independant of time period and scenario
        nd_num  = unique(nd)
        nd_num  = nd_num[1,1]
        idx_nd_nw_buses = findall(x->x==nd_num,rdata_buses[:,1])
        idx_nd_nw_buses = idx_nd_nw_buses[1,1]
        for s in 1:nSc
            for t in 1:nTP
                pinj_expr = []
                qinj_expr = []
                pinj_expr_nl = []
                qinj_expr_nl = []
                for j in 1:size(nd,1)                                                # length of each node vector in 'node_data' variable
                    gij_line     = node_data[i,1].node_gij_sr[j,1]
                    bij_line     = node_data[i,1].node_bij_sr[j,1]
                    gij_line_sh  = node_data[i,1].node_gij_sh[j,1]
                    bij_line_sh  = node_data[i,1].node_bij_sh[j,1]
                    cnctd_nd     = node_data[i,1].node_cnode[j,1]
                    idx_cnctd_nd = findall(x->x==cnctd_nd,rdata_buses[:,1])
                    idx_cnctd_nd = idx_cnctd_nd[1,1]
                    if  (node_data[i,1].node_tap_ratio[j,1]!=0.0)  &&  (node_data[i,1].node_tap_ratio_min[j,1]!=0.0 && node_data[i,1].node_tap_ratio_max[j,1]!=0.0)                          # Transformer branch

                        trsf_br = [nd_num cnctd_nd]
                        line_ft  = trsf_br.-ft_line
                        line_tf  = trsf_br.-tf_line
                        idx_ft = findall(all(line_ft.==idx_zero,dims=2))
                        idx_tf = findall(all(line_tf.==idx_zero,dims=2))

                        if !isempty(idx_ft) idx_ft = idx_ft[1][1] end
                        if !isempty(idx_tf) idx_tf = idx_tf[1][1] end

                        if !isempty(idx_ft) && isempty(idx_tf)
                            idx_trsf = idx_ft
                        elseif isempty(idx_ft) && !isempty(idx_tf)
                            idx_trsf = idx_tf
                        end

                        status_oltc = rdata_lines[idx_trsf,end]

                        if  nd_num==node_data[i,1].node_from[j,1]
                            if flex_oltc == 1 && status_oltc == 1
                                # tap =  tratio[s,t,node_data[i,1].node_idx_trsf[j,1]]   # Node is connected to the sending end of the transformer (variable reference for the transformer)
                                idx_trsf = node_data[i,1].node_idx_trsf[j,1]
                                # tap = 1.0
                                tap_nl = ratio_sql[s,t,idx_trsf]
                                tap    = ratio_sql[s,t,idx_trsf]
                            elseif (flex_oltc==0 || flex_oltc == 1) && status_oltc == 0
                                tap    = node_data[i,1].node_tap_ratio[j,1]      # Constant value of Transformer tap ratio
                                tap_nl = node_data[i,1].node_tap_ratio[j,1]
                            end
                            (pij,qij,pij_nl,qij_nl) = power_flow_linear_theta_sdg_trsf_recovery_milp(tap,gij_line_sh,bij_line_sh,gij_line,bij_line,vm,va,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,flex_oltc,trsf_fict_nodes,tap_nl,status_oltc)
                        else
                            tap1 = 1.00
                            if flex_oltc == 1 && status_oltc == 1
                                idx_trsf = node_data[i,1].node_idx_trsf[j,1]
                                tap_nl   = ratio_sql[s,t,idx_trsf]

                                # tap2 =  tratio[s,t,node_data[i,1].node_idx_trsf[j,1]]   # Node is connected to the receiving end of the transformer (variable reference for the transformer)
                                # tap2 = 1.00
                                tap2     = ratio_sql[s,t,idx_trsf]
                            elseif (flex_oltc==0 || flex_oltc == 1) && status_oltc == 0
                                tap2   =  node_data[i,1].node_tap_ratio[j,1]                # Constant value of Transformer tap ratio
                                tap_nl =  node_data[i,1].node_tap_ratio[j,1]
                            end
                            (pij,qij,pij_nl,qij_nl) = power_flow_linear_theta_rcv_trsf_recovery_milp(tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,vm,va,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,flex_oltc,trsf_fict_nodes,tap_nl,status_oltc)
                        end
                    elseif (node_data[i,1].node_tap_ratio[j,1]!=0.0)  &&  (node_data[i,1].node_tap_ratio_min[j,1]==0.0 && node_data[i,1].node_tap_ratio_max[j,1]==0.0)
                            (pij,qij,pij_nl,qij_nl) = power_flow_linear_branch_theta_recovery_milp(gij_line_sh,bij_line_sh,gij_line,bij_line,vm,va,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero)
                    else                                                                      # Transmission/Distribution line
                            (pij,qij,pij_nl,qij_nl) = power_flow_linear_branch_theta_recovery_milp(gij_line_sh,bij_line_sh,gij_line,bij_line,vm,va,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero)
                    end
                    push!(pinj_expr,pij)
                    push!(qinj_expr,qij)
                    push!(pinj_expr_nl,pij_nl)
                    push!(qinj_expr_nl,qij_nl)
                end

                p_inj_tmp = sum(pinj_expr[j,1] for j in 1:size(pinj_expr,1))
                q_inj_tmp = sum(qinj_expr[j,1] for j in 1:size(qinj_expr,1))

                p_inj_nl_tmp = sum(pinj_expr_nl[j,1] for j in 1:size(pinj_expr_nl,1))
                q_inj_nl_tmp = sum(qinj_expr_nl[j,1] for j in 1:size(qinj_expr_nl,1))

                if abs(p_inj_tmp) <= p_tol
                    p_inj[s,t,i] = 0.0
                else
                    p_inj[s,t,i] = sum(pinj_expr[j,1] for j in 1:size(pinj_expr,1))
                end
                if abs(q_inj_tmp) <=q_tol
                    q_inj[s,t,i] = 0.0
                else
                    q_inj[s,t,i] = sum(qinj_expr[j,1] for j in 1:size(qinj_expr,1))
                end

                if abs(p_inj_nl_tmp) <= p_tol
                    p_inj_nl[s,t,i] = 0.0
                else
                    p_inj_nl[s,t,i] = sum(pinj_expr_nl[j,1] for j in 1:size(pinj_expr_nl,1))
                end
                if abs(q_inj_nl_tmp) <=q_tol
                    q_inj_nl[s,t,i] = 0.0
                else
                    q_inj_nl[s,t,i] = sum(qinj_expr_nl[j,1] for j in 1:size(qinj_expr_nl,1))
                end

            end
        end
    end
    return p_inj,q_inj,p_inj_nl,q_inj_nl
end

################################################################################
##---- Import/Export Active and Reactive Power Ranges at TSO-DSO Interface------
################################################################################
function recovery_tso_dso_flow(p_imp_pu,q_imp_pu,prob_scs,sbase,p_imp_td,q_imp_td)
    p_imp_si = prob_scs.*p_imp_pu.*sbase;
    q_imp_si = prob_scs.*q_imp_pu.*sbase;

    p_imp_td = sum(p_imp_si,dims=1)
    q_imp_td = sum(q_imp_si,dims=1)
return p_imp_td,q_imp_td
end
################################################################################
################### Recovering Cost of Flexibility  ############################
################################################################################
function recovery_cost_flexibility(nSc,nTP,p_strg_ec,p_strg_nd,fl_inc_nd,fl_dec_nd,gen_Curt_P_gen,trsf_ratio,nStr_active,nFl,nCurt_gen,nTrsf,bus_data_Ssheet,cost_b_str,bus_data_lsheet,nd_fl,cost_load_inc,cost_load_dec,pg_max,i_curt_gens,tdata,trsf_ratio_init,prob_scs,flex_apc,flex_oltc,flex_str,flex_fl)

##------------------- Cost of Storage Flexibility Provision -----------------###
    if flex_str == 1
        cost_str = zeros(Float64,(nSc,nTP,nStr_active))
            for i in 1:nStr_active
                strbus = Int64(bus_data_Ssheet[i,1])
                idx_strbus_strcost = findall(x->x==strbus,bus_data_Ssheet)
                pstr_ec   = p_strg_ec[idx_strbus_strcost,2:end]                  # Storages net power from energy clearnace market
                pstr_ec = pstr_ec'
                    for s in 1:nSc
                        pstr_opf  = p_strg_nd[s,:,idx_strbus_strcost]            # Storages net power obtained during OPF program
                        lcost     = cost_b_str[s,:,idx_strbus_strcost]           # Cost of operation of storages
                        cost_str[s,:,i]  = prob_scs[s,1].*lcost.*abs.(pstr_ec.-pstr_opf)
                    end
            end
            cost_flex_str = sum(cost_str,dims=1)
            cost_flex_str = reshape(cost_flex_str,(nTP,nStr_active))
            cost_flex_str = sum(cost_flex_str,dims=2)
    else
        cost_flex_str = 0.0
    end
##------------- Cost of Flexible Load Flexibility Provision -----------------###
    if flex_fl == 1
        cost_fl = zeros(Float64,(nSc,nTP,nFl))
            for i in 1:nFl
                nd = nd_fl[i,1]
                idx_flbus_flcost = findall(x->x==nd,bus_data_lsheet)
                    for s in 1:nSc
                        cost_fl[s,:,i] = prob_scs[s,1].*cost_load_inc[s,:,i].*abs.(fl_inc_nd[s,:,i].-fl_dec_nd[s,:,i])
                    end
            end
            cost_flex_fl = sum(cost_fl,dims=1)
            cost_flex_fl = reshape(cost_flex_fl,(nTP,nFl))
            cost_flex_fl = sum(cost_flex_fl,dims=2)
    else
        cost_flex_fl = 0.0
    end
##------ Cost of Active Power Curtailement Flexibility Provision ------------###
    if flex_apc == 1
        cost_apc =  zeros(Float64,(nSc,nTP,nCurt_gen))
        cost_curt   = 10.0
            for i in 1:nCurt_gen
                for s in 1:nSc
                    for t in 1:nTP
                        cost_apc[s,:,i] = prob_scs[s,1].*cost_curt.*abs.(gen_Curt_P_gen[s,:,i])
                    end
                end
            end
            cost_flex_apc = sum(cost_apc,dims=1)
            cost_flex_apc = reshape(cost_flex_apc,(nTP,nCurt_gen))
            cost_flex_apc = sum(cost_flex_apc,dims=2)
    else
        cost_flex_apc = 0.0
    end
##----------------- Cost of OLTC Flexibility Provision ----------------------###
    if flex_oltc == 1
        cost_oltc = zeros(Float64,(nSc,nTP,nTrsf))
        cost_tap = 0.001
            for i in 1:nTrsf
                for s in 1:nSc
                        if (tdata[i,1]!=0.0) && (tdata[i,2]!=0.0) && (tdata[i,3]!=0.0)
                        idx  = Int64(tdata[i,1])
                        cost_oltc[s,:,i] = prob_scs[s,1].*cost_tap.*abs.(trsf_ratio_init[s,:,idx].-trsf_ratio[:,s,idx])
                    end
                end
            end
            cost_flex_oltc = sum(cost_oltc,dims=1)
            cost_flex_oltc = reshape(cost_flex_oltc,(nTP,nTrsf))
            cost_flex_oltc = sum(cost_flex_oltc,dims=2)
    else
        cost_flex_oltc = 0.0
    end
##----------------- Total cost of flexibility Provision ----------------------##
cost_flexibility = cost_flex_str.+cost_flex_fl.+cost_flex_apc.+cost_flex_oltc

return cost_flexibility
end

################################################################################
################### Recovering Cost of Flexibility  ############################
################################################################################
function recovery_cost_flexibility_no_str_dgs(nSc,nTP,p_strg_nd,fl_inc_nd,fl_dec_nd,gen_Curt_P_gen,trsf_ratio,nStr_active,nFl,nCurt_gen,nTrsf,bus_data_Ssheet,cost_b_str,bus_data_lsheet,nd_fl,cost_load_inc,cost_load_dec,pg_max,i_curt_gens,tdata,trsf_ratio_init,prob_scs,flex_apc,flex_oltc,flex_str,flex_fl)

##------------------- Cost of Storage Flexibility Provision -----------------###
    if flex_str == 1
        cost_str = zeros(Float64,(nSc,nTP,nStr_active))
            for i in 1:nStr_active
                strbus = Int64(bus_data_Ssheet[i,1])
                idx_strbus_strcost = findall(x->x==strbus,bus_data_Ssheet)
                    for s in 1:nSc
                        pstr_opf  = p_strg_nd[s,:,idx_strbus_strcost]            # Storages net power obtained during OPF program
                        lcost     = cost_b_str[s,:,idx_strbus_strcost]           # Cost of operation of storages
                        cost_str[s,:,i]  = prob_scs[s,1].*lcost.*abs.(-pstr_opf)
                    end
            end
            cost_flex_str = sum(cost_str,dims=1)
            cost_flex_str = reshape(cost_flex_str,(nTP,nStr_active))
            cost_flex_str = sum(cost_flex_str,dims=2)
    else
        cost_flex_str = 0.0
    end
##------------- Cost of Flexible Load Flexibility Provision -----------------###
    if flex_fl == 1
        cost_fl = zeros(Float64,(nSc,nTP,nFl))
            for i in 1:nFl
                nd = nd_fl[i,1]
                idx_flbus_flcost = findall(x->x==nd,bus_data_lsheet)
                    for s in 1:nSc
                        cost_fl[s,:,i] = prob_scs[s,1].*cost_load_inc[s,:,i].*abs.(fl_inc_nd[s,:,i].-fl_dec_nd[s,:,i])
                    end
            end
            cost_flex_fl = sum(cost_fl,dims=1)
            cost_flex_fl = reshape(cost_flex_fl,(nTP,nFl))
            cost_flex_fl = sum(cost_flex_fl,dims=2)
    else
        cost_flex_fl = 0.0
    end
##------ Cost of Active Power Curtailement Flexibility Provision ------------###
    if flex_apc == 1
        cost_apc =  zeros(Float64,(nSc,nTP,nCurt_gen))
        cost_curt   = 10.0
            for i in 1:nCurt_gen
                for s in 1:nSc
                    for t in 1:nTP
                        cost_apc[s,:,i] = prob_scs[s,1].*cost_curt.*abs.(-(pg_max[s,:,i_curt_gens[i][1]].-gen_Curt_P_gen[s,:,i]))
                    end
                end
            end
            cost_flex_apc = sum(cost_apc,dims=1)
            cost_flex_apc = reshape(cost_flex_apc,(nTP,nCurt_gen))
            cost_flex_apc = sum(cost_flex_apc,dims=2)
    else
        cost_flex_apc = 0.0
    end
##----------------- Cost of OLTC Flexibility Provision ----------------------###
    if flex_oltc == 1
        cost_oltc = zeros(Float64,(nSc,nTP,nTrsf))
        cost_tap = 0.001
            for i in 1:nTrsf
                for s in 1:nSc
                        if (tdata[i,1]!=0.0) && (tdata[i,2]!=0.0) && (tdata[i,3]!=0.0)
                        idx  = Int64(tdata[i,1])
                        cost_oltc[s,:,i] = prob_scs[s,1].*cost_tap.*abs.(trsf_ratio_init[s,:,idx].-trsf_ratio[:,s,idx])
                    end
                end
            end
            cost_flex_oltc = sum(cost_oltc,dims=1)
            cost_flex_oltc = reshape(cost_flex_oltc,(nTP,nTrsf))
            cost_flex_oltc = sum(cost_flex_oltc,dims=2)
    else
        cost_flex_oltc = 0.0
    end
##----------------- Total cost of flexibility Provision ----------------------##
cost_flexibility = cost_flex_str.+cost_flex_fl.+cost_flex_apc.+cost_flex_oltc

return cost_flexibility
end

################################################################################
################### Recovering Cost of Flexibility  ############################
################################################################################
function recovery_cost_flexibility_no_str(nSc,nTP,p_strg_nd,fl_inc_nd,fl_dec_nd,gen_Curt_P_gen,trsf_ratio,nStr_active,nFl,nCurt_gen,nTrsf,bus_data_Ssheet,cost_b_str,bus_data_lsheet,nd_fl,cost_load_inc,cost_load_dec,pg_max,i_curt_gens,tdata,trsf_ratio_init,prob_scs,flex_apc,flex_oltc,flex_str,flex_fl)

##------------------- Cost of Storage Flexibility Provision -----------------###
    if flex_str == 1
        cost_str = zeros(Float64,(nSc,nTP,nStr_active))
            for i in 1:nStr_active
                strbus = Int64(bus_data_Ssheet[i,1])
                idx_strbus_strcost = findall(x->x==strbus,bus_data_Ssheet)
                    for s in 1:nSc
                        pstr_opf  = p_strg_nd[s,:,idx_strbus_strcost]            # Storages net power obtained during OPF program
                        lcost     = cost_b_str[s,:,idx_strbus_strcost]           # Cost of operation of storages
                        cost_str[s,:,i]  = prob_scs[s,1].*lcost.*abs.(-pstr_opf)
                    end
            end
            cost_flex_str = sum(cost_str,dims=1)
            cost_flex_str = reshape(cost_flex_str,(nTP,nStr_active))
            cost_flex_str = sum(cost_flex_str,dims=2)
    else
        cost_flex_str = 0.0
    end
##------------- Cost of Flexible Load Flexibility Provision -----------------###
    if flex_fl == 1
        cost_fl = zeros(Float64,(nSc,nTP,nFl))
            for i in 1:nFl
                nd = nd_fl[i,1]
                idx_flbus_flcost = findall(x->x==nd,bus_data_lsheet)
                    for s in 1:nSc
                        cost_fl[s,:,i] = prob_scs[s,1].*cost_load_inc[s,:,i].*abs.(fl_inc_nd[s,:,i].-fl_dec_nd[s,:,i])
                    end
            end
            cost_flex_fl = sum(cost_fl,dims=1)
            cost_flex_fl = reshape(cost_flex_fl,(nTP,nFl))
            cost_flex_fl = sum(cost_flex_fl,dims=2)
    else
        cost_flex_fl = 0.0
    end
##------ Cost of Active Power Curtailement Flexibility Provision ------------###
    if flex_apc == 1
        cost_apc =  zeros(Float64,(nSc,nTP,nCurt_gen))
        cost_curt   = 10.0
            for i in 1:nCurt_gen
                for s in 1:nSc
                    for t in 1:nTP
                        cost_apc[s,:,i] = prob_scs[s,1].*cost_curt.*abs.(gen_Curt_P_gen[s,:,i])
                    end
                end
            end
            cost_flex_apc = sum(cost_apc,dims=1)
            cost_flex_apc = reshape(cost_flex_apc,(nTP,nCurt_gen))
            cost_flex_apc = sum(cost_flex_apc,dims=2)
    else
        cost_flex_apc = 0.0
    end
##----------------- Cost of OLTC Flexibility Provision ----------------------###
    if flex_oltc == 1
        cost_oltc = zeros(Float64,(nSc,nTP,nTrsf))
        cost_tap = 0.001
            for i in 1:nTrsf
                for s in 1:nSc
                        if (tdata[i,1]!=0.0) && (tdata[i,2]!=0.0) && (tdata[i,3]!=0.0)
                        idx  = Int64(tdata[i,1])
                        cost_oltc[s,:,i] = prob_scs[s,1].*cost_tap.*abs.(trsf_ratio_init[s,:,idx].-trsf_ratio[:,s,idx])
                    end
                end
            end
            cost_flex_oltc = sum(cost_oltc,dims=1)
            cost_flex_oltc = reshape(cost_flex_oltc,(nTP,nTrsf))
            cost_flex_oltc = sum(cost_flex_oltc,dims=2)
    else
        cost_flex_oltc = 0.0
    end
##----------------- Total cost of flexibility Provision ----------------------##
cost_flexibility = cost_flex_str.+cost_flex_fl.+cost_flex_apc.+cost_flex_oltc

return cost_flexibility
end
