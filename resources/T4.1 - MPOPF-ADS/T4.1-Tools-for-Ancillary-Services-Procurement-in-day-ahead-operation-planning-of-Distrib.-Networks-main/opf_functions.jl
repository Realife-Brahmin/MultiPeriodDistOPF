function  opf_model_initialization(v_mag_int,v_ang_int,nBus,nNcurt_gen,nCurt_gen,nTrsf,nStr_active,nFl,nSc,nTP,flex_oltc,flex_adpf,flex_str,flex_fl,str_bin,fl_bin,rdata_buses,nLines,trsf_fict_nodes,tap_ratio_range,oltc_bin,solver_ipopt,solver_cbc,solver_bonmin,idx_slack_from,idx_slack_to,num_slack_cnctd_nodes,solver_cplex,int_tol,opt_tol,num_thread,nTrsf_s1)
#-------------------------- Model Initialization -------------------------------
#---------------------------Decalration of OPF Variables------------------------
if solver_ipopt == 1
acopf = Model()
set_optimizer(acopf,Ipopt.Optimizer)                                             # Here the solver is being set
set_optimizer_attributes(acopf,"tol"=>1e-6)  # "max_iter"=>100
## --------------------------- Cbc Solver --------------------------------------#
elseif solver_cbc==1
acopf = Model()
set_optimizer(acopf,Cbc.Optimizer)                                               # Here the solver is being set
# set_optimizer_attributes(acopf,"logLevel"=>1)  # "max_iter"=>100
set_optimizer_attributes(acopf,"logLevel"=>2,"allowableGap"=>1e-1,"ratioGap"=>0.05)  # "max_iter"=>100
# set_optimizer_attributes(acopf,"logLevel"=>2,"allowableGap"=>1e-3,"ratioGap"=>0.001)
#-------------------- AmplNLWriter (Bonmin) Solver -----------------------------
elseif solver_bonmin == 1
optimizer  = AmplNLWriter.Optimizer
acopf = Model(() -> AmplNLWriter.Optimizer("C:/bonmin-win64/bonmin.exe",["print_level=1","bonmin.file_solution=yes","file_print_level=5"]))
acopf = Model(() -> AmplNLWriter.Optimizer("C:/bonmin-win64/bonmin.exe"))

elseif solver_cplex == 1
    acopf = direct_model(CPLEX.Optimizer())
    # set_optimizer_attribute(acopf, "CPX_PARAM_EPINT", int_tol) # Disabled on May 5th.
    # set_optimizer_attribute(acopf, "CPX_PARAM_EPOPT", opt_tol) # Disabled on May 5th.
    # set_optimizer_attribute(acopf, "CPX_PARAM_BARALG", 3) # Disabled on May 5th.
    # set_optimizer_attribute(acopf, "CPX_PARAM_THREADS", num_thread) # Disabled on May 5th.
    # set_optimizer_attribute(acopf, "CPX_PARAM_DATACHECK", CPX_ON) # Disabled on May 5th.
    set_optimizer_attribute(acopf, "CPX_PARAM_SOLUTIONTYPE", 2)
    set_optimizer_attribute(acopf, "CPX_PARAM_LPMETHOD", 4) # ADDED BY IMAN ON MARCH 31ST
end
################## Transformer Ratio variables #################################
if flex_oltc == 0
    @variables(acopf,begin
        (v_mag[i=1:nSc, j=1:nTP,k=1:nBus])
        (v_ang[i=1:nSc, j=1:nTP,k=1:nLines])                                     # This is theta_ij variable
    end
    )
        aux_tratio = 0.0
        tratio = 0.0
        tratio_aux =0.0

elseif flex_oltc == 1
        @variables(acopf,begin
        (v_mag[i=1:nSc, j=1:nTP,k=1:(nBus+nTrsf_s1)])                            # vol_mag variable is extended for fictitious nodes
        (v_ang[i=1:nSc, j=1:nTP,k=1:nLines])
        aux_tratio[i=1:nSc, j=1:nTP,k=1:nTrsf_s1]                                # This variabkle is used to set-up the cost of tap ratio in the OF
        tratio[i=1:nSc, j=1:nTP,k=1:nTrsf_s1]
        tratio_aux[i=1:nSc, j=1:nTP,k=1:nTrsf_s1]
    end
    )
end
###-------------- Linearizing transformer tap constraint --------------------###
if flex_oltc == 1 && oltc_bin == 1
@variables(acopf,begin
    (aux_tap[i=1:nSc,j=1:nTP,k=1:size(tap_ratio_range,1),l=1:nTrsf_s1])
    (bin_tap[i=1:nSc,j=1:nTP,k=1:size(tap_ratio_range,1),l=1:nTrsf_s1],Bin)
    end
    )
else
    aux_tap = 0.0
    bin_tap = 0.0
end
############### -------- Slack bus Angle Variables -----------------------#####
@variable(acopf,v_ang_node[i=1:nSc, j=1:nTP,k=1:size(slack_cnctd_nodes,1)+1])

#---------- nNcurt_gen and nCurt_gen now forms the complete nGens --------------
@variables(acopf,begin
    Pg[i=1:nSc, j=1:nTP,k=1:nNcurt_gen]                                          # Pg = non-curtailable Gen
    Qg[i=1:nSc, j=1:nTP,k=1:nNcurt_gen]
    Pg_curt[i=1:nSc, j=1:nTP,k=1:nCurt_gen]
    Qg_curt[i=1:nSc, j=1:nTP,k=1:nCurt_gen]
    end
    )
#--------------------- Adaptive power factor -----------------------------------
# if flex_adpf == 1
#     # @variable(acopf,adpf[i=1:nSc,j=1:nTP,k=1:nCurt_gen])
#     @variable(acopf,Qg_curt[i=1:nSc,j=1:nTP,k=1:nCurt_gen])
# elseif flex_adpf == 0
#     Qg_curt = 0
# end
#----------------- Variables for Energy Storage Model --------------------------
if flex_str == 1
# @timeit md_init "md_init_4"
@variables(acopf,begin
        p_ch[i=1:nSc,j=1:nTP,k=1:nStr_active]
        p_dis[i=1:nSc,j=1:nTP,k=1:nStr_active]
        soc[i=1:nSc,j=1:nTP,k=1:nStr_active]
        p_strg[i=1:nSc,j=1:nTP,k=1:nStr_active]
    end)
    if str_bin == 1
        @variable(acopf,bin_ch[i=1:nSc,j=1:nTP,k=1:nStr_active],Bin)              # Use these variables for mixed-integer model of a battery
    elseif str_bin == 0
        bin_ch  = 0
    end
elseif flex_str == 0
    @variables(acopf,begin
        p_ch[i=1:nSc,j=1:nTP,k=1:nStr_active]
        p_dis[i=1:nSc,j=1:nTP,k=1:nStr_active]
        soc[i=1:nSc,j=1:nTP,k=1:nStr_active]
        p_strg[i=1:nSc,j=1:nTP,k=1:nStr_active]
    end)
    bin_ch  = 0
end
## --------------- Variables for Felxible Loads --------------------------------
if flex_fl == 1
@variables(acopf,begin
        p_fl_inc[i=1:nSc,j=1:nTP,k=1:nFl]
        p_fl_dec[i=1:nSc,j=1:nTP,k=1:nFl]
        end )
        if fl_bin == 1
            @variable(acopf,bin_fl[i=1:nSc,j=1:nTP,k=1:nFl],Bin)
        elseif fl_bin == 0
            bin_fl = 0
        end
elseif flex_fl == 0
    @variables(acopf,begin
        p_fl_inc[i=1:nSc,j=1:nTP,k=1:nFl]
        p_fl_dec[i=1:nSc,j=1:nTP,k=1:nFl]
        end )
        bin_fl = 0
end
################################################################################
############# ---------- Fixing Slack Bus Angle --------------------############
################################################################################
fix.(v_ang_node[:,:,1],0.0)
obj_exp = Dict()                                                                 # Use to save the expression of objective
con_exp = Dict{Symbol,Any}(:pb_p=>Dict{String,Any}(),:pb_q=>Dict{String,Any}(),:vol => Dict{String,Any}(),:Icrnt=>Dict{String,Any}(),:S_from=>Dict{String,Any}(),:S_to=>Dict{String,Any}(),:pft=>Dict{String,Any}(),:qft=>Dict{String,Any}(),:ptf=>Dict{String,Any}(),:qtf=>Dict{String,Any}())    # pb = power balance, vol = voltage, Icrnt = Longitudnal current

return acopf,v_mag,v_ang,Pg,Qg,Pg_curt,aux_tratio,tratio,Qg_curt,p_ch,p_dis,soc,p_strg,bin_ch,p_fl_inc,p_fl_dec,bin_fl,obj_exp,con_exp,aux_tap,bin_tap,v_ang_node,tratio_aux
end
#----------------------Script for Objective Functions--------------------------
function opf_model_objective(acopf,flex_oltc,flex_str,flex_fl,flex_apc,nTrsf_s1,nStr_active,nCurt_gen,nFl,nd_fl,nSc,nTP,prob_scs,time_step,tdata,bus_data_Ssheet,bus_data_lsheet,cost_a_str,cost_b_str,cost_c_str,cost_load_inc,cost_load_dec,aux_tratio,p_ch,p_dis,p_fl_inc,p_fl_dec,Pg_curt,sbase,stoch_model)
################################################################################
##-------------------------- Cost of transformer tap ---------------------------
################################################################################
cost_trsf_tap = []
cost_tap = 0.001                                                                 # $dollar/tap
if flex_oltc == 1
    for i in 1:nTrsf_s1
        for s in 1:nSc
            for t in 1:nTP
                # if (tdata[i,1]!=0.0) && (tdata[i,2]!=0.0) && (tdata[i,3]!=0.0)
                # idx  = Int64(tdata[i,1])
                    if stoch_model == 0
                        cost = cost_tap*(aux_tratio[s,t,i])
                    elseif stoch_model == 1
                        cost = prob_scs[s,1]*cost_tap*(aux_tratio[s,t,i])
                    end
                cost_trsf = @expression(acopf,cost)
                push!(cost_trsf_tap,cost_trsf)
                # end
            end
        end
    end
elseif isempty(cost_trsf_tap)
    cost_trsf_tap = 0.0
end

# cost_trsf_tap = []
# cost_tap = 0.001                                                                 # $dollar/tap
# if flex_oltc == 1
#     for i in 1:nTrsf
#         for s in 1:nSc
#             for t in 1:nTP
#                 if (tdata[i,1]!=0.0) && (tdata[i,2]!=0.0) && (tdata[i,3]!=0.0)
#                 idx  = Int64(tdata[i,1])
#                     if stoch_model == 0
#                         cost = cost_tap*(aux_tratio[s,t,idx])
#                     elseif stoch_model == 1
#                         cost = prob_scs[s,1]*cost_tap*(aux_tratio[s,t,idx])
#                     end
#                 cost_trsf = @expression(acopf,cost)
#                 push!(cost_trsf_tap,cost_trsf)
#                 end
#             end
#         end
#     end
# elseif isempty(cost_trsf_tap)
#     cost_trsf_tap = 0.0
# end
################################################################################
##-------------------------- Cost of Storages ----------------------------------
################################################################################
cost_pstr = []
rep_Str_nodes = []
if flex_str == 1
    for i in 1:nStr_active
        for s in 1:nSc
            for t in 1:nTP
            strbus = Int64(bus_data_Ssheet[i,1])
            idx_strbus_strcost = findall(x->x==strbus,bus_data_Ssheet)           # Index of generator bus in the gen cost sheet

                if isempty(findall(x->x==strbus,rep_Str_nodes))
                    pstr  = -p_ch[s,t,idx_strbus_strcost]+p_dis[s,t,idx_strbus_strcost]
                    qcost = cost_a_str[s,t,idx_strbus_strcost]
                    lcost = cost_b_str[s,t,idx_strbus_strcost]
                    ccost = cost_c_str[s,t,idx_strbus_strcost]
                    aux_cost  = @expression(acopf,sum(lcost[j,1]*pstr[j,1] for j in 1:size(lcost,1))+sum(ccost[j,1] for j in 1:size(ccost,1)))
                    if stoch_model == 0
                        cost  = aux_cost*time_step                                   # Sbase should not be multiplied here because cost has already been multiplied by Sbase in DataLoad function
                    elseif stoch_model == 1
                        cost  = prob_scs[s,1]*aux_cost*time_step                     # Sbase should not be multiplied here because cost has already been multiplied by Sbase in DataLoad function
                    end
                    push!(cost_pstr,cost)
                end
                if size(idx_strbus_strcost,1)!=1                                 # More than 1 gen are conncetd at a node
                push!(rep_Str_nodes,strbus)
                end
            end
        end
    end
elseif isempty(cost_pstr)
    cost_pstr = 0.0
end
################################################################################
##------------------------- Cost of flexible loads -----------------------------
################################################################################
cost_pfl = []
rep_fl_nodes = []
if flex_fl == 1
    for i in 1:nFl
        for s in 1:nSc
            for t in 1:nTP
            nd = nd_fl[i,1]
            idx_flbus_flcost = findall(x->x==nd,bus_data_lsheet)                  # Index of generator bus in the gen cost sheet

                if isempty(findall(x->x==nd,rep_fl_nodes))
                    pfl_inc  = p_fl_inc[s,t,i]
                    pfl_dec  = p_fl_dec[s,t,i]
                    c_inc    = cost_load_inc[s,t,i]
                    c_dec    = cost_load_dec[s,t,i]
                    cost     = @expression(acopf,c_inc*pfl_inc+c_dec*pfl_dec)
                        if stoch_model == 0
                            cost  = cost*time_step                               # Sbase should not be multiplied here because cost has already been multiplied by Sbase in DataLoad function
                        elseif stoch_model == 1
                            cost  = prob_scs[s,1]*cost*time_step                 # Sbase should not be multiplied here because cost has already been multiplied by Sbase in DataLoad function
                        end
                    push!(cost_pfl,cost)
                end

                if size(idx_flbus_flcost,1)!=1                                   # More than 1 gen are conncetd at a node
                push!(rep_fl_nodes,nd)
                end
            end
        end
    end
elseif isempty(cost_pfl)
    cost_pfl = 0.0
end
################################################################################
##################### Cost of Active Power Curtailment #########################
##-------------- Minimizing the active power curtailment cost ------------------
################################################################################
curt_energy = []
cost_curt   = 10.0   # euro/MWh
if flex_apc == 1
    for i in 1:nCurt_gen
        for s in 1:nSc
            for t in 1:nTP
                if stoch_model == 0
                    crt_energy = cost_curt*time_step*sbase*((Pg_curt[s,t,i]))             # So, it has to be multiplied here!
                elseif stoch_model == 1
                    crt_energy = cost_curt*time_step*sbase*(prob_scs[s,1]*(Pg_curt[s,t,i])) # Sbase has to be multiplied here because cost of APC has not been defined anywhere else in the program.
                end
            push!(curt_energy,crt_energy)
            end
        end
    end
elseif isempty(curt_energy)
    curt_energy = 0.0
end

a = sum(curt_energy)
b = sum(cost_pfl)
c = sum(cost_pstr)
d = sum(cost_trsf_tap)
################################################################################
##------------------- Total cost of all components -----------------------------
################################################################################
# cost_expr = sum(curt_energy[j,1] for j in 1:size(curt_energy,1))+sum(cost_pfl[j,1] for j in 1:size(cost_pfl,1))+sum(cost_pstr[j,1] for j in 1:size(cost_pstr,1))+sum(cost_trsf_tap[j,1] for j in 1:size(cost_trsf_tap,1))
cost_expr = a+b+c+d
return cost_expr
end
################################################################################
##---------------------Script for Power Balance Constraint-----------------------
################################################################################
function opf_model_power_balance_cons(acopf,Pg_curt,Pg,Qg,p_strg,p_fl_inc,p_fl_dec,tratio,v_mag,v_ang,nw_buses,rdata_buses,rdata_loads,bus_data_lsheet,bus_data_Ssheet,node_data,nd_fl,nd_curt_gen,nd_ncurt_gen,p_load,q_load,idx_Gs_lsheet,idx_Bs_lsheet,pg_max,yii_sh,i_curt_gens,sbase,dg_pf,nSc,nTP,nBus,flex_oltc,flex_adpf,v_mag_int,v_ang_int,con_exp,trsf_fict_nodes,v_initial_pf,Qg_curt,min_vol_limit,max_vol_limit,lin_itr,p_dis,p_ch,rdata_trsf,nTrsf_s1,nw_trsf)
    ft_line = [rdata_lines[:,1] rdata_lines[:,2]]
    tf_line = [rdata_lines[:,2] rdata_lines[:,1]]
    idx_zero = [0.0 0.0]
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
        gii_sh  = check_shunt_milp(Pii_sh,idx_nd_nw_buses,nw_buses)
        bii_sh  = check_shunt_milp(Qii_sh,idx_nd_nw_buses,nw_buses)
        push!(yii_sh,gii_sh+(bii_sh)im)

        for s in 1:nSc
            for t in 1:nTP
##------------- If gens are both curtailable and  non-curtailable --------------
                (pg,qg)   = nodal_components_values_milp(Pg_curt,Pg,Qg,s,t,idx_ncurt_gen,idx_curt_gen,Qg_curt)
##-------------------- If gens are only non-curtailable ------------------------
            # pg = nodal_components_values(Pg,s,t,idx_bus_gsheet)
            # qg = nodal_components_values(Qg,s,t,idx_bus_gsheet)
##--------- Saving flexible load data, fixed load and storage data -----------##
                p_flx_inc = nodal_components_values_milp(p_fl_inc,s,t,idx_fl)
                p_flx_dec = nodal_components_values_milp(p_fl_dec,s,t,idx_fl)
                pd        = nodal_components_values_milp(p_load,s,t,idx_bus_lsheet)
                qd        = nodal_components_values_milp(q_load,s,t,idx_bus_lsheet)
                pstr      = nodal_components_values_milp(p_strg,s,t,idx_bus_Stsheet)
                p_str_ch  = nodal_components_values_milp(p_ch,s,t,idx_bus_Stsheet)
                p_str_dis = nodal_components_values_milp(p_dis,s,t,idx_bus_Stsheet)
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
                    if  (node_data[i,1].node_tap_ratio[j,1]!=0.0)  &&  (node_data[i,1].node_tap_ratio_min[j,1]!=0.0 && node_data[i,1].node_tap_ratio_max[j,1]!=0.0) # Transformer branch

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

                            if  nd_num==node_data[i,1].node_from[j,1]            # Sending ending of the transformer branch
                                if flex_oltc == 1 && status_oltc == 1            # Transformer being used as a flexibilty device
                                    tap = 1.0
                                elseif (flex_oltc == 1 || flex_oltc==0) && status_oltc == 0 # Tap ratio is fixed and transformer does not participate in the flexibility
                                    tap =  node_data[i,1].node_tap_ratio[j,1]         # The hard value is used to validate the code against Carleton model!
                                end
                                (pij,qij) = power_flow_linear_theta_sdg_trsf_milp(acopf,tap,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,trsf_fict_nodes,flex_oltc,status_oltc)
                            else
                                tap1 = 1.0
                                if flex_oltc == 1 && status_oltc == 1
                                    tap2 = 1.0
                                elseif (flex_oltc == 1||flex_oltc==0) && status_oltc == 0
                                    tap2 =  node_data[i,1].node_tap_ratio[j,1]                # The hard value is used to validate the code against Carleton model!
                                end
                                (pij,qij) = power_flow_linear_theta_rcv_trsf_milp(acopf,tap1,tap2,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero,trsf_fict_nodes,flex_oltc,status_oltc)
                            end
                    elseif (node_data[i,1].node_tap_ratio[j,1]!=0.0)  &&  (node_data[i,1].node_tap_ratio_min[j,1]==0.0 && node_data[i,1].node_tap_ratio_max[j,1]==0.0)
                            println("WARNING! The transformer tap ratio is non-zero but both minimum and maximum tap ratios are set to zero")
                            println("This transformer branch is treated as a simple line while setting the power injection constraint!")
                            tap = 1.0
                            (pij,qij) = power_flow_linear_branch_theta_milp(acopf,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero)
                    else                                                                      # Transmission/Distribution line
                            tap = 1.0
                            (pij,qij) = power_flow_linear_branch_theta_milp(acopf,gij_line_sh,bij_line_sh,gij_line,bij_line,v_mag,v_ang,v_mag_int,v_ang_int,s,t,idx_nd_nw_buses,idx_cnctd_nd,nd_num,cnctd_nd,ft_line,tf_line,idx_zero)
                    end
                    push!(pinj_expr,pij)
                    push!(qinj_expr,qij)
                end
                p_inj = @expression(acopf,sum(pinj_expr[j,1] for j in 1:size(pinj_expr,1)))
                q_inj = @expression(acopf,sum(qinj_expr[j,1] for j in 1:size(qinj_expr,1)))

## ------ Power balance constraints when DG curtailement, Storage and ----------
##------------- load flexibility variables are considered-----------------------

                if isempty(idx_ncurt_gen) && ~isempty(idx_curt_gen)
################################################################################
##----------Constant power factor of DGs is considered!-------------------------
################################################################################
                    if flex_adpf == 0                                                # The power factor can be 1 (No Reactive power) or other than 1 (Reactive power is injected)
                        if flex_apc == 1
                            p_avl    = pg_max[s,t,i_curt_gens[idx_curt_gen][1]]
                            @constraint(acopf,(sum(p_avl[j,1] for j in 1:size(p_avl,1))-sum(pg[j,1] for j in 1:size(pg,1)))+(sum(p_str_dis[j,1] for j in 1:size(p_str_dis,1))+sum(p_str_ch[j,1] for j in 1:size(p_str_ch,1)))-(sum(pd[j,1] for j in 1:size(pd,1))+sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(gii_sh[j,1] for j in 1:size(gii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==p_inj)
                            @constraint(acopf,tan(acos(dg_pf))*(sum(p_avl[j,1] for j in 1:size(p_avl,1))-sum(pg[j,1] for j in 1:size(pg,1)))-(sum(qd[j,1] for j in 1:size(qd,1))+load_theta*sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-load_theta*sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(bii_sh[j,1] for j in 1:size(bii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==q_inj) # Q_inj of DG DOES NOT depend upon P_curt

                        elseif flex_apc == 0
                            p_avl    = pg_max[s,t,i_curt_gens[idx_curt_gen][1]]
                            @constraint(acopf,sum(p_avl[j,1] for j in 1:size(p_avl,1))+(sum(p_str_dis[j,1] for j in 1:size(p_str_dis,1))+sum(p_str_ch[j,1] for j in 1:size(p_str_ch,1)))-(sum(pd[j,1] for j in 1:size(pd,1))+sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(gii_sh[j,1] for j in 1:size(gii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==p_inj)
                            @constraint(acopf,-(sum(qd[j,1] for j in 1:size(qd,1))+load_theta*sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-load_theta*sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(bii_sh[j,1] for j in 1:size(bii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==q_inj) # Q_inj of DG DOES NOT depend upon P_curt
                        end
################################################################################
##-------------- Adaptive power factor of DGs is considered --------------------
################################################################################
                    elseif flex_adpf == 1
                        p_avl = pg_max[s,t,i_curt_gens[idx_curt_gen][1]]
                        idx   = idx_curt_gen[1]
                        @constraint(acopf,sum(p_avl[j,1] for j in 1:size(p_avl,1))-sum(pg[j,1] for j in 1:size(pg,1))+(sum(p_str_dis[j,1] for j in 1:size(p_str_dis,1))+sum(p_str_ch[j,1] for j in 1:size(p_str_ch,1)))-(sum(pd[j,1] for j in 1:size(pd,1))+sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(gii_sh[j,1] for j in 1:size(gii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==p_inj)
                        @constraint(acopf,sum(qg[j,1] for j in 1:size(qg,1))-(sum(qd[j,1] for j in 1:size(qd,1))+load_theta*sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-load_theta*sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(bii_sh[j,1] for j in 1:size(bii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==q_inj)
                    end
##------------------------------------------------------------------------------
                else
                    @constraint(acopf,sum(pg[j,1] for j in 1:size(pg,1))+(sum(p_str_dis[j,1] for j in 1:size(p_str_dis,1))+sum(p_str_ch[j,1] for j in 1:size(p_str_ch,1)))-(sum(pd[j,1] for j in 1:size(pd,1))+sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(gii_sh[j,1] for j in 1:size(gii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==p_inj)
                    @constraint(acopf,sum(qg[j,1] for j in 1:size(qg,1))-(sum(qd[j,1] for j in 1:size(qd,1))+load_theta*sum(p_flx_inc[j,1] for j in 1:size(p_flx_inc,1))-load_theta*sum(p_flx_dec[j,1] for j in 1:size(p_flx_dec,1)))+(sum(bii_sh[j,1] for j in 1:size(bii_sh,1))*(v_mag[s,t,idx_nd_nw_buses]))==q_inj)
                end
################################################################################
############## Setting voltage constraints on network buses ####################
################################################################################
                     # @constraint(acopf,(nw_buses[idx_nd_nw_buses].bus_vmin)^2<=(v_mag[s,t,idx_nd_nw_buses])<=(nw_buses[idx_nd_nw_buses].bus_vmax)^2)
                if nw_buses[i].bus_type==3                                       # Fixing the voltage at that bus whose type is 3 (Assuming this to be the Point of Common Coupling!)
                    # @constraint(acopf,(v_mag[s,t,idx_nd_nw_buses]<=(nw_buses[idx_nd_nw_buses].bus_vmax)^2))
                    # @constraint(acopf,(nw_buses[idx_nd_nw_buses].bus_vmin)^2<=(v_mag[s,t,idx_nd_nw_buses]))
                    if nTP == 24
                        @constraint(acopf,(v_mag[s,t,idx_nd_nw_buses])==v_initial_pf[1,t].^2)
                    else
                        @constraint(acopf,(v_mag[s,t,idx_nd_nw_buses])==v_initial_pf.^2)
                    end
                else                                                             # Constraint on the remaining buses
                    @constraint(acopf,(v_mag[s,t,idx_nd_nw_buses])<=(nw_buses[idx_nd_nw_buses].bus_vmax)^2)
                    @constraint(acopf,(nw_buses[idx_nd_nw_buses].bus_vmin)^2<=(v_mag[s,t,idx_nd_nw_buses]))

                    # @constraint(acopf,(min_vol_limit[s,t,idx_nd_nw_buses,lin_itr])^2<=(v_mag[s,t,idx_nd_nw_buses])<=(max_vol_limit[s,t,idx_nd_nw_buses,lin_itr])^2)

                end
################################################################################
################################################################################
            end                # End Time Loop
        end                    # End Scneario Loop
    end                        # End Outer Node Loop
return acopf
end


################################################################################
##------------------------- Script for Network Constraint-----------------------
################################################################################
################################################################################
function network_constraints_no_function_new_milp(acopf,node_data,nw_buses,nBus,nNcurt_gen,nCurt_gen,nFl,nStr_active,bus_data_gsheet,Pg,Qg,Pg_curt,p_strg,p_ch,p_dis,soc,p_fl_inc,p_fl_dec,v_mag,v_ang,v_mag_int,v_ang_int,pg_min,pg_max,qg_min,qg_max,p_load,nw_storage,iStr_active,yij_line,dLines,rdata_buses,trsf_data,tratio,aux_tratio,aux_tap,bin_tap,bin_fl,bin_ch,flex_apc,flex_oltc,flex_adpf,flex_str,flex_fl,fl_bin,str_bin,i_ncurt_gens,i_curt_gens,trsf_fict_nodes,tap_ratio_range,oltc_bin,angle_lb,angle_ub,cycleList,v_ang_node,slack_cnctd_nodes,slack_nd,tratio_aux,oltc_ratio,Qg_curt,max_crnt_limit,rdata_trsf,nTrsf_s1,nw_trsf)
################################################################################
############## Constraint on Transformer Fictitious Node #######################
################################################################################

if flex_oltc == 1 && oltc_bin == 0
    for i in 1:size(trsf_fict_nodes,1)
        nd_num = trsf_fict_nodes[i,1]
        nd_num_fict     = trsf_fict_nodes[i,2]

        idx_nd_nw_buses = findall(x->x==nd_num,rdata_buses[:,1])
        idx_nd_nw_buses = idx_nd_nw_buses[1,1]

        idx_fict_nd_nw_buses = findall(x->x==nd_num_fict,rdata_buses[:,1])
        idx_fict_nd_nw_buses = idx_fict_nd_nw_buses[1,1]

        for s in 1:nSc
            for t in 1:nTP
                v_fbus      = v_mag[s,t,idx_fict_nd_nw_buses]
                v_fbus_fict = v_mag[s,t,nd_num_fict]
                @constraint(acopf,(nw_buses[idx_nd_nw_buses].bus_vmin)^2<=v_fbus_fict<=(nw_buses[idx_nd_nw_buses].bus_vmax).^2)
                @constraint(acopf,v_fbus_fict == tratio[s,t,i]*v_fbus)           # This is a non-linear constraint
            end
        end
    end
elseif flex_oltc == 1 && oltc_bin == 1
    trsf_ratio_range = sqrt.(tap_ratio_range)
    for i in 1:size(trsf_fict_nodes,1)    # Total number of OLTC containing Transformers
        nd_num          = trsf_fict_nodes[i,1]   # Actual From Node
        # nd_num_fict     = trsf_fict_nodes[i,2]
        nd_num_fict     = trsf_fict_nodes[i,4]    # Index of fictitious node

        idx_nd_nw_buses = findall(x->x==nd_num,rdata_buses[:,1]) # Index of Actual From node
        idx_nd_nw_buses = idx_nd_nw_buses[1,1]

        for s in 1:nSc
            for t in 1:nTP
                prod_cb = []    # Product of continuous and binary
                prod_rb = []    # Product of ratio and binary
                sum_bin = []
                v_fbus      = v_mag[s,t,idx_nd_nw_buses]
                v_fbus_fict = v_mag[s,t,nd_num_fict]
                @constraint(acopf,v_fbus_fict<=(nw_buses[idx_nd_nw_buses].bus_vmax).^2)
                @constraint(acopf,(nw_buses[idx_nd_nw_buses].bus_vmin)^2<=v_fbus_fict)

                iter_inner_loop = size(findall(x->x!=0,tap_ratio_range[:,i]),1)  # For total number of non-zero Tap ratio range values
                    # for j = 1:size(tap_ratio_range,1)
                    for j = 1:iter_inner_loop
                        # push!(prod_cb,tap_ratio_range[j,1]*aux_tap[s,t,j,i])
                        push!(prod_cb,tap_ratio_range[j,i]*aux_tap[s,t,j,i])
                        push!(sum_bin,bin_tap[s,t,j,i])
                        # push!(prod_rb,trsf_ratio_range[j,1]*bin_tap[s,t,j,i])
                        push!(prod_rb,trsf_ratio_range[j,i]*bin_tap[s,t,j,i])

                        @constraint(acopf,aux_tap[s,t,j,i]<=(nw_buses[idx_nd_nw_buses].bus_vmax).^2*bin_tap[s,t,j,i])
                        @constraint(acopf,aux_tap[s,t,j,i]>=(nw_buses[idx_nd_nw_buses].bus_vmin).^2*bin_tap[s,t,j,i])
                        @constraint(acopf,v_fbus-aux_tap[s,t,j,i]<=(nw_buses[idx_nd_nw_buses].bus_vmax).^2*(1-bin_tap[s,t,j,i]))
                        @constraint(acopf,v_fbus-aux_tap[s,t,j,i]>=(nw_buses[idx_nd_nw_buses].bus_vmin).^2*(1-bin_tap[s,t,j,i]))
                    end

                    if iter_inner_loop < n_tap
                        fix.(bin_tap[:,:,iter_inner_loop+1:n_tap,i],0)
                        fix.(aux_tap[:,:,iter_inner_loop+1:n_tap,i],0)
                    end
                @constraint(acopf,(v_fbus_fict-sum(prod_cb[j,1] for j in 1:size(prod_cb,1)))==0)
                @constraint(acopf,(sum(sum_bin[j,1] for j in 1:size(sum_bin,1)))==1)
                @constraint(acopf,+((sum(prod_rb[j,1] for j in 1:size(prod_rb,1)))-ratio_init)<=aux_tratio[s,t,i])    # This constraint links the OF variable to the other constraints variables
                @constraint(acopf,-((sum(prod_rb[j,1] for j in 1:size(prod_rb,1)))-ratio_init)<=aux_tratio[s,t,i])    # This constraint links the OF variable to the other constraints variables
            end
        end
    end
end

################################################################################
#################### Transformer Tap Ratio Constraint ##########################
################################################################################
if flex_oltc == 1 && oltc_bin == 0
    for i in 1:size(trsf_data,1)
        for s in 1:nSc
            for t in 1:nTP
                if (trsf_data[i,1]!=0.0) && (trsf_data[i,2]!=0.0) && (trsf_data[i,3]!=0.0)
                    idx = Int64(trsf_data[i,1])
                    @constraint(acopf,(trsf_data[idx,2]).^2<=tratio[s,t,idx]<=(trsf_data[idx,3]).^2)
                    @constraint(acopf,+(tratio_aux[s,t,idx]-ratio_init)<=aux_tratio[s,t,idx])
                    @constraint(acopf,-(tratio_aux[s,t,idx]-ratio_init)<=aux_tratio[s,t,idx])
                    @constraint(acopf,trsf_data[idx,2]<=tratio_aux[s,t,idx]<=trsf_data[idx,3])
                end
            end
        end
    end
end
################################################################################
####################### Non-Curtailable Gen Constraints ########################
################################################################################
for i in 1:nNcurt_gen
    for s in 1:nSc
        for t in 1:nTP
        @constraint(acopf,Pg[s,t,i]<=pg_max[s,t,i_ncurt_gens[i][1]])
        @constraint(acopf,pg_min[s,t,i_ncurt_gens[i][1]]<=Pg[s,t,i])
        @constraint(acopf,Qg[s,t,i]<=qg_max[s,t,i_ncurt_gens[i][1]])
        @constraint(acopf,qg_min[s,t,i_ncurt_gens[i][1]]<=Qg[s,t,i])
        end
    end
end
################################################################################
##### Curtailable Gen Constraints with/without adaptive power factor ###########
################################################################################
if flex_apc == 1
    for i in 1:nCurt_gen
        for s in 1:nSc
            for t in 1:nTP
                if flex_adpf == 0
                    @constraint(acopf,Pg_curt[s,t,i]<=pg_max[s,t,i_curt_gens[i][1]])
                    @constraint(acopf,0<=Pg_curt[s,t,i])
                elseif flex_adpf == 1
                    @constraint(acopf,Pg_curt[s,t,i]<=pg_max[s,t,i_curt_gens[i][1]])
                    @constraint(acopf,0<=Pg_curt[s,t,i])
                    @constraint(acopf,(Qg_curt[s,t,i]+tan(acos(dg_pf))*Pg_curt[s,t,i])<=+tan(acos(dg_pf))*(pg_max[s,t,i_curt_gens[i][1]]))
                    @constraint(acopf,(Qg_curt[s,t,i]-tan(acos(dg_pf))*Pg_curt[s,t,i])>=-tan(acos(dg_pf))*(pg_max[s,t,i_curt_gens[i][1]]))
                    @constraint(acopf,(Qg_curt[s,t,i]>=0))
                end
            end
        end
    end
end
################################################################################
##################### Flexible Loads Constraints Model #########################
################################################################################
if flex_fl == 1
    for i in 1:nFl
        for s in 1:nSc
            for t in 1:nTP
                if fl_bin == 0                                                   # Continuous Model
                    @constraint(acopf,0<=p_fl_inc[s,t,i]<=load_inc_prct*p_load[s,t,iFl[i,1]])
                    @constraint(acopf,0<=p_fl_dec[s,t,i]<=load_dec_prct*p_load[s,t,iFl[i,1]])
                elseif fl_bin == 1                                               # Binary Model
                    @constraint(acopf,0<=p_fl_inc[s,t,i])
                    @constraint(acopf,0<=p_fl_dec[s,t,i])
                    @constraint(acopf,p_fl_inc[s,t,i]<=load_inc_prct*p_load[s,t,iFl[i,1]]*bin_fl[s,t,i])
                    @constraint(acopf,p_fl_dec[s,t,i]<=load_dec_prct*p_load[s,t,iFl[i,1]]*(1-bin_fl[s,t,i]))
                end
            end
            @constraint(acopf,(sum(p_fl_inc[s,t,i] for t in 1:nTP)-sum(p_fl_dec[s,t,i] for t in 1:nTP))==0)
        end
    end
end
################################################################################
################## Storage Constraints (Continuous) ############################
################################################################################
if flex_str == 1
    for i in 1:nStr_active
        const_ch  = time_step*nw_storage[iStr_active[i,1]].storage_ch_eff/(nw_storage[iStr_active[i,1]].storage_e_rat/sbase)
        const_dis = time_step/(nw_storage[iStr_active[i,1]].storage_dis_eff*(nw_storage[iStr_active[i,1]].storage_e_rat/sbase))
        for s in 1:nSc
            for t in 1:nTP
                if str_bin == 0                                                  # Continuous Model
                    @constraint(acopf,-1*(nw_storage[iStr_active[i,1]].storage_ch_rat/sbase)<=p_ch[s,t,i]<=0)
                    @constraint(acopf,0<=p_dis[s,t,i]<=(nw_storage[iStr_active[i,1]].storage_dis_rat/sbase)*1)
                elseif str_bin == 1                                              # Binary Model
                    @constraint(acopf, p_ch[s,t,i]<=0)
                    @constraint(acopf,p_dis[s,t,i]>=0)
                    @constraint(acopf, -bin_ch[s,t,i]<=p_ch[s,t,i]/(nw_storage[i].storage_ch_rat/sbase))
                    @constraint(acopf,p_dis[s,t,i]/(nw_storage[i].storage_dis_rat/sbase)<=(1-bin_ch[s,t,i]))
                end

                if t==1
                    @constraint(acopf,soc[s,t,i]-nw_storage[iStr_active[i,1]].storage_soc_initial==(-const_ch*p_ch[s,t,i])+(-const_dis*p_dis[s,t,i]))
                else
                    @constraint(acopf,soc[s,t,i]-soc[s,t-1,i]==(-const_ch*p_ch[s,t,i])+(-const_dis*p_dis[s,t,i]))
                end

                if t==nTP
                    @constraint(acopf,soc[s,t,i]==nw_storage[iStr_active[i,1]].storage_soc_initial)
                else
                    @constraint(acopf,nw_storage[iStr_active[i,1]].storage_soc_min<=soc[s,t,i])
                    @constraint(acopf,soc[s,t,i]<=nw_storage[iStr_active[i,1]].storage_soc_max)
                end
            end
        end
    end
end
################################################################################
####################### Logitudinal Current Constraint #########################
################################################################################

for i in 1:nLines
    gij_line  = real(yij_line[i,1])
    bij_line  = imag(yij_line[i,1])
    tap       = nw_lines[i].line_tap_ratio
    idx_tap   = Int64(dLines[i,end])
    status_oltc = Int64(rdata_lines[i,end])
    f_bus     = nw_lines[i].line_from
    t_bus     = nw_lines[i].line_to
    idx_t_bus = findall(x->x==t_bus,rdata_buses[:,1])
    idx_f_bus = findall(x->x==f_bus,rdata_buses[:,1])
    idx_f_bus = idx_f_bus[1,1]
    idx_t_bus = idx_t_bus[1,1]

    Imax_from = (nw_lines[i].line_Smax_A*1000)/vbase[idx_f_bus,1]
    Imax_to   = (nw_lines[i].line_Smax_A*1000)/vbase[idx_t_bus,1]
    Imax_lng  = (nw_lines[i].line_Smax_A*1000)/vbase[idx_f_bus,1]

    Imax_fpu  = Imax_from/Ibase[idx_f_bus,1]
    Imax_tpu  = Imax_to/Ibase[idx_t_bus,1]
    Imax_lpu  = Imax_lng/Ibase[idx_f_bus,1]

    Imax_tpu_new = Imax_tpu.^2/(gij_line^2+bij_line^2)
    Imax_lpu_new = Imax_lpu.^2/(gij_line^2+bij_line^2)

    idx_line  = i
    # eta_ub = (nw_buses[idx_f_bus].bus_vmax+nw_buses[idx_t_bus].bus_vmax)/2
    # eta_lb = (nw_buses[idx_f_bus].bus_vmin+nw_buses[idx_t_bus].bus_vmin)/2
    #
    # beta_lb  = 0
    # beta_ub  = ((nw_buses[idx_f_bus].bus_vmax-nw_buses[idx_t_bus].bus_vmin)/2).^2
    for s in 1:nSc
        for t in 1:nTP
            theta_diff_int = v_ang_int[s,t,idx_f_bus]-v_ang_int[s,t,idx_t_bus]
            eta_ij_int     = (v_mag_int[s,t,idx_f_bus]+v_mag_int[s,t,idx_t_bus])/2

            c_int = cos(theta_diff_int)                                          # Cosine argument must be in Radians!

            vf = v_mag[s,t,idx_f_bus]
            vf_int = v_mag_int[s,t,idx_f_bus]
            vt = v_mag[s,t,idx_t_bus]
            vt_int = v_mag_int[s,t,idx_t_bus]

            coeff_vf = (((1+c_int)*vf_int)+((1-3*c_int)*vt_int))/(vf_int+vt_int)
            coeff_vt = (((1+c_int)*vt_int)+((1-3*c_int)*vf_int))/(vf_int+vt_int)
            cons_ij  = c_int*(vf_int-vt_int)^2

            if tap!=0.0
# This is to compare the longitiudinal current constraint with the Imax constraint
                if flex_oltc == 1 && status_oltc == 1
                    from_nd_trsf = findall(x->x==f_bus,trsf_fict_nodes[:,1])
                    to_nd_trsf   = findall(x->x==t_bus,trsf_fict_nodes[:,3])

                    idx_fict_nd  = intersect(from_nd_trsf,to_nd_trsf)
                    idx_fict_nd  = idx_fict_nd[1,1]
                    # fict_nd      = trsf_fict_nodes[idx_fict_nd,2]
                    fict_nd      = trsf_fict_nodes[idx_fict_nd,4]                # Index of fictitious node
                    vf           = v_mag[s,t,fict_nd]
                    vf_int       = v_mag_int[s,t,fict_nd]

                    coeff_vf = (((1+c_int)*vf_int)+((1-3*c_int)*vt_int))/(vf_int+vt_int)
                    coeff_vt = (((1+c_int)*vt_int)+((1-3*c_int)*vf_int))/(vf_int+vt_int)
                    cons_ij  = c_int*(vf_int-vt_int)^2

                elseif (flex_oltc == 0 || flex_oltc == 1) && status_oltc == 0
                    t_tap = tap
                    # Sending End Equation
                    coeff_vf = (((t_tap^4+t_tap^3*c_int)*vf_int)+((t_tap^4-3*t_tap^3*c_int)*vt_int))/(vf_int+vt_int)
                    coeff_vt = (((t_tap^2+t_tap^3*c_int)*vt_int)+((t_tap^2-3*t_tap^3*c_int)*vf_int))/(vf_int+vt_int)
                    cons_ij  = t_tap^3*c_int*(vf_int-vt_int)^2

                    # Receivng End Equation
                    # coeff_vf = (((t_tap^2+t_tap*c_int)*vf_int)+((t_tap^2-3*t_tap*c_int)*vt_int))/(vf_int+vt_int)
                    # coeff_vt = (((1+t_tap*c_int)*vt_int)+((1-3*t_tap*c_int)*vf_int))/(vf_int+vt_int)
                    # cons_ij  = t_tap*c_int*(vf_int-vt_int)^2
                end
####--------------------- Linear constraint ------------------------------######
                # @constraint(acopf,(gij_line^2+bij_line^2)*((coeff_vf*vf)+(coeff_vt*vt)-cons_ij)<=max_crnt_limit[s,t,i]^2) # OLTC Constraints, disabled on May 5th
                # @constraint(acopf,-max_crnt_limit[s,t,i]^2<=(gij_line^2+bij_line^2)*((coeff_vf*vf)+(coeff_vt*vt)-cons_ij))# OLTC Constraints, disabled on May 5th
                # max_crnt_limit is calculated in julia. it is not imported from excel. trace it back to find out the current limits of the OLTC.
            else
                @constraint(acopf,(gij_line^2+bij_line^2)*((coeff_vf*vf)+(coeff_vt*vt)-cons_ij)<=max_crnt_limit[s,t,i]^2)
                @constraint(acopf,-max_crnt_limit[s,t,i]^2<=(gij_line^2+bij_line^2)*((coeff_vf*vf)+(coeff_vt*vt)-cons_ij))
            end
        end
    end
end
################################################################################
############### Constraints for the positivity of Vij ##########################
################################################################################
for i in 1:nBus
    nd      = node_data[i,1].node_num                                            # Node num is independant of time period and scenario
    nd_num  = unique(nd)
    nd_num  = nd_num[1,1]
    idx_nd_nw_buses = findall(x->x==nd_num,rdata_buses[:,1])
    idx_nd_nw_buses = idx_nd_nw_buses[1,1]
    for s in 1:nSc
        for t in 1:nTP
            for j in 1:size(nd,1)
                cnctd_nd     = node_data[i,1].node_cnode[j,1]
                idx_cnctd_nd = findall(x->x==cnctd_nd,rdata_buses[:,1])
                idx_cnctd_nd = idx_cnctd_nd[1,1]

                var_term  = (v_mag[s,t,idx_nd_nw_buses]+v_mag[s,t,idx_cnctd_nd])./2
                var_coeff = (v_mag_int[s,t,idx_nd_nw_buses]-v_mag_int[s,t,idx_cnctd_nd])./(v_mag_int[s,t,idx_nd_nw_buses]+v_mag_int[s,t,idx_cnctd_nd])
                offset    = ((v_mag_int[s,t,idx_nd_nw_buses]-v_mag_int[s,t,idx_cnctd_nd]).^2)./2

                @constraint(acopf,var_term-((var_coeff)*(v_mag[s,t,idx_nd_nw_buses]-v_mag[s,t,idx_cnctd_nd]))+offset>=0)
            end
        end
    end
end
################################################################################
######## Constraints for the line angles to remain in 2pi Radian Limits ########
################################################################################
for i in 1:nLines
    for s in 1:nSc
        for t in 1:nTP
            @constraint(acopf,v_ang[s,t,i]<=angle_ub)
            @constraint(acopf,angle_lb<=v_ang[s,t,i])
        end
    end
end
################################################################################
###### Constraints for sum of the line angles to be zero remain in a loop ######
################################################################################
ft_line = [rdata_lines[:,1] rdata_lines[:,2]]
tf_line = [rdata_lines[:,2] rdata_lines[:,1]]
idx_zero = [0.0 0.0]
for i in 1:size(cycleList,1)                                                     # For all cycles
    list = cycleList[i][1]
        for s in 1:nSc
                for t in 1:nTP
                    cycle_tmp = []
                    for j in 1:size(list,1)
                        if j != length(list)
                            node1 = list[j,1]
                            node2 = list[j+1,1]
                        elseif j == length(list)
                            node1 = list[j,1]
                            node2 = list[j-length(list)+1,1]
                        end
                        line = [node1 node2]
                        id_ft   = line.-ft_line
                        id_tf   = line.-tf_line
                        line_ft = findall(all(id_ft.==idx_zero,dims=2))
                        line_tf = findall(all(id_tf.==idx_zero,dims=2))

                        if !isempty(line_ft) line_ft = line_ft[1][1] end
                        if !isempty(line_tf) line_tf = line_tf[1][1] end

                        if !isempty(line_ft)&&isempty(line_tf)
                            line_angle = @expression(acopf,+1*v_ang[s,t,line_ft])
                        elseif isempty(line_ft)&& !isempty(line_tf)
                            line_angle = @expression(acopf,-1*v_ang[s,t,line_tf])
                        end
                        push!(cycle_tmp,line_angle)
                    end
                    @constraint(acopf,sum(cycle_tmp[j,1] for j in 1:size(cycle_tmp,1))==0)
                end
        end
end
################################################################################
###################### Slack Bus Angle Constraint ##############################
################################################################################
ft_line = [rdata_lines[:,1] rdata_lines[:,2]]
tf_line = [rdata_lines[:,2] rdata_lines[:,1]]
idx_zero = [0.0 0.0]

for i in 2:size(slack_cnctd_nodes,1)
    for s in 1:nSc
        for t in 1:nTP
            line = [slack_nd slack_cnctd_nodes[i,1]]
            id_ft   = line.-ft_line
            id_tf   = line.-tf_line
            line_ft = findall(all(id_ft.==idx_zero,dims=2))
            line_tf = findall(all(id_tf.==idx_zero,dims=2))

            if !isempty(line_ft) line_ft = line_ft[1][1] end
            if !isempty(line_tf) line_tf = line_tf[1][1] end

            if !isempty(line_ft)&&isempty(line_tf)
                line_angle = @expression(acopf,+1*v_ang[s,t,line_ft])
            elseif isempty(line_ft)&& !isempty(line_tf)
                line_angle = @expression(acopf,-1*v_ang[s,t,line_tf])
            end
            @constraint(acopf,v_ang_node[s,t,1]-v_ang_node[s,t,i] == line_angle)

        end
    end
end
####-----------------------------------------------------------------------#####
return acopf
end
################################################################################
############################ Solve OPF Model ###################################
################################################################################
function opf_model_solve(acopf,error_msg,sql_itr)
    if isempty(error_msg)
        # print(acopf)
        status = optimize!(acopf)
    else
        println("The model is not formulated Successfully!")
    end

    term_status = JuMP.termination_status(acopf)                                 # Status of the resulting problem
    if isequal(termination_status(acopf),MOI.OPTIMAL) || isequal(termination_status(acopf),MOI.LOCALLY_SOLVED)
        println("Sql Itr: $sql_itr")
        println("Objective value: ", JuMP.objective_value(acopf))
        println("Solver Time ", JuMP.solve_time(acopf))
        sql_obj  =  JuMP.objective_value(acopf)
        sql_time =  JuMP.solve_time(acopf)

        JuMP.primal_status(acopf)
        JuMP.dual_status(acopf)
        file  = string("results_",filename)
        term_status = 1
    elseif termination_status(acopf) == MOI.TIME_LIMIT && has_values(acopf)
        suboptimal_objective = objective_value(acopf)
        term_status = 2
        sql_obj = 0
        sql_time = 0
    else
        println("PRINTLN: The model has not been solved Successfully!")
        # show("SHOW: The model has not been solved Successfully!")
        term_status = 3
        sql_obj = 0
        sql_time = 0
    end
    return term_status,sql_obj,sql_time
end
################################################################################
##################### Recovering the OPF Solution ##############################
################################################################################
function opf_model_solution(term_statusterm_status,
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
    if term_status == 1 || term_status == 2
################################################################################
############# Recovering Voltage Values (Magnitude & Angle) ####################
################################################################################
   @time vm = JuMP.value.(v_mag)[:,:,:]
   @time va = JuMP.value.(v_ang)[:,:,:]                                          # va recovered should be in Radians!

        vm = sqrt.(vm)                                                           # In linearization program, V^2 is taken as a variable!
        if flex_oltc==1 && (oltc_bin==1 || oltc_bin == 0)
            vm = vm[:,:,1:nBus+nTrsf_s1]
            # vm = vm[:,:,1:nBus+nTrsf]
            va = va[:,:,1:nLines]
        else
            vm = vm[:,:,1:nBus]
            va = va[:,:,1:nLines]
        end

        v_angle_node = angle_recovery(nw_buses,nw_lines,node_data,va,nSc,nTP,nBus,nTrsf,idx_slack_bus_pf,rdata_buses)
        v_angle_fict = zeros(Float64,nSc,nTP,size(trsf_fict_nodes,1))           # Recover the angle of Fictitious Node

        if flex_oltc==1 && (oltc_bin==1 || oltc_bin == 0)
            for i in 1:size(trsf_fict_nodes,1)
                nd_fbus = trsf_fict_nodes[i,1]
                idx_nd_fbus = findall(x->x==nd_fbus,rdata_buses[:,1])
                idx_nd_fbus = idx_nd_fbus[1,1]                                   # Index of Transfomer fictitious node
                    for s in 1:nSc
                        for t in 1:nTP
                            # v_angle_fict[s,t,i] = v_angle_node[s,t,nd_fbus]
                            v_angle_fict[s,t,i] = v_angle_node[s,t,idx_nd_fbus]
                        end
                    end
            end
            v_angle_node = cat(v_angle_node,v_angle_fict,dims=3)
        end

        v_acr = vm.*(cos.(v_angle_node)+sin.(v_angle_node)im)
        vr  = real(v_acr)
        vim = imag(v_acr)
        vm_plr = abs.(v_acr)
        va_plr = angle.(v_acr)

        Pgen_nCurt_pu = JuMP.value.(Pg[:,:,:])
        Qgen_nCurt_pu = JuMP.value.(Qg[:,:,:])
        Pgen_Curt_pu  = JuMP.value.(Pg_curt[:,:,:])

        if flex_str == 1
            p_chr         = JuMP.value.(p_ch[:,:,:]).*sbase
            p_dchr        = JuMP.value.(p_dis[:,:,:]).*sbase
            st_soc        = JuMP.value.(soc[:,:,:])
            p_str         = JuMP.value.(p_strg[:,:,:]).*sbase
            if str_bin == 1
                str_bv = JuMP.value.(bin_ch[:,:,:])
            else
                str_bv = 0.0
            end
        end

        if flex_fl == 1
            pfl_inc   = JuMP.value.(p_fl_inc[:,:,:]).*sbase
            pfl_dec   = JuMP.value.(p_fl_dec[:,:,:]).*sbase
            if fl_bin == 1
                fl_bv = JuMP.value.(bin_fl[:,:,:])
            else
                fl_bv = 0.0
            end
        end

        if flex_adpf == 1
            Qgen_Curt_pu   = JuMP.value.(Qg_curt[:,:,:])
        end

        if flex_oltc == 1
            trsf_tap = JuMP.value.(tratio[:,:,:])                                # Only for those transformers whose status is set to 1.00 and is used when oltc_bin is set to 0
            if oltc_bin == 1
                # ratio_sql   = zeros(Float64,(nSc,nTP,n_tap))
                # oltc_bv     = JuMP.value.(bin_tap[:,:,:])
                # oltc_aux_bv = JuMP.value.(aux_tap[:,:,:])
                ratio_sql   = zeros(Float64,(nSc,nTP,n_tap,nTrsf_s1))            # Everytime redefinition of Ratio_SQL is required as its dimension are changed in the code below.
                oltc_bv     = JuMP.value.(bin_tap[:,:,:,:])
                oltc_aux_bv = JuMP.value.(aux_tap[:,:,:,:])

                for i in 1:size(trsf_fict_nodes,1)
                    for j in 1:size(tap_ratio_range,1)    # Tap raio range is defined only for those transformers whose status is set to 1
                        # ratio_sql[:,:,i] = tap_ratio_range[i,1].*oltc_bv[:,:,i]
                        ratio_sql[:,:,j,i] = tap_ratio_range[j,i].*oltc_bv[:,:,j,i]
                    end
                end
                        ratio_sql = sum(ratio_sql,dims=3)
                        ratio_sql = reshape(ratio_sql,(nSc,nTP,nTrsf_s1))
                        ratio_sql = sqrt.(ratio_sql)

            elseif oltc_bin == 0
                ratio_sql   = sqrt.(trsf_tap)
                oltc_bv     = 0.0
                oltc_aux_bv = 0.0
            else
                oltc_bv     = 0.0
                oltc_aux_bv = 0.0
                ratio_sql   = 0.0
            end
        else
            trsf_tap = 0
        end

        oltc_ratio_tmp = zeros(nSc,nTP,nTrsf)
        idx_trsf_s1    = findall(x->x==1,rdata_trsf[:,end])
        oltc_ratio_tmp = oltc_ratio   #oltc_ratio_tmp is set to the values of oltc_ratio before solving the current SLA
        oltc_ratio_tmp[:,:,idx_trsf_s1] = ratio_sql #oltc_ratio_tmp is updated with the new ratio values only for those transformers whose status is set to 1

################################################################################
######################### Recovering Voltages ##################################
################################################################################
        (vm_ac,va_ac,v_rect,vm_min,vm_max) = recovery_voltages(nTP,nSc,nBus,vr,vim,vm_ac,va_ac,vm_min,vm_max,flex_oltc)
################################################################################
########### Recovering Generator Injection Values (Non-Curtailable) ############
################################################################################
        (gen_Ncurt_P,gen_Ncurt_Q) = recovery_ncurt_gen(nTP,nSc,nNcurt_gen,Pgen_nCurt_pu,Qgen_nCurt_pu,pgen_tol,qgen_tol,gen_Ncurt_P,gen_Ncurt_Q)
        gen_Ncurt_P = gen_Ncurt_P.*sbase
        gen_Ncurt_Q = gen_Ncurt_Q.*sbase
################################################################################
########### Recovering Generator Injection Values (Curtailable) ################
################################################################################
        (gen_Curt_P_sc,gen_Curt_P_gen,total_p_curt,avg_gen_curt) = recovery_curt_gen(nTP,nSc,nCurt_gen,Pgen_Curt_pu,pgen_tol,qgen_tol,gen_Curt_P_sc,gen_Curt_P_gen,total_p_curt,avg_gen_curt)
        gen_Curt_P_sc   = gen_Curt_P_sc*sbase
        gen_Curt_P_gen  = gen_Curt_P_gen*sbase
        total_p_curt    = total_p_curt*sbase
################################################################################
###################### Recovering Storage Values ###############################
################################################################################
        if flex_str == 1
            (p_ch_nd,p_dis_nd,soc_nd,p_strg_nd) = recovery_storage(nTP,nSc,nStr_active,p_chr,p_dchr,st_soc,p_str,p_ch_nd,p_dis_nd,soc_nd,p_strg_nd)# Write code here to store the values of storage values
            #added by Baraa upon request of Croatia
            # @infiltrate
            # any(i-> i!=0, p_ch_nd, dims=2)
            # any(i-> i!=0, p_dis_nd, dims=2)
        else
            p_ch_nd   = zeros(nSc,nTP,nStr_active)
            p_dis_nd  = zeros(nSc,nTP,nStr_active)
            # soc_nd    = zeros(nSc,nTP,nStr_active)
            soc_nd    = 0
            p_strg_nd = 0
            str_bv    = 0
        end
################################################################################
#################### Recovering Flexible Load Values ###########################
################################################################################
        if flex_fl == 1
            (fl_inc_sc,fl_dec_sc,fl_inc_nd,fl_dec_nd) = recovery_flexible(nTP,nSc,nFl,pfl_inc,pfl_dec,fl_inc_sc,fl_dec_sc,fl_inc_nd,fl_dec_nd)
        else
            fl_inc_sc = 0
            fl_dec_sc = 0
            fl_inc_nd = zeros(nSc,nTP,nFl)
            fl_dec_nd = zeros(nSc,nTP,nFl)
            fl_bv     = 0
        end
################################################################################
#################### Recovering Transfomer Tap Values ##########################
################################################################################
        if flex_oltc == 1
            trsf_ratio = recovery_oltc(nTP,nSc,nTrsf,trsf_tap,trsf_ratio)
        else
            trsf_ratio = 0.0
            oltc_bv    = 0.0
            oltc_aux_bv = 0.0
        end
################################################################################
################ Recovering Adaptive Power Factor Values #######################
################################################################################
        if flex_adpf == 1
            gen_Curt_Q = recovery_react_dg(nTP,nSc,nCurt_gen,Qgen_Curt_pu,sbase)
        else
            gen_Curt_Q = recovery_reactive_dg(nTP,nSc,nCurt_gen,dg_pf,pg_max,i_curt_gens,sbase,Pgen_Curt_pu)
        end
################################################################################
#################### Recovering Branch Current Values #########################
################################################################################
        if rcvr_branch == 1
            # line_alpha = JuMP.value.(alpha_ij[:,:,:])
            # line_beta  = JuMP.value.(beta_ij[:,:,:])
            # line_eta  = JuMP.value.(eta_ij[:,:,:])
            (br_crnt_pu,br_crnt_si,br_pwr_ij_pu,br_pwr_ij_SI,br_crnt_ij_pu_ex,br_crnt_ij_SI_ex,br_pwr_ij_pu_ex,br_pwr_ij_SI_ex,app_err) = recovery_branch_current(nSc,nTP,nLines,nw_lines,rdata_buses,vm,v_mag_int,va_plr,v_ang_int,yij_line,dLines,Ibase,br_crnt_pu,br_crnt_si,oltc_ratio,
            br_pwr_ij_pu,br_pwr_ij_SI,br_crnt_ij_pu_ex,br_crnt_ij_SI_ex,br_pwr_ij_pu_ex,br_pwr_ij_SI_ex,app_err)
        else
            br_crnt_pu = 0
            br_crnt_si = 0
            br_pwr_ij_pu = 0
            br_pwr_ij_SI = 0
            br_crnt_ij_pu_ex = 0
            br_crnt_ij_SI_ex = 0
            br_pwr_ij_pu_ex = 0
            br_pwr_ij_SI_ex = 0
            app_err
        end
################################################################################
#################### Recovering Power Injection Values #########################
################################################################################
        if rcvr_inj == 1
            # (p_inj,q_inj) = recovery_injection(nTP,nSc,nBus,p_tol,q_tol,node_data,rdata_buses,flex_oltc,trsf_tap,v_rect)
            # (p_inj,q_inj,p_inj_nl,q_inj_nl) = recovery_injection(rdata_lines,rdata_buses,nSc,nTP,nBus,node_data,flex_oltc,tratio,vm,va,v_mag_int,v_ang_int,trsf_fict_nodes,ratio_sql,p_inj,q_inj,p_inj_nl,q_inj_nl)
            (p_inj,q_inj,p_inj_nl,q_inj_nl) = recovery_injection(rdata_lines,rdata_buses,nSc,nTP,nBus,node_data,flex_oltc,tratio,vm,va,v_mag_int,v_ang_int,trsf_fict_nodes,oltc_ratio_tmp,p_inj,q_inj,p_inj_nl,q_inj_nl)

        else
            p_inj = 0
            q_inj = 0
            p_inj_nl = 0
            q_inj_nl = 0
        end
################################################################################
        if rcvr_tso_dso == 1
            (p_imp_td,q_imp_td) = recovery_tso_dso_flow(Pgen_nCurt_pu,Qgen_nCurt_pu,prob_scs,sbase,p_imp_td,q_imp_td)
        else
            p_imp_td = 0
            q_imp_td = 0
        end

    elseif term_status == 3
        # println("The solution cannot be retrieved as model has not been solved SUCCESSFULLY!")
        println("The solution cannot be retrieved because the model was not solved SUCCESSFULLY!")
        # show("The solution cannot be retrieved because model was not solved SUCCESSFULLY!")
    end

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
        oltc_ratio_tmp

end
