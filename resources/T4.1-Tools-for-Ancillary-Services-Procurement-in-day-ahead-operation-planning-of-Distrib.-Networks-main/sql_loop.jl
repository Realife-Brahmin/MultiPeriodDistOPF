# function sql_loop(vol_viol_max,vol_viol_tol,crnt_viol_max,crnt_viol_tol,lin_itr,lin_itr_max,vol_nodes_mag_pwr_flow,vol_nodes_theta_pf,nSc,nTP,nBus,nNcurt_gen,nCurt_gen,nTrsf,nStr_active,nFl,nd_fl,flex_apc,flex_oltc,
# flex_adpf,flex_str,flex_fl,str_bin,fl_bin,rdata_buses,nLines,trsf_fict_nodes,tap_ratio_range,oltc_bin,prob_scs,time_step,tdata,bus_data_Ssheet,bus_data_lsheet,cost_a_str,cost_b_str,cost_c_str,cost_load_inc,
# cost_load_dec,nw_buses,rdata_loads,node_data,nd_curt_gen,nd_ncurt_gen,p_load,q_load,idx_Gs_lsheet,idx_Bs_lsheet,pg_max,yii_sh,i_curt_gens,sbase,dg_pf,iStr_active,yij_line,dLines,i_ncurt_gens,nw_lines,bus_data_gsheet,
# pg_min,qg_min,qg_max,pgen_tol,qgen_tol,p_tol,q_tol,rcvr_branch,rcvr_inj,rcvr_tso_dso,Ibase,idx_slack_bus_pf,solver_ipopt,solver_cbc,solver_bonmin,sql_itr,sql_itr_max,s_inj_error_rel,s_error_tol,angle_lb,angle_ub,
# cycleList,idx_slack_from,idx_slack_to,num_slack_cnctd_nodes,slack_cnctd_nodes,slack_nd,stoch_model,tratio_init,rdata_storage,load_theta,nd_Str_active,slack_bus_type,vol_ctrl_bus_type,load_bus_type,
# rdata_gens,nw_pPrf_data_load,nw_qPrf_data_load,nw_pPrf_data_gen_max,nw_qPrf_data_gen_max,scenario_data_p_min,scenario_data_p_max,scenario_data_q_min,scenario_data_q_max,nw_buses_pf,nw_lines_pf,nw_loads_pf,nw_gens_pf,
# nw_gcost_pf,nw_sbase_pf,v_initial_pf,v_mag_pf,v_angle_pf,max_mismatch_pf,epsilon,iteration,itr_max,ordata_buses_pf,ybase,vbase,I_rat,vol_viol_nw,vol_viol_nw_max,crnt_viol_nw,crnt_viol_nw_max,min_vol_limit,max_vol_limit,max_crnt_limit,term_status,vol_cstr_tol,br_crnt_si_pf_cnv,solver_cplex,int_tol,opt_tol,num_thread,rdata_trsf,nTrsf_s1,nw_trsf,n_tap,br_pwr_si_pf_cnv)

function sql_loop(InputDict)
vol_viol_max            = InputDict["vol_viol_max"]
vol_viol_tol            = InputDict["vol_viol_tol"]
crnt_viol_max           = InputDict["crnt_viol_max"]
crnt_viol_tol           = InputDict["crnt_viol_tol"]
lin_itr                 = InputDict["lin_itr"]
lin_itr_max             = InputDict["lin_itr_max"]
vol_nodes_mag_pwr_flow  = InputDict["vol_nodes_mag_pf"]
vol_nodes_theta_pf      = InputDict["vol_nodes_theta_pf"]
nSc                     = InputDict["nSc"]
nTP                     = InputDict["nTP"]
nBus                    = InputDict["nBus"]
nNcurt_gen              = InputDict["nNcurt_gen"]
nCurt_gen               = InputDict["nCurt_gen"]
nTrsf                   = InputDict["nTrsf"]
nStr_active             = InputDict["nStr_active"]
nFl                     = InputDict["nFl"]
nd_fl                   = InputDict["nd_fl"]
flex_apc                = InputDict["flex_apc"]
flex_oltc               = InputDict["flex_oltc"]
flex_adpf               = InputDict["flex_adpf"]
flex_str                = InputDict["flex_str"]
flex_fl                 = InputDict["flex_fl"]
str_bin                 = InputDict["str_bin"]
fl_bin                  = InputDict["fl_bin"]
rdata_buses             = InputDict["rdata_buses"]
nLines                  = InputDict["nLines"]
trsf_fict_nodes         = InputDict["trsf_fict_nodes"]
tap_ratio_range         = InputDict["tap_ratio_range"]
oltc_bin                = InputDict["oltc_bin"]
prob_scs                = InputDict["prob_scs"]
time_step               = InputDict["time_step"]
tdata                   = InputDict["tdata"]
bus_data_Ssheet         = InputDict["bus_data_Ssheet"]
bus_data_lsheet         = InputDict["bus_data_lsheet"]
cost_a_str              = InputDict["cost_a_str"]
cost_b_str              = InputDict["cost_b_str"]
cost_c_str              = InputDict["cost_c_str"]
cost_load_inc           = InputDict["cost_load_inc"]
cost_load_dec           = InputDict["cost_load_dec"]
nw_buses                = InputDict["nw_buses"]
rdata_loads             = InputDict["rdata_loads"]
node_data               = InputDict["node_data"]
nd_curt_gen             = InputDict["nd_curt_gen"]
nd_ncurt_gen            = InputDict["nd_ncurt_gen"]
p_load                  = InputDict["p_load"]
q_load                  = InputDict["q_load"]
idx_Gs_lsheet           = InputDict["idx_Gs_lsheet"]
idx_Bs_lsheet           = InputDict["idx_Bs_lsheet"]
pg_max                  = InputDict["pg_max"]
yii_sh                  = InputDict["yii_sh"]
i_curt_gens             = InputDict["i_curt_gens"]
sbase                   = InputDict["sbase"]
dg_pf                   = InputDict["dg_pf"]
iStr_active             = InputDict["iStr_active"]
yij_line                = InputDict["yij_line"]
dLines                  = InputDict["dLines"]
i_ncurt_gens            = InputDict["i_ncurt_gens"]
nw_lines                = InputDict["nw_lines"]
bus_data_gsheet         = InputDict["bus_data_gsheet"]
pg_min                  = InputDict["pg_min"]
qg_min                  = InputDict["qg_min"]
qg_max                  = InputDict["qg_max"]
pgen_tol                = InputDict["pgen_tol"]
qgen_tol                = InputDict["qgen_tol"]
p_tol                   = InputDict["p_tol"]
q_tol                   = InputDict["q_tol"]
rcvr_branch             = InputDict["rcvr_branch"]
rcvr_inj                = InputDict["rcvr_inj"]
rcvr_tso_dso            = InputDict["rcvr_tso_dso"]
Ibase                   = InputDict["Ibase"]
idx_slack_bus_pf        = InputDict["idx_slack_bus_pf"]
solver_ipopt            = InputDict["solver_ipopt"]
solver_cbc              = InputDict["solver_cbc"]
solver_bonmin           = InputDict["solver_bonmin"]
sql_itr                 = InputDict["sql_itr"]
sql_itr_max             = InputDict["sql_itr_max"]
s_inj_error_rel         = InputDict["s_inj_error_rel"]
s_error_tol             = InputDict["s_error_tol"]
angle_lb                = InputDict["angle_lb"]
angle_ub                = InputDict["angle_ub"]
cycleList               = InputDict["cycleList"]
idx_slack_from          = InputDict["idx_slack_from"]
idx_slack_to            = InputDict["idx_slack_to"]
num_slack_cnctd_nodes   = InputDict["num_slack_cnctd_nodes"]
slack_cnctd_nodes       = InputDict["slack_cnctd_nodes"]
slack_nd                = InputDict["slack_nd"]
stoch_model             = InputDict["stoch_model"]
tratio_init             = InputDict["tratio_init"]
rdata_storage           = InputDict["rdata_storage"]
load_theta              = InputDict["load_theta"]
nd_Str_active           = InputDict["nd_Str_active"]
slack_bus_type          = InputDict["slack_bus_type"]
vol_ctrl_bus_type       = InputDict["vol_ctrl_bus_type"]
load_bus_type           = InputDict["load_bus_type"]
rdata_gens              = InputDict["rdata_gens"]
nw_pPrf_data_load       = InputDict["nw_pPrf_data_load"]
nw_qPrf_data_load       = InputDict["nw_qPrf_data_load"]
nw_pPrf_data_gen_max    = InputDict["nw_pPrf_data_gen_max"]
nw_qPrf_data_gen_max    = InputDict["nw_qPrf_data_gen_max"]
scenario_data_p_min     = InputDict["scenario_data_p_min"]
scenario_data_p_max     = InputDict["scenario_data_p_max"]
scenario_data_q_min     = InputDict["scenario_data_q_min"]
scenario_data_q_max     = InputDict["scenario_data_q_max"]
nw_buses_pf             = InputDict["nw_buses_pf"]
nw_lines_pf             = InputDict["nw_lines_pf"]
nw_loads_pf             = InputDict["nw_loads_pf"]
nw_gens_pf              = InputDict["nw_gens_pf"]
nw_gcost_pf             = InputDict["nw_gcost_pf"]
nw_sbase_pf             = InputDict["nw_sbase_pf"]
v_initial_pf            = InputDict["v_initial_pf"]
v_mag_pf                = InputDict["v_mag_pf"]
v_angle_pf              = InputDict["v_angle_pf"]
max_mismatch_pf         = InputDict["max_mismatch_pf"]
epsilon                 = InputDict["epsilon"]
iteration               = InputDict["iteration"]
itr_max                 = InputDict["itr_max"]
ordata_buses_pf         = InputDict["ordata_buses_pf"]
ybase                   = InputDict["ybase"]
vbase                   = InputDict["vbase"]
I_rat                   = InputDict["I_rat"]
vol_viol_nw             = InputDict["vol_viol_nw"]
vol_viol_nw_max         = InputDict["vol_viol_nw_max"]
crnt_viol_nw            = InputDict["crnt_viol_nw"]
crnt_viol_nw_max        = InputDict["crnt_viol_nw_max"]
min_vol_limit           = InputDict["min_vol_limit"]
max_vol_limit           = InputDict["max_vol_limit"]
max_crnt_limit          = InputDict["max_crnt_limit"]
term_status             = InputDict["term_status"]
vol_cstr_tol            = InputDict["vol_cstr_tol"]
br_crnt_si_pf_cnv       = InputDict["br_crnt_si_pf_cnv"]
solver_cplex            = InputDict["solver_cplex"]
int_tol                 = InputDict["int_tol"]
opt_tol                 = InputDict["opt_tol"]
num_thread              = InputDict["num_thread"]
rdata_trsf              = InputDict["rdata_trsf"]
nTrsf_s1                = InputDict["nTrsf_s1"]
nw_trsf                 = InputDict["nw_trsf"]
n_tap                   = InputDict["n_tap"]
br_pwr_si_pf_cnv        = InputDict["br_pwr_si_pf_cnv"]
soc_0                   = InputDict["soc_0"]

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
    line_alpha,
    line_eta,
    line_beta,
    dg_pf_value) = opf_results_variables(nTP,
                                         nSc,
                                         nBus,
                                         nLines,
                                         nNcurt_gen,
                                         nCurt_gen,
                                         nFl,
                                         nStr_active,
                                         nTrsf,
                                         nTrsf_s1,
                                         n_tap)

    max_err_per_bus               =    zeros(nBus,sql_itr_max,lin_itr_max)
    norm_s_err_max_power          =    zeros(nBus,sql_itr_max,lin_itr_max)
    p_err_max_power               =    zeros(nBus,sql_itr_max,lin_itr_max)
    q_err_max_power               =    zeros(nBus,sql_itr_max,lin_itr_max)

    avg_err_per_bus               =    zeros(nBus,sql_itr_max,lin_itr_max)
    norm_s_err_avg                =    zeros(nBus,sql_itr_max,lin_itr_max)
    p_err_avg                     =    zeros(nBus,sql_itr_max,lin_itr_max)
    q_err_avg                     =    zeros(nBus,sql_itr_max,lin_itr_max)

    rel_max_error                 =    zeros(sql_itr_max,lin_itr_max)
    rel_avg_error                 =    zeros(sql_itr_max,lin_itr_max)

    sql_obj                       =    zeros(sql_itr_max,lin_itr_max)
    sql_time                      =    zeros(sql_itr_max,lin_itr_max)

    vol_nodes_mag_pwr_flow_cnv    =    zeros(nSc,nTP,nBus)
    vol_nodes_ang_pf_cnv          =    zeros(nSc,nTP,nBus)

    br_crnt_opf_app               =    zeros(nSc,nTP,nLines)
    br_crnt_opf_ex                =    zeros(nSc,nTP,nLines)

    p_curt_gen                    =    zeros(nSc,nTP,nCurt_gen)
    q_curt_gen                    =    zeros(nSc,nTP,nCurt_gen)
    total_p_curt                  =    0
    str_ch_lin                    =    zeros(nSc,nTP,nStr_active)
    str_dis_lin                   =    zeros(nSc,nTP,nStr_active)
    viol_nodes                    =    0
    vol_viol_cn_milp              =    0
    p_curt_lin                    =    zeros(nTP,lin_itr_max,nCurt_gen)
    oltc_ratio                    =    zeros(Float64,(nSc,nTP,nTrsf)) # For all Transformers irrespective of their status

    opf_term_status               =    0 # opf_term_status plays a role in terminating the while loop. its value is updated as: opf_term_status = infsb_model
    model_infs                    =    []

    p_res_curt_op                 =    zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    q_res_op                      =    zeros(Float64,(length(1:nSc),length(1:nTP),nCurt_gen))
    fl_inc_op                     =    zeros(Float64,(length(1:nSc),length(1:nTP),nFl))
    fl_dec_op                     =    zeros(Float64,(length(1:nSc),length(1:nTP),nFl))
    p_ch_op                       =    zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    p_dis_op                      =    zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    soc_op                        =    zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    p_strg_op                     =    zeros(Float64,(length(1:nSc),length(1:nTP),nStr_active))
    p_slack_pf_ip                 =    zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))
    q_slack_pf_ip                 =    zeros(Float64,(length(1:nSc),length(1:nTP),nNcurt_gen))

    v_rect_euler                  =    zeros(ComplexF64,nSc,nTP,nBus)
    v_rect_cmplx                  =    zeros(ComplexF64,nSc,nTP,nBus)

    if flex_oltc == 0
        v_mag_int_prv_itr    = zeros(nSc,nTP,nBus,lin_itr_max+1)
        v_ang_int_prv_itr    = zeros(nSc,nTP,nBus,lin_itr_max+1)
    elseif flex_oltc == 1
        # v_mag_int_prv_itr    = zeros(nSc,nTP,nBus+nTrsf,lin_itr_max+1)
        v_mag_int_prv_itr    = zeros(nSc,nTP,nBus+nTrsf_s1,lin_itr_max+1)
        v_ang_int_prv_itr    = zeros(nSc,nTP,nBus,lin_itr_max+1)
    end

    if lin_itr == 1
        v_mag_int_prv_itr[:,:,1:nBus,lin_itr] = vol_nodes_mag_pwr_flow
        v_ang_int_prv_itr[:,:,:,lin_itr] = vol_nodes_theta_pf
    end

    while (vol_viol_max > vol_viol_tol || crnt_viol_max > crnt_viol_tol) && (lin_itr <=lin_itr_max) && (opf_term_status == 0)
        println(" ")
        show("Linear Iteration: $lin_itr")
        println(" ")

        (vm_ac,
        va_ac,
        vm_plr,
        va_plr,
        vm_min,
        vm_max,
        gen_Ncurt_P,
        gen_Ncurt_Q,
        gen_Curt_P_sc,
        gen_Curt_P_lin,
        total_p_curt,
        fl_inc_sc,
        fl_dec_sc,
        fl_inc_nd_lin,
        fl_dec_nd_lin,
        p_ch_nd_lin,
        p_dis_nd_lin,
        soc_nd,
        p_strg_nd,
        trsf_ratio,
        gen_Curt_Q_lin,
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
        term_status) = sql_opf_model(vol_nodes_mag_pwr_flow,
                                     vol_nodes_theta_pf,
                                     nSc,
                                     nTP,
                                     nBus,
                                     nNcurt_gen,
                                     nCurt_gen,
                                     nTrsf,
                                     nStr_active,
                                     nFl,
                                     nd_fl,
                                     flex_apc,
                                     flex_oltc,
                                     flex_adpf,
                                     flex_str,
                                     flex_fl,
                                     str_bin,
                                     fl_bin,
                                     rdata_buses,
                                     nLines,
                                     trsf_fict_nodes,
                                     tap_ratio_range,
                                     oltc_bin,
                                     prob_scs,
                                     time_step,
                                     tdata,
                                     bus_data_Ssheet,
                                     bus_data_lsheet,
                                     cost_a_str,
                                     cost_b_str,
                                     cost_c_str,
                                     cost_load_inc,
                                     cost_load_dec,
                                     nw_buses,
                                     rdata_loads,
                                     node_data,
                                     nd_curt_gen,
                                     nd_ncurt_gen,
                                     p_load,
                                     q_load,
                                     idx_Gs_lsheet,
                                     idx_Bs_lsheet,
                                     pg_max,
                                     yii_sh,
                                     i_curt_gens,
                                     sbase,
                                     dg_pf,
                                     iStr_active,
                                     yij_line,
                                     dLines,
                                     i_ncurt_gens,
                                     nw_lines,
                                     bus_data_gsheet,
                                     pg_min,
                                     qg_min,
                                     qg_max,
                                     pgen_tol,
                                     qgen_tol,
                                     p_tol,
                                     q_tol,
                                     rcvr_branch,
                                     rcvr_inj,
                                     rcvr_tso_dso,
                                     Ibase,
                                     idx_slack_bus_pf,
                                     solver_ipopt,
                                     solver_cbc,
                                     solver_bonmin,
                                     sql_itr,
                                     sql_itr_max,
                                     s_inj_error_rel,
                                     s_error_tol,
                                     angle_lb,
                                     angle_ub,
                                     cycleList,
                                     idx_slack_from,
                                     idx_slack_to,
                                     num_slack_cnctd_nodes,
                                     slack_cnctd_nodes,
                                     slack_nd,
                                     stoch_model,
                                     tratio_init,
                                     rdata_storage,
                                     load_theta,
                                     nd_Str_active,
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
                                     sql_obj,
                                     sql_time,
                                     lin_itr,
                                     v_initial_pf,
                                     min_vol_limit,
                                     max_vol_limit,
                                     max_crnt_limit,
                                     oltc_ratio,
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
                                     line_alpha,
                                     line_eta,
                                     line_beta,
                                     dg_pf_value,
                                     v_mag_int_prv_itr,
                                     v_ang_int_prv_itr,
                                     term_status,
                                     solver_cplex,
                                     int_tol,
                                     opt_tol,
                                     num_thread,
                                     rdata_trsf,
                                     nTrsf_s1,
                                     nw_trsf,
                                     n_tap)

        # @infiltrate; throw(error("Debugging breakpoint"))
        # println("Inside WHILE LOOP, lin_itr #$lin_itr")
        # println("opf_term_status is [$opf_term_status]") #(inserted by baraa)
        # println("term_status     is [$term_status]") #(inserted by baraa)
        # println("infsb_model     is [$infsb_model]") #(inserted by baraa)

        p_curt_opf_lin = gen_Curt_P_lin./sbase    # pu
        p_dis_opf_lin  = p_dis_nd_lin./sbase
        p_ch_opf_lin   = abs.(p_ch_nd_lin)./sbase
        p_od_opf_lin   = fl_inc_nd_lin./sbase
        p_ud_opf_lin   = fl_dec_nd_lin./sbase
        str_ch_lin     = p_ch_nd_lin
        str_dis_lin    = p_dis_nd_lin
        q_dg_opf_lin   = gen_Curt_Q_lin./sbase
        trsf_ratio_opf_lin  = oltc_ratio
        p_strg_opf_lin = p_dis_opf_lin-p_ch_opf_lin

        p_curt_tmp = sum(prob_scs.*gen_Curt_P_lin,dims=1)
        p_curt_tmp = reshape(p_curt_tmp,(nTP,nCurt_gen))

        for i in 1:nCurt_gen
            p_curt_lin[:,lin_itr,i] = p_curt_tmp[:,i]
        end

        br_crnt_opf_app = br_crnt_si
        br_crnt_opf_ex  = br_crnt_ij_SI_ex

################################################################################
#################### Power Flow through IPOPT Optimization #####################
################################################################################
        # if term_status == 1 || term_status == 2
        println("AC-Power Flow at SQL Loop Iteration: $lin_itr")
        println("")
        pg_max_pf = pg_max[pf_scenario,:,:];
        pg_min_pf = pg_min[pf_scenario,:,:];
        qg_max_pf = qg_max[pf_scenario,:,:];
        qg_min_pf = pg_min[pf_scenario,:,:];

        pg_max_pf = reshape(pg_max_pf,(nSc_pf_new,nTP_pf_new,nGens));
        pg_min_pf = reshape(pg_min_pf,(nSc_pf_new,nTP_pf_new,nGens));
        qg_max_pf = reshape(qg_max_pf,(nSc_pf_new,nTP_pf_new,nGens));
        qg_min_pf = reshape(qg_min_pf,(nSc_pf_new,nTP_pf_new,nGens));
        prob_scs_pf = prob_scs[pf_scenario,1];

        (vol_nodes_mag_pwr_flow_lin,
        vol_nodes_theta_pf_lin,
        vol_rect_pf_ip,
        p_slack_pf_ip,
        q_slack_pf_ip,
        v_rect_check) = ac_power_flow_opf(nBus,
                                         nLines,
                                         nNcurt_gen,
                                         rdata_buses,
                                         nSc_pf,
                                         nTP_pf,
                                         prob_scs_pf,
                                         time_step,
                                         stoch_model,
                                         p_curt_opf_lin,
                                         p_od_opf_lin,
                                         p_ud_opf_lin,
                                         p_strg_opf_lin,
                                         q_dg_opf_lin,
                                         nw_buses,
                                         nw_lines,
                                         rdata_loads,
                                         bus_data_lsheet,
                                         bus_data_Ssheet,
                                         node_data,
                                         nd_fl,
                                         nd_curt_gen,
                                         nd_ncurt_gen,
                                         p_load,
                                         q_load,
                                         idx_Gs_lsheet,
                                         idx_Bs_lsheet,
                                         yii_sh,
                                         i_curt_gens,
                                         sbase,
                                         flex_oltc,
                                         vbase,
                                         Ibase,
                                         pg_min_pf,
                                         pg_max_pf,
                                         qg_min_pf,
                                         qg_max_pf,
                                         yij_line,
                                         dLines,
                                         i_ncurt_gens,
                                         error_msg,
                                         pgen_tol,
                                         qgen_tol,
                                         trsf_ratio_opf_lin,
                                         vol_cstr_tol);

    (br_crnt_pu_pf_lin,
    br_crnt_si_pf_lin,
    br_pwr_pu_pf_lin,
    br_pwr_si_pf_lin) = recovery_branch_current_pf(nSc,
                                                    nTP,
                                                    nLines,
                                                    nw_lines_pf,
                                                    vol_rect_pf_ip,
                                                    yij_line,
                                                    Ibase,
                                                    sbase,
                                                    ybase,
                                                    vbase,
                                                    oltc_ratio,
                                                    dLines);

    (vol_viol,
    crnt_viol,
    max_vol_viol,
    max_crnt_viol,
    viol_nodes,
    nVol_viol_lin,
    avg_vol_viol_lin,
    above_avg_vol_viol_lin,
    below_avg_vol_viol_lin,
    nCrnt_viol_lin,
    avg_crnt_viol_lin,
    above_avg_crnt_viol_lin,
    below_avg_crnt_viol_lin,
    min_vol_limit,
    max_vol_limit,
    max_crnt_limit,
    vol_viol_cn_milp) = constraint_violation(nBus,
                                            nw_buses,
                                            nLines,
                                            nw_lines,
                                            nTP,
                                            nSc,
                                            vol_nodes_mag_pwr_flow_lin,
                                            br_crnt_si_pf_lin,
                                            I_rat,
                                            lin_itr,
                                            min_vol_limit,
                                            max_vol_limit,
                                            max_crnt_limit,
                                            vol_viol_tol,
                                            crnt_viol_tol,
                                            vol_cnst_fctr,
                                            crnt_cnst_fctr,
                                            lin_itr_max);

############## PoL update for next Successive Linear Iteration  ################

        vol_nodes_mag_pwr_flow   = deepcopy(vol_nodes_mag_pwr_flow_lin);
        vol_nodes_theta_pf       = deepcopy(vol_nodes_theta_pf_lin);

        if flex_oltc == 0
            v_mag_int_prv_itr[:,:,:,lin_itr+1] = deepcopy(vol_nodes_mag_pwr_flow_lin);
            v_ang_int_prv_itr[:,:,:,lin_itr+1] = deepcopy(vol_nodes_theta_pf_lin);

        elseif flex_oltc==1
            v_mag_int_prv_itr[:,:,1:nBus,lin_itr+1] = deepcopy(vol_nodes_mag_pwr_flow_lin);
            v_ang_int_prv_itr[:,:,1:nBus,lin_itr+1] = deepcopy(vol_nodes_theta_pf_lin);
            for i in 1:size(trsf_fict_nodes,1)
                trsf_nd_from = trsf_fict_nodes[i,1];
                idx_trsf_from = findall(x->x==trsf_nd_from,rdata_buses[:,1]);
                idx_trsf_from = idx_trsf_from[1,1];

                # vol_from  = vol_nodes_mag_pwr_flow_lin[:,:,trsf_fict_nodes[i,1]]
                # v_mag_int_prv_itr[:,:,trsf_fict_nodes[i,2],lin_itr+1] = vol_from

                vol_from  = vol_nodes_mag_pwr_flow_lin[:,:,idx_trsf_from];
                v_mag_int_prv_itr[:,:,trsf_fict_nodes[i,4],lin_itr+1] = vol_from;
            end
        end

        vol_nodes_mag_pwr_flow_cnv = vol_nodes_mag_pwr_flow_lin;
        vol_nodes_ang_pf_cnv = vol_nodes_theta_pf_lin;
        br_crnt_si_pf_cnv    = br_crnt_si_pf_lin;
        br_pwr_si_pf_cnv     = br_pwr_si_pf_lin;

        p_curt_gen = gen_Curt_P_lin;
        q_curt_gen = gen_Curt_Q_lin;

        min_vol_limit  = min_vol_limit;
        max_vol_limit  = max_vol_limit;
        max_crnt_limit = max_crnt_limit;

        v_rect_euler = vol_rect_pf_ip;
        v_rect_cmplx = v_rect_check;

        if !isempty(max_vol_viol)
            vol_viol_max  = maximum(max_vol_viol);
            vol_viol_less_tol = size(findall(x->x<=1,max_vol_viol),1);
        else
            vol_viol_max = 0.001;
            vol_viol_less_tol = 0.0;
        end

        if !isempty(max_crnt_viol)
            crnt_viol_max      = maximum(max_crnt_viol);
            crnt_viol_less_tol = size(findall(x->x<=1,max_crnt_viol),1);
        else
            crnt_viol_max = 0.001;
            crnt_viol_less_tol = 0.0;
        end

        push!(vol_viol_nw,vol_viol);
        push!(crnt_viol_nw,crnt_viol);
        vol_viol_data  = [nVol_viol_lin  vol_viol_max  avg_vol_viol_lin  above_avg_vol_viol_lin  below_avg_vol_viol_lin vol_viol_less_tol];
        crnt_viol_data = [nCrnt_viol_lin crnt_viol_max avg_crnt_viol_lin above_avg_crnt_viol_lin below_avg_crnt_viol_lin crnt_viol_less_tol];
        vol_viol_nw_max[lin_itr,:] = vol_viol_data;
        crnt_viol_nw_max[lin_itr,:] = crnt_viol_data;
        lin_itr = lin_itr+1;
        opf_term_status = infsb_model;
        push!(model_infs,(lin_itr-1,infsb_model));

        p_res_curt_op  = p_curt_opf_lin.*sbase;
        q_res_op       = q_dg_opf_lin.*sbase;
        fl_inc_op      = p_od_opf_lin.*sbase;
        fl_dec_op      = p_ud_opf_lin.*sbase;
        p_ch_op        = str_ch_lin.*sbase;
        p_dis_op       = str_dis_lin.*sbase;
        p_strg_op      = p_strg_opf_lin.*sbase;

        # println("INSIDE WHILE LOOP (SQL_LOOP), lin_itr #$lin_itr")
        # println("term_status     is [$term_status]") #(inserted by baraa)
    end    # End While Loop
    # println("OUTSIDE WHILE LOOP (SQL_LOOP)")
    # println("opf_term_status is [$opf_term_status]") #(inserted by baraa)
    # println("term_status     is [$term_status]") #(inserted by baraa)

    oltc_ratio = oltc_ratio.^-1;

    OutDict = Dict("vol_viol_nw"		        => vol_viol_nw	,
                   "vol_viol_nw_max"		    => vol_viol_nw_max	,
                   "crnt_viol_nw"		        => crnt_viol_nw	,
                   "crnt_viol_nw_max"		    => crnt_viol_nw_max	,
                   "vol_nodes_mag_pwr_flow_cnv"	=> vol_nodes_mag_pwr_flow_cnv	,
                   "vol_nodes_ang_pf_cnv"		=> vol_nodes_ang_pf_cnv	,
                   "br_crnt_si_pf_cnv"		    => br_crnt_si_pf_cnv	,
                   "sql_obj"		            => sql_obj	,
                   "sql_time"		            => sql_time	,
                   "rel_max_error"		        => rel_max_error	,
                   "rel_avg_error"		        => rel_avg_error	,
                   "norm_s_err_max_power"		=> norm_s_err_max_power	,
                   "norm_s_err_avg"		        => norm_s_err_avg	,
                   "p_curt_gen"		            => p_curt_gen	,
                   "q_curt_gen"		            => q_curt_gen	,
                   "total_p_curt"		        => total_p_curt	,
                   "str_dis_lin"		        => str_dis_lin	,
                   "str_ch_lin"		            => str_ch_lin	,
                   "str_soc"		            => soc_nd	    , # added by baraa
                   "viol_nodes"		            => viol_nodes	,
                   "p_curt_lin"		            => p_curt_lin	,
                   "min_vol_limit"		        => min_vol_limit	,
                   "max_vol_limit"		        => max_vol_limit	,
                   "oltc_ratio"		            => oltc_ratio	,
                   "p_err_max_power"		    => p_err_max_power	,
                   "q_err_max_power"		    => q_err_max_power	,
                   "vm_ac"		                => vm_ac	,
                   "va_ac"		                => va_ac	,
                   "vol_viol_cn_milp"		    => vol_viol_cn_milp	,
                   "model_infs"		            => model_infs	,
                   "br_crnt_opf_app"		    => br_crnt_opf_app	,
                   "br_crnt_opf_ex"		        => br_crnt_opf_ex	,
                   "p_res_curt_op"		        => p_res_curt_op	,
                   "q_res_op"		            => q_res_op	,
                   "fl_inc_op"		            => fl_inc_op	,
                   "fl_dec_op"		            => fl_dec_op	,
                   "p_ch_op"		            => p_ch_op	,
                   "p_dis_op"		            => p_dis_op	,
                   "p_strg_op"		            => p_strg_op	,
                   "p_slack_pf_ip"		        => p_slack_pf_ip	,
                   "q_slack_pf_ip"		        => q_slack_pf_ip	,
                   "br_pwr_si_pf_cnv"		    => br_pwr_si_pf_cnv	,
                   "v_rect_euler"		        => v_rect_euler	,
                   "v_rect_cmplx"		        => v_rect_cmplx	,
                   "term_status"                => term_status,
                #    "opf_term_status"            => opf_term_status, # see above: opf_term_status = infsb_model
                   )
    return OutDict

    # return vol_viol_nw,vol_viol_nw_max,crnt_viol_nw,crnt_viol_nw_max,vol_nodes_mag_pwr_flow_cnv,vol_nodes_ang_pf_cnv,br_crnt_si_pf_cnv,sql_obj,
    # sql_time,rel_max_error,rel_avg_error,norm_s_err_max_power,norm_s_err_avg,p_curt_gen,q_curt_gen,total_p_curt,str_dis_lin,str_ch_lin,
    # viol_nodes,p_curt_lin,min_vol_limit,max_vol_limit,oltc_ratio,p_err_max_power,q_err_max_power,vm_ac,va_ac,vol_viol_cn_milp,model_infs,br_crnt_opf_app,br_crnt_opf_ex,
    # p_res_curt_op,q_res_op,fl_inc_op,fl_dec_op,p_ch_op,p_dis_op,p_strg_op,p_slack_pf_ip,q_slack_pf_ip,br_pwr_si_pf_cnv,v_rect_euler,v_rect_cmplx
end
