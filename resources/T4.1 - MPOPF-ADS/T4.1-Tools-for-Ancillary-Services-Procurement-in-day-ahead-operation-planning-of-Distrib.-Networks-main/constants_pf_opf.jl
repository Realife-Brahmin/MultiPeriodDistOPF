################################################################################
##---------- Constants to be used in the AC Power Flow Model -------------------
################################################################################
nw_buses_pf = array_bus
nw_lines_pf = array_lines
nw_loads_pf = array_loads
nw_gens_pf  = array_gens
nw_gcost_pf = array_gcost
nw_sbase_pf = array_sbase

obus_type_pf    = rdata_buses[:,2]
ordata_buses_pf = obus_type_pf

nTP_pf = nTP
nSc_pf = nSc

load_bus_type     = 1
vol_ctrl_bus_type = 2
slack_bus_type    = 3

idx_slack_bus_pf = findall(x->x==slack_bus_type,rdata_buses[:,2])
nd_num     = rdata_buses[idx_slack_bus_pf[1,1],1]
idx_nd_gen = findall(x->x==nd_num,rdata_gens[:,1])
sp_vol = nw_gens_pf[idx_nd_gen[1,1]].gen_V_set

if nTP == 24
    v_initial_pf = rdata_vProfile[:,2:end]
else
    v_initial_pf = sp_vol   # Setup constant voltage at slack bus in AC-OPF and AC-PF through IPOPT in all scenarios and time-periods
end

v_initial = 1.00
v_sp_pf        = v_initial.*exp(1.0im*0*pi/180)
node_vol_pf    = v_sp_pf
node_vol_pf    = fill(v_sp_pf,(size(nw_buses_pf,1),nTP_pf,nSc_pf))
v_mag_pf       = abs.(node_vol_pf) # Setting up initial voltage at all buses in AC PF through NR method. The voltage at slack bus is then updated in ac_power_flow_model_check by voltage mentioned in Gens sheet
v_angle_pf     = angle.(node_vol_pf)

################################################################################
##------------- Constants to be used in AC MILP-OPF program --------------------
################################################################################
nw_buses   = array_bus
nw_lines   = array_lines
nw_loads   = array_loads
nw_gens    = array_gens
nw_gcost   = array_gcost
nw_sbase   = array_sbase
nw_storage = array_storage
nw_Strcost = array_Strcost
nw_trsf    = array_trsf

#--------------- Active and reactive power load profiles -----------------------
nw_pPrf_header_load     = rheader_pProfile_load
nw_qPrf_header_load     = rheader_qProfile_load
nw_pPrf_data_load       = rdata_pProfile_load
nw_qPrf_data_load       = rdata_qProfile_load
#--------- Active and reactive power generation profiles for F-MP-OPF ----------
nw_pPrf_header_gen_min  = rheader_pProfile_gen_min                               # All the following data is in MW/MVAr
nw_pPrf_header_gen_max  = rheader_pProfile_gen_max
nw_qPrf_header_gen_min  = rheader_qProfile_gen_min
nw_qPrf_header_gen_max  = rheader_qProfile_gen_max

nw_pPrf_data_gen_min    = rdata_pProfile_gen_min                                 # All the following data is in MW/MVAr
nw_pPrf_data_gen_max    = rdata_pProfile_gen_max
nw_qPrf_data_gen_min    = rdata_qProfile_gen_min
nw_qPrf_data_gen_max    = rdata_qProfile_gen_max
#------------------------- Network Constants -----------------------------------
sbase = nw_sbase[1].sbase                                                        # sbase is in MV
function vol_base(rdata_buses,nw_buses)
    vbase = []
    for i in 1:size(rdata_buses,1)
        push!(vbase,nw_buses[i].bus_vnom)
    end
    return vbase
end
vbase = vol_base(rdata_buses,nw_buses)                                           # vbase is in kVA
Ibase = (sbase*1000)./vbase                                                      # base_current = base kVA/base kV. Current base is in Amperes
zbase = (vbase).^2/sbase
ybase = zbase.^-1

nSc_pf_new = nSc
nTP_pf_new = nTP
pf_scenario = collect(1:nSc_pf_new)

##################### Bound Tightening Factors #######################
vol_cnst_fctr = 0.05/100       # 0.5%
crnt_cnst_fctr = 0
#####################################################################
term_status     = 0   # To start the outer while-loop of SLA
vol_viol_max    = 5   # To start the outer while-loop of SLA
crnt_viol_max   = 5   # To start the outer while-loop of SLA
s_inj_error_rel = 10  # To start the inner while-loop of SLA
s_error_tol     = 1e-5    # tol is in pu: New value of tolerance
time_step = 24/size(nw_pPrf_data_load[:,2:end],2)
load_theta = tan(acos(pf))                                                       # Used in Flexible Load
dg_ac_pf = dg_pf


cp_viol  = 100
pgen_tol = 1e-3
qgen_tol = 1e-3
p_tol    = 1e-4
q_tol    = 1e-4
eps_cmp  = 1e-2
angle_degree = 30
angle_radian = deg2rad(30)
angle_ub  = angle_radian
angle_lb  = -angle_radian
numCycles = 0
cycles    = []
#------------------------- Scenario profiles -----------------------------------
nw_pPrf_header_wind_sc = rheader_pProfile_wind_sc
nw_pPrf_data_wind_sc   = rdata_pProfile_wind_sc

nw_pPrf_header_pv_sc   = rheader_pProfile_pv_sc
nw_pPrf_data_pv_sc     = rdata_pProfile_pv_sc

nw_prob_header_sc      = rheader_prob_sc
nw_prob_data_sc        = rdata_prob_sc

##----------------- Creation of a scenario matrix------------------------------
sbase = nw_sbase[1].sbase                                                        # sbase is in MVA
function scenario_data(nw_gens,nw_pPrf_data_pv_sc,nw_pPrf_data_wind_sc,sbase,nSc,nTP,nGens)
    scenario_data_p_min  = zeros(nSc,nTP,nGens)                                  # nGens = nCurt_gen + nNon_curt_gen
    scenario_data_p_max  = zeros(nSc,nTP,nGens)
    scenario_data_q_min  = zeros(nSc,nTP,nGens)
    scenario_data_q_max  = zeros(nSc,nTP,nGens)
    for i in 1:nGens
        dg_type = nw_gens[i].gen_type
        for j in 1:nSc
            if dg_type==0       # Ref bus
                scenario_data_p_min[j,:,i] .= (nw_gens[i].gen_P_min)/sbase       # Scenarios are arranged in NODE ORDER which appears in 'Gens' sheet
                scenario_data_p_max[j,:,i] .= (nw_gens[i].gen_P_max)/sbase
                scenario_data_q_min[j,:,i] .= (nw_gens[i].gen_Qg_min)/sbase
                scenario_data_q_max[j,:,i] .= (nw_gens[i].gen_Qg_max)/sbase
            elseif dg_type==1   # PV DG
                scenario_data_p_min[j,:,i] .= (nw_gens[i].gen_P_min*0)/sbase
                scenario_data_p_max[j,:,i]  = 1.0*(nw_gens[i].gen_P_max*nw_pPrf_data_pv_sc[j,1:nTP])/sbase     # Scenarios given in 'Scenarios sheet' are in MW of 1MW installed capacity
                scenario_data_q_min[j,:,i] .= (nw_gens[i].gen_Qg_min*0)/sbase    # So, for solar DG bigger or smaller than 1 MW, it is needed to scale up/down the values by multiplying their rated capacity.
                scenario_data_q_max[j,:,i] .= (nw_gens[i].gen_P_max*tan(acos(dg_pf))*nw_pPrf_data_pv_sc[j,1:nTP])/sbase
            elseif dg_type==2  # WIND DG
                scenario_data_p_min[j,:,i] .= (nw_gens[i].gen_P_min*0)/sbase
                scenario_data_p_max[j,:,i]  = 1.5*(nw_gens[i].gen_P_max*nw_pPrf_data_wind_sc[j,1:nTP])/sbase
                scenario_data_q_min[j,:,i] .= (nw_gens[i].gen_Qg_min*0)/sbase
                scenario_data_q_max[j,:,i] .= (nw_gens[i].gen_P_max*tan(acos(dg_pf))*nw_pPrf_data_wind_sc[j,1:nTP])/sbase
            # else
                # Further DGs can be added here!
            end
        end
    end
return scenario_data_p_min, scenario_data_p_max, scenario_data_q_min, scenario_data_q_max
end
(scenario_data_p_min, scenario_data_p_max, scenario_data_q_min, scenario_data_q_max) = scenario_data(nw_gens,nw_pPrf_data_pv_sc,nw_pPrf_data_wind_sc,sbase,nSc,nTP,nGens)



## ------------------------- Scenario probabilities ----------------------------
prob_scs = rdata_prob_sc[1:nSc,1]

#----------------------------- ACOPF Model--------------------------------------
##------------------------------------------------------------------------------
node_data      = []
error_msg      = []
yii_sh         = []

##----------------------------- Rated Bracnh Current ------------------------###
I_rat = []
I_rat_pu = []
max_crnt_limit = zeros(nSc,nTP,nLines)
max_crnt_limit_si = zeros(nSc,nTP,nLines)
for i in 1:nLines
    f_bus     = nw_lines[i].line_from
    t_bus     = nw_lines[i].line_to
    tap       = nw_lines[i].line_tap_ratio
    idx_f_bus = findall(x->x==f_bus,rdata_buses[:,1])
    idx_t_bus = findall(x->x==t_bus,rdata_buses[:,1])

    idx_f_bus = idx_f_bus[1,1]
    idx_t_bus = idx_t_bus[1,1]

    Imax_from = (nw_lines[i].line_Smax_A*1000)/vbase[idx_f_bus,1]
    Imax_to   = (nw_lines[i].line_Smax_A*1000)/vbase[idx_t_bus,1]
    Imax_lng  = (nw_lines[i].line_Smax_A*1000)/vbase[idx_f_bus,1]

    Imax_from_pu = Imax_to/Ibase[idx_f_bus,1]
    Imax_to_pu = Imax_to/Ibase[idx_t_bus,1]
    Imax_lng_pu = Imax_lng/Ibase[idx_f_bus,1]

    if tap !=0 # current limits of the OLTC. comment added by baraa on May 5th
        push!(I_rat,Imax_from)
        push!(I_rat_pu,Imax_from_pu)
        max_crnt_limit[:,:,i] .= Imax_from_pu
        max_crnt_limit_si[:,:,i] .= Imax_from

        # push!(I_rat,Imax_to)
        # push!(I_rat_pu,Imax_to_pu)
        # max_crnt_limit[:,:,i] .= Imax_to_pu
        # max_crnt_limit_si[:,:,i] .= Imax_to
    else # current limits for a passive branch. comment added by baraa on May 5th
        push!(I_rat,Imax_lng)
        push!(I_rat_pu,Imax_lng_pu)
        max_crnt_limit[:,:,i] .= Imax_lng_pu
        max_crnt_limit_si[:,:,i] .= Imax_lng
    end
end

min_vol_limit = zeros(nSc,nTP,nBus,lin_itr_max)
max_vol_limit = zeros(nSc,nTP,nBus,lin_itr_max)
for i in 1:nBus
    min_vol_limit[:,:,i,:] .= nw_buses[i].bus_vmin
    max_vol_limit[:,:,i,:] .= nw_buses[i].bus_vmax
end

##------------------------------------------------------------------------------                                                        # Shunt elements connected to a node
################################################################################
################################################################################
vol_viol_nw = []
crnt_viol_nw = []
vol_viol_nw_max = zeros(lin_itr_max,6)
crnt_viol_nw_max = zeros(lin_itr_max,6)

vol_viol  = []
crnt_viol = []
max_vol_viol = []
max_crnt_viol = []


rcvr_branch   = 1
rcvr_inj      = 1
rcvr_tso_dso  = 0

if scenario == 0
    stoch_model   = 0
elseif scenario == 1
    stoch_model = 1
end
