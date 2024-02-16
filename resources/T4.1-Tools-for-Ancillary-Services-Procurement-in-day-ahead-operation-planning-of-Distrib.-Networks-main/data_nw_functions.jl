#------------------- Formatting of Lines Data ----------------------------------
#-------------------------------------------------------------------------------
function data_lines(array_lines,nw_lines,nLines)
    idx_from_line  = Float64[]
    idx_to_line    = Float64[]
    yij_line       = ComplexF64[]
    yij_line_sh    = ComplexF64[]
    tap_ratio      = Float64[]
    tap_ratio_max  = Float64[]
    tap_ratio_min  = Float64[]
    for i in 1:nLines                                                                # nLines = nTrsf + nTransmission lines
        push!(idx_from_line,array_lines[i].line_from)                                # Saving 'from' end of lines in a vector
        push!(idx_to_line,array_lines[i].line_to)                                    # Saving 'to' end of lines in a vector
        push!(yij_line,inv(nw_lines[i].line_r+(nw_lines[i].line_x)im))               # Line admittance calculated from the given r and x values
        push!(yij_line_sh,nw_lines[i].line_g_shunt+(nw_lines[i].line_b_shunt)im)     # Shunt line admittance
        push!(tap_ratio,nw_lines[i].line_tap_ratio)
        push!(tap_ratio_min,nw_lines[i].line_tap_ratio_min)
        push!(tap_ratio_max,nw_lines[i].line_tap_ratio_max)
    end
    nTrsf  = size(findall(x->x!=0.0,tap_ratio),1)
    dLines = [idx_from_line idx_to_line 1:nLines tap_ratio]                          # dLines = data_Lines
return  idx_from_line, idx_to_line, yij_line,yij_line_sh,tap_ratio,tap_ratio_min,tap_ratio_max,nTrsf,dLines
end
###############################################################################
function data_trsf(dLines)
    idx_tap = Float64[]
    ctr = 0
    for i in 1:size(dLines,1)
        tap = dLines[i,4]
        if tap !=0.0                                                             # Transformer branch
            ctr = ctr+1
            push!(idx_tap,ctr)
        elseif tap ==0.0                                                         # Transmission/Distribution line
            push!(idx_tap,0)
        end
    end
return idx_tap
end
# ################# Old Transformer Fictitious Nodes Program ###################
# function fictitious_trsf_node(tap_ratio,rdata_buses,nw_lines,n_tap,tdata)
#     idx_trsf_lines  = findall(x->x!=0.0,tap_ratio)
#     trsf_fict_nodes = zeros(Int64,(size(idx_trsf_lines,1),3))
#     max_nd_num      = Int64(maximum(rdata_buses[:,1]))
#     tap_ratio_range = zeros(n_tap,size(tdata,1))
#
#     for i in 1:size(idx_trsf_lines,1)
#         trsf_fict_nodes[i,1] = nw_lines[idx_trsf_lines[i,1]].line_from               # From_Old
#         trsf_fict_nodes[i,2] = max_nd_num+1                                          # From_New
#         trsf_fict_nodes[i,3] = nw_lines[idx_trsf_lines[i,1]].line_to                 # To
#         max_nd_num = max_nd_num+1
#
#     end
#     # for i in 1:size(tdata,1)
#     #     step_ratio = (tdata[i,3].^2-tdata[i,2].^2)/n_tap
#     #     for j in 1:n_tap+1
#     #         tap_ratio_range[j,i] = tdata[i,2].^2+(j-1)*step_ratio
#     #     end
#     # end
#
#     for i in 1:size(tdata,1)
#         step_ratio = (tdata[i,3]-tdata[i,2])/(n_tap-1)
#         ratio = collect(tdata[i,2]:step_ratio:tdata[i,3])
#         tap_ratio_range[:,i] = ratio.^2
#     end
#
#     tratio_init = zeros(Float64,(nSc,nTP,nTrsf))
#     for i in 1:size(tdata,1)
#         for s in 1:nSc
#             for t in 1:nTP
#                 tratio_init[s,t,i] = tdata[i,4]
#             end
#         end
#     end
# return trsf_fict_nodes,tap_ratio_range,tratio_init
# end

# ################# New Transformer Fictitious Nodes Program ###################
function fictitious_trsf_node_new(tap_ratio,rdata_buses,tdata,nw_trsf,rdata_trsf)

    # idx_trsf_lines  = findall(x->x!=0.0,tap_ratio)
    # trsf_fict_nodes = zeros(Int64,(size(idx_trsf_lines,1),3))
    # max_nd_num      = Int64(maximum(rdata_buses[:,1]))
    # tap_ratio_range = zeros(n_tap,size(tdata,1))
if flex_oltc == 1 && oltc_bin == 1

    idx_trsf_lines  = findall(x->x!=0.0,rdata_trsf[:,end])
    # trsf_fict_nodes = zeros(Int64,(size(idx_trsf_lines,1),3))
    trsf_fict_nodes = zeros(Int64,(size(idx_trsf_lines,1),4))
    max_nd_num      = Int64(maximum(rdata_buses[:,1]))
    n_tap           = Int64(maximum(rdata_trsf[idx_trsf_lines,6]))
    tap_ratio_range = zeros(n_tap,size(idx_trsf_lines,1))

    for i in 1:size(idx_trsf_lines,1)
        trsf_fict_nodes[i,1] = nw_trsf[idx_trsf_lines[i,1]].trsf_from               # From_Old
        trsf_fict_nodes[i,2] = max_nd_num+1                                          # From_New
        trsf_fict_nodes[i,3] = nw_trsf[idx_trsf_lines[i,1]].trsf_to                 # To
        trsf_fict_nodes[i,4] = size(rdata_buses,1)+i
        max_nd_num = max_nd_num+1
        # global max_nd_num = max_nd_num+1
    end

    # for i in 1:size(tdata,1)
    #     step_ratio = (tdata[i,3].^2-tdata[i,2].^2)/n_tap
    #     for j in 1:n_tap+1
    #         tap_ratio_range[j,i] = tdata[i,2].^2+(j-1)*step_ratio
    #     end
    # end

    for i in 1:size(idx_trsf_lines,1)
        # step_ratio = (tdata[i,3]-tdata[i,2])/(n_tap-1)
        step_ratio = rdata_trsf[idx_trsf_lines[i,1],3]
        ratio = collect(rdata_trsf[idx_trsf_lines[i,1],17]:step_ratio:rdata_trsf[idx_trsf_lines[i,1],18])
        tap_ratio_range[1:size(ratio,1),i] = ratio.^2
    end
else
    trsf_fict_nodes = 0
    tap_ratio_range = 0
    n_tap = 1
end

    tratio_init = zeros(Float64,(nSc,nTP,nTrsf))
    for i in 1:size(tdata,1)
        for s in 1:nSc
            for t in 1:nTP
                # @printf("i=%3i/%i,\ts=%3i/%i,\t t=%3i/%i\n", i, size(tdata,1), s, nSc, t, nTP) # added by baraa for debugging
                tratio_init[s,t,i] = tdata[i,4]
                # comment on error by baraa: tdata has 2 rows (size(tdata,1)=2), while nTrsf=1, so tratio_init has 1 row.
            end
        end
    end

return trsf_fict_nodes,tap_ratio_range,tratio_init,n_tap
end
###############################################################################
function data_nodes(nBus,nw_buses,idx_from_line,idx_to_line,yij_line,yij_line_sh,tap_ratio,tap_ratio_min,tap_ratio_max,idx_tap,node)
    node_data  = node[]
    # tdata = []
    for i in 1:nBus
        bus_num = nw_buses[i].bus_num
        ft_line = findall(x->x==bus_num,idx_from_line)                               # index of all lines which are connected to bus i (bus i is present in 'from' column)
        tf_line = findall(x->x==bus_num,idx_to_line)                                 # index of all lines which are connected to bus i (bus i is present in 'to' column)
        telem   = size(ft_line,1)+size(tf_line,1)

        ft_bus  = idx_to_line[ft_line]                                                # buses connected to bus i (bus i is present in 'from' column)
        tf_bus  = idx_from_line[tf_line]                                              # buses connected to bus i (bus i is present in 'to' column)

        b_lines  = union(ft_line,tf_line)                                            # indexes of all lines connected to bus (b) i
        b_cbuses = union(ft_bus,tf_bus)                                              # buses connected to bus i

        gij_line = real(yij_line[b_lines,1])
        bij_line = imag(yij_line[b_lines,1])

        gij_line_sh = real(yij_line_sh[b_lines,1])
        bij_line_sh = imag(yij_line_sh[b_lines,1])

        tp_rt     = vcat(tap_ratio[ft_line],tap_ratio[tf_line])                      # Tap ratio
        tp_rt_min = vcat(tap_ratio_min[ft_line],tap_ratio_min[tf_line])              # Tap ratio minimum
        tp_rt_max = vcat(tap_ratio_max[ft_line],tap_ratio_max[tf_line])              # Tap ratio maximum
        from_col  = idx_from_line[b_lines,1]
        to_col    = idx_to_line[b_lines,1]
        bus       = repeat([bus_num],telem)
        iTap      = idx_tap[b_lines,end]
        trsf_info = hcat(iTap,tp_rt_min,tp_rt_max)
        push!(node_data,node(bus,b_cbuses,b_lines,gij_line_sh,gij_line,bij_line_sh,bij_line,tp_rt,tp_rt_min,tp_rt_max,from_col,to_col,iTap))
        # push!(tdata,trsf_info)
    end
return node_data
end
###############################################################################
function data_trsf_tap(node_data)
# tdata = Array{Float64}(undef,(1,3))

tdata = []
    for i in 1:nBus
        tdata_tap       = node_data[i].node_idx_trsf
        tdata_tap_min   = node_data[i].node_tap_ratio_min
        tdata_tap_max   = node_data[i].node_tap_ratio_max
        tdata_tap_ratio = node_data[i].node_tap_ratio
        aa  = hcat(tdata_tap,tdata_tap_min,tdata_tap_max,tdata_tap_ratio)
        push!(tdata,aa)
    end
    tdata   = unique(vcat(tdata...),dims=1)
    itap_tr = findall(x->x!=0,tdata[:,1])                                            # Index of transformer not equal to zero
    itap_0  = findall(x->x==0,tdata[:,1])                                            # Index of transformer equal to zero
return tdata,itap_tr,itap_0
end
###############################################################################
function data_trsf_tap_index(tdata,itap_tr,itap_0)                                                  # Rows that has to be removed based on the idx_r values
itap_r = Float64[]
    for i in 1:size(itap_0,1)
        if iszero(tdata[itap_0[i,1],:])                                              # Complete row is zero
            push!(itap_r,itap_0[i,1])                                                # Saving the index of zero row in order to delete it
        elseif iszero(tdata[itap_0[i,1],1]) && (!iszero(tdata[itap_0[i,1],2]) || !iszero(tdata[itap_0[i,1],3]))
            # println("")
            # global error_msg = "ERROR! The tap ratio is set to 0 but min and max tap ratio values are not set to zero. Is this a transformer branch or a distribution line?"
            # println(error_msg)
            # println("Check the transformer data!")      # Transformer branch constraint is not set for this condition
        elseif !iszero(tdata[itap_0[i,1],1]) && (iszero(tdata[itap_0[i,1],2]) && iszero(tdata[itap_0[i,1],3]))
            # global error_msg = "ERROR! The transformer min and max tap ratios are both set to 0"
            # println(error_msg)
            # println("Check the transformer data!")                                      # Transformer branch power injection constraint is replaced by the transmission line constraint
        end
    end
return itap_r
end
################################################################################
function data_load(nw_loads,rheader_loads,rdata_loads,nLoads,nTP,nSc,sbase)
#--------------------------- Load Data -----------------------------------------
    bus_data_lsheet  = [nw_loads[j].load_bus_num for j in 1:nLoads]                  # Buses order in load sheet
    idx_Gs_lsheet = findall(x->x=="Gs",rheader_loads)
    idx_Bs_lsheet = findall(x->x=="Bs",rheader_loads)
    idx_St_lsheet = findall(x->x=="Status",rheader_loads)                            # Index of status in load sheet!

    iFl     =  findall(x->x==1,rdata_loads[:,idx_St_lsheet])                         # Index of flexible loads
    nFl     =  length(iFl)                                                           # Number of flexible loads in a system
    nd_fl   =  convert.(Int64,rdata_loads[iFl])                                      # Nodes to which flexible loads are connected. Node data is taken from the power profile sheet.

# ------------------- Formatting of Load Cost Data------------------------------
#-------------------------------------------------------------------------------
    dim_Fl           = (length(1:nFl),length(1:nTP))
    cost_load_inc    = [nw_loads[iFl[j,1]].load_cost_inc*sbase for j in 1:nFl]       # Here, if the time step is not hour and cost is given in $/hr, then do we have to multiply with the step or not?
    cost_load_dec    = [nw_loads[iFl[j,1]].load_cost_dec*sbase for j in 1:nFl]
    cfl_inc   = zeros(Float64,dim_Fl)                                            # cfl = cost of flexible loads
    cfl_dec   = zeros(Float64,dim_Fl)
        for i in 1:nFl
            nfl  = nd_fl[i,1]                                                        # Node to which the flexible load is connected
            cfl_inc[i,:]     = transpose(repeat([cost_load_inc[i,1]],nTP))           # The same cost is associated across all time periods
            cfl_dec[i,:]     = transpose(repeat([cost_load_dec[i,1]],nTP))           # The same cost is associated across all time periods
        end
#-------------------------------------------------------------------------------
# Assuming that flexible load cost will be different in each time period and for each scenario
    dim_Fl = (length(1:nSc),length(1:nTP),length(1:nFl))
    prof_cost_inc_fl    = zeros(Float64,dim_Fl)                                      # Profile of cost increment  for each scenario in each time step
    prof_cost_dec_fl    = zeros(Float64,dim_Fl)                                      # Profile of cost decrement  for each scenario in each time step
    for i in 1:nFl
        nd      = nd_fl[i,1]                                                         # Node num is independant of time period and scenario
        ifl     = iFl[i,1]                                                           # Index of active flexible load in a load sheet
        for s in 1:nSc
            load_cost_inc = cfl_inc[i,1:end]
            load_cost_dec = cfl_dec[i,1:end]

            prof_cost_inc_fl[s,:,i]    = load_cost_inc
            prof_cost_dec_fl[s,:,i]    = load_cost_dec
        end
    end
return cfl_inc,cfl_dec,prof_cost_inc_fl,prof_cost_dec_fl,bus_data_lsheet,idx_Gs_lsheet,idx_Bs_lsheet,idx_St_lsheet,nFl,iFl,nd_fl
end
################################################################################
function data_gen(nw_gens,nw_gcost,rheader_gens,rdata_gens,nGens,nTP)
    bus_data_gsheet = [nw_gens[j].gen_bus_num for j in 1:nGens]                      # Buses order in gen sheet
    bus_data_gcost = [nw_gcost[j].gcost_bus for j in 1:nGens]
    dim_Gen    = (length(1:nGens),length(1:nTP))
    Pg_max     = zeros(Float64,dim_Gen)
    Pg_min     = zeros(Float64,dim_Gen)
    Qg_max     = zeros(Float64,dim_Gen)
    Qg_min     = zeros(Float64,dim_Gen)
    cA_gen     = zeros(Float64,dim_Gen)                                              # Quadratic cost of generator
    cB_gen     = zeros(Float64,dim_Gen)                                              # Linear cost of generator
    cC_gen     = zeros(Float64,dim_Gen)                                              # Constant cost of generator
    pg_max = [nw_gens[j].gen_P_max for j in 1:nGens]                                 # Maximum and minimum values are in SI units
    pg_min = [nw_gens[j].gen_P_min for j in 1:nGens]
    qg_max = [nw_gens[j].gen_Qg_max for j in 1:nGens]
    qg_min = [nw_gens[j].gen_Qg_min for j in 1:nGens]

    cost_a_gen = [nw_gcost[j].gcost_a for j in 1:nGens]                              # Cost of generation (quadratic and linear) is in SI units
    cost_b_gen = [nw_gcost[j].gcost_b for j in 1:nGens]
    cost_c_gen = [nw_gcost[j].gcost_c for j in 1:nGens]

for i in 1:nGens
    ngen  = bus_data_gsheet[i,1]                                                 # Node to which the gen is connected
    Pg_max[i,:] = transpose(repeat([pg_max[i,1]],nTP))                           # In future, the value of pg_max has to be changed at this location
    Pg_min[i,:] = transpose(repeat([pg_min[i,1]],nTP))
    Qg_max[i,:] = transpose(repeat([qg_max[i,1]],nTP))
    Qg_min[i,:] = transpose(repeat([qg_min[i,1]],nTP))
    cA_gen[i,:] = transpose(repeat([cost_a_gen[i,1]],nTP))
    cB_gen[i,:] = transpose(repeat([cost_b_gen[i,1]],nTP))
    cC_gen[i,:] = transpose(repeat([cost_c_gen[i,1]],nTP))
end
    idx_curt_status = findall(x->x=="CurtStatus",rheader_gens)
    i_ncurt_gens    = findall(x->x==0,rdata_gens[:,idx_curt_status])                 # Index of generator whose power cannot be curtailed!
    i_curt_gens     = findall(x->x==1,rdata_gens[:,idx_curt_status])                 # Index of generator whose power can be curtailed!
    nNcurt_gen      = size(i_ncurt_gens,1)                                           # Number of non-curtailable DGs (Dispatchable DGs such as nuclear, coal, gas turbine are not allowed to curtail!)
    nCurt_gen       = size(i_curt_gens,1)                                            # Number of curtailable DGs (RERs power can be curtailed!)
    nd_ncurt_gen    = convert.(Int64,rdata_gens[i_ncurt_gens,1])                     # Nodes to which those DGs are connected whose power cannot be curtailed!
    nd_curt_gen     = convert.(Int64,rdata_gens[i_curt_gens,1])                      # Nodes to which those DGs are connected whose power can be curtailed!

return bus_data_gsheet,i_ncurt_gens,i_curt_gens,nNcurt_gen,nCurt_gen,nd_ncurt_gen,nd_curt_gen,Pg_max,Pg_min,Qg_max,Qg_min,cA_gen,cB_gen,cC_gen
end
################################################################################
function data_storage(rheader_storage,rdata_storage,nw_storage,nw_Strcost,nTP,nSc,sbase)
    idx_St_Strsheet = findall(x->x=="status",rheader_storage)                        # Index of Status (St) in Storage (Str) Sheet
    iStr_active     = findall(x->x==1,rdata_storage[:,idx_St_Strsheet])
    nStr_active     = length(iStr_active)                                                  # Number of Active Storages in a system
    nd_Str_active   = convert.(Int64,rdata_storage[iStr_active])

    i_soc_0 = findall(i-> i == "soc_0", rheader_storage)
    # i_soc_0 = findfirst(i-> isequal(i, "soc_0"), rheader_storage)
    soc_0 = rdata_storage[:,i_soc_0]
# ------------------- Formatting of Storage Cost Data---------------------------
##------------------------------------------------------------------------------
    dim_Str    = (length(1:nStr_active),length(1:nTP))
    cA_str     = zeros(Float64,dim_Str)
    cB_str     = zeros(Float64,dim_Str)
    cC_str     = zeros(Float64,dim_Str)

    bus_data_Ssheet  = [nw_storage[iStr_active[j,1]].storage_bus for j in 1:nStr_active]                   # Buses order in storage sheet
    bus_data_Strcost = [nw_Strcost[iStr_active[j,1]].str_cost_bus for j in 1:nStr_active]
    cost_a_str = [nw_Strcost[iStr_active[j,1]].str_cost_a*sbase^2*1 for j in 1:nStr_active]                # Here, if the time step is not hour and cost is given in $/hr, then do we have to multiply with the step or not?
    cost_b_str = [nw_Strcost[iStr_active[j,1]].str_cost_b*sbase*1 for j in 1:nStr_active]
    cost_c_str = [nw_Strcost[iStr_active[j,1]].str_cost_c for j in 1:nStr_active]

    for i in 1:nStr_active
        nstr  = bus_data_Ssheet[i,1]                                                 # Node to which the storage is connected
        cA_str[i,:]     = transpose(repeat([cost_a_str[i,1]],nTP))
        cB_str[i,:]     = transpose(repeat([cost_b_str[i,1]],nTP))
        cC_str[i,:]     = transpose(repeat([cost_c_str[i,1]],nTP))
    end

##------------------------------------------------------------------------------
# Assuming that storage cost will be different in each time period and for each scenario
# To run the following code, first the file storage_data_mpopf has to be run
    dim_Str = (length(1:nSc),length(1:nTP),length(1:nStr_active))
    prof_costA_str    = zeros(Float64,dim_Str)
    prof_costB_str    = zeros(Float64,dim_Str)
    prof_costC_str    = zeros(Float64,dim_Str)
    for i in 1:nStr_active
        nd      = nw_storage[iStr_active[i,1]].storage_bus                           # Node num is independant of time period and scenario
        nd_num  = unique(nd)
        nd_num  = nd_num[1,1]
        iStr = iStr_active[i,1][1]                                                   # Index of storage in a Storage sheet
        for s in 1:nSc
            str_costA = cA_str[iStr,1:end]
            str_costB = cB_str[iStr,1:end]
            str_costC = cC_str[iStr,1:end]

            prof_costA_str[s,:,i]    = str_costA
            prof_costB_str[s,:,i]    = str_costB
            prof_costC_str[s,:,i]    = str_costC
        end
    end
return bus_data_Ssheet,bus_data_Strcost,idx_St_Strsheet,iStr_active,nStr_active,nd_Str_active,cA_str,cB_str,cC_str,prof_costA_str,prof_costB_str,prof_costC_str, soc_0
end
################################################################################
function data_load_gen(nLoads,nGens,nTP,nSc,bus_data_lsheet,bus_data_gsheet,Pg_max,Pg_min,Qg_max,Qg_min,nw_pPrf_data_load,nw_qPrf_data_load,nw_pPrf_data_gen_max,nw_pPrf_data_gen_min,nw_qPrf_data_gen_max,nw_qPrf_data_gen_min,nw_loads,sbase,scenario,scenario_data_p_max,scenario_data_p_min,scenario_data_q_max,scenario_data_q_min,cA_gen,cB_gen,cC_gen)
    dim_Load    = (length(1:nSc),length(1:nTP),length(1:nLoads))
    prof_Ploads = zeros(Float64,dim_Load)
    prof_Qloads = zeros(Float64,dim_Load)
    for i in 1:nLoads
        nd      = bus_data_lsheet[i,1]                                           # Node num is independant of time period and scenario
        nd_num  = unique(nd)
        nd_num  = nd_num[1,1]
        iload   = findall(x->x==nd_num,nw_pPrf_data_load[:,1])                   # Index of load in a Load Profile sheet (sheet containing load profiles over horizon)
        for s in 1:nSc
            # nd_pPrf = nw_pPrf[iload,2:end]./sbase                                  # Here, the data is being extracted from the profile sheet. However, later when the scenarios will be created, the data has to be extracted from scenario profiles
            # nd_qPrf = nw_qPrf[iload,2:end]./sbase
            if nTP == 1                                                              # For time period equal to 1, pick the data from Loads sheet
                nd_pPrf_load = nw_loads[i].load_P./sbase
                nd_qPrf_load = nw_loads[i].load_Q./sbase
                # nd_pPrf_load = nw_loads[iload[1,1]].load_P./sbase
                # nd_qPrf_load = nw_loads[iload[1,1]].load_Q./sbase
                prof_Ploads[s,1:nTP,i] .= nd_pPrf_load
                prof_Qloads[s,1:nTP,i] .= nd_qPrf_load
            else                                                                     # For time period, greater than 1, pick the data from profile sheet
                nd_pPrf_load = nw_pPrf_data_load[iload,2:nTP+1]./sbase               # Here, the data is being extracted from the profile sheet. However, later # when the scenarios will be created, the data has to be extracted from scenario profiles
                nd_qPrf_load = nw_qPrf_data_load[iload,2:nTP+1]./sbase
                prof_Ploads[s,1:nTP,i] = nd_pPrf_load
                prof_Qloads[s,1:nTP,i] = nd_qPrf_load
            end
        end
    end
    p_load = prof_Ploads
    q_load = prof_Qloads
#-------------------- Code for formatting the gen data -------------------------
# Assuming that generation limit will be different in each time period and for each scenario
# To run the following code, first the file gen_mpopf has to be run!!!
    dim_Gen = (length(1:nSc),length(1:nTP),length(1:nGens))
    prof_Pmaxgens = zeros(Float64,dim_Gen)
    prof_Pmingens = zeros(Float64,dim_Gen)
    prof_Qmaxgens = zeros(Float64,dim_Gen)
    prof_Qmingens = zeros(Float64,dim_Gen)
    prof_costA_gen    = zeros(Float64,dim_Gen)
    prof_costB_gen    = zeros(Float64,dim_Gen)
    prof_costC_gen    = zeros(Float64,dim_Gen)
    for i in 1:nGens
        nd      = bus_data_gsheet[i,1]                                               # Node num is independant of time period and scenario
        nd_num  = unique(nd)
        nd_num  = nd_num[1,1]
        # igen    = findall(x->x==nd_num,bus_data_gsheet)                               # This is WRONG!
        igen    = findall(x->x==nd_num,nw_pPrf_data_gen_max[:,1])                    # Index of load in a Gen sheet
        for s in 1:nSc
            if nTP==1
                g_pmax  = Pg_max[i,1:end]./sbase
                g_pmin  = Pg_min[i,1:end]./sbase
                g_qmax  = Qg_max[i,1:end]./sbase
                g_qmin  = Qg_min[i,1:end]./sbase

                g_costA = cA_gen[i,1:end].*sbase^2
                g_costB = cB_gen[i,1:end].*sbase
                g_costC = cC_gen[i,1:end]

                # g_pmax  = Pg_max[igen,1:end]./sbase
                # g_pmin  = Pg_min[igen,1:end]./sbase
                # g_qmax  = Qg_max[igen,1:end]./sbase
                # g_qmin  = Qg_min[igen,1:end]./sbase
            else
                g_pmax  = nw_pPrf_data_gen_max[igen,2:nTP+1]./sbase
                g_pmin  = nw_pPrf_data_gen_min[igen,2:nTP+1]./sbase
                g_qmax  = nw_qPrf_data_gen_max[igen,2:nTP+1]./sbase
                g_qmin  = nw_qPrf_data_gen_min[igen,2:nTP+1]./sbase

                g_costA = cA_gen[igen,1:end].*sbase^2
                g_costB = cB_gen[igen,1:end].*sbase
                g_costC = cC_gen[igen,1:end]
            end
            # g_costA = cA_gen[igen,1:end].*sbase^2
            # g_costB = cB_gen[igen,1:end].*sbase
            # g_costC = cC_gen[igen,1:end]

            prof_Pmaxgens[s,1:nTP,i] = g_pmax
            prof_Pmingens[s,1:nTP,i] = g_pmin
            prof_Qmaxgens[s,1:nTP,i] = g_qmax
            prof_Qmingens[s,1:nTP,i] = g_qmin

            prof_costA_gen[s,1:nTP,i] = g_costA
            prof_costB_gen[s,1:nTP,i] = g_costB
            prof_costC_gen[s,1:nTP,i] = g_costC
        end
    end
################################################################################
##-------------------------------------------------------------------------- ###
##------ Which data to be used for generation minimum and maximum values?? -----
    if scenario == 0
        pg_max = prof_Pmaxgens
        pg_min = prof_Pmingens
        qg_max = prof_Qmaxgens
        qg_min = prof_Qmingens
    else
        pg_max = scenario_data_p_max
        pg_min = scenario_data_p_min
        qg_max = scenario_data_q_max
        qg_min = scenario_data_q_min
    end
    cost_a_gen = prof_costA_gen
    cost_b_gen = prof_costB_gen
    cost_c_gen = prof_costC_gen
return p_load,q_load,pg_max,pg_min,qg_max,qg_min,cost_a_gen,cost_b_gen,cost_c_gen
end
