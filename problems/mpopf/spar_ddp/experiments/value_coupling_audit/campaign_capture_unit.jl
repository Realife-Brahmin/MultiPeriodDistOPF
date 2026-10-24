using Printf,Serialization,SharedArrays
include(joinpath(@__DIR__,"campaign_capture_helper.jl"))
ENV["FILTERDDP_CAMPAIGN_CAPTURE_TAGS"]="late"
ENV["FILTERDDP_CAMPAIGN_CAPTURE_STAGES"]="4"
state=Dict{Symbol,Any}(:tags=>Set{String}())
opts=(optimality_tolerance=1e-5,)
@assert isempty(fc_select!(state,(dual_inf=1.0,cs_inf_0=0.0,primal_inf=0.0),opts,nothing))
@assert fc_select!(state,(dual_inf=1e-6,cs_inf_0=0.0,primal_inf=0.0),opts,nothing)==["late"]
push!(state[:tags],"late")
@assert isempty(fc_select!(state,(dual_inf=1e-6,cs_inf_0=0.0,primal_inf=0.0),opts,nothing))
@assert campaign_stages(3)==[1,2,3]&&campaign_stages(24)==[1,2,3,4]
W=(ocp=(N=24,control_limits=nothing),nx=2,nu=3,nc=4);S=(x=[1.0],)
for t in 1:24;fc_cache_stage(t,S,W);end
@assert sort(collect(keys(_FC_STAGE)))==[1,2,3,4]
S.x[1]=2.0;@assert _FC_STAGE[1].S.x[1]==1.0
struct PintElem
    Maa::Matrix{Float64};Mab::Matrix{Float64};Mbb::Matrix{Float64};ma::Vector{Float64};mb::Vector{Float64}
end
n=2;sh=SharedArray{Float64}((3n*n+2n,));sh[:]=collect(1.0:length(sh))
a=fc_before((exits=1,shared=[sh]),n,0,0,0.2)
@assert vec(a.elems[1].Maa)==collect(1.0:4.0)&&a.elems[1].mb==[15.0,16.0]
@assert a.k==0&&a.j==0&&a.mu==0.2&&a.index==1
println("CAMPAIGN_CAPTURE_UNIT_PASS stages1to4=true late_only=true payload_layout=true deep_copy=true")
