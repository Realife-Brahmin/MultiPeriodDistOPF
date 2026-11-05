using Printf,Serialization,SharedArrays
include(joinpath(@__DIR__,"capture_helper.jl"))
@assert length(AUDIT_KEYS)==8
@assert audit_tag(0,0)=="early_before"&&audit_tag(1,0)=="early_after"
@assert audit_tag(7,3)=="middle_before"&&audit_tag(8,3)=="middle_after"
@assert audit_tag(12,7)=="late_before"&&audit_tag(13,7)=="late_after"
@assert audit_tag(4,0)=="barrier_before"&&audit_tag(4,1)=="barrier_after"
@assert audit_tag(14,8)==""&&audit_tag(7,0)==""
W=(ocp=(N=24,control_limits=nothing),nx=2,nu=3,nc=4);S=(x=[1.0],)
for t in 1:24;fc_cache_stage(t,S,W);end
@assert sort(collect(keys(_FC_STAGE)))==collect(1:24)
S.x[1]=2.0;@assert _FC_STAGE[1].S.x[1]==1.0
struct PintElem
    Maa::Matrix{Float64};Mab::Matrix{Float64};Mbb::Matrix{Float64};ma::Vector{Float64};mb::Vector{Float64}
end
n=2;sh=SharedArray{Float64}((3n*n+2n,));sh[:]=collect(1.0:length(sh))
a=fc_before((exits=1,shared=[sh]),n,0,0,1.0,0.125)
@assert vec(a.elems[1].Maa)==collect(1.0:4.0)&&a.elems[1].mb==[15.0,16.0]
@assert a.index==1&&a.reg_before==0.125&&a.tag=="early_before"
@assert isnothing(fc_before((exits=1,shared=[sh]),n,2,0,1.0,0.0))
b=fc_before((exits=1,shared=[sh]),n,4,0,1.0,0.0)
@assert b.index==3&&b.tag=="barrier_before"
println("INFLUENCE_CAPTURE_UNIT_PASS all24stages=true eightkeys=true layout=true reg_metadata=true deepcopy=true unselected_not_persisted=true")


no_parse_error(x)=!(x isa Expr)||(!(x.head in (:error,:incomplete))&&all(no_parse_error,x.args))
for name in ("capture_helper.jl","verify_capture.jl")
    @assert no_parse_error(Meta.parseall(read(joinpath(@__DIR__,name),String))) name
end
println("INFLUENCE_CAPTURE_PARSE_PASS")
