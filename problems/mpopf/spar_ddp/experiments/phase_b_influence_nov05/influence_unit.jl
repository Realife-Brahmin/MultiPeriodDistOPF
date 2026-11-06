using LinearAlgebra,Random,Test
include(joinpath(@__DIR__,"influence_core.jl"))
struct TinyElem
 Maa::Matrix{Float64};Mab::Matrix{Float64};Mbb::Matrix{Float64};ma::Vector{Float64};mb::Vector{Float64}
end
struct TinyWork
 K::Matrix{Float64};g::Vector{Float64};Y1::Matrix{Float64};Y2::Matrix{Float64}
end
Random.seed!(9271);n=5
x=TinyElem(Matrix(4.0I,n,n),randn(n,n),Matrix(5.0I,n,n),randn(n),randn(n))
y=TinyElem(Matrix(6.0I,n,n),randn(n,n),Matrix(7.0I,n,n),randn(n),randn(n))
z=merge_inputs(x,y);old=(H=z.H+Diagonal(ones(n)),BL=z.BL.+0.1,BR=z.BR.-0.2,h=z.h.+0.3)
for mode in INTERVENTIONS
 q=retained_inputs(z,old,mode)
 @test z.H==x.Mbb+y.Maa # no mutation of current inputs
 @test mode==:H ? q.H==old.H : mode==:Hdiag ? diag(q.H)==diag(old.H) : q.H==z.H
 @test mode==:BLoffdiag ? diag(q.BL)==diag(z.BL) : true
 a=deepcopy(x);w=TinyWork(zeros(n,n),zeros(n),zeros(n,n),zeros(n,n));input_combine!(a,y,w,q)
 K=Symmetric(q.H,:U)
 @test a.Maa≈x.Maa-q.BL*(K\q.BL')
 @test a.Mab≈-q.BL*(K\q.BR)
 @test a.Mbb≈y.Mbb-q.BR'*(K\q.BR)
 @test a.ma≈x.ma-q.BL*(K\q.h)
 @test a.mb≈y.mb-q.BR'*(K\q.h)
end
ops,changed=dependency_trace(24)
@test count(x->x.injected,ops)==1
@test findall(changed)==collect(1:5)
@test any(x->x.right==5&&x.right_dependency,ops)
@test length(scan_ops(24))==68
@test all(x->x in scan_ops(24),LOCATIONS)
@test passes_metric(metric(ones(4),ones(4)))
@test !passes_metric(metric(ones(4),zeros(4)))
println("INFLUENCE_SYNTHETIC_PASS substitutions=6 dependency_suffixes=1:5 consumers=1:4 full_forward=1:24")
