include(joinpath(@__DIR__,"compact_unit.jl"))
include(joinpath(@__DIR__,"compact_chain.jl"))
p0=unit_elem();p=construct_blocks(p0,tree,1e-12)
s=Dict{String,Int}();chain=cc_prepend(p,c;stats=s)
K1=Matrix(Symmetric(p0.Mbb+x.Maa,:U));K2=Matrix(Symmetric(x.Mbb+y.Maa,:U))
H=[K1 x.Mab;x.Mab' K2]
Z=randn(2n,3)
@assert norm(ck_apply(chain.H,Z)-H*Z)/norm(H*Z)<1e-10
@assert norm(ck_solve(chain.factor,Z)-H\Z)/norm(H\Z)<1e-9
B=[p0.Mab zeros(n,n)];D=[zeros(n,n);y.Mab]
refs=(Matrix(Symmetric(p0.Maa,:U))-B*(H\B'),-B*(H\D),Matrix(Symmetric(y.Mbb,:U))-D'*(H\D))
for (j,block) in enumerate((:Maa,:Mab,:Mbb))
    @assert norm(ck_action(chain,block,X)-refs[j]*X)/max(1,norm(refs[j]*X))<1e-9
end
T=5Matrix{Float64}(I,n,n);reports=Any[]
terminal=(rep=block_payload(T,tree,1e-12,reports,:terminal),v=randn(n))
closed=cc_close(chain,terminal;stats=s)
HH=[K1 x.Mab zeros(n,n);x.Mab' K2 y.Mab;zeros(n,n) y.Mab' Matrix(Symmetric(y.Mbb+T,:U))]
ZZ=randn(3n,3)
@assert norm(ck_apply(closed.H,ZZ)-HH*ZZ)/norm(HH*ZZ)<1e-10
BB=[p0.Mab zeros(n,2n)];V=Matrix(Symmetric(p0.Maa,:U))-BB*(HH\BB')
@assert norm(cc_value(closed,X)-V*X)/norm(V*X)<1e-9
g=vcat(p0.mb+x.ma,x.mb+y.ma,y.mb+terminal.v)
@assert norm(closed.v-(p0.ma-BB*(HH\g)))<1e-9
SPD=Matrix{Float64}(I,n,n)+X*X'
sol=cc_cg(z->SPD*z,X[:,1],diag(SPD))
@assert norm(sol.x-SPD\X[:,1])/norm(sol.x)<1e-9
rejected=try;cc_cg(z->-z,ones(n),ones(n));false;catch ex;occursin("NOT_SPD",sprint(showerror,ex));end
@assert rejected
@assert s["max_local_dense_dim"]<=CC_LOCAL_FALLBACK_CAP
println("COMPACT_CHAIN_UNIT_PASS two_merges=true terminal_actions=true raw_skew=true interleaved=true interior_growth=true")
