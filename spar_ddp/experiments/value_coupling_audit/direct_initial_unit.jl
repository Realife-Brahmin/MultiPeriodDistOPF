include(joinpath(@__DIR__,"compact_chain_unit.jl"))
include(joinpath(@__DIR__,"direct_initial.jl"))
using SparseArrays
P=sprandn(n,n,.03);U=randn(n,3);V=randn(n,3);stats=Dict{String,Int}()
a=di_payload(P,U,V,tree,stats);A=Matrix(P)+U*V'
@assert norm(expand_block(a,n)-A)/norm(A)<1e-10
@assert norm(ck_apply(a,X)-A*X)/norm(A*X)<1e-10
@assert norm(ck_apply(a,X;trans=true)-A'*X)/norm(A'*X)<1e-10
@assert norm(A-A')>1
println("DIRECT_INITIAL_UNIT_PASS sparse_plus_raw_separator=true interleaved=true transpose=true no_dense_input_discovery=true")
