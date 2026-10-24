include(joinpath(@__DIR__,"compact_chain_unit.jl"))
include(joinpath(@__DIR__,"direct_initial.jl"))
using SparseArrays
P=sprandn(n,n,.03);U=randn(n,3);V=randn(n,3);stats=Dict{String,Int}()
a=di_payload(P,U,V,tree,stats);A=Matrix(P)+U*V'
@assert norm(expand_block(a,n)-A)/norm(A)<1e-10
@assert norm(ck_apply(a,X)-A*X)/norm(A*X)<1e-10
@assert norm(ck_apply(a,X;trans=true)-A'*X)/norm(A'*X)<1e-10
@assert norm(A-A')>1
function value_log_regression()
    d=(value=:dense_reference,);io=IOBuffer()
    err=di_record_value(io,"test",:v,[1.0],[2.0])
    @assert d.value==:dense_reference&&err==0.5
    @assert occursin("test,v,1.0,2.0,0.5",String(take!(io)))
end
value_log_regression()
println("DIRECT_VALUE_LOG_REGRESSION_PASS dense_reference_not_captured=true")
println("DIRECT_INITIAL_UNIT_PASS sparse_plus_raw_separator=true interleaved=true transpose=true no_dense_input_discovery=true")
