# Read-only raw saved leaves; no model, factors, optimizer or representation.
using Serialization,LinearAlgebra,SparseArrays,SHA
BLAS.set_num_threads(1)
include(joinpath(@__DIR__,"../../examples/power_system/pint_scan.jl"))
length(ARGS)==2 || error("symmetry_audit.jl LEAVES OUTPUT_CSV")
payload=deserialize(ARGS[1]);@assert !ispath(ARGS[2])
open(ARGS[2],"w") do io
    println(io,"stage,block,rows,cols,fro,skew_fro,skew_relative,exact_nnz,nnz_rel_1e12")
    for (t,e) in enumerate(payload.elems),f in (:Maa,:Mab,:Mbb)
        A=getproperty(e,f);skew=norm(A-A');scale=norm(A);threshold=1e-12*maximum(abs,A)
        println(io,join((t,f,size(A,1),size(A,2),scale,skew,skew/max(scale,floatmin(Float64)),count(!iszero,A),count(x->abs(x)>threshold,A)),','))
    end
end
println("RAW_LEAF_SYMMETRY_AUDIT_COMPLETE source_sha256=$(bytes2hex(sha256(read(ARGS[1]))))")
