using Serialization, LinearAlgebra, SparseArrays, SHA
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"../../.."))
include(joinpath(ROOT,"spar_ddp/examples/power_system/pint_scan.jl"))
length(ARGS)==2 || error("Usage: julia inspect.jl INPUT_DIRECTORY OUTPUT_DIRECTORY")
const DIR=abspath(ARGS[1])
const OUT=abspath(ARGS[2]); mkpath(OUT)
d=deserialize(joinpath(DIR,"network_data_large10kC_1ph_T96_periodic.jls"))
d3=deserialize(joinpath(DIR,"network_data_large10kC_1ph_T3_periodic.jls"))
@assert d[:Bset]==d3[:Bset] && d[:parent]==d3[:parent]
function groupof(b)
    while d[:parent][b] != d[:substationBus]; b=d[:parent][b]; end
    b
end
groups=groupof.(d[:Bset]); keys=sort(unique(groups)); n=length(groups)
@assert n==1020
left=findall(in(Set(keys[1:length(keys)÷2])),groups)
right=setdiff(1:n,left)
largest=argmax(k->count(==(k),groups),keys)
one=findall(==(largest),groups)
println("METADATA n=",n," groups=",length(keys)," left=",length(left)," right=",length(right)," largest_group=",length(one)," battery_order_metadata_match=true")
es=[deserialize(joinpath(DIR,"stage_$t.jls")) for t in 1:4]
open(joinpath(OUT,"input.sha256"),"w") do io
    for f in ("stage_1.jls","stage_2.jls","stage_3.jls","stage_4.jls","network_data_large10kC_1ph_T96_periodic.jls","network_data_large10kC_1ph_T3_periodic.jls")
        println(io,bytes2hex(sha256(read(joinpath(DIR,f)))),"  ",f)
    end
end
pair=pint_combine(es[1],es[2]); four=pint_combine(pair,pint_combine(es[3],es[4]))
function spectrum(A)
    s=svdvals(A); scale=isempty(s) ? 0.0 : s[1]
    ranks=[scale==0 ? 0 : count(>(tol*scale),s) for tol in (1e-12,1e-8,1e-6,1e-4)]
    tail=scale==0 ? 0.0 : s[min(4,length(s))]/scale
    join((ranks...,tail),",")
end
open(joinpath(OUT,"matrix_summary.csv"),"w") do io
    println(io,"interval,block,n,exact_nonzeros,relative_symmetry_error,cross_rows,cross_cols,rank_1e12,rank_1e8,rank_1e6,rank_1e4,sigma4_over_sigma1,cross_frobenius_fraction,within_group_rank_1e12,within_group_rank_1e8,within_group_rank_1e6,within_group_rank_1e4,within_group_sigma4_over_sigma1")
    for (label,e) in (("0:1",es[1]),("1:2",es[2]),("0:2",pair),("0:4",four))
        for name in (:Maa,:Mab,:Mbb)
            A=getproperty(e,name); X=A[left,right]
            println(io,join((label,name,n,count(!iszero,A),norm(A-A')/max(norm(A),eps()),length(left),length(right)),","),",",spectrum(X),",",norm(X)/max(norm(A),eps()),",",spectrum(A[one,one]))
        end
    end
end
println(read(joinpath(OUT,"matrix_summary.csv"),String))
println("INSPECTION_COMPLETE first_pass_only=true local_diagnostic_not_benchmark=true")
