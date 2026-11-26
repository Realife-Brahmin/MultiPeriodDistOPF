# Compare fixed-node implementation against pinned pre-repair code, not an
# independently approximated expected answer. No optimizer or changed inputs.
include(joinpath(@__DIR__,"followup_unit.jl"))
BLAS.set_num_threads(1)
legacy_source=if haskey(ENV,"VALUE_FOLLOWUP_LEGACY")
    read(joinpath(ENV["VALUE_FOLLOWUP_LEGACY"],"followup_blocks.jl"),String)
else
    read(`git show f8f60d39:spar_ddp/experiments/value_coupling_audit/followup_blocks.jl`,String)
end
legacy=Module(:PinnedFollowupBlocks)
Core.eval(legacy,:(using LinearAlgebra))
Core.eval(legacy,:(const PintElem=Main.PintElem))
Base.include_string(legacy,legacy_source,"pinned_followup_blocks")
function compare_implementations(raw,hier;budgets=(0.0,1e-6,1e-9,1e-12))
    for tau in budgets
        a=construct_blocks(raw,hier,tau);b=legacy.construct_blocks(raw,hier,tau)
        x=materialize_blocks(a,raw,tau);y=legacy.materialize_blocks(b,raw,tau)
        @assert a.bytes==b.bytes && a.reports==b.reports
        @assert x.fallback==y.fallback && x.error==y.error && x.bytes==y.bytes
        @assert all(getproperty(x.elem,f)==getproperty(y.elem,f) for f in (:Maa,:Mab,:Mbb,:ma,:mb))
    end
end
compare_implementations(e,tree)
println("FIXED_NODE_EQUIVALENCE_PASS synthetic=true byte_identical_arrays=true ranks_budgets_payload_unchanged=true")
if length(ARGS)==2
    payload=deserialize(ARGS[1]);data=deserialize(ARGS[2]);hier=topology_hierarchy(data;maxleaf=8)
    compare_implementations(payload.elems[1],hier;budgets=(0.0,1e-6))
    println("FIXED_NODE_REAL_LEAF_EQUIVALENCE_PASS byte_identical_arrays=true")
end
