# A targeted algebra regression with asymmetric raw diagonal blocks and
# non-symmetric cross blocks. Exercise both Cholesky and unchanged fallbacks.
using LinearAlgebra, Random
include(joinpath(ENV["REPLAY_REPO"], "ddp/examples/power_system/ieee123c_filterddp.jl"))
include(joinpath(ENV["REPLAY_REPO"], "ddp/examples/power_system/procs_pass.jl"))
include(joinpath(@__DIR__, "right_kernel.jl"))
Random.seed!(20261105)
BLAS.set_num_threads(1)
for n in (1, 3, 17), terminal in (false, true), positive in (true, false)
    # Indefinite interface is invertible and should take the original fallback.
    aa=randn(n,n); bb=(positive ? 8.0 : -8.0)*Matrix{Float64}(I,n,n)+0.01randn(n,n)
    x=PintElem(aa,randn(n,n),bb,randn(n),randn(n))
    y=terminal ? PintVal(0.01randn(n,n),randn(n)) :
        PintElem(0.01randn(n,n),randn(n,n),randn(n,n),randn(n),randn(n))
    reference=_pp_combine!(deepcopy(x),deepcopy(y),PairWork(n))
    candidate=right_combine!(deepcopy(x),deepcopy(y),PairWork(n))
    fields=terminal ? (:V,:v) : (:Maa,:Mab,:Mbb,:ma,:mb)
    for field in fields
        a=getproperty(candidate,field); b=getproperty(reference,field)
        @assert norm(a-b,Inf)/max(1,norm(b,Inf))<1e-11
    end
end
println("RIGHT_KERNEL_SYNTHETIC_PASS")
