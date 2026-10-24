# Synthetic radial KKT correctness test; no exported cases or external solvers.
using SPARDDP, LinearAlgebra, SparseArrays, Random, Test, Printf
include(joinpath(@__DIR__, "control_layout.jl"))
include(joinpath(@__DIR__, "tree_kkt.jl"))
ENV["FILTERDDP_HALF_BLOCK"] = "1"
ENV["FILTERDDP_LEAN_EXACT"] = "1"
BLAS.set_num_threads(1)

function fixture()
    # A branching tree with batteries at both parent 2 and child 3; no special
    # independent-feeder assumption. Both arms see exactly the same matrix.
    data = Dict(:Nset=>[1,2,3,4], :Lset=>[(1,2),(2,3),(2,4)],
        :Nm1set=>[2,3,4], :substationBus=>1, :Bset=>[2,3], :Dset=>[4],
        :parent=>Dict(2=>1,3=>2,4=>2))
    idx, nu = control_layout(data)
    nc = 2length(data[:Nset]) + 2length(data[:Lset]) + 1 + 2
    n = nu + nc
    K = spdiagm(0=>fill(8.0, n))
    # Construct physical network groups using the production index convention.
    lay0 = tree_kkt_layout(data, idx, nu, K)
    rng = MersenneTwister(42)
    for j in 1:4
        own = lay0.own[j]
        for a in 1:length(own)-1
            i,k = own[a],own[a+1]; K[i,k]=K[k,i]=0.1randn(rng)
        end
        j == 1 && continue
        parent = lay0.parent[j]
        i,k = first(own),first(lay0.own[parent]); K[i,k]=K[k,i]=0.2
    end
    for (b,j) in enumerate([2,3])
        p,e = lay0.E[b],lay0.E[b+2]
        # Network coupling only on power rows; NN comes from energy auxiliaries.
        net = first(lay0.own[j]); K[p,net]=K[net,p]=0.3
        K[p,e]=K[e,p]=-1.0
        a=lay0.iso[b]; K[a,e]=K[e,a]=0.4
    end
    return K,tree_kkt_layout(data,idx,nu,K)
end

@testset "tree versus generic sparse network elimination" begin
    K,lay=fixture(); stat=tree_kkt_static(lay,K); nB=length(lay.E)÷2
    # Non-diagonal future curvature tests the actual exact-curvature path.
    Vc=[1.0 0.2;0.2 1.5]
    withenv("FILTERDDP_NETWORK_BACKEND"=>"tree") do
        global tree=tree_kkt_solver(lay,stat,K;Vc=Vc)
    end
    withenv("FILTERDDP_NETWORK_BACKEND"=>"umfpack") do
        global generic=tree_kkt_solver(lay,stat,K;Vc=Vc)
    end
    Kfull=copy(K); Kfull[lay.E[1:nB],lay.E[1:nB]] += sparse(Vc)
    rng=MersenneTwister(71)
    b=randn(rng,lay.n)
    xt=copy(b);xg=copy(b);ldiv!(tree,xt);ldiv!(generic,xg)
    xr=lu(Kfull)\b
    @test xt ≈ xr rtol=1e-10 atol=1e-11
    @test xg ≈ xr rtol=1e-10 atol=1e-11
    @test xg ≈ xt rtol=1e-11 atol=1e-12
    RE=randn(rng,length(lay.E),nB)
    full=zeros(lay.n,nB);full[lay.E,:]=RE
    @test tree_kkt_battery_rows(generic,RE) ≈ (lu(Kfull)\full)[lay.E,:] rtol=1e-10
    B=randn(rng,nB,nB);cx=sparse(-Matrix{Float64}(I,nB,nB))
    @test tree_kkt_feedback_rows(tree,B,cx) ≈ tree_kkt_feedback_rows(generic,B,cx) rtol=1e-10
    @test last(tree_kkt_value_halved(tree,B,cx)) ≈ last(tree_kkt_value_halved(generic,B,cx)) rtol=1e-10
    @test_throws ErrorException withenv("FILTERDDP_NETWORK_BACKEND"=>"invalid") do
        network_kkt_factor(lay,K)
    end
    mktempdir() do directory
        _NETWORK_CAPTURE_COUNT[]=0
        withenv("TREE_ABLATION_CAPTURE_DIR"=>directory,"FILTERDDP_NETWORK_BACKEND"=>"tree") do
            network_kkt_factor(lay,K)
        end
        file=joinpath(directory,"network_$(getpid()).jls")
        @test isfile(file)
        @test SPARDDP.Serialization.deserialize(file).K==K
        if !isempty(ARGS)
            cp(file,ARGS[1];force=false)
        end
    end
    @printf("TREE_ABLATION_SMOKE_PASS n=%d nE=%d nB=%d relative_solution_difference=%.3e\n",
        lay.n,length(lay.E),nB,norm(xt-xg)/norm(xt))
end
