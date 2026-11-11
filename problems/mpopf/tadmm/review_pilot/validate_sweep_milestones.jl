using Test,Printf,SHA,Serialization,TOML
include("sweep_milestones.jl")
@testset "Milestone definitions distinguish convergence and quality" begin
    early=sweep_flags(1e-7,1e-7,10.0,.06,1e-6)
    @test early["owned_raw_1e-5"] && !early["final_rms"] && !early["physical_NO"]
    diagonal=sweep_flags(1e-6,1e-3,.001,1e-4,.01)
    @test diagonal["owned_dual_1e-5"] && !diagonal["all_dual_1e-5"] && !diagonal["final_rms"]
    final=sweep_flags(1e-6,1e-5,.01,1e-4,1e-5)
    @test final["final_rms"] && final["physical_NO"] && final["current_cost_1e-5"]
    historical=sweep_flags(1e-3,.01,1.0,.004,.01)
    @test historical["paper_cost_1e-3"] && !historical["current_cost_1e-3"]
    groups=[[1,3],[2,4]];times=[1.0,2.0,3.0,4.0]
    @test maximum(sum(times[t] for t in group) for group in groups)==6.0
end
@testset "Saved first crossings are deduplicated and immutable" begin
    mktempdir(pwd()) do folder
        observe,marks=sweep_observer(folder,Dict(:T=>1),"test",1.0)
        metrics=(;objective=1.0,original_primal=1e-6,maxeq=1e-6,maxbound=0.0,maxcone=0.0,soc_equality_gap=0.0,fictitious_loss_kwh=0.0)
        x=(;k=1,metrics,owned_r=1e-6,rr=1e-6,owned_s=1e-3,rmax=1e-6,smax=.001,rho=10.0,elapsed=2.0,us=[[1.0]],bb=zeros(1,1))
        @test observe(x)
        @test observe(merge(x,(;k=2,elapsed=3.0)))
        @test all(v["iteration"]==1 for v in values(marks))
        @test length(filter(p->endswith(p,".jls"),readdir(folder)))==1
        point=joinpath(folder,"point_00001.jls")
        @test marks["final_rms"]["snapshot_sha256"]==bytes2hex(sha256(read(point)))
        @test deserialize(point)[:iteration]==1
        @test TOML.parsefile(joinpath(folder,"milestones.toml"))["final_rms"]["core_s"]==2.0
    end
end
