include("historical_stopping.jl")
using Test
@testset "Historical owned-period RMS stopping" begin
    windows=[[1,2],[1,2,3],[2,3]]
    z=zeros(2,3); old=fill(0.01,2,3)
    copies=[zeros(2,length(w)) for w in windows]
    copies[1][:,2].=10 # Non-owned copy must not enter the historical metric.
    copies[2][:,2].=sqrt(3)*1e-4
    r,s=historical_rms(copies,z,old,windows,2.0)
    @test r≈1e-4
    @test s≈0.02
    @test historical_stop(1e-4,1e-2,[100.0]).met
    @test !historical_stop(1.001e-4,0.0,fill(100.0,5)).met
    @test !historical_stop(0.0,0.02,fill(100.0,4)).met
    @test historical_stop(1e-4,0.02,[100.0,100.01,100.02,100.01,100.0]).met
    @test !historical_stop(0.0,0.02,[100.0,101,100,100,100]).met
    @test !historical_stop(0.0,0.02,zeros(5)).met
end
