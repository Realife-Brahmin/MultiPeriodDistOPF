using Test
include("sequential_dss.jl")
const S=SequentialDSS
const D=S.D
function toy(;initial=50.0,eff=100.0,dt=1.0)
    D.Text.Command("clear")
    D.Text.Command("new circuit.seq basekv=12.47 pu=1.0 phases=1 bus1=1")
    D.Text.Command("new storage.b phases=1 bus1=1 kv=12.47 kwrated=25 kva=25 kwhrated=100 kwhstored=$initial %reserve=30 %EffCharge=$eff %EffDischarge=$eff %IdlingkW=0 %cutin=0 %cutout=0 DispMode=External model=1")
    D.Text.Command("set mode=daily stepsize=$(3600dt)s number=1 controlmode=off hour=0 sec=0 tolerance=1e-10 maxiterations=100")
    D.Circuit.SetActiveElement("Storage.b")
    @test lowercase(S.prop("DispMode"))=="external"
end
function step(p;dt=1.0,eff=1.0)
    a=S.storage_step!(["b"],[p],dt,[eff],[eff])
    b=S.finish_storage_step!(["b"],a,dt,[eff],[eff])
    @test a.converged
    @test isapprox(b.h1,a.h0+dt;atol=1e-10)
    return a,b
end
@testset "Chronology, signs and independent energy conservation" begin
    toy();a,b=step(-10.0)
    @test a.before[1]≈50 atol=1e-7
    @test a.actual[1]≈-10 atol=1e-6
    @test b.after[1]≈60 atol=1e-6
    @test abs(b.balance[1])<1e-6
    a,b=step(20.0)
    @test a.before[1]≈60 atol=1e-6
    @test b.after[1]≈40 atol=1e-6
    @test b.h1==2
    @test abs(b.balance[1])<1e-6
end
@testset "Reserve clamping cannot masquerade as executable energy" begin
    toy(initial=35.0);a,b=step(10.0)
    @test a.actual[1]≈10 atol=1e-6
    @test b.after[1]≈30 atol=1e-6
    @test b.balance[1]≈5 atol=1e-6
    @test S.bound_error(25.0,30.0,95.0)==5
    a,b=step(10.0)
    @test a.before[1]≈30 atol=1e-6
    @test abs(a.actual[1])<1e-6
    @test abs(a.actual[1]-10)>9.99
end
@testset "Capacity clamping and operating upper bound remain distinct" begin
    toy(initial=95.0);a,b=step(-10.0)
    @test b.after[1]≈100 atol=1e-6
    @test b.balance[1]≈-5 atol=1e-6
    @test S.bound_error(b.after[1],30.0,95.0)≈5 atol=1e-6
    a,b=step(-10.0)
    @test abs(a.actual[1])<1e-6
end
@testset "Efficiencies and quarter-hour durations" begin
    toy(eff=95.0);a,b=step(-10.0;eff=.95)
    @test b.after[1]≈59.5 atol=1e-6
    @test abs(b.balance[1])<1e-6
    a,b=step(10.0;eff=.95)
    @test b.after[1]≈59.5-10/.95 atol=1e-6
    @test abs(b.balance[1])<1e-6
    toy(dt=.25);a,b=step(10.0;dt=.25)
    @test b.after[1]≈47.5 atol=1e-6
    @test b.h1==.25
    a,b=step(-4.0;dt=.25)
    @test a.before[1]≈47.5 atol=1e-6
    @test b.after[1]≈48.5 atol=1e-6
    @test b.h1==.5
end
@testset "Explicit unity shape does not double-apply load profile" begin
    toy()
    D.Text.Command("new LoadShape.seq_unity npts=1 interval=1 mult=(1)")
    D.Text.Command("new Load.test phases=1 bus1=1 kv=12.47 model=1 kw=10 kvar=2 daily=seq_unity yearly=seq_unity duty=seq_unity")
    for t in 1:3
        a=S.storage_step!(["b"],[0.0],1.0,[1.0],[1.0])
        D.Circuit.SetActiveElement("Load.test")
        @test real(first(D.CktElement.Powers()))≈10 atol=1e-6
        @test imag(first(D.CktElement.Powers()))≈2 atol=1e-6
        b=S.finish_storage_step!(["b"],a,1.0,[1.0],[1.0])
        @test b.h1==t
    end
end
@testset "Lossless idle periods do not leak stored energy" begin
    toy(initial=50.0,eff=100.0)
    a,b=step(0.0)
    @test abs(a.actual[1])<1e-7
    @test b.after[1]≈50 atol=1e-7
    @test abs(b.balance[1])<1e-7
    a,b=step(0.0)
    @test b.after[1]≈50 atol=1e-7
    @test b.h1==2
end
