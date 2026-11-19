include("validate_stage_reuse.jl")
using SHA, TOML

function cap_gate_trial(tolerance;expect_accept,screen="feasible")
    withenv("TADMM_CAP_ITERATIONS"=>"1","TADMM_CAP_END_ITERATION"=>"2",
            "TADMM_CAP_PRIMAL_TOL"=>string(tolerance),"TADMM_CAP_SCREEN"=>screen) do
        s=toy_stage("primal_dual","batched";policy="early_cap30")
        # A nonlinear equality exposes a genuine primal error after one step.
        @constraint(s.m,s.b[1,2]*s.b[2,3]==0.8)
        s.reuse.constraints=all_constraints(s.m;include_variable_in_set_constraints=true)
        d=Dict(:Bset=>[1,2]);u=zeros(2,3)
        first=stage_step!(s,d,fill(0.8,2,3),u,100.0,1)
        @test first.status in ("LOCALLY_SOLVED","OPTIMAL")
        @test !first.inexact_accepted && first.cap_attempt_iterations==0
        capped=stage_step!(s,d,[2.0 0.2 1.7;0.3 1.8 0.2],u,100.0,2)
        @test capped.cap_attempt_iterations==1
        @test capped.inexact_accepted==expect_accept
        @test capped.cap_recovery==!expect_accept
        @test isfinite(capped.cap_attempt_primal)
        @test capped.elapsed_s>=capped.local_primal_gate_s>=0
        if expect_accept
            @test capped.status=="ITERATION_LIMIT"
            @test capped.local_primal_evaluated && isfinite(capped.local_primal)
            @test screen=="finite" || capped.local_primal<=tolerance
            @test screen!="finite" || capped.local_primal>tolerance
            @test capped.iterations==1
        else
            @test capped.status in ("LOCALLY_SOLVED","OPTIMAL")
            @test capped.cap_attempt_primal>tolerance
            @test capped.iterations>1
        end
        final=stage_step!(s,d,fill(0.9,2,3),u,100.0,3)
        @test final.status in ("LOCALLY_SOLVED","OPTIMAL")
        @test !final.inexact_accepted && !final.cap_recovery && final.cap_attempt_iterations==0
        @test maximum(values(primal_feasibility_report(s.m;atol=0.0));init=0.0)<1e-7
    end
end

@testset "Persistent caps apply from first round with no exact recovery or switch" begin
    withenv("TADMM_CAP_ITERATIONS"=>"1","TADMM_CAP_SCREEN"=>"finite") do
        s=toy_stage("primal_dual","batched";policy="persistent_cap")
        d=Dict(:Bset=>[1,2]);limits=0
        for k in (1,2,31,1001)
            r=stage_step!(s,d,fill(isodd(k) ? 0.8 : 1.9,2,3),zeros(2,3),100.0,k)
            @test r.iterations<=1
            @test r.status in ("LOCALLY_SOLVED","OPTIMAL","ITERATION_LIMIT")
            @test !r.cap_recovery && !r.cold_retry
            @test !r.local_primal_evaluated
            @test r.inexact_accepted==r.iteration_limit_hit
            @test all(isfinite,value.(s.vars))
            @test MOI.get(s.m,MOI.RawOptimizerAttribute("max_iter"))==1
            limits+=r.iteration_limit_hit
        end
        @test limits>0
    end
end

@testset "Early cap accepts screened iterate, recovers rejection, restores exact solve" begin
    cap_gate_trial(1e3;expect_accept=true)
    cap_gate_trial(1e-20;expect_accept=false)
    cap_gate_trial(1e-20;expect_accept=true,screen="finite")
end

include("sweep_milestones.jl")
@testset "No final endpoint before exact phase" begin
    withenv("TADMM_LOCAL_SOLVE_POLICY"=>"early_cap30","TADMM_CAP_END_ITERATION"=>"30") do
        mktempdir() do folder
            observe,marks=sweep_observer(folder,Dict(:T=>24),"large10kC_1ph",1.0;save_points=false)
            metrics=(objective=1.0,original_primal=0.0,maxeq=0.0,maxbound=0.0,maxcone=0.0,
                soc_equality_gap=0.0,fictitious_loss_kwh=0.0)
            x=(metrics=metrics,owned_r=0.0,rr=0.0,owned_s=0.0,elapsed=0.0)
            @test !observe(merge(x,(k=30,)))
            @test !haskey(marks,"final_rms")
            @test observe(merge(x,(k=31,)))
            @test marks["final_rms"]["iteration"]==31
        end
    end
end
