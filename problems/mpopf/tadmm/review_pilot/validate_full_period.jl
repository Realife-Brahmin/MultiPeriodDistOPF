include("run.jl")
using Test
T=parse(Int,ARGS[1]);out=ARGS[2];mkpath(out);register_process!()
@assert parse(Int,ENV["TADMM_NWORKERS"])==T
d=deserialize(joinpath(SHARED,"..","..","results/network_filterddp/network_data_large10kC_1ph_T$(T)_periodic.jls"))
apply_review_cb!(d,"large10kC_1ph";profile="july")
ENV["TADMM_CHECKPOINT_INTERVAL"]="0"
try
    for (name,warm,update) in [("legacy","primal","legacy"),("batched","primal","batched"),("dual","primal_dual","batched")]
        ENV["TADMM_WARM_START"]=warm;ENV["TADMM_OBJECTIVE_UPDATE"]=update
        r=run_admm(d;maxiter=8,rho0=625.0T,out=joinpath(out,name*".csv"),phase="gate_"*name)
        @test length(POOL_LAYOUT)==T
        @test all(w->length(w["periods"])==1,POOL_LAYOUT)
        @test length(unique((w["hostname"],w["affinity"]) for w in POOL_LAYOUT))==T
        @test length(unique(w["hostname"] for w in POOL_LAYOUT))==parse(Int,ENV["SLURM_JOB_NUM_NODES"])
        open(io->TOML.print(io,Dict("workers"=>deepcopy(POOL_LAYOUT),"core_s"=>r.elapsed,"iterations"=>r.k)),joinpath(out,name*".toml"),"w")
        serialize(joinpath(out,name*"-dispatch.jls"),Dict(:system=>"large10kC_1ph",:T=>T,:objective=>r.metrics.objective,:u=>r.us,:B=>r.bb,:status=>"SHORT_EXECUTION_GATE"))
        check=`$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) $(joinpath(SHARED,"opendss_check.jl")) large10kC_1ph $T $(joinpath(out,name*"-dispatch.jls")) $name`
        open(joinpath(out,name*"-opendss.log"),"w") do io
            @test success(pipeline(check;stdout=io,stderr=io))
        end
    end
finally
    shutdown_pool!()
end
# Only the numerical no-op objective change should retain the same trajectory.
function columns(path)
    lines=readlines(path);header=split(first(lines),',')
    wanted=[findfirst(==(key),header) for key in ("objective","original_primal","owned_r_rms","r_rms","rho")]
    [parse.(Float64,split(line,',')[wanted]) for line in lines[2:end]]
end
a=columns(joinpath(out,"legacy.csv"));b=columns(joinpath(out,"batched.csv"))
@test length(a)==length(b)
@test all(isapprox(x,y;atol=1e-9,rtol=1e-10) for (x,y) in zip(a,b))
write(joinpath(out,"FULL_PERIOD_GATE_VALIDATED"),"One worker per period; unique host/CPU; exact local statuses checked separately; objective-update path equivalence. Short gate, not a convergence benchmark.\n")
