using Test
include("model.jl")

# A nonlinear cone, equalities and lower/upper bounds exercise the multiplier
# classes. Compare the same parameter sequence, including unchanged/changed rho.
function toy_stage(warm,update)
    m=Model(Ipopt.Optimizer);set_silent(m)
    set_optimizer_attribute(m,"linear_solver",get(ENV,"TADMM_LINEAR_SOLVER","mumps"))
    for kv in split(get(ENV,"TADMM_IPOPT_OPTIONS",""),';';keepempty=false)
        key,val=split(kv,'=';limit=2)
        set_optimizer_attribute(m,key,something(tryparse(Int,val),tryparse(Float64,val),String(val)))
    end
    set_optimizer_attribute(m,"tol",1e-10)
    set_optimizer_attribute(m,"bound_relax_factor",0.0)
    @variable(m,0.1<=b[1:2,1:3]<=3,start=1)
    @variable(m,0.1<=v<=3,start=1)
    @variable(m,0.1<=ell<=3,start=1)
    @constraint(m,b[1,1]==1)
    @constraint(m,b[2,1]+b[2,2]==2)
    @constraint(m,b[1,2]^2+b[2,2]^2<=v*ell)
    @objective(m,Min,50sum(b.^2)+v+ell)
    vars=all_variables(m)
    refs=all_constraints(m;include_variable_in_set_constraints=true)
    reuse=StageReuse(100.0,warm,update,VariableRef[b[j,t] for t in 1:3 for j in 1:2],zeros(6),refs,0)
    return (;m,vars,reuse,b,times=collect(1:3),owned=[2],local_policy="exact",local_primal_report=true)
end

function trial(warm,update)
    s=toy_stage(warm,update);history=[]
    for (k,rho) in enumerate([100.0,100.0,20.0,20.0,50.0])
        z=fill(0.8+0.01k,2,3)
        r=stage_step!(s,Dict(:Bset=>[1,2]),z,zeros(2,3),rho,k)
        @test r.status in ("LOCALLY_SOLVED","OPTIMAL")
        @test !r.iteration_limit_hit
        @test maximum(values(primal_feasibility_report(s.m;atol=0.0));init=0.0)<1e-7
        for a in 1:3,n in 1:2
            @test coefficient(objective_function(s.m),s.b[n,a])≈-rho*z[n,a]
            @test coefficient(objective_function(s.m),s.b[n,a],s.b[n,a])≈rho/2
        end
        push!(history,(value.(s.vars),objective_value(s.m)))
    end
    @test s.reuse.quadratic_updates==(update=="batched" ? 2 : 5)
    return history
end

@testset "Repeated exact local solves and batched objective updates" begin
    legacy=trial("primal","legacy")
    batched=trial("primal","batched")
    warm=trial("primal_dual","batched")
    for k in eachindex(legacy)
        @test legacy[k][1]≈batched[k][1] atol=1e-7 rtol=1e-7
        @test legacy[k][2]≈batched[k][2] atol=1e-7 rtol=1e-7
        @test legacy[k][1]≈warm[k][1] atol=1e-6 rtol=1e-6
        @test legacy[k][2]≈warm[k][2] atol=1e-7 rtol=1e-7
    end
end
