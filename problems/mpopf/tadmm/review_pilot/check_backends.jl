# A loading/correctness smoke check, not an OPF performance benchmark.
using JuMP, Ipopt, LinearAlgebra, TOML
BLAS.set_num_threads(1)
folder=get(ENV,"HSL_LIB_DIR","C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build")
ENV["PATH"]=folder*";"*ENV["PATH"]
results=Dict{String,Any}()
for backend in ("mumps","ma57","ma97")
    m=Model(Ipopt.Optimizer); set_silent(m)
    set_optimizer_attribute(m,"linear_solver",backend)
    set_optimizer_attribute(m,"linear_system_scaling","none")
    backend=="mumps" || set_optimizer_attribute(m,"hsllib",joinpath(folder,backend=="ma57" ? "libma57.dll" : "libhsl_ma97.dll"))
    @variable(m,x>=0)
    @variable(m,y>=0)
    @constraint(m,x+y==3)
    @objective(m,Min,(x-1)^2+(y-2)^2)
    optimize!(m)
    status=string(termination_status(m))
    ok=status=="LOCALLY_SOLVED" && abs(value(x)-1)<1e-6 && abs(value(y)-2)<1e-6
    results[backend]=Dict("status"=>status,"passed"=>ok)
    println("BACKEND_CHECK ",backend," ",results[backend]); flush(stdout)
    ok || error("backend preflight failed")
end
isempty(ARGS) || open(io->TOML.print(io,results),ARGS[1],"w")
