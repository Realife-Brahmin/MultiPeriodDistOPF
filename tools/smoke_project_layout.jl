# Small import/export regression plus a scalar Ipopt/MUMPS check; no MPOPF run.
import Pkg
Pkg.instantiate(; allow_autoprecomp=false)
using Test, Serialization, OpenDSSDirect, Parameters, SPARDDP, JuMP, Ipopt
const ROOT = normpath(joinpath(@__DIR__, ".."))
const DDP = joinpath(ROOT, "problems", "mpopf", "spar_ddp")
module SharedParser
    include(joinpath(@__DIR__, "..", "test_networks", "parse_opendss.jl"))
end
module LegacyParser
    # The unchanged pre-migration numerical parser, evaluated only for parity.
    const ROOT = normpath(joinpath(@__DIR__, ".."))
    const old = read(`git -C $ROOT show 5742a5d1:envs/tadmm/parse_opendss.jl`, String)
    include_string(@__MODULE__, old, "legacy_parse_opendss.jl")
end
@testset "Network parser parity" begin
    kwargs = (; rawDataFolderPath=joinpath(ROOT, "test_networks"),
              LoadShapeLoad=[0.9,1.0,0.8], LoadShapeCost=[0.1,0.2,0.12],
              LoadShapePV=[0.0,1.0,0.0], kV_B=2.4018)
    old = LegacyParser.parse_system_from_dss("ieee123C_1ph",3;kwargs...)
    new = SharedParser.parse_system_from_dss("ieee123C_1ph",3;kwargs...)
    @test isequal(old, new)
    @test new[:T] == 3
end
@testset "SPAR-DDP export using only its environment" begin
    driver=joinpath(DDP,"examples","power_system","export_ieee123c_data.jl")
    # Export locally generated smoke input only, never a frozen campaign file.
    mktempdir() do output
        withenv("PROFILE_PERIODIC"=>"1", "FILTERDDP_EXPORT_DIR"=>output) do
            run(`$(Base.julia_cmd()) --startup-file=no --project=$DDP $driver ieee123C_1ph 3`)
        end
        file=joinpath(output,"network_data_ieee123C_1ph_T3_periodic.jls")
        data=deserialize(file)
        @test data[:T]==3
        @test data[:profile_periodic]
        @test !isempty(data[:Bset])
    end
end
@testset "Centralized reference needs no proprietary solver" begin
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "linear_solver", "mumps")
    @variable(model, x >= 0, start=0.5)
    @objective(model, Min, (x-1)^2)
    optimize!(model)
    @test termination_status(model) == JuMP.MOI.LOCALLY_SOLVED
    @test isapprox(value(x),1.0;atol=1e-5)
end
println("PROJECT_LAYOUT_SMOKE_PASS (imports, parser parity, export and scalar MUMPS check; no MPOPF benchmark)")
