# Small import/export regression for the problem-first layout; no optimization.
import Pkg
Pkg.instantiate(; allow_autoprecomp=false)
using Test, Serialization, OpenDSSDirect, Parameters, DDP4OPF, JuMP, Ipopt
const ROOT = normpath(joinpath(@__DIR__, ".."))
const DDP = joinpath(ROOT, "problems", "mpopf", "spar_ddp")
module SharedParser
    include(joinpath(@__DIR__, "..", "test_networks", "parse_opendss.jl"))
end
module LegacyParser
    # The unchanged pre-migration numerical parser, evaluated only for parity.
    const ROOT = normpath(joinpath(@__DIR__, ".."))
    const old = read(`git -C $ROOT show 0efa8015:envs/tadmm/parse_opendss.jl`, String)
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
println("PROJECT_LAYOUT_SMOKE_PASS (imports, parser parity and export; no optimizer run)")
