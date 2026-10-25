# Regression: a primary at 1.06 pu is allowed when its customer stays in band.
using Test
empty!(ARGS); push!(ARGS, "none")
include(joinpath(@__DIR__, "..", "root_level", "multiperiod_study.jl"))
commands = ["Clear", "New Circuit.bandtest phases=1 bus1=sourcebus basekv=11.547 pu=1.06 R1=0 X1=1e-8 R0=0 X0=1e-8",
    "New Line.feed phases=1 bus1=sourcebus bus2=primary R1=0.001 X1=0.001 C1=0 C0=0",
    "New Line.svc_customer phases=1 bus1=primary bus2=customer R1=2 X1=0.01 C1=0 C0=0",
    "New Load.customer phases=1 bus1=customer kV=11.547 kW=1000 kvar=0 model=1 Vminpu=0.5 Vmaxpu=1.5",
    "New Storage.idle phases=1 bus1=primary kV=11.547 kWrated=1 kWhrated=4 kVA=1 %idlingkW=0 state=idling"]
compile_deck("small2poi_1ph"; extra=commands, snapshot=true, load_band=(0.5, 1.5))
net = read_network()
s = setup(net, []; T=1)
@test s.band_buses == Set(["customer"])
@test study_banded(s) == ["customer"]
@test problem(s, :cyclic).band_buses == s.band_buses
@test problem(s, :nobatt; hour=1).band_buses == s.band_buses
ang = equal_angles(net)
x = solve_case(s, ang, :cyclic)
@test x.ok
@test abs(x.V_pu["primary"][1]) > 1.05
@test 0.95 <= abs(x.V_pu["customer"][1]) <= 1.05
rp = replay_opendss(net, x, ang; load_factor=s.lam)
@test rp.converged && rp.pq_band
@test rp.dP < 1e-3 && rp.dV < 1e-7
# A margin solve exposes the legacy primary-band violation without a slow infeasibility proof.
sa = setup(net, []; T=1, band_buses=nothing)
y = solve_case(sa, ang, :cyclic; objective=:margin)
@test y.status in RAY_OK
@test y.max_violation > 0.01
rpa = replay_opendss(net, y, ang; load_factor=sa.lam)
@test rpa.converged && rpa.dP < 1e-3 && rpa.dV < 1e-7
@test Set(study_banded(sa)) == Set(["primary", "customer"])
# The existing primary-only deck keeps its historical behavior.
legacy, pv = load_sources(["subs3", "subs4"]; system="ieee123_5poi_1ph")
sl = setup(legacy, pv; T=1)
@test sl.band_buses === nothing
@test Set(study_banded(sl)) == Set(banded(legacy))
println("PASS: customer-only band, hourly propagation, margin diagnostics, legacy behavior; both solves replayed in OpenDSS")

# Converter variant commands and new stage names must resolve without launching a study.
@test haskey(MULTI, "four_ieee9500xf")
@test MULTI["four_ieee9500xf"].add == ["l2804256"]
for (key, off) in (("four_smartds_p1u_typeiii", "subs5"), ("four_smartds_p1u_typeii", "subs3"))
    mc = MULTI[key]
    @test length(mc.sources) == 4
    @test !(off in mc.sources)
    @test "Vsource.$off.enabled=false" in variant_switches(mc.system, mc.variant)
end
primary_lines(deck) = filter(l -> startswith(lowercase(l), "new line.") && !startswith(lowercase(l), "new line.svc_"),
    readlines(joinpath(REPO_ROOT, "test_networks", deck, "BranchData.dss")))
@test primary_lines("ieee9500_1ph") == primary_lines("ieee9500xf_1ph")
println("PASS: both P1U switching variants; IEEE9500 fourth-site transfer preserves all primary line impedances")
