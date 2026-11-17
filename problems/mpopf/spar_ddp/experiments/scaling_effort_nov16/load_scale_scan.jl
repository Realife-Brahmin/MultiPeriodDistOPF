# Lowest and highest bus voltage of a single-substation deck at load scale
# lambda (PV and batteries off, snapshot power flow), for choosing how much load
# one substation can carry within 0.95-1.05 pu.
#   julia --project=problems/mpopf/spar_ddp .../load_scale_scan.jl <deck> [lambda ...]
using OpenDSSDirect, Printf
const ODD = OpenDSSDirect
const REPO = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
deck = ARGS[1]
lams = length(ARGS) > 1 ? parse.(Float64, ARGS[2:end]) : [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]
ODD.Text.Command("clear")
ODD.Text.Command("redirect [$(joinpath(REPO, "test_networks", deck, "Master.dss"))]")
ODD.Text.Command("BatchEdit PVSystem..* enabled=false")
ODD.Text.Command("BatchEdit Storage..* enabled=false")
ODD.Text.Command("Set mode=snapshot")
ODD.Text.Command("Set maxiterations=200")
for lam in lams
    ODD.Text.Command("Set LoadMult=$lam")
    ODD.Text.Command("Solve")
    ok = ODD.Solution.Converged()
    v = ODD.Circuit.AllBusMagPu()
    @printf("%s lambda=%.2f converged=%s lowest %.4f highest %.4f pu import %.0f kW\n", deck, lam, ok,
            ok ? minimum(v) : NaN, ok ? maximum(v) : NaN, ok ? -real(first(ODD.Circuit.TotalPower())) : NaN)
end
