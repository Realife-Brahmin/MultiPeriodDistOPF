# Verify FilterDDP against the centralized JuMP + Ipopt solution.
#
# Both halves run from the unified environment (see envs/ddp2026/README.md).
# Produce the centralized reference first, then diff against it:
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/copper_plate_centralized.jl
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/verify_against_centralized.jl
#
# The FilterDDP models are reused verbatim from copper_plate_battery_bounds.jl
# rather than restated here, so there is exactly one definition of the DDP model.
# Ipopt solves eq:cp_all in its original variables with the true per-period data;
# FilterDDP solves the optimal-control form with Lagrange-interpolated data and a
# slack reformulation of the energy bound. Agreement therefore tests the
# reformulation as well as the solver.

include("copper_plate_battery_bounds.jl")

const CSV_PATH = joinpath(@__DIR__, "..", "..", "results", "copper_plate",
                          "centralized_reference.csv")

"""Read centralized_reference.csv into case => (psub, pb, B, J, feasible)."""
function load_centralized(path)
    isfile(path) || error("""
        Centralized reference not found at
          $path
        Run copper_plate_centralized.jl under --project=envs/ddp2026 first.""")
    acc = Dict{String,Any}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue                       # header
        f = split(strip(line), ',')
        length(f) == 7 || continue
        name = f[1]; T = parse(Int, f[2]); t = parse(Int, f[3])
        e = get!(acc, name, Dict("T" => T, "psub" => Dict{Int,Float64}(),
                                 "pb" => Dict{Int,Float64}(), "B" => Dict{Int,Float64}(),
                                 "J" => NaN, "feasible" => true))
        if t == -1                                # infeasibility marker
            e["feasible"] = false
            continue
        end
        f[4] == "NaN" || (e["psub"][t] = parse(Float64, f[4]))
        f[5] == "NaN" || (e["pb"][t]   = parse(Float64, f[5]))
        f[6] == "NaN" || (e["B"][t]    = parse(Float64, f[6]))
        f[7] == "NaN" || (e["J"]       = parse(Float64, f[7]))
    end
    out = Dict{String,Any}()
    for (name, e) in acc
        T = e["T"]
        out[name] = (T = T, feasible = e["feasible"], J = e["J"],
                     psub = e["feasible"] ? [e["psub"][t] for t = 1:T] : Float64[],
                     pb   = e["feasible"] ? [e["pb"][t]   for t = 1:T] : Float64[],
                     B    = e["feasible"] ? [e["B"][t]    for t = 1:T+1] : Float64[])
    end
    return out
end

# Case name => (solver function, kwargs). Must mirror CASES in the centralized script.
const VERIFY_CASES = [
    ("base_T6", solve_ddp_bounds, 6, NamedTuple()),
    ("6a_T6",   solve_ddp_bounds, 6, (ps_lo = 0.0,)),
    ("6b_T6",   solve_ddp_bounds, 6, (pb_lo = -0.45, pb_hi = 0.45)),
    ("6c_T6",   solve_ddp_energy, 6, (B_lo = 1.20, B_hi = 1.95)),
    ("6d_T6",   solve_ddp_energy, 6, (ps_lo = 0.0, pb_lo = -0.45, pb_hi = 0.45,
                                      B_lo = 1.20, B_hi = 1.95)),
    ("6g_T6",   solve_ddp_energy, 6, (ps_lo = 0.0, pb_lo = 0.0, pb_hi = 0.30,
                                      B_lo = 2.05, B_hi = 2.10)),
]

function verify()
    ref = load_centralized(CSV_PATH)

    println("="^94)
    println("FilterDDP  vs  centralized JuMP/Ipopt solution of eq:cp_all")
    println("="^94)
    @printf("%-10s %6s %7s %14s %14s %11s %11s %11s\n",
            "case", "iters", "status", "J (FilterDDP)", "J (Ipopt)", "|ΔJ|",
            "max|ΔPsub|", "max|ΔB|")
    println("-"^94)

    worst = 0.0
    for (name, solver_fn, T, cfg) in VERIFY_CASES
        haskey(ref, name) || (println("$name: absent from the CSV, skipped"); continue)
        r = ref[name]
        g = solver_fn(T; cfg...)
        sd = g.solver.data
        c = tadmm_cost(T)
        J = sum(c .* g.psub .* Δt) + C_B * sum(g.pb .^ 2) * Δt + W * (g.B[T+1] - B_0)^2

        if !r.feasible
            @printf("%-10s %6d %7d %14.8f %14s %11s %11s %11s\n",
                    name, sd.k, sd.status, J, "INFEASIBLE", "-", "-", "-")
            continue
        end
        dJ  = abs(J - r.J)
        dPs = maximum(abs.(g.psub .- r.psub))
        dB  = maximum(abs.(g.B .- r.B))
        worst = max(worst, dPs, dB)
        @printf("%-10s %6d %7d %14.8f %14.8f %11.3e %11.3e %11.3e\n",
                name, sd.k, sd.status, J, r.J, dJ, dPs, dB)
    end

    println("-"^94)
    @printf("worst trajectory disagreement across all feasible cases: %.3e\n", worst)
    println("""
Note on 6g_T6: Ipopt returns LOCALLY_INFEASIBLE and so certifies the bound set
empty. FilterDDP returns status 7, which it also returns for a merely hard
problem — the caller must inspect residuals to tell the two apart. This is the
missing feasibility-restoration phase described in the experiment README.""")
end

verify()
