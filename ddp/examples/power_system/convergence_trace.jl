# Stage 9: per-iteration convergence trace for FilterDDP on the copper plate.
#
# FilterDDP's `SolverData` holds only the CURRENT iterate -- there is no history
# buffer anywhere in the package (src/solver_data.jl:8-32). The only per-iteration
# record it produces is the verbose print stream, so this script runs each case
# with `verbose = true`, captures stdout, and parses the table back out. The
# upstream clone stays unmodified.
#
# Columns, and which half of the iteration produces them (src/solve.jl:20-48):
#
#   BACKWARD pass, evaluated at the current iterate before any step is taken
#     pr_inf   ∞-norm of the constraint violation      (primal infeasibility)
#     du_inf   ∞-norm of the Lagrangian gradient       (dual infeasibility)
#     cs_inf   ∞-norm of the complementarity error
#     lg(reg)  log10 of the inertia regularisation applied to the stage KKT system
#     lg(mu)   log10 of the barrier parameter for the current subproblem
#   FORWARD pass, the filter line search that follows
#     alpha    accepted step size
#     ls       data.l -- see the warning below, this is NOT the trial count
#
# WARNING -- `data.l` does not count line-search trials, despite the "ls" header.
# `forward_pass!` backtracks from three places (src/forward_pass.jl):
#     line 17  rollout! failed          -> step_size *= 0.5, NO increment of l
#     line 25  rejected by the filter   -> step_size *= 0.5, l += 1
#     line 37  insufficient decrease    -> step_size *= 0.5, l += 1
# so a backtrack caused by the rollout (the fraction-to-boundary rule) is
# invisible in `l`. Every trace below shows alpha = 0.5 alongside ls = 0, which
# is exactly that case and would read as a contradiction otherwise. The reliable
# quantity is alpha itself, so the table reports the true backtrack count as
# log2(1/alpha) and relabels l as what it actually is: filter/decrease
# rejections only.
#
# IMPORTANT -- the one-iteration stagger. `iteration_status` is called BEFORE
# `forward_pass!` (src/solve.jl:38 then :40), so in the raw print stream the
# alpha and ls on row k are the forward pass that PRECEDED it, i.e. iteration
# k-1's. This is noted in README_FILTERDDP_EXPERIMENT.md as an incidental
# observation; here it actually matters, because a table that pairs a backward
# pass with the wrong forward pass is simply wrong. The parser therefore shifts
# alpha and ls up by one row, so every row holds one iteration's backward pass
# and the forward pass it actually fed. The last row's forward pass does not
# exist (the solver converged or stopped instead of stepping) and is left blank.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/convergence_trace.jl

include("copper_plate_battery_bounds.jl")

# ---------------------------------------------------------------------------
# Capture and parse the verbose stream
# ---------------------------------------------------------------------------

"""Run `f` with stdout captured, returning what it printed. The reader runs as a
task so a full pipe buffer cannot deadlock the solve."""
function capture_stdout(f)
    old = stdout
    rd, wr = redirect_stdout()
    reader = @async read(rd, String)
    try
        f()
    finally
        redirect_stdout(old)
        close(wr)
    end
    return fetch(reader)
end

"""
Parse FilterDDP's verbose iteration table.

Row format is fixed by src/print.jl:24 --
    iter  objective  pr_inf  du_inf  cs_inf  lg(mu)  lg(reg)  alpha  ls
with lg(reg) printed as a bare "-" when no regularisation fired.
"""
function parse_trace(out::String)
    rows = NamedTuple[]
    for line in split(out, '\n')
        f = split(strip(line))
        length(f) == 9 || continue
        all(!isempty, f) || continue
        occursin(r"^\d+$", f[1]) || continue          # skip the header line
        k = parse(Int, f[1])
        vals = try
            (obj = parse(T_, f[2]), pr = parse(T_, f[3]), du = parse(T_, f[4]),
             cs = parse(T_, f[5]), lgmu = parse(T_, f[6]),
             lgreg = f[7] == "-" ? NaN : parse(T_, f[7]),
             alpha = parse(T_, f[8]), ls = parse(Int, f[9]))
        catch
            continue
        end
        push!(rows, merge((k = k,), vals))
    end
    # de-duplicate on k, keeping the last print for each (the post-loop call to
    # iteration_status re-emits the final iterate)
    seen = Dict{Int,Int}()
    for (i, r) in enumerate(rows)
        seen[r.k] = i
    end
    rows = [rows[seen[k]] for k in sort(collect(keys(seen)))]

    # Undo the stagger: alpha/ls printed on row k belong to iteration k-1.
    out_rows = NamedTuple[]
    for i = 1:length(rows)
        r = rows[i]
        nxt = i < length(rows) ? rows[i+1] : nothing
        push!(out_rows, (k = r.k, obj = r.obj, pr = r.pr, du = r.du, cs = r.cs,
                         lgmu = r.lgmu, lgreg = r.lgreg,
                         alpha = nxt === nothing ? NaN : nxt.alpha,
                         ls    = nxt === nothing ? -1  : nxt.ls))
    end
    return out_rows
end

"""Solve one case with verbose on and return its parsed trace."""
function trace_case(solver_fn, T::Int; cfg...)
    local rows
    out = capture_stdout() do
        solver_fn(T; cfg..., verbose = true)
    end
    return parse_trace(out)
end

# ---------------------------------------------------------------------------
# Cases to trace
# ---------------------------------------------------------------------------

const TRACE_CASES = [
    ("6e_T6", solve_ddp_energy, 6, (ps_lo = 0.0, pb_lo = -0.45, pb_hi = 0.45,
                                    B_lo = 1.20, B_hi = 1.95)),
    ("6d_T3", solve_ddp_energy, 3, (ps_lo = 0.0, pb_lo = -0.5, pb_hi = 0.004,
                                    B_lo = 1.9895, B_hi = 1.9965)),
]

function main()
    outdir = joinpath(@__DIR__, "..", "..", "results", "copper_plate")
    mkpath(outdir)
    csv = joinpath(outdir, "convergence.csv")
    lines = ["case,k,objective,pr_inf,du_inf,cs_inf,lg_mu,lg_reg,alpha,ls"]

    figdir = joinpath(@__DIR__, "..", "..", "paper", "figures")

    for (name, fn, T, cfg) in TRACE_CASES
        rows = trace_case(fn, T; cfg...)
        # True backtrack count, which data.l does not give (see header).
        backtracks(a) = isnan(a) ? -1 : round(Int, log2(1 / a))

        println("\n" * "="^100)
        println("$name — FilterDDP convergence, $(length(rows)) iterations (k = 0..$(rows[end].k))")
        println("="^100)
        println("            |------------------ backward pass ------------------|------ forward ------|")
        @printf("%5s %16s %11s %11s %11s %8s %8s %9s %6s %6s\n",
                "iter", "objective", "pr_inf", "du_inf", "cs_inf", "lg(mu)", "lg(reg)",
                "alpha", "bktrk", "filt")
        println("-"^100)
        for r in rows
            @printf("%5d %16.10f %11.3e %11.3e %11.3e %8.2f %8s %9s %6s %6s\n",
                    r.k, r.obj, r.pr, r.du, r.cs, r.lgmu,
                    isnan(r.lgreg) ? "-" : @sprintf("%.2f", r.lgreg),
                    isnan(r.alpha) ? "-" : @sprintf("%.3e", r.alpha),
                    backtracks(r.alpha) < 0 ? "-" : string(backtracks(r.alpha)),
                    r.ls < 0 ? "-" : string(r.ls))
            push!(lines, string(name, ",", r.k, ",", r.obj, ",", r.pr, ",", r.du, ",",
                                r.cs, ",", r.lgmu, ",", r.lgreg, ",", r.alpha, ",", r.ls))
        end
        nreg  = count(r -> !isnan(r.lgreg), rows)
        nback = count(r -> backtracks(r.alpha) > 0, rows)
        nfilt = count(r -> r.ls > 0, rows)
        @printf("\nregularisation fired on %d of %d iterations\n", nreg, length(rows))
        @printf("%d iterations backtracked (alpha < 1); of those, %d were filter/decrease rejections\n",
                nback, nfilt)
        println("  -> the rest were rollout rejections, which data.l does not count (forward_pass.jl:17)")
        @printf("pr_inf %.2e -> %.2e,  du_inf %.2e -> %.2e,  cs_inf %.2e -> %.2e\n",
                rows[1].pr, rows[end].pr, rows[1].du, rows[end].du, rows[1].cs, rows[end].cs)

        # Figure data: one file per case, log-plottable (no zeros).
        floor_(v) = max(v, 1e-17)
        open(joinpath(figdir, "convergence_$(name).csv"), "w") do io
            println(io, "k,pr,du,cs,mu,obj")
            for r in rows
                println(io, r.k, ",", floor_(r.pr), ",", floor_(r.du), ",",
                        floor_(r.cs), ",", 10.0^r.lgmu, ",", r.obj)
            end
        end
        println("wrote ", joinpath(figdir, "convergence_$(name).csv"))

        # The COMPLETE tabular, generated rather than hand-copied. It emits the
        # rules itself instead of leaving \bottomrule to the wrapper: a \noalign
        # (which is what \bottomrule expands to) placed immediately after an
        # \input boundary inside a tabular is a misplaced \noalign, because the
        # end of the included file starts a cell. Keeping the whole alignment in
        # one file sidesteps that entirely.
        open(joinpath(figdir, "convergence_table_$(name).tex"), "w") do io
            println(io, "% Generated by convergence_trace.jl -- do not edit by hand.")
            println(io, "\\begin{tabular}{rcccccrc}")
            println(io, "\\toprule")
            println(io, "& & \\multicolumn{4}{c}{\\textit{backward pass}} & \\multicolumn{2}{c}{\\textit{forward pass}} \\\\")
            println(io, "\\cmidrule(lr){3-6}\\cmidrule(lr){7-8}")
            println(io, "\$k\$ & objective & \$\\mathrm{pr}_\\infty\$ & \$\\mathrm{du}_\\infty\$ & \$\\mathrm{cs}_\\infty\$ & \$\\log_{10}\\mu\$ & \$\\alpha\$ & bk \\\\")
            println(io, "\\midrule")
            for r in rows
                bk = backtracks(r.alpha)
                @printf(io, "\$%d\$ & \$%.6f\$ & \$%.1e\$ & \$%.1e\$ & \$%.1e\$ & \$%.2f\$ & %s & %s \\\\\n",
                        r.k, r.obj, r.pr, r.du, r.cs, r.lgmu,
                        isnan(r.alpha) ? "---" : @sprintf("\$%.2f\$", r.alpha),
                        bk < 0 ? "---" : "\$$bk\$")
            end
            println(io, "\\bottomrule")
            println(io, "\\end{tabular}")
        end
        println("wrote ", joinpath(figdir, "convergence_table_$(name).tex"))
    end

    open(csv, "w") do io
        for l in lines
            println(io, l)
        end
    end
    println("\nwrote $csv")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
