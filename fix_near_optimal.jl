#!/usr/bin/env julia
# Recompute near-optimal iterations with corrected criterion: r <= eps_pri required

using Printf

const EPS_PRI = 2e-3
const GAP_TOL = 0.005

function parse_csv(filename)
    lines = readlines(filename)
    header = split(lines[1], ",")
    data = Dict(col => Float64[] for col in header)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, ",")
        for (i, col) in enumerate(header)
            try push!(data[col], parse(Float64, vals[i]))
            catch; push!(data[col], NaN); end
        end
    end
    return data
end

function find_bf_objective(dir)
    for depth in 0:3
        bf_file = joinpath(dir, "results_socp_bf.txt")
        if isfile(bf_file)
            for line in eachline(bf_file)
                m = match(r"Total Cost: \$([0-9.]+)", line)
                !isnothing(m) && return parse(Float64, m.captures[1])
            end
        end
        dir = dirname(dir)
    end
    return NaN
end

function compute_near_optimal(csv_path, bf_obj)
    data = parse_csv(csv_path)
    obj = data["objective"]
    r = data["r_norm"]

    eff_col = haskey(data, "eff_time_this_k") ? "eff_time_this_k" :
              haskey(data, "eff_time_s") ? "eff_time_s" : nothing
    ref = isnan(bf_obj) ? obj[end] : bf_obj

    cum_eff = 0.0
    for k in 1:length(obj)
        if !isnothing(eff_col)
            cum_eff += data[eff_col][k]
        end
        gap = abs(obj[k] - ref) / abs(ref)
        if gap <= GAP_TOL && r[k] <= EPS_PRI
            return (k=k, eff_time=cum_eff, gap=gap*100, r=r[k], s=data["s_norm"][k], obj=obj[k])
        end
    end
    return nothing
end

function rewrite_results(results_path, near_opt, bf_obj)
    lines = readlines(results_path)
    out = String[]
    bf_time = NaN

    for line in lines
        m_bf = match(r"Wall-clock time: ([0-9.]+) seconds", line)
        if !isnothing(m_bf) && isnan(bf_time)
            bf_time = parse(Float64, m_bf.captures[1])
        end
    end

    for (i, line) in enumerate(lines)
        if contains(line, "Effective time (near-optimal)")
            if isnothing(near_opt)
                push!(out, "  Effective time (near-optimal): N/A (0.5% gap with r<=eps_pri never reached)")
            else
                ref_label = isnan(bf_obj) ? "final tADMM" : "BF"
                push!(out, @sprintf("  Effective time (near-optimal): %.4f seconds (k=%d, gap=%.3f%% vs %s, r=%.2e, s=%.2e)",
                    near_opt.eff_time, near_opt.k, near_opt.gap, ref_label, near_opt.r, near_opt.s))
            end
        elseif contains(line, "tADMM eff(near-optimal)")
            if isnothing(near_opt) || isnan(bf_time) || near_opt.eff_time == 0
                push!(out, "  tADMM eff(near-optimal) vs BF: N/A")
            else
                push!(out, @sprintf("  tADMM eff(near-optimal) vs BF: %.2fx", bf_time / near_opt.eff_time))
            end
        else
            push!(out, line)
        end
    end

    open(results_path, "w") do io
        println(io, join(out, "\n"))
    end
end

function rewrite_near_optimal_csv(csv_path, near_opt, bf_obj)
    !isfile(csv_path) && return
    lines = readlines(csv_path)
    length(lines) < 2 && return

    header = lines[1]
    if isnothing(near_opt)
        open(csv_path, "w") do io
            println(io, header)
            cols = split(header, ",")
            vals = join(fill("N/A", length(cols)), ",")
            println(io, vals)
        end
    else
        ref_label = isnan(bf_obj) ? "final_tADMM" : "BF"
        open(csv_path, "w") do io
            println(io, header)
            @printf(io, "%d,%.4f,%.6f,%.3f,%s,%.2e,%.2e\n",
                near_opt.k, near_opt.eff_time, near_opt.obj, near_opt.gap, ref_label, near_opt.r, near_opt.s)
        end
    end
end

# Find all result directories
base = "envs/tadmm/processedData"
dirs = String[]
for d in readdir(base, join=true)
    isdir(d) || continue
    isfile(joinpath(d, "results_socp_tadmm.txt")) && push!(dirs, d)
    sweep = joinpath(d, "rho_sweep")
    if isdir(sweep)
        for sd in readdir(sweep, join=true)
            isdir(sd) && isfile(joinpath(sd, "results_socp_tadmm.txt")) && push!(dirs, sd)
        end
    end
end

println("Recomputing near-optimal for $(length(dirs)) result directories...\n")

for dir in sort(dirs)
    csv = joinpath(dir, "convergence_data.csv")
    results = joinpath(dir, "results_socp_tadmm.txt")
    near_csv = joinpath(dir, "near_optimal_summary.csv")

    !isfile(csv) && continue
    !isfile(results) && continue

    # Skip if results file doesn't have near-optimal line
    txt = read(results, String)
    !contains(txt, "near-optimal") && continue

    bf_obj = find_bf_objective(dir)
    near_opt = compute_near_optimal(csv, bf_obj)

    short = replace(dir, "envs/tadmm/processedData/" => "")
    if isnothing(near_opt)
        println("  $short: near-optimal → N/A (never reached)")
    else
        @printf("  %s: near-optimal → k=%d (gap=%.3f%%, r=%.2e)\n", short, near_opt.k, near_opt.gap, near_opt.r)
    end

    rewrite_results(results, near_opt, bf_obj)
    isfile(near_csv) && rewrite_near_optimal_csv(near_csv, near_opt, bf_obj)
end

println("\nDone.")
