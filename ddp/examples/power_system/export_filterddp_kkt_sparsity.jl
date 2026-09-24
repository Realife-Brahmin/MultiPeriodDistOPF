# Extract per-stage KKT sparsity from FilterDDP's backward pass.
#
# Runs FilterDDP on a network instance for a few iterations with
# FILTERDDP_PERIODIC_CAPTURE_DIR, then reads back the serialized K matrices
# and writes (row,col) CSV files suitable for spy-plot rendering.
#
# Usage:
#   PROFILE_PERIODIC=1 TERMINAL_SOC_SOFT=1 \
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/export_filterddp_kkt_sparsity.jl \
#         [system] [T]

using Printf
using Serialization
using SparseArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T      = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3

capture_dir = joinpath(REPO_ROOT, "ddp", "results", "sparsity",
                       "filterddp_$(system)_T$(T)")
mkpath(capture_dir)

ENV["FILTERDDP_PERIODIC_CAPTURE_DIR"] = capture_dir
ENV["FILTERDDP_PERIODIC_CAPTURE_STRIDE"] = "1"
ENV["FILTERDDP_MAX_ITERATIONS"] = "3"
ENV["FILTERDDP_SKIP_SOLUTION_WRITE"] = "1"
ENV["PROFILE_PERIODIC"] = get(ENV, "PROFILE_PERIODIC", "1")
ENV["TERMINAL_SOC_SOFT"] = get(ENV, "TERMINAL_SOC_SOFT", "1")

@printf("Running FilterDDP for %s T=%d to capture per-stage KKT...\n", system, T)
flush(stdout)

include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))
main([system, string(T), "solve", "quiet"])

@printf("\nExtracting sparsity from captures in %s ...\n", capture_dir)
flush(stdout)

outdir = joinpath(REPO_ROOT, "ddp", "results", "sparsity")

function write_spy_csv(filename, S)
    open(filename, "w") do io
        println(io, "row,col")
        I, J, _ = findnz(S)
        for k in eachindex(I)
            println(io, I[k], ",", J[k])
        end
    end
end

for t in 1:T
    found = false
    for iter in 0:3
        capfile = joinpath(capture_dir, @sprintf("iter%04d_stage%03d.jls", iter, t))
        isfile(capfile) || continue

        cap = deserialize(capfile)
        K = cap.K
        nu, nc = cap.nu, cap.nc
        dim = size(K, 1)

        @printf("  stage %d (iter %d): K is %d x %d, nnz = %d (nu=%d, nc=%d)\n",
                t, iter, dim, dim, nnz(K), nu, nc)

        prefix = "$(system)_T$(T)_filterddp_stage$(t)"
        write_spy_csv(joinpath(outdir, "$(prefix)_kkt.csv"), K)

        H_block = K[1:nu, 1:nu]
        J_block = K[nu+1:end, 1:nu]
        write_spy_csv(joinpath(outdir, "$(prefix)_hessian.csv"), sparse(H_block))
        write_spy_csv(joinpath(outdir, "$(prefix)_jacobian.csv"), sparse(J_block))

        open(joinpath(outdir, "$(prefix)_dims.txt"), "w") do io
            @printf(io, "system=%s T=%d stage=%d nvar=%d ncon=%d nnz_H=%d nnz_J=%d nnz_KKT=%d\n",
                    system, T, t, nu, nc, nnz(sparse(H_block)), nnz(sparse(J_block)), nnz(K))
        end
        found = true
        break
    end
    found || @printf("  stage %d: no capture found\n", t)
end

@printf("done.\n")
