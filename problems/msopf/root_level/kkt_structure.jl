# Structure of the KKT matrix Clarabel factors for the dispatched certificate (2026-11-01). DIAGNOSTIC: no solve.
# Builds the CD relaxation at one point (direction FB_DIR, FB_RATIO x the dispatched reach), lets Clarabel assemble
# its KKT matrix (max_iter = 0) and compares orderings of that same sparsity pattern:
#   amd    QDLDL's AMD ordering, as Clarabel uses it (one sequential LDL' factorization)
#   hours  the variables that couple hours (B, P_B, m) taken out; each remaining connected block (one per hour,
#          if the coupling is only through those) ordered by AMD on its own, the coupling variables last
# For each: nonzeros of L, factorization work (sum of squared column counts) and the critical path of the
# elimination tree (work on its heaviest leaf-to-root path). Work / critical path bounds the speedup of a
# factorization that runs independent subtrees in parallel; for hours, the share of the work left in the
# coupling block says how much stays sequential.
#   FB_SYSTEM=sfp1u_four FB_DIR=24 FB_RATIO=1.1 FB_IDLE=... MULTIPERIOD_OUT_DIR=<campaign folder> [FB_KKT_T=2] julia ... kkt_structure.jl
ENV["FB_PHASE"] = "none"
include(joinpath(@__DIR__, "farkas_dispatched.jl"))
using SparseArrays
const QD = Clarabel.QDLDL; const AMDm = QD.AMD
amdp(M) = AMDm.amd(SparseMatrixCSC{Float64,Int}(size(M)..., M.colptr, M.rowval, ones(nnz(M))))

i = parse(Int, ENV["FB_DIR"]); d = DIRS[i]; ratio = parse(Float64, get(ENV, "FB_RATIO", "1.1"))
rs = parse(Float64, split(readlines(joinpath(OUT, "fb_AD_$(d.id).csv"))[2], ',')[4]); r = ratio * rs
amp = Dict{String,Float64}(); senv = Dict{String,Float64}()
for l in readlines(joinpath(OUT, "fb_envelope_disp.csv"))[2:end]
    c = split(l, ','); c[2] == "line" ? (amp[String(c[1])] = parse(Float64, c[4])) : (senv[String(c[1])] = parse(Float64, c[4]))
end
TK = parse(Int, get(ENV, "FB_KKT_T", string(kwb.T)))     # fewer hours (the first TK of the day) for a quick test
kwL = merge(kwb, (; T = TK, ampacity_pu = amp, source_ampacity_pu = senv))
say(@sprintf("%s: %.4f x reach %.5f = %.5f deg; cuts %s; T = %d", d.label, ratio, rs, r, join(CUTS, ","), TK))

t0 = time()
R = relax_margin(net, angles(d, r), kwL; k = KV, cuts = CUTS, solver_attributes = [CLARABEL_DEFAULTS; ("max_iter", 0)])
nv = num_variables(R.model)
# coupling variables by name; Clarabel's columns follow the model's variable order (checked below by the split)
sepvars = [v for v in all_variables(R.model) if (nm = name(v); startswith(nm, "B[") || startswith(nm, "P_B[") || nm == "m")]
sepidx = [index(v).value for v in sepvars]
optimize!(R.model)
S = JuMP.unsafe_backend(R.model).solver; ks = S.kktsystem.kktsolver
K = ks.KKT; F = ks.ldlsolver.factors
N = size(K, 1); nnzL_qdldl = nnz(F.L); pamd = copy(F.perm)
say(@sprintf("  KKT: %d x %d (n = %d variables, model has %d), %d nonzeros in its upper triangle; QDLDL L %d nonzeros; setup %.0f s, max RSS %.1f GB",
             N, N, ks.n, nv, nnz(K), nnzL_qdldl, time() - t0, Sys.maxrss() / 2^30))
ks.n == nv || say("  WARNING: Clarabel has a different number of variables than the model; the hour split below may be wrong")
Kp = SparseMatrixCSC{Bool,Int}(N, N, K.colptr, K.rowval, fill(true, nnz(K)))
R = nothing; S = nothing; ks = nothing; K = nothing; F = nothing; GC.gc(true)
Ksym = Kp .| permutedims(Kp); Kp = nothing; GC.gc(true)

"Column counts of L, elimination tree, total work and critical path for ordering p."
function analyse(Ksym, p)
    n = size(Ksym, 1)
    T = triu(Ksym[p, p]); T = SparseMatrixCSC{Float64,Int}(n, n, T.colptr, T.rowval, ones(nnz(T)))
    work = zeros(Int, n); Lnz = zeros(Int, n); et = zeros(Int, n)
    s = QD.QDLDL_etree!(n, T.colptr, T.rowval, work, Lnz, et)
    s < 0 && error("etree failed")
    w = Float64.(Lnz) .^ 2; best = zeros(n); cp = 0.0
    for j in 1:n
        c = w[j] + best[j]; cp = max(cp, c)
        et[j] > 0 && (best[et[j]] = max(best[et[j]], c))
    end
    return (; nnzL = s, work = sum(w), crit = cp, Lnz)
end
A = analyse(Ksym, pamd)
say(@sprintf("  amd:   L %d nonzeros (QDLDL's own count %d), work %.3e, critical path %.3e, parallel bound %.1fx",
             A.nnzL, nnzL_qdldl, A.work, A.crit, A.work / A.crit))

# hours: connected components of the pattern once the coupling columns are removed
function components(Ksym, sep)
    n = size(Ksym, 1); comp = zeros(Int, n); nc = 0; q = Int[]
    for s0 in 1:n
        (sep[s0] || comp[s0] > 0) && continue
        nc += 1; comp[s0] = nc; push!(q, s0)
        while !isempty(q)
            u = pop!(q)
            for k in Ksym.colptr[u]:(Ksym.colptr[u+1]-1)
                v = Ksym.rowval[k]; (sep[v] || comp[v] > 0) && continue
                comp[v] = nc; push!(q, v)
            end
        end
    end
    return comp, nc
end
sep = falses(N); sep[sepidx] .= true
comp, nc = components(Ksym, sep)
sizes = zeros(Int, nc); for c in comp; c > 0 && (sizes[c] += 1); end
big = sort(findall(>(1000), sizes); by = c -> -sizes[c])
say(@sprintf("  hours: %d coupling columns (B, P_B, m) removed -> %d blocks, %d of them over 1000 columns (sizes %s)",
             count(sep), nc, length(big), join(sizes[big[1:min(5, end)]], ", ") * (length(big) > 5 ? ", ..." : "")))
members = [Int[] for _ in 1:nc]; for j in 1:N; comp[j] > 0 && push!(members[comp[j]], j); end
ph = Int[]
for c in sortperm(sizes)
    mc = members[c]; length(mc) == 1 ? push!(ph, mc[1]) : append!(ph, mc[amdp(Ksym[mc, mc])])
end
sv = findall(sep); append!(ph, sv[amdp(Ksym[sv, sv])])
H = analyse(Ksym, ph)
tail = N - length(sv) + 1
sepwork = sum(Float64.(H.Lnz[tail:end]) .^ 2); bigw = maximum((sum(Float64.(H.Lnz[findall(==(c), comp[ph])]) .^ 2) for c in big); init = 0.0)
say(@sprintf("  hours: L %d nonzeros, work %.3e, critical path %.3e, parallel bound %.1fx; coupling block %d columns, %.1f%% of the work; largest hour block %.1f%% of the work",
             H.nnzL, H.work, H.crit, H.work / H.crit, length(sv), 100sepwork / H.work, 100bigw / H.work))
say(@sprintf("  done in %.0f s, max RSS %.1f GB", time() - t0, Sys.maxrss() / 2^30))
