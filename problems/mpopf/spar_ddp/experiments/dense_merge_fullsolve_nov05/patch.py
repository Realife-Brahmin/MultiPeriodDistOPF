import sys
from pathlib import Path

root = Path(sys.argv[1])
ps = root / "spar_ddp/examples/power_system"
p = ps / "procs_pass.jl"
s = p.read_text()
needle = "    BLAS.set_num_threads(1); _BLAS_DEFAULT[] = 1"
assert s.count(needle) == 1
s = s.replace(needle, "    BLAS.set_num_threads(STUDY_BLAS_TARGET); _BLAS_DEFAULT[] = STUDY_BLAS_TARGET")
s += '\ninclude(joinpath(@__DIR__, "study_control.jl"))\n'
a = s.index("    n = N + 1; h = 1\n")
b = s.index("        res = Vector{Any}", a)
s = s[:a] + '''    for pairs in study_scan_rounds(N, get(ENV,"STUDY_SCAN","balanced")=="grouped")
''' + s[b:]
needle = "        h *= 2\n    end\n    t_2 = time_ns()"
assert s.count(needle) == 1
s = s.replace(needle, "    end\n    t_2 = time_ns()")
s = '''include(joinpath(@__DIR__, "study_schedule.jl"))
const STUDY_ROUNDS = Dict{Tuple{Int,Bool},Vector{Vector{Tuple{Int,Int}}}}()
study_scan_rounds(N,grouped) = get!(() -> replay_rounds(N,grouped ? 2 : 1), STUDY_ROUNDS,(N,grouped))
''' + s
start = s.index("function pw_combine(")
end = s.index("\nend", start) + len("\nend")
part = s[start:end]
assert part.count("_pp_combine!(") == 3
part = part.replace("_pp_combine!(", "STUDY_MERGE(")
part = part.replace("    t0 = time_ns(); g0 = _pw_gc()", "    study_thread_check()\n    t0 = time_ns(); g0 = _pw_gc()")
s = s[:start] + part + s[end:]
for name in ("pw_pass_start", "pw_step1", "pw_step3"):
    start = s.index(f"function {name}(")
    end = s.index("\nend", start) + len("\nend")
    part = s[start:end]
    first = part.index("\n")
    part = part[:first+1] + "    study_thread_check()\n" + part[first+1:]
    marker = "    return ("
    where = part.index(marker)
    part = part[:where] + "    study_thread_check()\n" + part[where:]
    if name in ("pw_step1", "pw_step3"):
        part = part.replace(marker, marker + "study_threads=BLAS.get_num_threads(), ", 1)
    s = s[:start] + part + s[end:]
marker = "    t_2 = time_ns()"
assert s.count(marker) == 1
s = s.replace(marker, marker + '''
    println("FULL_STUDY_PHASE_B_THREADS iteration=$(d.k) mu=$μ threads=$(STUDY_BLAS_TARGET) checked_every_merge=1")
''')
p.write_text(s)

p = ps / "procs_main.jl"
s = p.read_text()
i = s.index("let system = ARGS[1]")
s = s[:i] + '''maps=readlines(ENV["STUDY_CPU_MAP"])
@assert length(maps)==length(workers())+1
run(pipeline(Cmd(["taskset","-apc",maps[1],string(getpid())]);stdout=devnull))
BLAS.set_num_threads(1)
@sync for (j,p) in enumerate(workers())
    @async remotecall_fetch(study_pin,p,maps[j+1])
end
''' + s[i:]
needle = "main()\nrmprocs(workers())"
assert s.count(needle) == 1
s = s.replace(needle, '''function audit_workers(tag)
    info=[remotecall_fetch(study_audit,p) for p in workers()]
    open(joinpath(ENV["STUDY_OUT"],"worker-audit-$tag.txt"),"w") do io
        for row in info;println(io,row);end
    end
    println("FULL_STUDY_AUDIT tag=$tag workers=$(length(info)) threads_min=$(minimum(x.blas for x in info)) threads_max=$(maximum(x.blas for x in info)) coordinator_threads=$(BLAS.get_num_threads())")
end
audit_workers("after-setup")
main()
audit_workers("after-solve")
rmprocs(workers())''')
p.write_text(s)

p = ps / "ieee123c_filterddp.jl"
s = p.read_text()
needle = 'warmup_spec = get(ENV, "FILTERDDP_WARMUP", "0")'
assert s.count(needle) == 1
s = "using SHA\n" + s
s = s.replace(needle, '''ctx=SHA.SHA2_256_CTX()
for u in ubar; SHA.update!(ctx,reinterpret(UInt8,u)); end
println("FULL_STUDY_START x0_sha256=$(bytes2hex(sha256(reinterpret(UInt8,x0)))) ubar_sha256=$(bytes2hex(SHA.digest!(ctx)))")
warmup_spec = get(ENV, "FILTERDDP_WARMUP", "0")''')
p.write_text(s)

p = root / "spar_ddp/experiments/battery_count/run.sh"
s = p.read_text()
needle = "FILTERDDP_PROCS_WORKERS=$((T<63?T:63))"
assert s.count(needle) == 1
s = s.replace(needle, "FILTERDDP_PROCS_WORKERS=$((T<STUDY_WORKERS?T:STUDY_WORKERS))")
p.write_text(s)
