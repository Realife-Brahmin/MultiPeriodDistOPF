# Loaded only by the isolated full-solve snapshot.
include(joinpath(@__DIR__, "study_right_kernel.jl"))
const STUDY_BLAS_TARGET = parse(Int, get(ENV, "STUDY_BLAS_THREADS", "1"))
const STUDY_USE_RIGHT = get(ENV, "STUDY_KERNEL", "current") == "right"
const STUDY_MERGE = STUDY_USE_RIGHT ? right_combine! : _pp_combine!

function study_thread_check()
    actual = BLAS.get_num_threads()
    @assert actual == STUDY_BLAS_TARGET "Worker BLAS setting drifted"
    @assert _BLAS_DEFAULT[] == STUDY_BLAS_TARGET "Worker BLAS default drifted"
    actual
end

function study_audit()
    n = study_thread_check()
    tasks = sort(readdir("/proc/self/task"))
    masks = [(tid=tid, allowed=strip(match(r"(?m)^Cpus_allowed_list:\s*(.+)$",
        read("/proc/self/task/$tid/status", String)).captures[1])) for tid in tasks]
    (pid=myid(), os_pid=getpid(), blas=n, blas_default=_BLAS_DEFAULT[],
     backend=string(BLAS.get_config()), affinity=masks, maxrss=Sys.maxrss())
end

function study_pin(cpus)
    run(pipeline(Cmd(["taskset", "-apc", cpus, string(getpid())]); stdout=devnull))
    nothing
end
