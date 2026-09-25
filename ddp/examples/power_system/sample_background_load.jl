# Measure CPU used by processes OTHER than this pipeline's own solver runs, so
# a timing taken while the machine was shared is flagged instead of silently
# mixed in with clean ones. (Added 2026-09-19: the user runs their own MSOPF
# analysis on the same machine, which inflates whichever arm it overlaps.)
#
# "Ours" = a Julia run whose command line names one of the pipeline scripts
# below. Under juliaup the arguments live on the julialauncher.exe PARENT, not on
# julia.exe itself, so ownership is read from the launcher and inherited by its
# julia.exe child (and by any julia.exe that child spawns). The sampler's own
# PowerShell query counts as ours. Everything else -- including the user's own
# Julia jobs -- is "other". Only julia/julialauncher command lines are queried,
# which keeps the WMI cost of sampling itself small.
#
#   sample mode:  julia sample_background_load.jl sample <stopfile> <out.csv> [interval_s=15]
#                 appends "time,other_cores,our_cores" rows until <stopfile> is removed
#   wait mode:    julia sample_background_load.jl wait <max_minutes> <threshold_cores>
#                 blocks until one interval shows other_cores < threshold, or timeout;
#                 prints "QUIET_WAIT waited_s=... other_cores=... timed_out=..."

using Dates, Printf

const OURS = "ieee123c_filterddp.jl|centralized_ipopt_matched.jl|hold_awake.jl|" *
             "sample_background_load.jl|export_ieee123c_data.jl|summarize_|" *
             "kkt_ordering_benchmark.jl|blocked_multirhs_solve_benchmark.jl"

const PS = """
\$w = @{}
\$w[[int]\$PID] = 1
Get-CimInstance Win32_Process -Filter "Name='julialauncher.exe'" | ForEach-Object { \$w[[int]\$_.ProcessId] = [int]([bool](\$_.CommandLine -match '$OURS')) }
\$jl = @(Get-CimInstance Win32_Process -Filter "Name='julia.exe'")
foreach (\$pass in 1..3) { foreach (\$p in \$jl) { if (\$w[[int]\$p.ParentProcessId] -eq 1) { \$w[[int]\$p.ProcessId] = 1 } } }
Get-Process | Where-Object { \$_.CPU -ne \$null -and \$_.Id -ne 0 } | ForEach-Object { '{0},{1},{2}' -f \$_.Id, \$_.CPU, [int]\$w[[int]\$_.Id] }
"""

function snapshot()
    out = read(`powershell -NoProfile -NonInteractive -Command $PS`, String)
    snap = Dict{Int,Tuple{Float64,Bool}}()
    for line in split(out, '\n')
        f = split(strip(line), ',')
        length(f) == 3 || continue
        pid = tryparse(Int, f[1]); cpu = tryparse(Float64, f[2])
        (pid === nothing || cpu === nothing) && continue
        snap[pid] = (cpu, strip(f[3]) == "1")
    end
    return snap, time()
end

"""Cores busy (other, ours) between two snapshots. A process that appears only
in the second snapshot started inside the interval, so all its CPU counts."""
function cores(a, ta, b, tb)
    other = 0.0; ours = 0.0
    for (pid, (cpu, mine)) in b
        d = haskey(a, pid) ? cpu - a[pid][1] : cpu
        d < 0 && continue
        mine ? (ours += d) : (other += d)
    end
    dt = max(tb - ta, 1e-3)
    return other / dt, ours / dt
end

function sample_mode(stopfile, outfile, interval)
    isfile(outfile) || write(outfile, "time,other_cores,our_cores\n")
    a, ta = snapshot()
    while isfile(stopfile)
        sleep(interval)
        b, tb = snapshot()
        o, m = cores(a, ta, b, tb)
        open(outfile, "a") do io
            @printf(io, "%s,%.3f,%.3f\n", Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"), o, m)
        end
        a, ta = b, tb
    end
end

function wait_mode(maxmin, thresh)
    t0 = time(); a, ta = snapshot(); o = NaN
    while true
        sleep(15)
        b, tb = snapshot()
        o, _ = cores(a, ta, b, tb)
        a, ta = b, tb
        (o < thresh || time() - t0 > 60 * maxmin) && break
    end
    @printf("QUIET_WAIT waited_s=%.0f other_cores=%.2f timed_out=%d\n",
            time() - t0, o, o < thresh ? 0 : 1)
end

if ARGS[1] == "sample"
    sample_mode(ARGS[2], ARGS[3], length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 15.0)
elseif ARGS[1] == "wait"
    wait_mode(parse(Float64, ARGS[2]), parse(Float64, ARGS[3]))
else
    error("mode must be sample or wait")
end
