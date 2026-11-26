# Keep the machine from idle-sleeping while a long unattended run is in
# progress, so no single timed solve straddles a sleep and reports inflated wall
# time.
#
# This is a per-process power request (SetThreadExecutionState, the call video
# players make). It changes no power setting and is released automatically when
# this process exits. The process exits when the lock file disappears -- the
# driver removes it on exit -- or after max_hours as a backstop, in case the
# driver was killed too hard to run its exit trap.
#
#   julia --startup-file=no hold_awake.jl <lockfile> [max_hours=16]

lockfile  = ARGS[1]
max_hours = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 16.0

const ES_CONTINUOUS      = 0x80000000
const ES_SYSTEM_REQUIRED = 0x00000001

if Sys.iswindows()
    prev = ccall((:SetThreadExecutionState, "kernel32"), stdcall, UInt32, (UInt32,),
                 ES_CONTINUOUS | ES_SYSTEM_REQUIRED)
    println("hold_awake: request ", prev == 0 ? "FAILED" : "active", " (max ", max_hours, " h)")
else
    println("hold_awake: not Windows, nothing to do")
end
flush(stdout)

deadline = time() + 3600 * max_hours
while isfile(lockfile) && time() < deadline
    sleep(60)
end
println("hold_awake: released (", isfile(lockfile) ? "deadline" : "lock removed", ")")
