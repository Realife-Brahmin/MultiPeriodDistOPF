using Profile
@everywhere using Profile
let output = ENV["PROFILE_RUN_DIR"]
    sampled = workers()[round.(Int, range(1, length(workers()); length=8))]
    @sync for pid in sampled
        @async remotecall_fetch(Core.eval, pid, Main, quote
            Profile.init(n=8_000_000, delay=0.02)
            Profile.clear()
            Profile.start_timer()
            nothing
        end)
    end
    Profile.init(n=8_000_000, delay=0.02)
    Profile.clear()
    try
        Profile.@profile main()
    finally
        @sync for pid in sampled
            @async remotecall_fetch(Core.eval, pid, Main, quote
                Profile.stop_timer()
                open(joinpath($output, "profile-worker-$(myid()).txt"), "w") do io
                    Profile.print(io; C=true, format=:tree)
                end
                nothing
            end)
        end
        open(joinpath(output, "profile-coordinator.txt"), "w") do io
            Profile.print(io; C=true, format=:tree)
        end
    end
end
