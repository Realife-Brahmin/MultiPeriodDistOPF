include("schedule.jl")
for (T,b,merges,rounds) in ((96,1,368,7),(96,2,256,8),(288,1,1392,9),(288,4,708,11))
    result=replay_rounds(T,b)
    @assert (sum(length,result),length(result)) == (merges,rounds)
    println("SCHEDULE_CHECK T=$T b=$b merges=$merges rounds=$rounds")
end
