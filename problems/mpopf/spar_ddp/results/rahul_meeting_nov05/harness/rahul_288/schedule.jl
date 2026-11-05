function replay_rounds(N, b)
    @assert N > 0 && ispow2(b)
    rounds = Vector{Tuple{Int,Int}}[]
    h = 1
    while h <= N
        pairs = [(t,i+h) for i in 1:2h:N+1-h for t in i:b:(h < b ? i : i+h-1)]
        isempty(pairs) || push!(rounds,pairs)
        h *= 2
    end
    h = b ÷ 2
    while h >= 1
        pairs = [(i,i+h) for i in 1+h:2h:N+1-h]
        isempty(pairs) || push!(rounds,pairs)
        h ÷= 2
    end
    ends = min.(collect(2:N+2),N+1)
    for pairs in rounds
        lefts = Set(first.(pairs))
        @assert length(lefts) == length(pairs)
        @assert isempty(intersect(lefts,Set(last.(pairs))))
        @assert all(ends[i] == j for (i,j) in pairs)
        for (i,j) in pairs
            ends[i] = ends[j]
        end
    end
    @assert all(==(N+1),ends)
    rounds
end
