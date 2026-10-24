# Frozen reference matching pw_combine's in-place arithmetic and round pattern.
# This is NOT a production solver hook. Do not replace it with allocating
# pint_combine when checking saved policies near ill-conditioned late iterates.
function production_values(elems,terminal)
    xs=Any[deepcopy(e) for e in elems];push!(xs,deepcopy(terminal))
    w=PairWork(length(first(elems).ma));h=1
    while h<length(xs)
        for i in 1:2h:length(xs)
            mid=i+h;mid>length(xs) && continue
            right=xs[mid]
            for t in i:mid-1;xs[t]=_pp_combine!(xs[t],right,w);end
        end
        h*=2
    end
    @assert all(x->x isa PintVal,xs)
    xs
end
