# Position of every stage variable in the FilterDDP control vector. Shared by
# the FilterDDP driver, the Ipopt driver (when it saves its solution in the
# same layout) and the OpenDSS check, so that none of them can drift.
function control_layout(data)
    N, L, B, D = length(data[:Nset]), length(data[:Lset]),
                  length(data[:Bset]), length(data[:Dset])
    k = 0
    take(n) = (r = (k+1):(k+n); k += n; r)
    idx = (ps=first(take(1)), qs=first(take(1)), P=take(L), Q=take(L),
           v=take(N), ell=take(L), pb=take(B), qnorm=take(D),
           soc_slack=take(L), energy_slack=take(B))
    return idx, k
end
