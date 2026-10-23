# Posthoc capture validation and old-state equivalence. Never optimizer timing.
ENV["VALUE_REPLAY_NUMERICAL_ROOT"]=joinpath(@__DIR__,"source")
include("/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v3/replay.jl")
out,old=ARGS
function find_tag(root,tag)
    rows=split.(readlines(joinpath(root,"captures","coverage.csv")),',')
    hits=[r[2] for r in rows[2:end] if r[1]==tag];@assert length(hits)==1
    joinpath(root,"captures",only(hits))
end
cap=find_tag(out,"late");meta=TOML.parsefile(joinpath(cap,"metadata.toml"))
@assert meta["campaign_capture"]&&meta["tags"]==["late"]&&meta["barrier_iteration"]>0
leaf=joinpath(cap,"leaves.jls");@assert bytes2hex(sha256(read(leaf)))==meta["leaves_sha256"]
p=deserialize(leaf);N=length(p.elems);n=length(first(p.elems).ma)
@assert N==meta["stages"]&&p.k==meta["iteration"]&&p.j==meta["barrier_iteration"]&&p.mu==meta["mu"]
expected=collect(1:min(4,N));packets=Any[]
for row in meta["stage_files"]
    f=joinpath(cap,row["file"]);@assert filesize(f)==row["bytes"]&&bytes2hex(sha256(read(f)))==row["sha256"]
    push!(packets,deserialize(f))
end
@assert sort([x.stage for x in packets])==expected
full=production_values(p.elems,PintVal(zeros(n,n),zeros(n)))
ord=ordinary(p.elems,PintVal(zeros(n,n),zeros(n)))
oe=maximum(max(scaled_error(full[t].V,ord[t].V),scaled_error(full[t].v,ord[t].v)) for t in 1:N);@assert oe<1e-7
probe=0.0
for packet in packets
    rule=policy(packet,full[packet.stage+1],p.mu)
    for pr in packet.probes,field in (:du,:dphi)
        probe=max(probe,scaled_error(getproperty(directions(rule,pr.dx),field),getproperty(pr,field)))
    end
end
@assert probe<1e-7 "unchanged dense saved-direction gate"
oldcap=find_tag(old,"late");om=TOML.parsefile(joinpath(oldcap,"metadata.toml"));op=deserialize(joinpath(oldcap,"leaves.jls"))
leaferror=0.0;exact=true
for t in 1:N,field in (:Maa,:Mab,:Mbb,:ma,:mb)
    A=getproperty(p.elems[t],field);B=getproperty(op.elems[t],field)
    leaferror=max(leaferror,scaled_error(A,B));exact&=A==B
end
nominalerror=0.0;nominalexact=true
for row in om["stage_files"]
    t=parse(Int,match(r"stage_(\d+)",row["file"]).captures[1]);t in expected||continue
    previous=deserialize(joinpath(oldcap,row["file"]));now=only(filter(x->x.stage==t,packets))
    for field in (:x,:u,:phi,:zl,:zu)
        A=getproperty(now.nominal,field);B=getproperty(previous.nominal,field)
        nominalerror=max(nominalerror,scaled_error(A,B));nominalexact&=A==B
    end
end
scalars=p.k==op.k&&p.j==op.j&&p.mu==op.mu
equiv=scalars&&leaferror<1e-7&&nominalerror<1e-7
classification=exact&&nominalexact&&scalars ? "bit_exact_original_Newton_state" : equiv ? "equivalent_within_reference_gates_new_capture_used" : "new_separately_identified_Newton_state"
report=Dict("validated"=>true,"classification"=>classification,"capture"=>cap,"old_capture"=>oldcap,
 "N"=>N,"nB"=>n,"selected_stages"=>expected,"k"=>p.k,"j"=>p.j,"mu"=>p.mu,
 "leaf_bit_exact"=>exact,"nominal_bit_exact"=>nominalexact,"scalar_state_matches"=>scalars,
 "old_leaf_max_scaled"=>leaferror,"old_nominal_max_scaled"=>nominalerror,
 "ordinary_error"=>oe,"saved_probe_error"=>probe,"leaves_sha256"=>meta["leaves_sha256"])
open(joinpath(out,"capture_validation.toml"),"w") do io;TOML.print(io,report);end
println("CAMPAIGN_CAPTURE_PASS classification=$classification stages=$expected leaf_error=$leaferror nominal_error=$nominalerror probe=$probe ordinary=$oe")
