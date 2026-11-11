#!/usr/bin/env bash
set -euo pipefail
root=/home/aryan.jha/kamiak_jobs
reference="$root/tadmm_review_nov09_med2522_single_node"
replay_source="$root/tadmm_review_nov09_med2522_caps_sequence"
tools="$root/tadmm_review_nov11_large10k_interim/campaign"
C="$root/tadmm_review_nov11_med2522_july_inputs"
P=problems/mpopf/tadmm/review_pilot
network=ieee2522C_1ph
inputs=problems/mpopf/spar_ddp/results/network_filterddp
data="$inputs/network_data_${network}_T96_periodic.jls"
tag=50b8901c324ba9eefd0d672c6e484e0ededc7a57
case $1 in
setup)
    mkdir "$C"
    mkdir "$C"/{campaign,logs,results}
    cp "$0" "$2" "$C/campaign/"
    cp -a "$reference/source" "$C/source"
    cp "$reference/FULL_PERIOD_READY" "$reference/ANCHOR_GATE_VALIDATED" "$C/"
    ln -s "$reference/julia-depot" "$C/julia-depot"
    cd "$C/source"
    sha256sum --check --quiet "$reference/execution-source.sha256"
    mkdir "$C/campaign/release"
    tar -xf "$C/campaign/med2522-july-release-nov11.tar" -C "$C/campaign/release"
    release="$C/campaign/release/rawData"
    mv "test_networks/$network" "$C/campaign/current-feeder"
    cp -a "$release/ieee2552C_1ph" "test_networks/$network"
    cp "$release"/LoadShape*Default.dss test_networks/
    diff -rq "$release/ieee2552C_1ph" "test_networks/$network"
    for f in "$release"/LoadShape*Default.dss; do cmp "$f" "test_networks/${f##*/}"; done
    printf '%s\n' "$tag" > "$C/campaign/release-commit.txt"
    mv "$data" "$C/campaign/current-T96.jls"
    # The frozen benchmark resolves the data through its existing periodic filename.
    # Store the actual release sampling flag in the benchmark metadata.
    python3 - "$P/run_horizon_sweep.jl" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1])
s=p.read_text()
old='"periodic"=>true'
assert s.count(old)==1
p.write_text(s.replace(old,'"periodic"=>get(d,:profile_periodic,true)'))
PY
    find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
    ;;
export)
    cd "$C/source"
    sha256sum --check --quiet "$C/execution-source.sha256"
    source "$root/tadmm_review_nov09_sequential_dss_v2/campaign/env.sh"
    export JULIA_DEPOT_PATH="$C/julia-depot:$JULIA_DEPOT_PATH"
    julia --startup-file=no --project="$P" - "$network" "$data" <<'JULIA'
using Serialization
include(abspath("test_networks/parse_opendss.jl"))
network, path = ARGS
T = 96
phase = range(0, 2pi; length=T)
load = 0.8 .+ 0.2 .* (sin.(phase .- 0.8) .+ 1) ./ 2
price = 0.08 .+ 0.12 .* (sin.(phase) .+ 1) ./ 2
pv = zeros(T)
pv[24:72] .= sin.(range(0, pi; length=49))
d = parse_system_from_dss(network, T; kV_B=4.1569, delta_t_h=0.25,
    LoadShapeLoad=load, LoadShapeCost=price, LoadShapePV=pv, C_B=1e-6*minimum(price))
d[:profile_periodic] = false
root = d[:substationBus]
# Same inactive substation-device policy as the current-input benchmark.
d[:Bset] = filter(!=(root), d[:Bset])
d[:Dset] = filter(!=(root), d[:Dset])
d[:drop_root_ders] = true
@assert length(d[:Nset]) == 2522 && length(d[:Lset]) == 2521
@assert length(d[:Bset]) == 249 && all(isfinite, load) && all(isfinite, price)
@assert isapprox(first(load),last(load)) && isapprox(first(price),last(price))
OpenDSSDirect.Circuit.SetActiveElement("Vsource.source")
@assert isapprox(parse(Float64,OpenDSSDirect.Properties.Value("basekv")),4.1569)
serialize(path,d)
println("July export: T=$T, kV=",d[:kV_B],", batteries=",length(d[:Bset]),", PV=",length(d[:Dset]))
JULIA
    sha256sum "$data" "test_networks/$network/"*.dss test_networks/LoadShape*Default.dss > "$C/input.sha256"
    date -u +%FT%TZ > "$C/JULY_INPUTS_VALIDATED"
    ;;
solve)
    test -s "$C/JULY_INPUTS_VALIDATED"
    bash "$reference/campaign/med2522-legacy-controller-nov11.sh" "$C"
    ;;
audit)
    j=$2
    run="$C/results/run_T96_primal_dual_$j"
    report="$C/results/chronological-$j-$SLURM_JOB_ID"
    test -s "$run/CASE_RECORDED"
    mkdir "$report"
    cd "$C/source"
    sha256sum --check --quiet "$C/input.sha256"
    (cd "$replay_source/source"; sha256sum --check --quiet "$replay_source/execution-source.sha256")
    cp "$C/input.sha256" "$C/execution-source.sha256" "$report/"
    sha256sum "$run/instance.jls" "$run/central_dispatch.jls" "$run/rho20000/endpoint.jls" > "$report/replay-inputs.sha256"
    source "$root/tadmm_review_nov09_sequential_dss_v2/campaign/env.sh"
    export PYTHONPATH="$root/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}"
    for name in central endpoint; do
        if [[ $name == central ]]; then dispatch="$run/central_dispatch.jls"; else dispatch="$run/rho20000/endpoint.jls"; fi
        julia --startup-file=no --project="$P" - "$replay_source/source/$P/sequential_dss.jl" "$run/instance.jls" "$dispatch" "$report/$name" "$C/source" <<'JULIA'
code, instance, dispatch, out, repo = ARGS
include(code)
SequentialDSS.replay(instance, dispatch, out; efficiency="model", repo=repo)
for file in readdir(out; join=true)
    endswith(file,".csv") && run(`gzip $file`)
end
JULIA
        python3 "$tools/audit-socp-storage-nov09.py" "$report/$name" "test_networks/$network/Storage.dss"
    done
    mkdir "$report/rounds"
    cp "$run/rho20000/trace.csv" "$run/rho20000/trace_subproblems.csv" "$report/rounds/"
    python3 "$tools/summarize-tadmm-round-timing-nov10.py" "$report/rounds"
    sha256sum --check --quiet "$report/replay-inputs.sha256"
    sacct -j "$j,$SLURM_JOB_ID" --format=JobID,State,Elapsed,ExitCode,MaxRSS,MaxVMSize,ReqMem,AllocCPUS,NodeList -P > "$report/slurm.txt"
    date -u +%FT%TZ > "$report/REPLAY_RECORDED"
    ;;
*) exit 2 ;;
esac
