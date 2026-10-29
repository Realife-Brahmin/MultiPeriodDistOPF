#!/usr/bin/env bash
# Local, serial sensitivity preparation. Each solve is independently replayed in OpenDSS.
# Run from the repository root; optional arguments select configured four-source stages.
set -euo pipefail
unset NO_COLOR FORCE_COLOR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
export MULTIPERIOD_LEAN=1 MULTI_SENS_ONLY=1 MULTI_SAMPLES=0
export MULTI_DIRSET=fib78 MULTI_VMARGIN=20000 LARGE_T=24
export MULTI_CASES=nobatt,cyclic MULTI_SENS_H=0.1
export MULTIPERIOD_MAX_WALL=600
root="${MSOPF_SENS_OUT:-problems/msopf/results/test_networks/sensitivity_oct25}"
mkdir -p "$root"
if [[ $# -eq 0 ]]; then
    set -- four_ieee9500xf four_smartds_p1u_typeiii
fi
for stage in "$@"; do
    out="$root/$stage"
    if [[ -e "$out/metadata.txt" ]]; then
        echo "Refusing to overwrite $out; select a fresh MSOPF_SENS_OUT" >&2
        exit 2
    fi
    mkdir -p "$out"
    {
        date -u '+start_utc=%Y-%m-%dT%H:%M:%SZ'
        echo "stage=$stage"
        git rev-parse HEAD
        git status --short
        hostname
        julia --version
        env | sort | grep -E '^(OPENBLAS|OMP_|JULIA_|MULTI|LARGE_|IPOPT_)'
        sha256sum problems/msopf/{Project.toml,Manifest.toml,full_angle_pf.jl,multiperiod_full_angle.jl} problems/msopf/root_level/{multiperiod_study.jl,run_test_edition_sensitivities.sh,check_test_edition_sensitivities.py}
        for deck in ieee9500xf_1ph smartds_p1u_xf_1ph; do
            sha256sum test_networks/"$deck"/*.dss
        done
    } > "$out/metadata.txt"
    git diff -- problems/msopf > "$out/source.patch"
    export MULTIPERIOD_OUT_DIR="$out"
    started=$SECONDS
    status=0
    julia --startup-file=no --project=problems/msopf problems/msopf/root_level/multiperiod_study.jl "$stage" > "$out/run.log" 2> "$out/run.err" || status=$?
    {
        echo "exit_code=$status"
        echo "wall_seconds=$((SECONDS-started))"
        date -u '+end_utc=%Y-%m-%dT%H:%M:%SZ'
    } >> "$out/metadata.txt"
    [[ $status -eq 0 ]] || exit "$status"
    if python problems/msopf/root_level/check_test_edition_sensitivities.py "$out"; then
        echo "validation=passed" >> "$out/metadata.txt"
    else
        echo "validation=failed" >> "$out/metadata.txt"
        exit 1
    fi
    echo "$stage finished and validated: $((SECONDS-started)) seconds; see $out"
done
