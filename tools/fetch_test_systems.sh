#!/usr/bin/env bash
# Download the public distribution test systems into test_systems/, each in its own folder,
# exactly as published (no edits). Safe to re-run: a folder that exists is left alone.
#   bash tools/fetch_test_systems.sh
# Needs curl and the GitHub CLI (gh) for listing repository trees; no SVN or AWS CLI.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$ROOT/test_systems"
mkdir -p "$DEST"

# Every file under <prefix> of a GitHub repo at <ref>, into <dir>, keeping the relative paths.
github_folder() {   # repo ref prefix dir
    local repo=$1 ref=$2 prefix=$3 dir=$4
    [ -d "$dir" ] && { echo "skip $dir (exists)"; return; }
    mkdir -p "$dir"
    gh api "repos/$repo/git/trees/$ref?recursive=1" --jq ".tree[] | select(.type==\"blob\") | .path | select(startswith(\"$prefix/\"))" |
    while read -r p; do
        local rel=${p#"$prefix"/}
        mkdir -p "$dir/$(dirname "$rel")"
        curl -sfL "https://raw.githubusercontent.com/$repo/$ref/$(printf '%s' "$p" | sed 's/ /%20/g')" -o "$dir/$rel" || echo "FAILED $p"
    done
    echo "got $dir ($(find "$dir" -type f | wc -l) files)"
}

# 1. EPRI test circuits (Ckt5, Ckt7, Ckt24; J1, K1, M1 under epri_dpv). OpenDSS's SVN trunk,
#    Version8/Distrib/EPRITestCircuits, through its GitHub mirror dss-extensions/electricdss-tst.
github_folder dss-extensions/electricdss-tst master Version8/Distrib/EPRITestCircuits "$DEST/epri_test_circuits"

# 2. IEEE 8500-node (OpenDSS) from GRIDAPPSD/Powergrid-Models, branch cimhub_refactor (the
#    default branch keeps only its CIM XML), and IEEE8500_3subs, published only as CIM XML, from
#    GRIDAPPSD/CIMHub model_output_tests.
github_folder GRIDAPPSD/Powergrid-Models cimhub_refactor models/feeders/OpenDSS/IEEE/IEEE8500 "$DEST/ieee8500"
if [ ! -d "$DEST/ieee8500_3subs" ]; then
    mkdir -p "$DEST/ieee8500_3subs"
    curl -sfL https://raw.githubusercontent.com/GRIDAPPSD/CIMHub/master/model_output_tests/IEEE8500_3subs.xml -o "$DEST/ieee8500_3subs/IEEE8500_3subs.xml"
    echo "got $DEST/ieee8500_3subs"
fi

# 3. NREL SMART-DS v1.0, San Francisco sub-region P1U, base_timeseries OpenDSS, from the public
#    OEDI bucket over plain HTTPS (the listing has fewer than 1000 keys, so one page).
S3=https://oedi-data-lake.s3.amazonaws.com
PFX=SMART-DS/v1.0/2018/SFO/P1U/scenarios/base_timeseries/opendss
if [ ! -d "$DEST/smartds_sfo_p1u" ]; then
    mkdir -p "$DEST/smartds_sfo_p1u"
    curl -sf "$S3/?list-type=2&prefix=$PFX/&max-keys=1000" | grep -o "<Key>[^<]*</Key>" | sed 's/<Key>//; s/<\/Key>//' |
    while read -r k; do
        rel=${k#"$PFX"/}; [ -z "$rel" ] && continue
        mkdir -p "$DEST/smartds_sfo_p1u/$(dirname "$rel")"
        curl -sfL "$S3/$k" -o "$DEST/smartds_sfo_p1u/$rel" || echo "FAILED $k"
    done
    echo "got $DEST/smartds_sfo_p1u ($(find "$DEST/smartds_sfo_p1u" -type f | wc -l) files)"
fi
