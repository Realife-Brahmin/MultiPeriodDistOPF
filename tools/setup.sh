#!/usr/bin/env bash
# Run from any directory after cloning. Julia installation is an explicit prior step.
set -euo pipefail
setup_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
setup_julia="${MPOPF_JULIA:-julia}"
if ! command -v "$setup_julia" >/dev/null 2>&1; then
    echo "Julia not found. Install Julia/Juliaup: https://julialang.org/downloads/" >&2
    exit 1
fi
setup_args=(--startup-file=no)
if [[ -n "${MPOPF_JULIA_CHANNEL:-}" ]]; then
    setup_args=("+${MPOPF_JULIA_CHANNEL}" "${setup_args[@]}")
fi
exec "$setup_julia" "${setup_args[@]}" "$setup_root/tools/setup.jl" "$@"
