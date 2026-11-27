#!/usr/bin/env bash
set -euo pipefail
[[ ${TADMM_PIN_CPU:?} =~ ^[0-9]+$ ]]
exec taskset -c "$TADMM_PIN_CPU" "${TADMM_JULIA_EXE:?}" "$@"
