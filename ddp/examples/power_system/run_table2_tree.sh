#!/usr/bin/env bash
# The Table II/III FilterDDP cells in two further configurations, on top of
# the typed residuals (run_table2_typed.sh):
#   sd    FILTERDDP_STRUCTURED_DYNAMICS=1: sparse f_x and f_u products and a
#         structured terminal stage, still UMFPACK + the blocked solve;
#   tree2 the same plus FILTERDDP_TREE_KKT=1: the radial tree solver
#         (tree_kkt.jl) replaces the stage's sparse LU and multi-column solve.
# ieee123 on one thread, med2522 and large10k on eight, as in Table II. Cells
# whose log is complete are skipped. Starts only when no julia.exe is running.
#
#   bash ddp/examples/power_system/run_table2_tree.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
until [ "$(tasklist //FI "IMAGENAME eq julia.exe" 2>/dev/null | grep -c julia.exe)" = "0" ]; do sleep 30; done
echo "[$(date '+%H:%M:%S')] queue start"
S=ddp/examples/power_system/run_blocked_solve_fullrun.sh
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 CB=system
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1

run_cells() {  # $1 = tag, rest = "system T" cells
  local tag=$1; shift
  for cell in "$@"; do
    set -- $cell
    case $1 in ieee123C_1ph) export JULIA_NUM_THREADS=1 RUN_TAG_SUFFIX=_tableII_cbsys_typed_$tag;;
               *)            export JULIA_NUM_THREADS=8 RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_$tag;; esac
    bash "$S" "$1" "$2" diag 16 1
  done
}

export FILTERDDP_TREE_KKT=1
run_cells tree2 "large10kC_1ph 6" "large10kC_1ph 24" "large10kC_1ph 48"
unset FILTERDDP_TREE_KKT
run_cells sd "ieee123C_1ph 6" "ieee123C_1ph 24" "ieee123C_1ph 96" \
             "ieee2522C_1ph 6" "ieee2522C_1ph 24" "ieee2522C_1ph 96" "large10kC_1ph 6"
export FILTERDDP_TREE_KKT=1
run_cells tree2 "ieee123C_1ph 6" "ieee123C_1ph 24" "ieee123C_1ph 96" \
                "ieee2522C_1ph 6" "ieee2522C_1ph 24" "ieee2522C_1ph 96"
echo "[$(date '+%H:%M:%S')] TABLE2_TREE_DONE"
