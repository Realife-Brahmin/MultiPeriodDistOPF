#!/usr/bin/env bash
# Build HSL MA57 and HSL_MA97 as Windows DLLs callable from Julia, from the
# user's own licensed HSL source packages. No HSL code lives in this repository:
# pass the directory that holds the unpacked packages (ma57-3.11.3/,
# hsl_ma97-2.8.1/) and a build directory outside the repo.
#
# Toolchain: MinGW gfortran/gcc (tested with the one bundled with Strawberry
# Perl, GCC 8.3). -march=haswell rather than native: this GCC's assembler cannot
# encode AVX-512 registers in Windows unwind data. BLAS/LAPACK come from the
# LP64 OpenBLAS32 artifact that Julia's MUMPS already uses (Julia's own OpenBLAS
# is ILP64 and cannot be linked by 32-bit-integer Fortran). METIS is replaced by
# HSL's placeholder, so METIS orderings are unavailable.
#
#   bash ddp/examples/power_system/hsl/build_hsl_windows.sh <hsl_src_dir> <build_dir> <libopenblas.dll>

set -eu
SRC=$(cd "$1" && pwd); OUTDIR=$2; OPENBLAS=$3
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUTDIR/o97"
cd "$OUTDIR"
cp "$OPENBLAS" ./libopenblas.dll
FLAGS="-O3 -march=haswell"

S57="$SRC/ma57-3.11.3/src"
gfortran $FLAGS -shared -o libma57.dll "$S57/ddeps.f" "$S57/ma57d.f" "$S57/fakemetis.f" libopenblas.dll

S97="$SRC/hsl_ma97-2.8.1/src"
(cd o97 && gfortran $FLAGS -fopenmp -c "$S97/common.f" "$S97/common90.f90" "$S97/ddeps90.f90" \
    "$S97/hsl_ma97d.f90" "$S97/hsl_ma97d_ciface.f90" "$S97/fakemetis.f" &&
 gcc -O2 -c "$HERE/s97.c" -I "$SRC/hsl_ma97-2.8.1/include")
gfortran -shared -fopenmp -o libhsl_ma97.dll o97/*.o libopenblas.dll

objdump -p libma57.dll | grep -qE "ma57cd_" && objdump -p libhsl_ma97.dll | grep -qE "s97_solve" &&
  echo "built libma57.dll and libhsl_ma97.dll in $OUTDIR"
