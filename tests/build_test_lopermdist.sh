#!/usr/bin/env bash
# Build the lpermdist test harness against your REAL project modules.
#
# Edit PROJECT_SRC and MODULES below to point at your actual source tree.
# List modules in dependency order (a module must be compiled before
# anything that USEs it). The order shown is what lpermdist.f90 itself
# needs; add any further modules that YOUR versions of these pull in.
#
# Two builds are produced:
#   test_lpermdist_debug   -- with -fcheck=bounds -fbacktrace: turns the
#                              known standard_orient natoms/nats bug into
#                              an immediate, pinpointed runtime crash
#                              instead of silent out-of-bounds memory access.
#   test_lpermdist_release -- without bounds checking, closer to a normal
#                              optimized build; may run to completion
#                              despite the bug, but a corrupted jmax2
#                              selection could still cost you a wrong
#                              orientation on a larger/different test case.
#
# Requires: gfortran, LAPACK/BLAS (liblapack, libblas) for DSYEV.
set -euo pipefail

# --- EDIT THIS: point at your project's real dependency modules -----------
PROJECT_SRC="../QCI-library"

MODULES=(
   "$PROJECT_SRC/source/utils/qciprec.f90"
   "$PROJECT_SRC/source/qci_keys.f90"
   "$PROJECT_SRC/source/perm/perm_defs.f90"
   "$PROJECT_SRC/source/perm/minperm.f90"

)

LOCAL_SRC=(
   "$PROJECT_SRC/source/perm/lopermdist.f90"
   "$PROJECT_SRC/tests/Test_lopermdist.F90"
)


# ----------------------------------------------------------------------------

for f in "${MODULES[@]}"; do
   if [[ ! -f "$f" ]]; then
      echo "ERROR: module not found: $f" >&2
      echo "Edit PROJECT_SRC / MODULES at the top of build.sh to match your tree." >&2
      exit 1
   fi
done

for f in "${MODULES[@]}" "${LOCAL_SRC[@]}"; do
   if [[ ! -f "$f" ]]; then
      echo "ERROR: source file not found: $f" >&2
      echo "Edit PROJECT_SRC / MODULES at the top of build.sh to match your tree," >&2
      exit 1
   fi
done
 
SRC=( "${MODULES[@]}" "${LOCAL_SRC[@]}" )


echo "== Debug build (bounds-checked) =="
gfortran -O0 -g -fcheck=bounds -fbacktrace -Wall "${SRC[@]}" -llapack -lblas -o test_lpermdist_debug

echo "== Release build (no bounds check) =="
gfortran -O2 "${SRC[@]}" -llapack -lblas -o test_lpermdist_release

echo
echo "Run with:"
echo "  ./test_lpermdist_debug"
echo "  ./test_lpermdist_release"