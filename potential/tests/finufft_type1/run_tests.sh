#!/usr/bin/env bash
# Build and run the FINUFFT Type 1 GetRho tests (see README.md).
# All artifacts go to a fresh mktemp directory; nothing is written to the source tree.
#
# Variants:
#   direct   : no FINUFFT macro, include path, or library (proves the default build needs none)
#   finufft  : -DCCPOTENTIAL_USE_FINUFFT
#   fixed2   : -DCCPOTENTIAL_USE_FINUFFT -DHAVETWOFIXEDPARTICLES
#
# Environment overrides:
#   FINUFFT_ROOT  FINUFFT tree with source/finufft-2.5.1/include, install/include, install/lib
#                 (default: <task>/control/native-finufft, resolved relative to this script)
#   GSL_PREFIX    GSL prefix (default: $HOME/install)
#   CXX           compiler (default: g++)
#   KEEP_WORK=1   keep the temporary directory for inspection
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODS="$(cd "$SCRIPT_DIR/../../.." && pwd)"
SRC="$SCRIPT_DIR/test_finufft_type1.cpp"
FINUFFT_ROOT="${FINUFFT_ROOT:-$SCRIPT_DIR/../../../../../control/native-finufft}"
GSL_PREFIX="${GSL_PREFIX:-$HOME/install}"
CXX="${CXX:-g++}"

FINUFFT_CFLAGS=(-DCCPOTENTIAL_USE_FINUFFT
	-I"$FINUFFT_ROOT/source/finufft-2.5.1/include" -I"$FINUFFT_ROOT/install/include")
FINUFFT_LFLAGS=("$FINUFFT_ROOT/install/lib/libfinufft.so" -Wl,-rpath,"$FINUFFT_ROOT/install/lib")

# No -DNDEBUG: GeometryVector's layout depends on NDEBUG and must match the prebuilt objects
# (config.mk builds them without it).
CXXFLAGS=(-std=c++17 -O2 -fopenmp -Wall -I"$MODS/core" -I"$MODS/potential" -I"$GSL_PREFIX/include")
LIBS=(-L"$GSL_PREFIX/lib" -Wl,-rpath,"$GSL_PREFIX/lib" -lgsl -lgslcblas -lm)

# Potential.o and the core objects do not include CollectiveCoordinatePotential.h, so the
# FINUFFT macro cannot create ODR mismatches with them.
MIN_OBJS=("$MODS/core/obj/GeometryVector.o" "$MODS/core/obj/etc.o" "$MODS/core/obj/RandomGenerator.o"
	"$MODS/potential/obj/Potential.o")
FULL_OBJS=("${MIN_OBJS[@]}" "$MODS/core/obj/PeriodicCellList.o" "$MODS/core/obj/SpherePacking.o"
	"$MODS/core/obj/ConfigurationIO.o")

for f in "$SRC" "${FULL_OBJS[@]}"; do
	if [[ ! -f "$f" ]]; then
		echo "missing: $f" >&2
		echo "build objects with: make -C \"$MODS/core\" objs && make -C \"$MODS/potential\" objs" >&2
		exit 2
	fi
done
for f in "$FINUFFT_ROOT/source/finufft-2.5.1/include/finufft.h" "$FINUFFT_ROOT/install/lib/libfinufft.so"; do
	if [[ ! -e "$f" ]]; then
		echo "missing FINUFFT file: $f (set FINUFFT_ROOT)" >&2
		exit 2
	fi
done

WORK="$(mktemp -d "${TMPDIR:-/tmp}/finufft_type1.XXXXXX")"
ARTIFACTS=()
cleanup() {
	if [[ "${KEEP_WORK:-0}" == 1 ]]; then
		echo "kept work directory: $WORK"
		return
	fi
	local f
	for f in "${ARTIFACTS[@]+"${ARTIFACTS[@]}"}"; do
		rm -f -- "$f"
	done
	rmdir -- "$WORK" 2>/dev/null || echo "note: $WORK is not empty; left in place"
}
trap cleanup EXIT

# build NAME [compile flags...] -- [link flags...]
build() {
	local name=$1
	shift
	local cflags=() lflags=()
	while [[ $# -gt 0 && $1 != -- ]]; do cflags+=("$1"); shift; done
	[[ $# -gt 0 ]] && shift
	lflags=("$@")
	local obj="$WORK/$name.o" exe="$WORK/$name" log="$WORK/$name.link.log"
	ARTIFACTS+=("$obj" "$exe" "$log")
	echo "== compile $name"
	"$CXX" "${CXXFLAGS[@]}" "${cflags[@]+"${cflags[@]}"}" -c "$SRC" -o "$obj"
	if "$CXX" -fopenmp "$obj" "${MIN_OBJS[@]}" "${lflags[@]+"${lflags[@]}"}" "${LIBS[@]}" -o "$exe" 2>"$log"; then
		echo "   linked with minimal object set"
	else
		echo "   minimal object set did not link; retrying with the baseline set. Linker output:"
		cat "$log"
		"$CXX" -fopenmp "$obj" "${FULL_OBJS[@]}" "${lflags[@]+"${lflags[@]}"}" "${LIBS[@]}" -o "$exe"
		echo "   linked with baseline object set"
	fi
}

build direct --
build finufft "${FINUFFT_CFLAGS[@]}" -- "${FINUFFT_LFLAGS[@]}"
build fixed2 "${FINUFFT_CFLAGS[@]}" -DHAVETWOFIXEDPARTICLES -- "${FINUFFT_LFLAGS[@]}"

status=0
if command -v ldd >/dev/null 2>&1; then
	if ldd "$WORK/direct" | grep -q libfinufft; then
		echo "FAIL: direct-only binary depends on libfinufft"; status=1
	else
		echo "ok: direct-only binary has no libfinufft dependency"
	fi
	if ldd "$WORK/finufft" | grep -q 'libfinufft.*not found'; then
		echo "FAIL: libfinufft not resolvable at run time"; status=1
	fi
fi

for name in direct finufft fixed2; do
	echo "== run $name (OMP_NUM_THREADS=2)"
	if OMP_NUM_THREADS=2 OMP_DYNAMIC=false "$WORK/$name"; then
		echo "== $name: exit 0"
	else
		echo "== $name: FAILED"; status=1
	fi
done

if [[ $status -eq 0 ]]; then echo "ALL VARIANTS PASSED"; else echo "SOME VARIANTS FAILED"; fi
exit $status
