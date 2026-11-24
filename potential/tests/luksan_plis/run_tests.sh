#!/usr/bin/env bash
# Build and run the LuksanPLIS tests (see README.md). Artifacts go to a fresh mktemp directory.
# Environment overrides: NLOPT_PREFIX (default $HOME/install), CXX (default g++), KEEP_WORK=1.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POTENTIAL="$(cd "$SCRIPT_DIR/../.." && pwd)"
NLOPT_PREFIX="${NLOPT_PREFIX:-$HOME/install}"
CXX="${CXX:-g++}"
WORK="$(mktemp -d)"
[[ "${KEEP_WORK:-0}" == 1 ]] || trap 'rm -rf "$WORK"' EXIT

# -ffp-contract=off as in potential/makefile: strict mode promises NLopt's arithmetic bit for bit.
"$CXX" -std=c++17 -O3 -g -fopenmp -Wall -ffp-contract=off -I"$POTENTIAL" -I"$NLOPT_PREFIX/include" \
	-o "$WORK/test_luksan_plis" "$SCRIPT_DIR/test_luksan_plis.cpp" "$POTENTIAL/LuksanPLIS.cpp" \
	-L"$NLOPT_PREFIX/lib" -Wl,-rpath,"$NLOPT_PREFIX/lib" -lnlopt -lm
OMP_NUM_THREADS=4 "$WORK/test_luksan_plis"
OMP_NUM_THREADS=4 OMP_WAIT_POLICY=PASSIVE OMP_THREAD_LIMIT=2 "$WORK/test_luksan_plis" | tail -1
