#!/usr/bin/env bash
# Run from any directory; build artifacts stay outside the repository.
set -euo pipefail
test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mods="$(cd "$test_dir/../../.." && pwd)"
gsl_prefix="${GSL_PREFIX:-$HOME/install}"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/multistealthy-export.XXXXXX")"
"${CXX:-g++}" -std=c++17 -O2 -fopenmp \
  -I"$test_dir/.." -I"$mods/core" -I"$mods/potential" -I"$gsl_prefix/include" \
  "$test_dir/test_export.cpp" "$test_dir/../MultiStealthyPotential.cpp" \
  "$mods/core/obj/GeometryVector.o" "$mods/core/obj/etc.o" \
  "$mods/core/obj/RandomGenerator.o" "$mods/potential/obj/Potential.o" \
  "$mods/core/obj/PeriodicCellList.o" "$mods/core/obj/SpherePacking.o" \
  "$mods/core/obj/ConfigurationIO.o" \
  -L"$gsl_prefix/lib" -Wl,-rpath,"$gsl_prefix/lib" -lgsl -lgslcblas \
  -o "$build_dir/test_export"
"$build_dir/test_export"
echo "Test binary: $build_dir/test_export"
