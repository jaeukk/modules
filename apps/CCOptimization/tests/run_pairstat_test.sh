#!/usr/bin/env bash
set -euo pipefail
test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mods="$(cd "$test_dir/../../.." && pwd)"
gsl_prefix="${GSL_PREFIX:-$HOME/install}"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/multistealthy-pairstat.XXXXXX")"
"${CXX:-g++}" -std=c++17 -O2 -fopenmp \
  -I"$test_dir/.." -I"$mods/core" -I"$mods/potential" -I"$mods/pairstat" -I"$gsl_prefix/include" \
  "$test_dir/test_pairstat.cpp" "$test_dir/../MultiStealthyPotential.cpp" \
  "$test_dir/../MultiComponentPairStat.cpp" \
  "$mods/core/obj/GeometryVector.o" "$mods/core/obj/etc.o" \
  "$mods/core/obj/RandomGenerator.o" "$mods/potential/obj/Potential.o" \
  "$mods/core/obj/PeriodicCellList.o" "$mods/core/obj/SpherePacking.o" \
  "$mods/core/obj/ConfigurationIO.o" "$mods/pairstat/obj/StructureFactor.o" \
  "$mods/pairstat/obj/NumberVariance.o" "$mods/pairstat/obj/VolumeFractionVariance.o" \
  -L"$gsl_prefix/lib" -Wl,-rpath,"$gsl_prefix/lib" -lgsl -lgslcblas \
  -o "$build_dir/test_pairstat"
"$build_dir/test_pairstat"
echo "Test binary: $build_dir/test_pairstat"
