#!/usr/bin/env bash
#==============================================================================
# run_example.sh — minimal end-to-end example for the CC-inversion app.
#
# Reconstructs small 1D weakly-perturbed-lattice targets from a subset of their
# collective coordinates and tabulates how many distinct solutions each
# (N, M = number of constraints) pair admits — i.e. a tiny instance of the
# uniqueness study in Kim, Zhang, Stillinger & Torquato, J. Stat. Mech. (2018)
# 113302.
#
# Usage:   ./run_example.sh
# (builds CCInversion.out via the makefile if it is missing, then runs it)
#==============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# 1. Build the executable if needed.
if [ ! -x ./CCInversion.out ]; then
    echo ">> building CCInversion.out ..."
    make
fi

# 2. Run directory.  CCInversionProblem1 writes to a CWD-relative path
#    results/Reconst/dE_<EnergyTolerance>/Nk_skip_<NKSkip>/  and does NOT create
#    it, so we make it here.  The dE_ token is printf("%.0e") of EnergyTolerance.
RUNDIR="./example_run"
ENERGY_TOL="1e-17"     # -> dE_1e-17
NKSKIP=0               # -> Nk_skip_0
OUTDIR="$RUNDIR/results/Reconst/dE_${ENERGY_TOL}/Nk_skip_${NKSKIP}"
mkdir -p "$OUTDIR"

# 3. Feed the CCInversion command + its parameter block on stdin.
#    Order (see CCInversionProblem1 in CCInversionProblem.cpp):
#      CCInversion <printAllSolutions 0|1>
#      NKSkip dim minN maxN dN numConfigs Nc systType1 systType2
#      methodNum1 methodNum2 EnergyTolerance delta1 delta2 d_delta M1 M2
#      numThreads doProxMetric
#    systType 0 = weakly perturbed lattice (delta = perturbation amplitude),
#             1 = Poisson, 2 = non-Bravais (1D only).
#    methodNum 0 = BFGS+MINOP, 1 = SteepestDescent, 2 = CG, 3 = LocalGradDescent.
echo ">> running CC inversion (1D perturbed lattice, N = 4) ..."
( cd "$RUNDIR" && "$OLDPWD/CCInversion.out" . ) <<'EOF'
CCInversion
0
0
1
4
4
1
5
10
0
0
0
0
1e-17
0.05
0.05
0.05
-1
1
1
0
Exit
EOF

# 4. Report.
echo
echo ">> done.  output files:"
find "$OUTDIR" -type f -printf '   %p\n'
echo
echo ">> example (first result table):"
sed -n '1,12p' "$(find "$OUTDIR" -type f | head -1)"
