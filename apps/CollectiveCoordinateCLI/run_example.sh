#!/usr/bin/env bash
# run_example.sh -- minimal end-to-end example for ./CollectiveCoordinateCLI.out
#
# Generates a small 2D stealthy-hyperuniform ground state (N=40, chi=0.3), relaxes
# it to (near) zero energy, writes the relaxed configuration, then does a 3-sample
# MultiRun that produces a ConfigurationPack.  All outputs land next to this script
# with the clearly-named prefix "stealthy_2D_N40".
#
# The interpreter reads whitespace-separated command tokens from stdin.  The token
# order used below (see README.md for the full command table):
#   Dimension <d> / NumParticles <N> / BoxSize <L> / Chi <chi>   set parameters
#   NumThreads <n> / RandomSeed <s>                              set run controls
#   GenerateConfiguration                                        random init in a cubic box
#   GeneratePotential                                            Type-0 (stealthy) CC potential
#   Run                                                          relax via NLopt/L-BFGS (prints E)
#   Output <prefix>                                              write <prefix>.pos (+ stub plot msg)
#   Prefix <p> / SampleNumber <k> / MultiRun / <Algorithm>       k-sample ground-state search
#   Exit
set -euo pipefail
cd "$(dirname "$0")"

make >/dev/null

echo "== 2D stealthy ground state: N=40, chi=0.3, seed=12345 =="
./CollectiveCoordinateCLI.out <<'EOF'
Dimension 2
NumParticles 40
BoxSize 10
Chi 0.3
NumThreads 4
RandomSeed 12345
GenerateConfiguration
GeneratePotential
Run
Output stealthy_2D_N40
Prefix stealthy_2D_N40
SampleNumber 3
MultiRun
LBFGS
Exit
EOF

echo
echo "Outputs:"
ls -1 stealthy_2D_N40* 2>/dev/null || echo "  (none)"
