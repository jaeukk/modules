#! /bin/bash

# Examples on a 2D RSA point pattern (5 configurations); weights in column 1.
# Each run below reproduces one set of the bundled output files
# (test1_*, test2_*, test3_*), exercising all three computations.
#
# Input data (MarkedPointPattern = ConfigurationPack + weights txt):
#   <centers>.ConfigPack / .ConfigPackIndex   — the point configurations
#   <weights> = "<printf-prefix> <column>"    — %04d is the configuration index

centers="/home/jaeukk/Projects/markedPtsProcess/configs/RSA/test2"
weights="/home/jaeukk/Projects/markedPtsProcess/configs/RSA/2D_RSA_0.10_10000_weights__%04d.txt 1"
config="GetConfigsFunction ReadConfigPack $centers $weights"

echo "$config"

# test1 — local variance of the weights; 5000 random window centers, seed 123
#         -> test1_varf.txt, test1_means.txt
./weighted.out <<< "NumThreads 2
NumConfig 5
$config
LocalVariance 10 0.05 5000 123
OutputPrefix test1
Calculation
Exit
"

# test2 — same, but use the particle positions as window centers (UseParticleCenters)
#         -> test2_varf_self.txt, test2_means_self.txt
./weighted.out <<< "NumThreads 2
NumConfig 5
$config
LocalVariance 10 0.05 5000 123
UseParticleCenters
OutputPrefix test2
Calculation
Exit
"

# test3 — weighted g2(r) up to r=10 (dr=0.05): two-sided, then single-sided
#         -> test3_s2.txt, test3_single-sided-s2.txt
./weighted.out <<< "NumThreads 2
NumConfig 5
$config
g2Computation 10.0 0.05
OutputPrefix test3
Calculation
Exit
"
./weighted.out <<< "NumThreads 2
NumConfig 5
$config
g2Computation 10.0 0.05
SingleSided
OutputPrefix test3
Calculation
Exit
"
