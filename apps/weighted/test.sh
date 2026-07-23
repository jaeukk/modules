#! /bin/bash

# Example: weighted structure factor S(k) of a 3D MRJ sphere packing.
# Weights = Voronoi-cell volumes, read from column 1 of the per-configuration
# weight files. 50 configurations, 5 threads. Writes <prefix>_chif.txt.
#
# Input data (MarkedPointPattern = ConfigurationPack + weights txt):
#   <centers>.ConfigPack / .ConfigPackIndex   — the point configurations
#   <weights> = "<printf-prefix> <column>"    — %04d is the configuration index

centers="/home/jaeukk/Projects/markedPtsProcess/VoronoiCellVolume/MRJ/3D_MRJ_5000_unitRho"
weights="/home/jaeukk/Projects/markedPtsProcess/VoronoiCellVolume/MRJ/3D_MRJ_5000_unitRho_weights__%04d.txt 1"
config="GetConfigsFunction ReadConfigPack $centers $weights"

echo "$config"

./weighted.out <<< "NumThreads 5
NumConfig 50
$config
SkComputation 3 50 0.1
OutputPrefix $centers
Calculation
Exit
"
