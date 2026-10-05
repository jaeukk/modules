#ifndef MULTICOMPONENT_PAIRSTAT_H
#define MULTICOMPONENT_PAIRSTAT_H

#include "MultiStealthyPotential.h"

struct MultiComponentPairStatOptions {
    double kMax = 2.0, kBinWidth = 0.05;
    double rMax = 2.0, rStep = 0.1;
    size_t windowSamples = 1000;
    int seed = 0;
    void Validate() const;
};

struct ComponentPairStatistics {
    int species = -1; // -1: combined packing; otherwise the original species label.
    size_t particles = 0;
    double packingFraction = 0.0;
    bool overlaps = false;
    // Native PairStat tables, unchanged. Empty species have empty tables.
    std::vector<GeometryVector> structureFactor, numberVariance, volumeFractionVariance;
};

// Reuse PairStat on in-memory snapshots; no configuration serialization or disk I/O.
// Combined statistics are computed directly, never by averaging species statistics.
// Shared sampling centers across all species and both variance estimators.
std::vector<ComponentPairStatistics> ComputeMultiComponentPairStatistics(
    const MultiStealthyPotential &potential, const MultiComponentPairStatOptions &options);

// Optional final report only; callers may instead consume the returned vectors.
void WriteMultiComponentPairStatistics(std::ostream &out,
    const std::vector<ComponentPairStatistics> &results,
    const MultiComponentPairStatOptions &options);

#endif
