#include "MultiComponentPairStat.h"
#include <StructureFactor.h>
#include <NumberVariance.h>
#include <VolumeFractionVariance.h>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <omp.h>

size_t Verbosity = 0;
static void check(bool ok) { if (!ok) throw std::runtime_error("PairStat adapter test failed"); }
template<class F> static void mustThrow(F f) {
    bool threw = false;
    try { f(); } catch (const std::exception &) { threw = true; }
    check(threw);
}
static void same(const std::vector<GeometryVector> &a, const std::vector<GeometryVector> &b) {
    check(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i)
        for (int d = 0; d < 4; ++d)
            check(std::isfinite(a[i].x[d]) && std::abs(a[i].x[d]-b[i].x[d]) < 1e-12);
}

int main() {
    omp_set_num_threads(2);
    MultiComponentPairStatOptions options;
    options.kMax = 3.0; options.kBinWidth = 0.2;
    options.rMax = 0.8; options.rStep = 0.2;
    options.windowSamples = 32; options.seed = 13;
    for (DimensionType dim : {2, 3}) for (bool skewed : {false, true}) {
        GeometryVector basis[3] = {GeometryVector(dim), GeometryVector(dim), GeometryVector(dim)};
        for (DimensionType d = 0; d < dim; ++d) basis[d].x[d] = 4.0 + d;
        if (skewed) basis[1].x[0] = 0.3;
        SpherePacking original(dim, basis, 1.0);
        for (size_t i = 0; i < 10; ++i) {
            GeometryVector s(dim);
            s.x[0] = (i % 5 + 0.5)/5.0;
            s.x[1] = (i / 5 + 0.5)/2.0;
            if (dim == 3) s.x[2] = (i + 0.5)/10.0;
            original.Insert(0.07 + 0.001*i, s);
        }
        MultiStealthyPotential potential(dim, 1.0);
        potential.SetPacking(original);
        potential.SetNumSpecies(6); // Five occupied species, one empty.
        mustThrow([&] { ComputeMultiComponentPairStatistics(potential, options); });
        for (size_t i = 0; i < 10; ++i) potential.SetSpeciesOfParticle(i, i % 5);
        const auto snapshot = potential.GetPackingSnapshot();
        for (size_t i = 0; i < 10; ++i) {
            check(snapshot.GetCharacteristics(i) == original.GetCharacteristics(i));
            check((snapshot.GetCartesianCoordinates(i)-original.GetCartesianCoordinates(i)).Modulus2() == 0);
        }
        const auto results = ComputeMultiComponentPairStatistics(potential, options);
        check(results.size() == 7 && results[0].particles == 10);
        check(results[6].particles == 0 && results[6].structureFactor.empty());
        check(results[6].numberVariance.empty() && results[6].volumeFractionVariance.empty());
        RandomGenerator rng(options.seed);
        std::vector<GeometryVector> centers;
        GetSamplingPoints(dim, options.windowSamples, centers, rng);
        for (int species = -1; species < 5; ++species) {
            SpherePacking expected(original);
            if (species >= 0) {
                expected.RemoveParticles();
                for (size_t i = 0; i < 10; ++i)
                    if (int(i % 5) == species)
                        expected.Insert(original.GetCharacteristics(i), original.GetRelativeCoordinates(i));
            }
            Configuration points(expected, "a");
            const auto &result = results[species + 1];
            check(result.species == species && result.particles == expected.NumParticle());
            check(result.packingFraction == expected.PackingFraction() && !result.overlaps);
            std::vector<GeometryVector> sk, nv, vfv;
            IsotropicStructureFactor([&](size_t) { return points; }, 1, options.kMax,
                options.kMax, sk, options.kBinWidth, 1.0);
            MC_NV_single(points, nv, options.rMax, options.rStep, centers);
            MC_vfv([&](size_t) { return expected; }, 1, vfv,
                options.rMax, options.rStep, centers, 1);
            check(!sk.empty() && !nv.empty() && !vfv.empty());
            same(sk, result.structureFactor); same(nv, result.numberVariance);
            same(vfv, result.volumeFractionVariance);
        }
        std::ostringstream report;
        WriteMultiComponentPairStatistics(report, results, options);
        check(report.str().find("# component -1 particles 10") != std::string::npos);
        check(report.str().find("# component 5 particles 0") != std::string::npos);
        report.setstate(std::ios::badbit);
        mustThrow([&] { WriteMultiComponentPairStatistics(report, results, options); });
        potential.species_constraints[0].idx.push_back(0);
        mustThrow([&] { ComputeMultiComponentPairStatistics(potential, options); });
    }
    options.windowSamples = 1;
    mustThrow([&] { options.Validate(); });
    options.windowSamples = 2; options.rStep = 0;
    mustThrow([&] { options.Validate(); });
    std::cout << "PairStat adapter matches direct library calls: total + five species, 2D/3D, rectangular/skewed.\n";
}
