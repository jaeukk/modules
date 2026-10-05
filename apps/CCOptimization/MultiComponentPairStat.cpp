#include "MultiComponentPairStat.h"
#include <StructureFactor.h>
#include <NumberVariance.h>
#include <VolumeFractionVariance.h>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

void MultiComponentPairStatOptions::Validate() const {
    for (double value : {kMax, kBinWidth, rMax, rStep})
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument("PairStat ranges and increments must be finite and positive.");
    if (windowSamples < 2 || windowSamples > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("PairStat requires 2..INT_MAX sampling centers.");
    if (!std::isfinite(rMax/rStep) || std::ceil(rMax/rStep) > std::numeric_limits<int>::max() ||
        !std::isfinite(kMax/kBinWidth) || kMax/kBinWidth >= std::numeric_limits<int>::max())
        throw std::invalid_argument("PairStat bin count exceeds supported range.");
}

std::vector<ComponentPairStatistics> ComputeMultiComponentPairStatistics(
    const MultiStealthyPotential &potential, const MultiComponentPairStatOptions &options) {
    options.Validate();
    const SpherePacking total = potential.GetPackingSnapshot();
    const auto dim = total.GetDimension();
    if (dim != 2 && dim != 3)
        throw std::invalid_argument("This packing analysis adapter supports 2D and 3D.");
    if (!total.NumParticle()) throw std::invalid_argument("Cannot analyze an empty total packing.");
    // Validate a disjoint, complete partition before calling any statistics routines.
    std::vector<bool> seen(total.NumParticle(), false);
    for (const auto &component : potential.species_constraints)
        for (size_t id : component.idx) {
            if (id >= seen.size() || seen[id])
                throw std::invalid_argument("Invalid or duplicate component particle index.");
            seen[id] = true;
        }
    for (bool assigned : seen)
        if (!assigned) throw std::invalid_argument("Component partition omits particles.");

    RandomGenerator rng(options.seed);
    std::vector<GeometryVector> centers;
    GetSamplingPoints(dim, options.windowSamples, centers, rng);
    std::vector<ComponentPairStatistics> results;
    auto analyze = [&](const SpherePacking &packing, int species) {
        ComponentPairStatistics result;
        result.species = species;
        result.particles = packing.NumParticle();
        result.packingFraction = packing.PackingFraction();
        if (result.particles) {
            result.overlaps = packing.CheckOverlap();
            const Configuration points(packing, "a");
            IsotropicStructureFactor([&](size_t) { return points; }, 1,
                options.kMax, options.kMax, result.structureFactor, options.kBinWidth, 1.0);
            MC_NV_single(points, result.numberVariance, options.rMax, options.rStep, centers);
            // One configuration: no ensemble parallelism needed. PairStat owns all estimators.
            MC_vfv([&](size_t) { return packing; }, 1, result.volumeFractionVariance,
                options.rMax, options.rStep, centers, 1);
        }
        results.push_back(std::move(result));
    };
    analyze(total, -1);
    for (size_t s = 0; s < potential.species_constraints.size(); ++s) {
        SpherePacking component(total);
        component.RemoveParticles();
        for (size_t id : potential.species_constraints[s].idx)
            component.Insert(total.GetCharacteristics(id), total.GetRelativeCoordinates(id));
        analyze(component, static_cast<int>(s));
    }
    return results;
}

void WriteMultiComponentPairStatistics(std::ostream &out,
    const std::vector<ComponentPairStatistics> &results,
    const MultiComponentPairStatOptions &options) {
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::setprecision(std::numeric_limits<double>::max_digits10);
    text << "# PairStat: species -1 is the combined packing\n"
         << "# k_max " << options.kMax << " k_bin_width " << options.kBinWidth
         << " r_max " << options.rMax << " r_step " << options.rStep
         << " window_samples " << options.windowSamples << " seed " << options.seed << '\n'
         << "# S: x=k, value=S, aux1=k standard error, aux2=S standard error\n"
         << "# NV: x=R, value=number variance, aux1/aux2=unused\n"
         << "# VFV: x=R, value=volume-fraction variance, aux1=ensemble error (unavailable for one config), aux2=running average\n"
         << "# VFV assumes nonoverlapping spheres; overlaps flag warns of additive-volume interpretation\n"
         << "# quantity species particles phi overlaps x value aux1 aux2\n";
    for (const auto &result : results) {
        text << "# component " << result.species << " particles " << result.particles
             << " phi " << result.packingFraction << " overlaps " << result.overlaps << '\n';
        auto table = [&](const char *name, const std::vector<GeometryVector> &rows) {
            for (const auto &row : rows) {
                text << name << '\t' << result.species << '\t' << result.particles << '\t'
                     << result.packingFraction << '\t' << result.overlaps;
                for (int d = 0; d < 4; ++d) text << '\t' << row.x[d];
                text << '\n';
            }
        };
        table("S", result.structureFactor);
        table("NV", result.numberVariance);
        table("VFV", result.volumeFractionVariance);
    }
    const std::string data = text.str();
    out.write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!out) throw std::runtime_error("Could not write PairStat report.");
}
