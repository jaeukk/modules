#ifndef CCO_THREADING_FIXTURES_H
#define CCO_THREADING_FIXTURES_H

#include "MultiStealthyPotential.h"
#include <RandomGenerator.h>
#include <etc.h>
#include <sstream>
#include <stdexcept>

namespace threading_test {
using Backend = ShiftedCCPotential::RhoBackend;

inline void Require(bool ok, const char * message)
{
    if (!ok) throw std::runtime_error(message);
}

// Scope console suppression to serial fixture construction, before workers start.
class QuietOutput {
    std::ostringstream sink;
    std::streambuf * previous = std::cout.rdbuf(sink.rdbuf());
public:
    ~QuietOutput() { std::cout.rdbuf(previous); }
    QuietOutput() = default;
    QuietOutput(const QuietOutput &) = delete;
    QuietOutput & operator=(const QuietOutput &) = delete;
    QuietOutput(QuietOutput &&) = delete;
    QuietOutput & operator=(QuietOutput &&) = delete;
};

inline SpherePacking Packing(int dim, const std::vector<size_t> & sizes, unsigned seed = 17, bool skew = false)
{
    size_t n = 0;
    for (auto size : sizes) n += size;
    const double length = std::pow(static_cast<double>(n), 1.0 / dim);
    GeometryVector basis[3] = {GeometryVector(dim), GeometryVector(dim), GeometryVector(dim)};
    for (int d = 0; d < dim; ++d) basis[d].x[d] = length;
    if (skew && dim > 1) basis[1].x[0] = 0.21 * length;
    SpherePacking packing(dim, basis, 1.0);
    RandomGenerator random(seed);
    const double radius = std::pow(0.70 / HyperSphere_Volume(dim, 1.0), 1.0 / dim);
    for (size_t i = 0; i < n; ++i) {
        GeometryVector position(dim);
        for (int d = 0; d < dim; ++d) position.x[d] = random.RandomDouble();
        packing.Insert(radius, position);
    }
    return packing;
}

template<class T = MultiStealthyPotential>
inline std::unique_ptr<T> MakePotential(const SpherePacking & packing,
    const std::vector<size_t> & sizes, int threads, bool sequential,
    Backend backend = Backend::FinufftType1, double contact = 100.0)
{
    QuietOutput quiet;
    auto potential = std::make_unique<T>(packing.GetDimension(), contact);
    potential->ParallelNumThread = threads;
    potential->SetPacking(packing);
    potential->SetNumSpecies(sizes.size());
    size_t first = 0;
    for (size_t s = 0; s < sizes.size(); ++s) {
        for (size_t j = 0; j < sizes[s]; ++j) potential->SetSpeciesOfParticle(first++, s);
        if (sizes[s] == 0) continue;
        const double k = 2 * pi * std::pow(2 * packing.GetDimension() * double(sizes[s]) * .3 /
            (packing.PeriodicVolume() * HyperSphere_Volume(packing.GetDimension(), 1)),
            1.0 / packing.GetDimension());
        for (const auto & mode : GetKs(packing, k, k, 1)) potential->AddConstraint(s, mode, {1.0});
    }
    potential->SetRhoBackend(backend);
    potential->SetSequentialFinufftSpecies(sequential);
    potential->SetPacking(packing);
    return potential;
}

struct Result {
    double energy;
    std::vector<GeometryVector> forces;
};
inline Result Evaluate(MultiStealthyPotential & potential, const Configuration & config)
{
    potential.SetConfiguration(config);
    Result result;
    result.energy = potential.EnergyAndAllForce(result.forces);
    return result;
}
inline void Compare(const Result & actual, const Result & expected, double tolerance = 2e-10)
{
    Require(actual.forces.size() == expected.forces.size(), "Force counts differ");
    Require(std::isfinite(actual.energy) && std::abs(actual.energy - expected.energy) <=
        tolerance * std::max(1.0, std::abs(expected.energy)), "Energy mismatch");
    double scale = 1;
    for (const auto & force : expected.forces)
        for (double value : force.x) scale = std::max(scale, std::abs(value));
    for (size_t i = 0; i < actual.forces.size(); ++i)
        for (size_t d = 0; d < ::MaxDimension; ++d)
            Require(std::isfinite(actual.forces[i].x[d]) &&
                std::abs(actual.forces[i].x[d] - expected.forces[i].x[d]) <= tolerance * scale,
                "Force component mismatch");
}
}
#endif
