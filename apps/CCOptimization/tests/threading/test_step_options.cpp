/** \file test_step_options.cpp
 * \brief The options that shorten an optimizer step, each against the path it replaces.
 * - Paired Type 2 forces (two components per transform) agree with unpaired Type 2 forces and with direct sums.
 * - Reusing the configuration (the optimizer's workspace and the potential's stored copy) changes no bit of an
 *   evaluation or of an optimization end point.
 * - The in-repository L-BFGS in strict mode ends at the bits of the NLopt library; its fast mode ends at the same
 *   bits on one thread and on all threads.
 * The contact grid has its own test (test_contact_grid.cpp), the L-BFGS its own in potential/tests/luksan_plis. */
#include "fixtures.h"
#include <StructureOptimization.h>
#include <cstring>
#include <ctime>
#include <limits>

size_t Verbosity = 0;
using namespace threading_test;

namespace {
int failures = 0;

/** \brief Report one check.
 *  \param ok     Whether the check passed.
 *  \param label  What was checked. */
void Check(bool ok, const std::string & label)
{
    failures += !ok;
    std::cout << (ok ? "ok   " : "FAIL ") << label << '\n';
}

/** \brief Largest force difference between two results, relative to the largest reference force component (at least 1).
 *  \param a  Result.
 *  \param b  Reference result.
 *  \return The scaled maximum difference. */
double ForceDifference(const Result & a, const Result & b)
{
    double scale = 1.0, worst = 0.0;
    for (const auto & f : b.forces)
        for (double v : f.x) scale = std::max(scale, std::abs(v));
    for (size_t i = 0; i < a.forces.size(); ++i)
        for (size_t d = 0; d < ::MaxDimension; ++d) worst = std::max(worst, std::abs(a.forces[i].x[d] - b.forces[i].x[d]));
    return worst / scale;
}

/** \brief Whether two results are bitwise equal in energy and every force component.
 *  \param a  First result.
 *  \param b  Second result.
 *  \return True if identical. */
bool Identical(const Result & a, const Result & b)
{
    if (a.forces.size() != b.forces.size() || std::memcmp(&a.energy, &b.energy, sizeof(double)) != 0) return false;
    for (size_t i = 0; i < a.forces.size(); ++i)
        if (std::memcmp(a.forces[i].x, b.forces[i].x, sizeof(a.forces[i].x)) != 0) return false;
    return true;
}

/** \brief Whether two configurations hold bitwise the same relative coordinates.
 *  \param a  First configuration.
 *  \param b  Second configuration.
 *  \return True if identical. */
bool SameCoordinates(const Configuration & a, const Configuration & b)
{
    if (a.NumParticle() != b.NumParticle()) return false;
    for (size_t i = 0; i < a.NumParticle(); ++i)
        if (std::memcmp(a.GetRelativeCoordinates(i).x, b.GetRelativeCoordinates(i).x, sizeof(a.GetRelativeCoordinates(i).x)) != 0) return false;
    return true;
}

/** \brief End point of a bounded optimization from the fixture's initial state.
 *  \param packing         Initial packing.
 *  \param sizes           Particles per species.
 *  \param implementation  StructureOptimization_LBFGSImplementation.
 *  \param reuse           Reuse the configuration in the optimizer and in the potential.
 *  \return The configuration after RelaxStructure_NLOPT with 60 evaluations per scaling phase. */
Configuration EndPoint(const SpherePacking & packing, const std::vector<size_t> & sizes, int implementation, bool reuse)
{
    Configuration config(packing, "a");
    auto potential = MakePotential(packing, sizes, 4, false);
    potential->ReuseConfigurationStorage = reuse;
    StructureOptimization_ReuseConfiguration = reuse;
    StructureOptimization_LBFGSImplementation = implementation;
    QuietOutput quiet;
    RelaxStructure_NLOPT(config, *potential, 0.0, 0, 0.0, 60);
    return config;
}
}

/** \brief Run every check.
 *  \return Number of failed checks. */
int main()
{
    try {
        omp_set_num_threads(4);
        ProgramStart = std::time(nullptr);
        TimeLimit = std::numeric_limits<time_t>::max();
        // Paired forces.
        for (int dim = 1; dim <= 3; ++dim)
            for (bool skew : {false, true}) {
                const std::vector<size_t> sizes = {1000, 2000};
                const auto packing = Packing(dim, sizes, 21 + dim, skew);
                Configuration config(packing, "a");
                const std::string label = "dim " + std::to_string(dim) + (skew ? " skew" : " orthogonal");
                auto paired = MakePotential(packing, sizes, 4, false), plain = MakePotential(packing, sizes, 4, false);
                auto direct = MakePotential(packing, sizes, 4, false, Backend::Direct);
                paired->SetFinufftPairedForces(true);
                const Result a = Evaluate(*paired, config), b = Evaluate(*plain, config), c = Evaluate(*direct, config);
                size_t type2 = 0;
                for (auto & species : paired->species_constraints)
                    type2 += species.CCPotential->GetLastForceBackendUsed() == ShiftedCCPotential::ForceBackend::FinufftType2;
                Check(type2 == sizes.size(), label + ": both species use Type 2 forces");
                Check(std::memcmp(&a.energy, &b.energy, sizeof(double)) == 0 && ForceDifference(a, b) <= 1e-11 && ForceDifference(a, c) <= 1e-10,
                    label + ": paired forces agree with unpaired Type 2 and with direct sums");
                Check(Identical(Evaluate(*paired, config), a), label + ": paired forces repeat bitwise");
                paired->SetFinufftPairedForces(false);
                Check(Identical(Evaluate(*paired, config), b), label + ": switching pairing off restores the unpaired forces bitwise");
            }
        // Sequential schedule with pairing.
        {
            const std::vector<size_t> sizes = {800, 800};
            const auto packing = Packing(2, sizes, 5);
            Configuration config(packing, "a");
            auto sequential = MakePotential(packing, sizes, 4, true), direct = MakePotential(packing, sizes, 4, false, Backend::Direct);
            sequential->SetFinufftPairedForces(true);
            Check(ForceDifference(Evaluate(*sequential, config), Evaluate(*direct, config)) <= 1e-10, "sequential schedule: paired forces agree with direct sums");
        }
        // Reuse and L-BFGS implementations: end points of a bounded optimization.
        {
            const std::vector<size_t> sizes = {400, 400, 400};
            const auto packing = Packing(2, sizes, 17);
            Configuration config(packing, "a");
            auto kept = MakePotential(packing, sizes, 4, false), copied = MakePotential(packing, sizes, 4, false);
            copied->ReuseConfigurationStorage = false;
            Configuration moved(config);
            moved.MoveParticle(3, GeometryVector(0.25, 0.75));
            bool same = true;
            for (const Configuration * state : {&config, &moved, &config})
                same = same && Identical(Evaluate(*kept, *state), Evaluate(*copied, *state));
            Check(same, "stored configuration reused: evaluations bitwise equal to fresh copies");
            const Configuration nlopt = EndPoint(packing, sizes, 0, false);
            Check(!SameCoordinates(nlopt, config), "the bounded optimization moves the particles");
            Check(SameCoordinates(EndPoint(packing, sizes, 0, true), nlopt), "configuration reuse: end point bitwise unchanged");
            Check(SameCoordinates(EndPoint(packing, sizes, 1, true), nlopt), "L-BFGS ring: end point bitwise that of the NLopt library");
            const Configuration fast = EndPoint(packing, sizes, 2, true);
            Check(!SameCoordinates(fast, nlopt), "L-BFGS fast: end point differs from NLopt's (another summation order)");
            Check(SameCoordinates(EndPoint(packing, sizes, 2, false), fast), "L-BFGS fast: configuration reuse changes no bit");
            // Enough variables (12,600) for the thread team of the fast mode to run.
            const std::vector<size_t> large = {2100, 2100, 2100};
            const auto big = Packing(2, large, 18);
            Check(SameCoordinates(EndPoint(big, large, 3, true), EndPoint(big, large, 2, true)),
                "L-BFGS fast: one thread and the thread team end at the same bits");
            StructureOptimization_LBFGSImplementation = 0;
            StructureOptimization_ReuseConfiguration = true;
        }
    }
    catch (const std::exception & error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << (failures == 0 ? "ALL PASSED\n" : "FAILURES\n");
    return failures;
}
