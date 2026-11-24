/** \file bench_step.cpp
 * \brief Wall time per optimizer step and per warm evaluation of GetMultiCCO's potential, for one variant per process.
 * Variants share the fixture of fixtures.h (three equal species, 2D, contact strength 100) and the same random
 * initial state per seed; see README.md in this directory and bench_step.py, which runs the comparison. */
#include "fixtures.h"
#include <StructureOptimization.h>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <limits>

size_t Verbosity = 0;
using namespace threading_test;
using Clock = std::chrono::steady_clock;

namespace {
/** \brief Seconds since a time point.
 *  \param start  Earlier time point.
 *  \return Elapsed wall time [s]. */
double Since(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

/** \brief Potential that counts the optimizer-visible calls and accumulates their wall time. */
struct TimedPotential : MultiStealthyPotential {
    using MultiStealthyPotential::MultiStealthyPotential;
    size_t evaluations = 0;	//< Energy and EnergyAndAllForce calls
    double inside = 0.0;	//< wall time in SetConfiguration, Energy and EnergyAndAllForce [s]
    /** \brief Timed SetConfiguration.
     *  \param c  Configuration to bind. */
    void SetConfiguration(const Configuration & c) override
    { const auto t = Clock::now(); MultiStealthyPotential::SetConfiguration(c); inside += Since(t); }
    /** \brief Timed and counted Energy.
     *  \return Total energy. */
    double Energy() override
    { const auto t = Clock::now(); const double e = MultiStealthyPotential::Energy(); inside += Since(t); ++evaluations; return e; }
    /** \brief Timed and counted EnergyAndAllForce.
     *  \param force  Receives the forces.
     *  \return Total energy. */
    double EnergyAndAllForce(std::vector<GeometryVector> & force) override
    { const auto t = Clock::now(); const double e = MultiStealthyPotential::EnergyAndAllForce(force); inside += Since(t); ++evaluations; return e; }
};

/** \brief FNV-1a hash of the relative coordinates of a configuration, to compare states across processes.
 *  \param config  Configuration.
 *  \return 64-bit hash of the coordinate bits in particle order. */
std::uint64_t Hash(const Configuration & config)
{
    std::uint64_t hash = 1469598103934665603ull;
    for (size_t i = 0; i < config.NumParticle(); ++i)
        for (DimensionType d = 0; d < config.GetDimension(); ++d) {
            std::uint64_t bits;
            const double value = config.GetRelativeCoordinates(i).x[d];
            std::memcpy(&bits, &value, sizeof(bits));
            hash = (hash ^ bits) * 1099511628211ull;
        }
    return hash;
}

/** \brief Whether a '+'-separated variant string contains a token.
 *  \param variant  Variant, e.g. "A+grid+fast".
 *  \param token    Token, e.g. "grid".
 *  \return True if the token is one of the parts. */
bool Has(const std::string & variant, const std::string & token)
{
    return ("+" + variant + "+").find("+" + token + "+") != std::string::npos;
}
}

/** \brief Run one measurement and print one JSON line.
 *  \param argc  Argument count.
 *  \param argv  eval|opt VARIANT N_PER_SPECIES THREADS SEED COUNT REPETITION [HISTORY [CONTACT [PHI]]].
 *  \return 0, or 1 after printing an error. */
int main(int argc, char ** argv)
{
    try {
        if (argc < 8 || argc > 11)
            throw std::invalid_argument("Usage: bench_step.out eval|opt VARIANT N_PER_SPECIES THREADS SEED COUNT REPETITION [HISTORY [CONTACT [PHI]]]\n"
                "CONTACT: relative contact strength (default 100); PHI: packing fraction (default 0.70)\n"
                "VARIANT: direct | A | B, then any of +grid (contact grid, bitwise the generic walk), +plis (in-repository L-BFGS, NLopt's iterates), "
                "+fast (in-repository L-BFGS, reordered dot products), +team (the same on all threads), +reuse (no configuration copies per evaluation), +pair (two force components per Type 2 transform); "
                "other suffixes, such as +active, only label a run");
        const std::string mode = argv[1], variant = argv[2];
        const int particles = std::stoi(argv[3]), threads = std::stoi(argv[4]);
        const int seed = std::stoi(argv[5]), count = std::stoi(argv[6]), repetition = std::stoi(argv[7]);
        if (argc >= 9) {
            const int history = std::stoi(argv[8]);
            Require(history > 0, "History must be positive");
            StructureOptimization_LBFGSHistory = static_cast<size_t>(history);
        }
        const double contact = argc >= 10 ? std::stod(argv[9]) : 100.0, phi = argc >= 11 ? std::stod(argv[10]) : 0.70;
        Require(contact > 0.0 && phi > 0.0 && phi < 1.0, "Contact strength and packing fraction must be positive");
        const bool direct = Has(variant, "direct"), sequential = Has(variant, "B");
        Require((mode == "eval" || mode == "opt") && (direct || sequential || Has(variant, "A")) && particles > 1 && threads > 0 && count > 0,
            "Invalid benchmark arguments");
        const int lbfgs = Has(variant, "team") ? 3 : Has(variant, "fast") ? 2 : Has(variant, "plis") ? 1 : 0;
        omp_set_num_threads(threads);
        ProgramStart = std::time(nullptr);
        TimeLimit = std::numeric_limits<time_t>::max();
        const std::vector<size_t> sizes(3, static_cast<size_t>(particles));
        auto packing = Packing(2, sizes, seed);
        // The fixture's radii give packing fraction 0.70 in 2D; rescale them for another one.
        for (size_t i = 0; i < packing.NumParticle(); ++i) packing.GetCharacteristics(i) *= std::sqrt(phi / 0.70);
        Configuration config(packing, "a");
        auto potential = MakePotential<TimedPotential>(packing, sizes, threads, sequential, direct ? Backend::Direct : Backend::FinufftType1, contact);
        potential->SetContactGridEnabled(Has(variant, "grid"));
        potential->ReuseConfigurationStorage = Has(variant, "reuse");
        potential->SetFinufftPairedForces(Has(variant, "pair"));
        // Reference for the correctness checks: direct sums and the generic contact walk.
        auto reference = MakePotential(packing, sizes, threads, false, Backend::Direct, contact);
        reference->SetContactGridEnabled(false);
        Compare(Evaluate(*potential, config), Evaluate(*reference, config));
        for (int warm = 0; warm < 10; ++warm) Evaluate(*potential, config);

        double wall = 0.0, cpu = 0.0;
        size_t evaluations = 0;
        std::uint64_t stateHash = Hash(config);
        std::cout << std::setprecision(17);
        if (mode == "eval") {
            // Three states with different positions, the same in every process: the initial one and two reached by the
            // unmodified schedule (Option A, generic contact walk, NLopt), which repeats bitwise.
            std::vector<Configuration> states(1, config);
            {
                auto control = MakePotential(packing, sizes, threads, false, Backend::FinufftType1, contact);
                control->SetContactGridEnabled(false);
                StructureOptimization_LBFGSImplementation = 0;
                Configuration work(config);
                QuietOutput quiet;
                for (int steps : {40, 160}) {
                    RelaxStructure_NLOPT(work, *control, 0.0, 0, 0.0, steps);
                    states.push_back(work);
                }
            }
            for (const auto & state : states) {
                stateHash = stateHash * 31 + Hash(state);
                Compare(Evaluate(*potential, state), Evaluate(*reference, state));
            }
            potential->evaluations = 0;
            const auto start = Clock::now();
            const auto cpuStart = std::clock();
            for (int i = 0; i < count; ++i) Evaluate(*potential, states[i % states.size()]);
            cpu = double(std::clock() - cpuStart) / CLOCKS_PER_SEC;
            wall = Since(start);
            evaluations = potential->evaluations;
        }
        else {
            StructureOptimization_LBFGSImplementation = lbfgs;
            StructureOptimization_ReuseConfiguration = Has(variant, "reuse");
            potential->evaluations = 0;
            potential->inside = 0.0;
            QuietOutput quiet;
            const auto start = Clock::now();
            const auto cpuStart = std::clock();
            RelaxStructure_NLOPT(config, *potential, 0.0, 0, 0.0, count);
            cpu = double(std::clock() - cpuStart) / CLOCKS_PER_SEC;
            wall = Since(start);
            evaluations = potential->evaluations;
        }
        const double inside = potential->inside;
        const auto result = Evaluate(*potential, config);
        const auto check = Evaluate(*reference, config);
        Compare(result, check);
        double maxForce = 0.0;
        for (const auto & force : check.forces)
            for (double value : force.x) maxForce = std::max(maxForce, std::abs(value));
        std::cout << "{\"mode\":\"" << mode << "\",\"variant\":\"" << variant << "\",\"particles_per_species\":" << particles
                  << ",\"species\":3,\"threads\":" << threads << ",\"seed\":" << seed << ",\"repetition\":" << repetition
                  << ",\"history\":" << StructureOptimization_LBFGSHistory << ",\"contact_strength\":" << contact << ",\"phi\":" << phi
                  << ",\"count\":" << count
                  << ",\"wall_seconds\":" << wall << ",\"cpu_seconds\":" << cpu << ",\"evaluations\":" << evaluations
                  << ",\"ms_per_step\":" << 1e3 * wall / double(evaluations)
                  << ",\"potential_ms_per_step\":" << (mode == "opt" ? 1e3 * inside / double(evaluations) : 1e3 * wall / double(evaluations))
                  << ",\"optimizer_ms_per_step\":" << (mode == "opt" ? 1e3 * (wall - inside) / double(evaluations) : 0.0)
                  << ",\"utilization_percent\":" << 100 * cpu / (threads * wall)
                  << ",\"energy\":" << result.energy << ",\"direct_energy\":" << check.energy
                  << ",\"accepted\":" << (check.energy < 1e-14 ? "true" : "false") << ",\"max_force\":" << maxForce
                  << ",\"contact_grid\":" << (potential->ContactGridInUse() ? "true" : "false") << ",\"lbfgs\":" << lbfgs
                  << ",\"state_hash\":\"" << std::hex << (mode == "eval" ? stateHash : Hash(config)) << std::dec << "\"}\n";
    }
    catch (const std::exception & error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
