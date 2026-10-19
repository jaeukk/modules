#include "fixtures.h"
#include <array>
#include <fstream>
#include <nlopt.h>
#include <StructureOptimization.h>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <limits>

size_t Verbosity = 0;
using namespace threading_test;
using Clock = std::chrono::steady_clock;

namespace {
struct CountingPotential : MultiStealthyPotential {
    using MultiStealthyPotential::MultiStealthyPotential;
    size_t evaluations = 0;
    double Energy() override { ++evaluations; return MultiStealthyPotential::Energy(); }
    double EnergyAndAllForce(std::vector<GeometryVector> & force) override
    { ++evaluations; return MultiStealthyPotential::EnergyAndAllForce(force); }
};
}
namespace {
std::array<std::pair<int, int>, 3> phaseReports{};
size_t phaseCount = 0;
}
extern "C" nlopt_result __real_nlopt_optimize(nlopt_opt, double *, double *);
extern "C" nlopt_result __wrap_nlopt_optimize(nlopt_opt optimizer, double * coordinates, double * value)
{
    const auto status = __real_nlopt_optimize(optimizer, coordinates, value);
    if (phaseCount < phaseReports.size())
        phaseReports[phaseCount] = {static_cast<int>(status), nlopt_get_numevals(optimizer)};
    ++phaseCount;
    return status;
}
int main(int argc, char ** argv)
{
    try {
        if (argc != 9) throw std::invalid_argument("Usage: bench_threading.out eval|opt A|B N_PER_SPECIES THREADS SEED EVALUATIONS REPETITION HISTORY");
        const std::string mode = argv[1], option = argv[2];
        const int particles = std::stoi(argv[3]), threads = std::stoi(argv[4]);
        const int seed = std::stoi(argv[5]), iterations = std::stoi(argv[6]);
        const int repetition = std::stoi(argv[7]);
        const int history = std::stoi(argv[8]);
        Require((mode == "eval" || mode == "opt") && (option == "A" || option == "B") &&
            particles > 1 && threads > 0 && iterations > 0 && history > 0, "Invalid benchmark arguments");
        omp_set_num_threads(threads);
        ProgramStart = std::time(nullptr);
        TimeLimit = std::numeric_limits<time_t>::max();
        const std::vector<size_t> sizes(3, static_cast<size_t>(particles));
        const auto packing = Packing(2, sizes, seed);
        Configuration config(packing, "a");
        auto potential = MakePotential<CountingPotential>(packing, sizes, threads, option == "B");
        auto reference = MakePotential(packing, sizes, threads, false, Backend::Direct);
        Compare(Evaluate(*potential, config), Evaluate(*reference, config));
        for (int warm = 0; warm < 10; ++warm) Evaluate(*potential, config);
        potential->evaluations = 0;
        {
            std::ofstream snapshot("initial.mspb", std::ios::binary);
            potential->ExportConfiguration(snapshot, MultiStealthyPotential::ExportFormat::Binary);
            snapshot.close();
            Require(bool(snapshot), "Initial snapshot write failed");
        }
        const auto start = Clock::now();
        const auto cpuStart = std::clock();
        if (mode == "eval") {
            for (int i = 0; i < iterations; ++i) Evaluate(*potential, config);
        }
        else {
            StructureOptimization_LBFGSHistory = static_cast<size_t>(history);
            QuietOutput quiet;
            RelaxStructure_NLOPT(config, *potential, 0.0, 0, 0.0, iterations);
        }
        const double cpu = double(std::clock() - cpuStart) / CLOCKS_PER_SEC;
        const double wall = std::chrono::duration<double>(Clock::now() - start).count();
        const size_t evaluations = potential->evaluations;
        const auto result = Evaluate(*potential, config);
        const auto direct = Evaluate(*reference, config);
        Compare(result, direct);
        if (mode == "opt") Require(phaseCount == 3, "Expected three NLopt phases");
        {
            std::ofstream snapshot("endpoint.mspb", std::ios::binary);
            potential->ExportConfiguration(snapshot, MultiStealthyPotential::ExportFormat::Binary);
            snapshot.close();
            Require(bool(snapshot), "Endpoint snapshot write failed");
        }
        double directForceMax = 0, directForceSum2 = 0;
        for (const auto & force : direct.forces)
            for (int d = 0; d < config.GetDimension(); ++d) {
                directForceMax = std::max(directForceMax, std::abs(force.x[d]));
                directForceSum2 += force.x[d] * force.x[d];
            }

        std::cout << std::setprecision(17) << "{\"mode\":\"" << mode << "\",\"option\":\"" << option
                  << "\",\"particles_per_species\":" << particles << ",\"species\":3,\"threads\":" << threads
                  << ",\"seed\":" << seed << ",\"repetition\":" << repetition
                  << ",\"wall_seconds\":" << wall << ",\"cpu_seconds\":" << cpu
                  << ",\"active_cpus\":" << cpu / wall
                  << ",\"utilization_percent\":" << 100 * cpu / (threads * wall)
                  << ",\"evaluations\":" << evaluations << ",\"energy\":" << result.energy
                  << ",\"direct_energy\":" << direct.energy << ",\"accepted\":" << (direct.energy < 1e-14 ? "true" : "false")
                  << ",\"rho_backends\":[";
        for (size_t s = 0; s < sizes.size(); ++s) {
            if (s) std::cout << ',';
            std::cout << static_cast<int>(potential->species_constraints[s].CCPotential->GetLastRhoBackendUsed());
        }
        std::cout << "],\"force_backends\":[";
        for (size_t s = 0; s < sizes.size(); ++s) {
            if (s) std::cout << ',';
            std::cout << static_cast<int>(potential->species_constraints[s].CCPotential->GetLastForceBackendUsed());
        }
        std::cout << "],\"internal_threads\":[";
        for (size_t s = 0; s < sizes.size(); ++s) {
            if (s) std::cout << ',';
            std::cout << potential->species_constraints[s].CCPotential->GetFinufftNumThreads();
        }
        std::cout << "],\"history\":" << history
                  << ",\"direct_force_max_abs\":" << directForceMax
                  << ",\"direct_force_rms_component\":" << std::sqrt(directForceSum2 / (config.NumParticle() * config.GetDimension()))
                  << ",\"nlopt_phases\":[";
        for (size_t phase = 0; phase < phaseCount; ++phase) {
            if (phase) std::cout << ',';
            std::cout << "{\"status\":" << phaseReports[phase].first
                      << ",\"evaluations\":" << phaseReports[phase].second << '}';
        }
        std::cout << "]}\n";
    }
    catch (const std::exception & error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
