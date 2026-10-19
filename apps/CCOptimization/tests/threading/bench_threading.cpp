#include "fixtures.h"
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
int main(int argc, char ** argv)
{
    try {
        if (argc != 8) throw std::invalid_argument("Usage: bench_threading.out eval|opt A|B N_PER_SPECIES THREADS SEED EVALUATIONS REPETITION");
        const std::string mode = argv[1], option = argv[2];
        const int particles = std::stoi(argv[3]), threads = std::stoi(argv[4]);
        const int seed = std::stoi(argv[5]), iterations = std::stoi(argv[6]);
        const int repetition = std::stoi(argv[7]);
        Require((mode == "eval" || mode == "opt") && (option == "A" || option == "B") &&
            particles > 1 && threads > 0 && iterations > 0, "Invalid benchmark arguments");
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
        const auto start = Clock::now();
        const auto cpuStart = std::clock();
        if (mode == "eval") {
            for (int i = 0; i < iterations; ++i) Evaluate(*potential, config);
        }
        else {
            StructureOptimization_LBFGSHistory = 100;
            QuietOutput quiet;
            RelaxStructure_NLOPT(config, *potential, 0.0, 0, 0.0, iterations);
        }
        const double cpu = double(std::clock() - cpuStart) / CLOCKS_PER_SEC;
        const double wall = std::chrono::duration<double>(Clock::now() - start).count();
        const size_t evaluations = potential->evaluations;
        const auto result = Evaluate(*potential, config);
        const auto direct = Evaluate(*reference, config);
        Compare(result, direct);
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
        std::cout << "]}\n";
    }
    catch (const std::exception & error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
