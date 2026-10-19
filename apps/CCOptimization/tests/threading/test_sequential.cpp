#include "fixtures.h"
#include <atomic>
#include <map>

size_t Verbosity = 0;
using namespace threading_test;

namespace {
struct Plan { int type; int threads; int transforms; };
std::map<finufft_plan, Plan> plans;
std::atomic<int> active{0}, peak{0}, failExecute{0};
bool failPlan = false;
size_t made = 0;
int expectedThreads = 1;
bool sequential = false;
}
extern "C" {
int __real_finufft_makeplan(int, int, const int64_t *, int, int, double, finufft_plan *, const finufft_opts *);
int __real_finufft_destroy(finufft_plan);
int __real_finufft_execute(finufft_plan, std::complex<double> *, std::complex<double> *);
}
extern "C" int __wrap_finufft_makeplan(int type, int dim, const int64_t * modes, int sign,
    int transforms, double tolerance, finufft_plan * plan, const finufft_opts * options)
{
    Require(!omp_in_parallel(), "Plan created inside an OpenMP region");
    ++made;
    const int status = __real_finufft_makeplan(type, dim, modes, sign, transforms, tolerance, plan, options);
    if (*plan) plans.emplace(*plan, Plan{type, options->nthreads, transforms});
    if (type == 2 && failPlan) { failPlan = false; return 1; }
    return status;
}
extern "C" int __wrap_finufft_destroy(finufft_plan plan)
{
    Require(plans.erase(plan) == 1, "Unknown or double-destroyed FINUFFT plan");
    return __real_finufft_destroy(plan);
}
extern "C" int __wrap_finufft_execute(finufft_plan plan, std::complex<double> * values, std::complex<double> * modes)
{
    const auto info = plans.at(plan);
    if (sequential) {
        Require(!omp_in_parallel(), "Sequential transform called inside an outer OpenMP region");
        Require(info.threads == expectedThreads, "Wrong internal FINUFFT thread count");
        Require(info.type != 2 || info.transforms >= 1, "Empty force batch");
    }
    if (failExecute == info.type && failExecute.exchange(0) == info.type) return 1;
    const int running = ++active;
    int previous = peak;
    while (running > previous && !peak.compare_exchange_weak(previous, running)) {}
    const int status = __real_finufft_execute(plan, values, modes);
    --active;
    return status;
}

namespace {
struct ThrowingCC : ShiftedCCPotential {
    using ShiftedCCPotential::ShiftedCCPotential;
    bool SupportsIndexedView() const override { return true; }
    double Energy() override { throw std::runtime_error("Injected species exception"); }
};
void CheckCase(int dim, const std::vector<size_t> & sizes, int threads, bool skew, double contact)
{
    const auto packing = Packing(dim, sizes, 17, skew);
    Configuration config(packing, "a");
    auto reference = MakePotential(packing, sizes, threads, false, Backend::Direct, contact);
    auto candidate = MakePotential(packing, sizes, threads, true, Backend::FinufftType1, contact);
    const auto direct = Evaluate(*reference, config);
    sequential = true;
    expectedThreads = std::min(threads, omp_get_thread_limit());
    peak = 0;
    Compare(Evaluate(*candidate, config), direct);
    Require(peak <= 1, "Species transforms overlap in sequential mode");
    for (const auto & entry : plans)
        if (entry.second.type == 2)
            Require(entry.second.transforms == dim, "Force components are not batched");
    const size_t warmPlans = made;
    Compare(Evaluate(*candidate, config), direct);
    Require(made == warmPlans, "Warm evaluation rebuilt plans");

    candidate->SetConfiguration(config);
    Result separate;
    candidate->AllForce(separate.forces);
    separate.energy = candidate->Energy();
    Compare(separate, direct);
    candidate->SetConfiguration(config);
    separate.energy = candidate->Energy();
    candidate->AllForce(separate.forces);
    Compare(separate, direct);

    // Change positions: reusing plans must still use the new points.
    GeometryVector moved(config.GetRelativeCoordinates(0));
    moved.x[0] += .013;
    config.MoveParticle(0, moved);
    Compare(Evaluate(*candidate, config), Evaluate(*reference, config));
    Require(made == warmPlans, "Particle motion rebuilt plans");
    for (auto & species : candidate->species_constraints) species.CCPotential->SetFinufftMaxModes(1);
    Compare(Evaluate(*candidate, config), Evaluate(*reference, config));
    sequential = false;
    candidate->SetSequentialFinufftSpecies(false);
    for (const auto & species : candidate->species_constraints)
        Require(species.CCPotential->GetFinufftNumThreads() == 1, "Option A did not restore internal threads");
    Compare(Evaluate(*candidate, config), Evaluate(*reference, config));
    candidate->SetSequentialFinufftSpecies(true);
    sequential = true;
    candidate->SetRhoBackend(Backend::Direct);
    Compare(Evaluate(*candidate, config), Evaluate(*reference, config));
    sequential = false;
}
void Failures()
{
    const std::vector<size_t> sizes(3, 100);
    const auto packing = Packing(2, sizes);
    Configuration config(packing, "a");
    auto reference = MakePotential(packing, sizes, 4, false, Backend::Direct);
    const auto direct = Evaluate(*reference, config);
    expectedThreads = std::min(4, omp_get_thread_limit());
    sequential = true;
    for (int kind : {1, 2}) {
        auto candidate = MakePotential(packing, sizes, 4, true);
        failExecute = kind;
        Compare(Evaluate(*candidate, config), direct);
        Require(failExecute == 0, "Execution failure was not exercised");
        Compare(Evaluate(*candidate, config), direct);
    }
    {
        auto candidate = MakePotential(packing, sizes, 4, true);
        failPlan = true;
        Compare(Evaluate(*candidate, config), direct);
        Require(!failPlan, "Plan failure was not exercised");
    }
    {
        auto candidate = MakePotential(packing, sizes, 4, true);
        candidate->species_constraints[0].CCPotential = std::make_unique<ThrowingCC>(2);
        candidate->SetConfiguration(config);
        bool threw = false;
        try { Evaluate(*candidate, config); }
        catch (const std::runtime_error & error) { threw = std::string(error.what()) == "Injected species exception"; }
        Require(threw, "Species exception was swallowed");
        Require(candidate->species_constraints[0].CCPotential->ParallelNumThread == 4,
            "Exception did not restore the species thread setting");
        candidate->species_constraints[0].CCPotential = std::make_unique<ShiftedCCPotential>(2);
        candidate->SetConfiguration(config);
        Evaluate(*candidate, config); // Reuse after the exception; async worker must be finished.
    }
    sequential = false;
}
void FiniteDifference()
{
    const std::vector<size_t> sizes(3, 30);
    const auto packing = Packing(2, sizes);
    Configuration config(packing, "a");
    auto candidate = MakePotential(packing, sizes, 4, true);
    const auto result = Evaluate(*candidate, config);
    const double step = 1e-6;
    const auto original = config.GetRelativeCoordinates(0);
    for (int d = 0; d < 2; ++d) {
        auto moved = original;
        moved.x[d] += step;
        config.MoveParticle(0, moved);
        candidate->SetConfiguration(config);
        const double plus = candidate->Energy();
        moved.x[d] -= 2 * step;
        config.MoveParticle(0, moved);
        candidate->SetConfiguration(config);
        const double minus = candidate->Energy();
        config.MoveParticle(0, original);
        const double expected = -result.forces[0].Dot(config.GetBasisVector(d));
        Require(std::abs((plus - minus) / (2 * step) - expected) < 2e-5 * std::max(1.0, std::abs(expected)),
            "Finite-difference gradient mismatch");
    }
}
}
int main()
{
    try {
        for (int dim : {1, 2, 3})
            for (int threads : {1, 4})
                CheckCase(dim, {100, 100, 100}, threads, dim == 3, 100);
        CheckCase(2, {30, 0, 50}, 4, true, 0);
        Failures();
        FiniteDifference();
        Require(plans.empty(), "FINUFFT plans leaked");
        std::cout << "PASS: sequential species, batched components, internal thread budget, direct reference,\n"
                     "cache reuse, force/energy orders, fallbacks, exceptions, and finite differences.\n";
    }
    catch (const std::exception & error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
