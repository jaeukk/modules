/** \file test_species_concurrent.cpp
 * \brief Checks for species-concurrent evaluation in MultiStealthyPotential (FINUFFT backend).
 *
 * For S = 1, 2, 5, 10 species (below and above the thread count), dimensions 1-3, orthogonal and skew
 * cells, contact on and off, and both backends at 4 and 8 threads:
 *  - energy and every force component match a single-threaded Direct reference within 1e-10
 *    (relative for energy, times max(1, max|F|) for forces);
 *  - AllForce without a prior Energy (MD order) equals the warm AllForce bitwise;
 *  - repeated evaluations on one potential and on a fresh potential are bitwise identical;
 *  - FINUFFT results are bitwise identical at 4 and 8 threads (no thread-count-dependent sums);
 *  - each species' rho is filled exactly once per evaluation (never lazily inside the force loop);
 *  - every FINUFFT species used the backend FinufftEligible predicts; species forced to fall back
 *    (grid cap, or a sheared cell tripping the sparse-grid guard) report Direct, and when all fall back
 *    the results equal the Direct backend bitwise; mixed eligible/fallback cases stay repeatable;
 *  - while the potentials of one thread count run, the sampled peak OS thread count never exceeds that
 *    team (no oversized outer team, no nested team, no extra library thread).
 * The makefile target runs it with default waiting, OMP_WAIT_POLICY=passive and nesting enabled.
 * Build and run: make test-species-concurrent */

#include "MultiStealthyPotential.h"
#include <RandomGenerator.h>
#include <etc.h>

#include <atomic>
#include <chrono>
#include <cstring>
#include <fstream>
#include <sstream>
#include <thread>

size_t Verbosity = 0;

namespace {

using Backend = ShiftedCCPotential::RhoBackend;

/** \brief One test case. */
struct Case {
	DimensionType dim;          //< spatial dimension
	bool skew;                  //< non-orthogonal cell
	std::vector<size_t> sizes;  //< particles per species
	double chi;                 //< per-species chi
	double e0;                  //< contact strength (0 disables contact)
	unsigned fallback = 0;      //< bit s: species s falls back to Direct (FinufftMaxModes = 1)
	double shear = 0;           //< if nonzero, second basis vector (shear L, L): sparse-grid fallback
};

/** \brief Species potential that counts rho cache misses, to catch a rho filled lazily inside the
 *  parallel force loop (a race whose racing threads write identical values, invisible in results). */
struct CountingCC : ShiftedCCPotential {
	std::atomic<int> fills{0};  //< GetRho calls that found the cache invalid

	/** \brief Construct an empty species potential.
	 *  \param d  Dimension. */
	explicit CountingCC(DimensionType d) : ShiftedCCPotential(d) {}

	/** \brief Opt in to MultiStealthyPotential's indexed species views: particles are read only through the base class.
	 *  \return true. */
	bool SupportsIndexedView() const override { return true; }

	/** \brief Count a cache miss, then fill rho as the base class does. */
	void GetRho() override
	{
		if (!this->RhoValid) ++this->fills;
		ShiftedCCPotential::GetRho();
	}
};

/** \brief Energy and forces of one evaluation. */
struct Result {
	double energy;
	std::vector<GeometryVector> force;
};

int failures = 0, checks = 0, maxPeak = 0;

/** \brief Record a failed check without stopping the suite.
 *  \param ok    Check result.
 *  \param what  Description printed on failure. */
void Check(bool ok, const std::string & what)
{
	++checks;
	if (!ok) {
		++failures;
		std::cout << "FAIL: " << what << '\n';
	}
}

/** \brief Bitwise equality of two evaluations.
 *  \param a    First result.
 *  \param b    Second result.
 *  \param dim  Dimension.
 *  \return     true if the energy and every force component have identical bits. */
bool Same(const Result & a, const Result & b, DimensionType dim)
{
	if (std::memcmp(&a.energy, &b.energy, sizeof(double)) != 0 || a.force.size() != b.force.size()) return false;
	for (size_t i = 0; i < a.force.size(); i++)
		if (std::memcmp(a.force[i].x, b.force[i].x, sizeof(double) * dim) != 0) return false;
	return true;
}

/** \brief Random configuration in a cell of volume N, optionally skewed.
 *  \param tc  Test case.
 *  \return    Configuration with species blocks contiguous in particle order. */
Configuration Positions(const Case & tc)
{
	size_t n = 0;
	for (size_t s : tc.sizes) n += s;
	const DimensionType d = tc.dim;
	const double L = std::pow(double(n), 1.0 / d);
	GeometryVector basis[3] = {GeometryVector(d), GeometryVector(d), GeometryVector(d)};
	for (DimensionType i = 0; i < d; i++)
		basis[i].x[i] = L;
	if (tc.skew && d >= 2) basis[1].x[0] = 0.27 * L;
	if (tc.shear != 0) basis[1].x[0] = tc.shear * L;
	if (tc.skew && d == 3) { basis[2].x[0] = -0.19 * L; basis[2].x[1] = 0.21 * L; }
	Configuration c(d, basis, 0.1 * L);
	RandomGenerator gen(77 + n + d);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x(d);
		for (DimensionType k = 0; k < d; k++) x.x[k] = gen.RandomDouble();
		c.Insert("a", x);
	}
	return c;
}

/** \brief Build a species-resolved potential bound to c.
 *  \param c        Configuration.
 *  \param tc       Test case.
 *  \param backend  Collective-density backend.
 *  \param threads  ParallelNumThread.
 *  \return         The potential. */
std::unique_ptr<MultiStealthyPotential> Make(const Configuration & c, const Case & tc, Backend backend, int threads)
{
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	auto p = std::make_unique<MultiStealthyPotential>(tc.dim, tc.e0);
	p->ParallelNumThread = threads;
	p->SetConfiguration(c);
	p->SetNumSpecies(tc.sizes.size());
	for (auto & sp : p->species_constraints)
		sp.CCPotential = std::make_unique<CountingCC>(tc.dim);
	const double radius = std::pow(0.2 / HyperSphere_Volume(tc.dim, 1.0), 1.0 / tc.dim);
	p->SetSecureRadii(std::vector<double>(c.NumParticle(), radius));
	size_t first = 0;
	for (size_t s = 0; s < tc.sizes.size(); s++) {
		for (size_t i = 0; i < tc.sizes[s]; i++)
			p->SetSpeciesOfParticle(first + i, s);
		first += tc.sizes[s];
		const double rhoS = double(tc.sizes[s]) / c.PeriodicVolume();
		const double K = 2 * pi * std::pow(2 * tc.dim * rhoS * tc.chi / HyperSphere_Volume(tc.dim, 1.0), 1.0 / tc.dim);
		for (const auto & k : GetKs(c, K, K, 1))
			p->AddConstraint(s, k, {1.0 + 0.1 * s});
	}
	p->SetRhoBackend(backend);
	for (size_t s = 0; s < tc.sizes.size(); s++)
		if (tc.fallback >> s & 1u)
			p->species_constraints[s].CCPotential->SetFinufftMaxModes(1);
	p->SetConfiguration(c);
	std::cout.rdbuf(old);
	return p;
}

/** \brief Evaluate in the MD order (AllForce, then Energy) or the optimizer order (Energy, then AllForce),
 *  and check that each species' rho was filled exactly once, before any force used it.
 *  \param p        Potential built by Make.
 *  \param c        Configuration to bind first (invalidates every cache).
 *  \param mdOrder  Call AllForce before Energy.
 *  \return         Energy and forces. */
Result Evaluate(MultiStealthyPotential & p, const Configuration & c, bool mdOrder)
{
	Result r;
	p.SetConfiguration(c);
	if (mdOrder) { p.AllForce(r.force); r.energy = p.Energy(); }
	else { r.energy = p.Energy(); p.AllForce(r.force); }
	for (auto & sp : p.species_constraints) {
		auto & cc = static_cast<CountingCC &>(*sp.CCPotential);
		Check(cc.fills == 1, "rho filled " + std::to_string(cc.fills) + " times in one evaluation" + (mdOrder ? " (MD order)" : ""));
		cc.fills = 0;
	}
	return r;
}

/** \brief Number of OS threads in this process.
 *  \return The Threads: field of /proc/self/status, or 0 if unavailable. */
int OsThreads()
{
	std::ifstream status("/proc/self/status");
	std::string key;
	while (status >> key)
		if (key == "Threads:") { int n = 0; status >> n; return n; }
	return 0;
}

/** \brief Run work at one team size while a helper thread samples the OS thread count (libgomp releases
 *  surplus threads, so a count taken afterwards misses an oversized team). When a smaller team follows a
 *  larger one, libgomp destroys the surplus pool threads, and when a larger one follows it creates new ones;
 *  the exiting and the new threads can coexist for a moment, so a process-wide peak over work with
 *  alternating team sizes can exceed the largest team without any oversubscription. So the pool is first set
 *  to the team size and allowed to settle (up to 2 s), and the peak is compared with that team only.
 *  Sampling can miss very short teams but never reports more threads than existed.
 *  \param team  Team size of every parallel region of work (regions of one thread add none).
 *  \param work  Work to observe.
 *  \return      Peak thread count seen during work, excluding the sampler. */
template <class F>
int PeakThreads(int team, F work)
{
	int size = 0;	// a region with a side effect: an empty one may be compiled away
#pragma omp parallel num_threads(team)
	{
#pragma omp single
		size = omp_get_num_threads();
	}
	Check(size == team, "settling team has " + std::to_string(size) + " threads, not " + std::to_string(team));
	for (int wait = 0; wait < 2000 && OsThreads() > team; wait++)
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	std::atomic<bool> done{false};
	int peak = 0;
	std::thread sampler([&] { while (!done) peak = std::max(peak, OsThreads() - 1); });
	work();
	done = true;
	sampler.join();
	return peak;
}

/** \brief Run every check for one case.
 *  \param tc  Test case. */
void Run(const Case & tc)
{
	std::ostringstream name;
	name << tc.dim << "D" << (tc.skew ? " skew" : "") << " S=" << tc.sizes.size() << " e0=" << tc.e0 << (tc.shear != 0 ? " sheared" : "") << (tc.fallback ? " fallback=" + std::to_string(tc.fallback) : std::string());
	const Configuration c = Positions(tc);
	const Result ref = Evaluate(*Make(c, tc, Backend::Direct, 1), c, false);
	double scale = 0;
	for (const auto & f : ref.force)
		for (DimensionType d = 0; d < tc.dim; d++) scale = std::max(scale, std::abs(f.x[d]));

	Result direct[2];  // direct[threads == 8]: Direct-backend result at that thread count
	for (Backend backend : {Backend::Direct, Backend::FinufftType1}) {
		const bool fast = backend == Backend::FinufftType1;
		Result firstT;
		for (int threads : {4, 8}) {
			const std::string label = name.str() + (fast ? " finufft" : " direct") + " T=" + std::to_string(threads);
			const int peak = PeakThreads(threads, [&] {
			auto p = Make(c, tc, backend, threads);
			const Result cold = Evaluate(*p, c, true), warm = Evaluate(*p, c, false);
			Check(Same(cold, warm, tc.dim), label + ": MD-order (cold AllForce) result differs from warm result");
			for (int rep = 0; rep < 5; rep++)
				Check(Same(Evaluate(*p, c, false), warm, tc.dim), label + ": repeated evaluation not bitwise identical");
			Check(Same(Evaluate(*Make(c, tc, backend, threads), c, false), warm, tc.dim), label + ": fresh potential not bitwise identical");

			double df = 0;
			for (size_t i = 0; i < ref.force.size(); i++)
				for (DimensionType d = 0; d < tc.dim; d++) df = std::max(df, std::abs(warm.force[i].x[d] - ref.force[i].x[d]));
			Check(std::abs(warm.energy - ref.energy) <= 1e-10 * std::abs(ref.energy),
			      label + ": energy rel diff " + std::to_string(std::abs(warm.energy - ref.energy) / std::abs(ref.energy)));
			Check(df <= 1e-10 * std::max(1.0, scale), label + ": force diff " + std::to_string(df));
			if (!fast)
				direct[threads == 8] = warm;
			else {
				size_t fallbacks = 0;
				for (size_t s = 0; s < tc.sizes.size(); s++) {
					const auto & cc = *p->species_constraints[s].CCPotential;
					const bool expectFallback = tc.shear != 0 || (tc.fallback >> s & 1u);
					fallbacks += expectFallback;
					Check(cc.FinufftEligible() == !expectFallback, label + ": FinufftEligible wrong for species " + std::to_string(s));
					Check(cc.GetLastRhoBackendUsed() == (expectFallback ? Backend::Direct : Backend::FinufftType1),
					      label + ": unexpected backend used by species " + std::to_string(s));
				}
				if (fallbacks == tc.sizes.size())
					Check(Same(warm, direct[threads == 8], tc.dim), label + ": all-fallback result differs from the Direct backend");
				else if (fallbacks == 0 && threads == 4) firstT = warm;
				else if (fallbacks == 0) Check(Same(warm, firstT, tc.dim), label + ": FINUFFT result depends on the thread count");
			}
			});
			Check(peak <= threads, label + ": peak of " + std::to_string(peak) + " OS threads exceeds the team");
			maxPeak = std::max(maxPeak, peak);
		}
	}
}

}  // namespace

/** \brief Run the suite.
 *  \return 0 if every check passed, 1 otherwise. */
int main()
{
	// FINUFFT's FFT (ducc0) starts an idle pool of hardware-1 threads on first use; with nthreads=1 it
	// never receives work. Size it to one so the OS thread count below measures OpenMP teams only.
	// The idle pool also exists at baseline (c352d19) and in the application, which leaves this unset.
	setenv("DUCC0_NUM_THREADS", "1", 1);
	if (!ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1)) {
		std::cout << "FINUFFT backend not compiled in; build with make USE_FINUFFT=1\n";
		return 1;
	}
	const std::vector<Case> cases = {
		{1, false, {400}, 0.3, 1.0},
		{2, true, {300, 500}, 0.25, 1.0},
		{2, false, std::vector<size_t>(5, 200), 0.4, 1.0},
		{2, false, {60, 80, 100, 120, 140, 160, 180, 200, 220, 240}, 0.3, 1.0},
		{3, false, std::vector<size_t>(5, 150), 0.2, 0.0},
		{3, true, {90, 100, 110, 120, 130, 140, 150, 160, 170, 180}, 0.2, 1.0},
		// Direct fallback of FINUFFT-selected species: all (sheared cell, sparse-grid guard; grid cap) and mixed.
		{2, false, {600}, 0.4, 0.0, 0, 100.0},
		{2, false, {60, 80, 100, 120, 140, 160, 180, 200, 220, 240}, 0.3, 1.0, 0x3ff},
		{2, true, {300, 500}, 0.25, 1.0, 0x2},
		{3, true, {90, 100, 110, 120, 130}, 0.2, 1.0, 0xa},
	};
	for (const auto & tc : cases)
		Run(tc);
	std::cout << "peak OS threads (excluding sampler, largest team 8): " << maxPeak << '\n';
	std::cout << (failures ? "FAILED: " + std::to_string(failures) + " of " : std::string("All ")) << checks << " species-concurrent checks" << (failures ? " failed\n" : " passed\n");
	return failures ? 1 : 0;
}
