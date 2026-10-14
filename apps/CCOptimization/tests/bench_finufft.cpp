/** \file bench_finufft.cpp
 * \brief Benchmark of MultiStealthyPotential backends on the GetMultiCCO ground-state setup.
 *
 * Usage: bench_finufft.out eval|opt DIM N SPECIES CHI THREADS SEED [direct|finufft [HISTORY]]
 *   eval  per-evaluation timings (Configuration build, SetConfiguration, Energy, AllForce) for
 *         Direct and, when compiled in, FINUFFT, plus energy/force agreement with Direct.
 *   opt   one RelaxStructure_NLOPT run (as CollectiveCoordinateMultiRun) with the given backend;
 *         prints wall time, evaluation count, final energy (also re-evaluated by a fresh Direct
 *         potential) and a hash of the final coordinates for repeatability checks.
 * The setup mirrors GetMultiCCO with "random" initial packings: unit number density, equal radii
 * at packing fraction 0.2, contact strength 1, contiguous equal-size species, flat v(k), and every
 * reciprocal-lattice vector with |k| < K(CHI) per species (CHI is per species). One JSON line is
 * printed per result. Build: make USE_FINUFFT=1 bench_finufft.out (or without USE_FINUFFT for
 * Direct only). Runner: tests/bench_finufft.py. */

#include "MultiStealthyPotential.h"
#include <StructureOptimization.h>
#include <RandomGenerator.h>
#include <etc.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

size_t Verbosity = 0;

namespace {

using Backend = ShiftedCCPotential::RhoBackend;
using Clock = std::chrono::steady_clock;

/** \brief Benchmark parameters read from the command line. */
struct Setup {
	DimensionType dim;  //< spatial dimension
	size_t n;           //< total particle number
	size_t species;     //< number of equal-size species
	double chi;         //< per-species stealthiness target
	int threads;        //< OpenMP threads (ParallelNumThread)
	int seed;           //< random seed of the initial packing
};

/** \brief Seconds elapsed since a time point.
 *  \param t  Start time.
 *  \return   Elapsed seconds. */
double Since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

/** \brief MultiStealthyPotential that counts energy evaluations (one per optimizer callback). */
struct CountingPotential : MultiStealthyPotential {
	size_t energyCalls = 0;  //< number of Energy() calls
	using MultiStealthyPotential::MultiStealthyPotential;
	/** \brief Count and forward to MultiStealthyPotential::Energy.
	 *  \return Total potential energy. */
	double Energy() override { ++energyCalls; return MultiStealthyPotential::Energy(); }
};

/** \brief Random packing as in GetMultiCCO ("random"): cubic box of unit number density.
 *  \param s  Benchmark parameters.
 *  \return   Packing with equal radii at packing fraction 0.2. */
SpherePacking RandomPacking(const Setup & s)
{
	RandomGenerator gen(s.seed);
	const double radius = std::pow(0.2 / HyperSphere_Volume(s.dim, 1.0), 1.0 / s.dim);
	SpherePacking packing(GetUnitCubicBox(s.dim, 0.1), 0);
	packing.Rescale(std::pow(double(s.n), 1.0 / s.dim));
	for (size_t i = 0; i < s.n; i++) {
		GeometryVector x;
		RandomVectorInCube(s.dim, gen, x);
		packing.Insert(radius, x);
	}
	return packing;
}

/** \brief Build the potential exactly as GetMultiCCO does for flat stealthy species.
 *  \param s        Benchmark parameters.
 *  \param packing  Initial packing.
 *  \param backend  Collective-density backend.
 *  \param[out] modes  Number of constrained wavevectors per species.
 *  \return         The configured potential. */
std::unique_ptr<CountingPotential> MakePotential(const Setup & s, const SpherePacking & packing, Backend backend, size_t & modes)
{
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	auto p = std::make_unique<CountingPotential>(s.dim, 1.0);
	p->ParallelNumThread = s.threads;
	p->SetPacking(packing);
	p->SetNumSpecies(s.species);
	const size_t ns = s.n / s.species;
	for (size_t i = 0; i < s.n; i++)
		p->SetSpeciesOfParticle(i, std::min(i / ns, s.species - 1));
	// K for a per-species chi at the species density (getK in main.cpp).
	const double rhoS = double(ns) / packing.PeriodicVolume();
	const double K = 2 * pi * std::pow(2 * s.dim * rhoS * s.chi / HyperSphere_Volume(s.dim, 1.0), 1.0 / s.dim);
	std::vector<GeometryVector> ks = GetKs(packing, K, K, 1);
	modes = ks.size();
	for (size_t sp = 0; sp < s.species; sp++)
		for (const auto & k : ks)
			p->AddConstraint(sp, k, {1.0});
	p->SetRhoBackend(backend);
	p->SetPacking(packing);
	std::cout.rdbuf(old);
	return p;
}

/** \brief 64-bit FNV-1a hash of the bit patterns of all relative coordinates.
 *  \param c  Configuration to hash.
 *  \return   Hash value. */
std::uint64_t Hash(const Configuration & c)
{
	std::uint64_t h = 1469598103934665603ULL;
	for (size_t i = 0; i < c.NumParticle(); i++) {
		const GeometryVector & x = c.GetRelativeCoordinates(i);
		for (DimensionType d = 0; d < c.GetDimension(); d++) {
			unsigned char bytes[sizeof(double)];
			std::memcpy(bytes, &x.x[d], sizeof(double));
			for (unsigned char b : bytes) { h ^= b; h *= 1099511628211ULL; }
		}
	}
	return h;
}

/** \brief Name of a backend for the JSON output.
 *  \param b  Backend.
 *  \return   "direct" or "finufft". */
const char * Name(Backend b) { return b == Backend::Direct ? "direct" : "finufft"; }

/** \brief Time the parts of one optimizer evaluation and compare energies/forces with Direct.
 *  \param s  Benchmark parameters. */
void Eval(const Setup & s)
{
	const SpherePacking packing = RandomPacking(s);
	const Configuration start(packing, "a");
	size_t modes = 0;
	auto direct = MakePotential(s, packing, Backend::Direct, modes);
	std::vector<Backend> backends{Backend::Direct};
	if (ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1))
		backends.push_back(Backend::FinufftType1);
	RandomGenerator jitter(s.seed + 1);
	for (Backend b : backends) {
		auto p = MakePotential(s, packing, b, modes);
		const int reps = 100;
		std::vector<double> total;
		double tBuild = 0, tSet = 0, tEnergy = 0, tForce = 0, dE = 0, dF = 0, fMax = 0;
		std::vector<GeometryVector> f, fd;
		for (int r = 0; r < reps; r++) {
			auto t0 = Clock::now();
			Configuration moved(start);  // the optimizer builds a fresh Configuration per callback
			for (size_t i = 0; i < s.n; i++) {
				GeometryVector x = moved.GetRelativeCoordinates(i);
				x.x[0] += 1e-6 * jitter.RandomDouble();
				moved.MoveParticle(i, x);
			}
			const double build = Since(t0);
			t0 = Clock::now(); p->SetConfiguration(moved); const double set = Since(t0);
			t0 = Clock::now(); const double e = p->Energy(); const double energy = Since(t0);
			t0 = Clock::now(); p->AllForce(f); const double force = Since(t0);
			tBuild += build; tSet += set; tEnergy += energy; tForce += force;
			total.push_back(build + set + energy + force);
			if (r == 0) {
				direct->SetConfiguration(moved);
				const double ed = direct->Energy();
				direct->AllForce(fd);
				dE = std::abs(e - ed) / std::max(std::abs(ed), 1e-300);
				for (size_t i = 0; i < s.n; i++)
					for (DimensionType d = 0; d < s.dim; d++) {
						dF = std::max(dF, std::abs(f[i].x[d] - fd[i].x[d]));
						fMax = std::max(fMax, std::abs(fd[i].x[d]));
					}
			}
		}
		std::sort(total.begin(), total.end());
		const double ms = 1e3 / reps;
		std::cout << std::setprecision(4) << "{\"mode\":\"eval\",\"backend\":\"" << Name(b) << "\",\"dim\":" << s.dim
			<< ",\"N\":" << s.n << ",\"species\":" << s.species << ",\"chi\":" << s.chi << ",\"modes\":" << modes
			<< ",\"threads\":" << s.threads << ",\"build_ms\":" << tBuild * ms << ",\"set_ms\":" << tSet * ms
			<< ",\"energy_ms\":" << tEnergy * ms << ",\"allforce_ms\":" << tForce * ms
			<< ",\"median_total_ms\":" << 1e3 * total[total.size() / 2]
			<< std::setprecision(3) << ",\"energy_rel_diff\":" << dE << ",\"force_rel_diff\":" << dF / std::max(fMax, 1e-300) << "}\n";
	}
}

/** \brief One ground-state optimization as in CollectiveCoordinateMultiRun.
 *  \param s  Benchmark parameters.
 *  \param b  Backend. */
void Opt(const Setup & s, Backend b)
{
	const SpherePacking packing = RandomPacking(s);
	Configuration c(packing, "a");
	size_t modes = 0;
	auto p = MakePotential(s, packing, b, modes);
	p->SetConfiguration(c);
	p->energyCalls = 0;
	const auto t0 = Clock::now();
	RelaxStructure_NLOPT(c, *p, 0.0, 0, 0.0, 10000);
	const double seconds = Since(t0);
	const size_t calls = p->energyCalls;
	p->SetConfiguration(c);
	const double e = p->Energy();
	auto direct = MakePotential(s, packing, Backend::Direct, modes);
	direct->SetConfiguration(c);
	const double ed = direct->Energy();
	std::cout << std::setprecision(6) << "{\"mode\":\"opt\",\"backend\":\"" << Name(b) << "\",\"dim\":" << s.dim
		<< ",\"N\":" << s.n << ",\"species\":" << s.species << ",\"chi\":" << s.chi << ",\"modes\":" << modes
		<< ",\"threads\":" << s.threads << ",\"seed\":" << s.seed << ",\"history\":" << StructureOptimization_LBFGSHistory
		<< ",\"seconds\":" << seconds
		<< ",\"energy_calls\":" << calls << ",\"ms_per_call\":" << 1e3 * seconds / std::max<size_t>(calls, 1)
		<< std::setprecision(17) << ",\"final_energy\":" << e << ",\"final_energy_direct\":" << ed
		<< ",\"hash\":\"" << std::hex << Hash(c) << std::dec << "\"}\n";
}

}  // namespace

/** \brief Entry point; see the file comment for usage.
 *  \param argc  Argument count.
 *  \param argv  Arguments: eval|opt DIM N SPECIES CHI THREADS SEED [direct|finufft].
 *  \return      0 on success, 1 on error. */
int main(int argc, char ** argv)
try {
	if (argc < 8)
		throw std::invalid_argument("usage: bench_finufft.out eval|opt DIM N SPECIES CHI THREADS SEED [direct|finufft]");
	const std::string mode = argv[1];
	Setup s{static_cast<DimensionType>(std::stoi(argv[2])), std::stoul(argv[3]), std::stoul(argv[4]),
	        std::stod(argv[5]), std::stoi(argv[6]), std::stoi(argv[7])};
	if (s.species < 1 || s.n < 2 * s.species || s.threads < 1)
		throw std::invalid_argument("need SPECIES >= 1, N >= 2*SPECIES, THREADS >= 1");
	omp_set_num_threads(s.threads);
	ProgramStart = std::time(nullptr);
	TimeLimit = std::numeric_limits<time_t>::max();
	if (mode == "eval")
		Eval(s);
	else if (mode == "opt") {
		const std::string backend = argc > 8 ? argv[8] : "direct";
		if (argc > 9)
			StructureOptimization_LBFGSHistory = std::stoul(argv[9]);  // L-BFGS correction pairs (default 100)
		if (backend != "direct" && backend != "finufft")
			throw std::invalid_argument("backend must be direct or finufft");
		Opt(s, backend == "direct" ? Backend::Direct : Backend::FinufftType1);
	}
	else
		throw std::invalid_argument("mode must be eval or opt");
	return 0;
}
catch (const std::exception & e) {
	std::cerr << e.what() << '\n';
	return 1;
}
