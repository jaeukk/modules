/** \file check_gate.cpp
 * \brief End-to-end check of the Type 2 force gate in MultiStealthyPotential::AllForce.
 *
 * For each case (S equal species in a cubic cell of unit number density, random positions, every reciprocal-lattice
 * vector with |k| < K(chi) per species at the species density, contact on, the FINUFFT backend) it times a cold
 * AllForce (after SetConfiguration, as the optimizer calls it) of the potential as built against the same potential
 * with Type 2 disabled (species of a derived class, which the gate never accepts), alternating the two over 60
 * evaluations and keeping the medians of the last 50. It prints the number of species the gate sent to Type 2
 * (jointly, as AllForce decides) and how many of them would be accepted alone (UseFinufftForces), the two median
 * times and their ratio. A sound calibration gives ratio < 1 whenever a species is accepted, and about 1 otherwise.
 * Results are machine-specific and noisy on a loaded machine; output goes to stdout only (see README.md).
 * Usage (from apps/CCOptimization, after make finufft-gate):
 *   ./check_gate.out                    the default list at 4 and 8 threads
 *   ./check_gate.out DIM N S CHI T      one case */

#ifndef CCPOTENTIAL_USE_FINUFFT
#error "check_gate.cpp needs the FINUFFT backend: build it with make finufft-gate"
#endif

#include "MultiStealthyPotential.h"
#include <RandomGenerator.h>
#include <etc.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <sstream>
#include <string>

size_t Verbosity = 0;

namespace {

using Clock = std::chrono::steady_clock;

/** \brief Species potential of a derived class: FINUFFT rho as usual, but the gate never accepts Type 2 forces. */
struct NoType2 : ShiftedCCPotential
{
	/** \brief Construct an empty species potential.
	 *  \param d  Dimension. */
	explicit NoType2(DimensionType d) : ShiftedCCPotential(d) {}

	/** \brief Allow MultiStealthyPotential's indexed species views (particles are read only through the base class).
	 *  \return true. */
	bool SupportsIndexedView() const override { return true; }
};

/** \brief Build a species-resolved FINUFFT potential bound to c.
 *  \param c        Configuration (species are contiguous blocks of N/S particles).
 *  \param S        Number of species.
 *  \param chi      Per-species chi at the species density.
 *  \param threads  ParallelNumThread.
 *  \param noType2  Use NoType2 species.
 *  \return         The potential. */
std::unique_ptr<MultiStealthyPotential> Make(const Configuration & c, size_t S, double chi, int threads, bool noType2)
{
	const DimensionType d = c.GetDimension();
	const size_t n = c.NumParticle(), ns = n / S;
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	auto p = std::make_unique<MultiStealthyPotential>(d, 1.0);
	p->ParallelNumThread = threads;
	p->SetConfiguration(c);
	p->SetNumSpecies(S);
	if (noType2)
		for (auto & sp : p->species_constraints)
			sp.CCPotential = std::make_unique<NoType2>(d);
	p->SetSecureRadii(std::vector<double>(n, std::pow(0.2 / HyperSphere_Volume(d, 1.0), 1.0 / d)));
	for (size_t i = 0; i < n; i++)
		p->SetSpeciesOfParticle(i, std::min(i / ns, S - 1));
	const double K = 2 * pi * std::pow(2 * d * (double(ns) / c.PeriodicVolume()) * chi / HyperSphere_Volume(d, 1.0), 1.0 / d);
	for (size_t s = 0; s < S; s++)
		for (const auto & k : GetKs(c, K, K, 1))
			p->AddConstraint(s, k, {1.0});
	p->SetRhoBackend(ShiftedCCPotential::RhoBackend::FinufftType1);
	p->SetConfiguration(c);
	std::cout.rdbuf(old);
	return p;
}

/** \brief Time one case and print its line.
 *  \param d        Dimension (1-3).
 *  \param n        Particle number.
 *  \param S        Number of species (n >= 2 S).
 *  \param chi      Per-species chi.
 *  \param threads  Team size. */
void Case(int d, size_t n, size_t S, double chi, int threads)
{
	const double L = std::pow(double(n), 1.0 / d);
	GeometryVector basis[3] = {GeometryVector(d), GeometryVector(d), GeometryVector(d)};
	for (int i = 0; i < d; i++)
		basis[i].x[i] = L;
	Configuration c(d, basis, 0.1 * L);
	RandomGenerator gen(11 + n);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x(d);
		for (int k = 0; k < d; k++) x.x[k] = gen.RandomDouble();
		c.Insert("a", x);
	}
	auto gated = Make(c, S, chi, threads, false), plain = Make(c, S, chi, threads, true);
	std::vector<GeometryVector> f;
	std::vector<double> tg, tp;
	RandomGenerator jitter(5);
	for (int r = 0; r < 60; r++) {
		Configuration moved(c);
		for (size_t i = 0; i < n; i++) {
			GeometryVector x = moved.GetRelativeCoordinates(i);
			x.x[0] += 1e-6 * jitter.RandomDouble();
			moved.MoveParticle(i, x);
		}
		for (int which = 0; which < 2; which++) {	// alternate the order so drifts in load hit both
			MultiStealthyPotential & p = (which ^ (r & 1)) ? *plain : *gated;
			p.SetConfiguration(moved);
			const auto t0 = Clock::now();
			p.AllForce(f);
			const double s = std::chrono::duration<double>(Clock::now() - t0).count();
			if (r >= 10)
				(&p == gated.get() ? tg : tp).push_back(s);
		}
	}
	std::sort(tg.begin(), tg.end());
	std::sort(tp.begin(), tp.end());
	size_t joint = 0, alone = 0;
	for (const auto & sp : gated->species_constraints) {
		joint += sp.CCPotential->GetLastForceBackendUsed() == ShiftedCCPotential::ForceBackend::FinufftType2;
		alone += sp.CCPotential->UseFinufftForces();
	}
	const double mg = tg[tg.size() / 2], mp = tp[tp.size() / 2];
	std::printf("d=%d N=%zu S=%zu chi=%.2f T=%d M_s=%zu type2_species=%zu alone=%zu gated_ms=%.3f notype2_ms=%.3f ratio=%.2f\n",
	            d, n, S, chi, threads, gated->species_constraints[0].CCPotential->constraints.size(), joint, alone,
	            1e3 * mg, 1e3 * mp, mg / mp);
	std::fflush(stdout);
}

}  // namespace

/** \brief Run the default case list or one case from the command line.
 *  \param argc  Argument count.
 *  \param argv  Nothing, or DIM N S CHI T.
 *  \return      0 on success, 1 on bad arguments. */
int main(int argc, char ** argv)
{
	if (argc == 6) {
		const int d = std::atoi(argv[1]), threads = std::atoi(argv[5]);
		const size_t n = std::stoul(argv[2]), S = std::stoul(argv[3]);
		if (d < 1 || d > 3 || S < 1 || n < 2 * S || threads < 1) {
			std::fprintf(stderr, "need DIM in 1-3, S >= 1, N >= 2 S, T >= 1\n");
			return 1;
		}
		Case(d, n, S, std::atof(argv[4]), threads);
		return 0;
	}
	if (argc != 1) {
		std::fprintf(stderr, "usage: check_gate.out [DIM N S CHI T]\n");
		return 1;
	}
	struct Row { int d; size_t n, S; double chi; };
	const Row rows[] = {{2, 2000, 5, 0.1}, {2, 1000, 5, 0.2}, {2, 500, 5, 0.4}, {2, 300, 5, 0.4}, {2, 200, 5, 0.4},
	                    {2, 2000, 10, 0.2}, {2, 1000, 10, 0.1}, {2, 1200, 3, 0.3}, {3, 2000, 5, 0.4}, {3, 2000, 1, 0.4},
	                    {3, 1500, 3, 0.8}, {1, 2000, 5, 0.3}, {1, 500, 5, 0.3}};
	for (int threads : {4, 8})
		for (const Row & r : rows)
			Case(r.d, r.n, r.S, r.chi, threads);
	return 0;
}
