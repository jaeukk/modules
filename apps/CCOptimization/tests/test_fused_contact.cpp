/** \file test_fused_contact.cpp
 * \brief Checks that MultiStealthyPotential::EnergyAndAllForce equals Energy() followed by AllForce() bitwise.
 *
 * For dimensions 1-3, orthogonal and skew cells, S = 1..10 species (interleaved in particle order), equal and
 * unequal radii, contact on and off, both backends, and 4 and 8 threads:
 *  - the fused call on a fresh potential (cold) and again after moving every particle (warm plans) gives the
 *    energy and every force component of the separate calls bitwise;
 *  - a separate evaluation after fused ones, a cold AllForce (MD order) and an energy-only call are unchanged;
 *  - energy and forces match a single-threaded Direct reference within 1e-10 (relative for energy, times
 *    max(1, max|F|) for forces).
 * Special cases: tiny periodic cells where a particle overlaps several images of one neighbor and its own
 * images, coincident particles (the overlap rule of ContactForce), FINUFFT species falling back to Direct,
 * and species whose Type 2 forces the gate accepts (the test requires that some evaluation used Type 2).
 * Also: the Potential default (a ShiftedCCPotential) equals Energy then AllForce bitwise, and
 * RelaxStructure_NLOPT, which now calls EnergyAndAllForce, ends at the same coordinates after the same number
 * of evaluations as with the separate calls.
 * The makefile target runs it with default waiting and OMP_WAIT_POLICY=passive.
 * Build and run: make test-fused-contact */

#include "MultiStealthyPotential.h"
#include <StructureOptimization.h>
#include <RandomGenerator.h>
#include <etc.h>

#include <cstring>
#include <limits>
#include <map>
#include <sstream>

size_t Verbosity = 0;

namespace {

using Backend = ShiftedCCPotential::RhoBackend;
using ForceBackend = ShiftedCCPotential::ForceBackend;

/** \brief One test case. */
struct Case {
	DimensionType dim;          //< spatial dimension
	bool skew;                  //< non-orthogonal cell
	std::vector<size_t> sizes;  //< particles per species
	double chi;                 //< per-species chi
	double e0;                  //< contact strength (0 disables contact)
	double phi = 0.2;           //< packing fraction at unit number density that sets the mean radius
	bool unequal = false;       //< radii spread over 0.7-1.3 times the mean
	size_t coincident = 0;      //< particle n-1-k is put exactly on particle k, for k < coincident
	unsigned fallback = 0;      //< bit s: species s falls back to Direct (FinufftMaxModes = 1)
	bool images = false;        //< require a pair in contact through two or more periodic images
};

/** \brief Energy and forces of one evaluation. */
struct Result {
	double energy;
	std::vector<GeometryVector> force;
};

int failures = 0, checks = 0, type2Evaluations = 0;

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

/** \brief Bitwise equality of two force lists.
 *  \param a    First list.
 *  \param b    Second list.
 *  \param dim  Dimension.
 *  \return     true if both have the same size and every component has identical bits. */
bool Same(const std::vector<GeometryVector> & a, const std::vector<GeometryVector> & b, DimensionType dim)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); i++)
		if (std::memcmp(a[i].x, b[i].x, sizeof(double) * dim) != 0) return false;
	return true;
}

/** \brief Bitwise equality of two evaluations.
 *  \param a    First result.
 *  \param b    Second result.
 *  \param dim  Dimension.
 *  \return     true if the energy and every force component have identical bits. */
bool Same(const Result & a, const Result & b, DimensionType dim)
{
	return std::memcmp(&a.energy, &b.energy, sizeof(double)) == 0 && Same(a.force, b.force, dim);
}

/** \brief Random configuration at unit number density, optionally skewed, with coincident particles.
 *  \param tc  Test case.
 *  \return    Configuration. */
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
	if (tc.skew && d == 3) { basis[2].x[0] = -0.19 * L; basis[2].x[1] = 0.21 * L; }
	Configuration c(d, basis, 0.1 * L);
	RandomGenerator gen(31 + n + d);
	std::vector<GeometryVector> x(n, GeometryVector(d));
	for (size_t i = 0; i < n; i++)
		for (DimensionType k = 0; k < d; k++) x[i].x[k] = gen.RandomDouble();
	for (size_t k = 0; k < tc.coincident; k++)
		x[n - 1 - k] = x[k];
	for (size_t i = 0; i < n; i++)
		c.Insert("a", x[i]);
	return c;
}

/** \brief Particle radii of a case.
 *  \param tc  Test case.
 *  \param n   Particle count.
 *  \return    Radii (equal, or spread deterministically over 0.7-1.3 times the mean). */
std::vector<double> Radii(const Case & tc, size_t n)
{
	const double mean = std::pow(tc.phi / HyperSphere_Volume(tc.dim, 1.0), 1.0 / tc.dim);
	std::vector<double> radii(n, mean);
	if (tc.unequal)
		for (size_t i = 0; i < n; i++)
			radii[i] = mean * (0.7 + 0.6 * std::fmod(0.6180339887 * double(i), 1.0));
	return radii;
}

/** \brief Build a species-resolved potential bound to c, species interleaved in particle order.
 *  \param c        Configuration.
 *  \param tc       Test case.
 *  \param backend  Collective-density backend.
 *  \param threads  ParallelNumThread.
 *  \return         The potential. */
template <class P = MultiStealthyPotential>
std::unique_ptr<P> Make(const Configuration & c, const Case & tc, Backend backend, int threads)
{
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	auto p = std::make_unique<P>(tc.dim, tc.e0);
	p->ParallelNumThread = threads;
	p->SetConfiguration(c);
	p->SetNumSpecies(tc.sizes.size());
	p->SetSecureRadii(Radii(tc, c.NumParticle()));
	std::vector<size_t> label;
	for (size_t s = 0; s < tc.sizes.size(); s++)
		label.insert(label.end(), tc.sizes[s], s);
	RandomGenerator gen(13 + label.size());
	for (size_t i = label.size(); i > 1; i--)
		std::swap(label[i - 1], label[size_t(gen.RandomDouble() * i) % i]);
	for (size_t i = 0; i < label.size(); i++)
		p->SetSpeciesOfParticle(i, label[i]);
	for (size_t s = 0; s < tc.sizes.size(); s++) {
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

/** \brief Energy() then AllForce() on c (the behavior the fused call must reproduce).
 *  \param p  Potential.
 *  \param c  Configuration to bind first.
 *  \return   Energy and forces. */
Result Separate(MultiStealthyPotential & p, const Configuration & c)
{
	Result r;
	p.SetConfiguration(c);
	r.energy = p.Energy();
	p.AllForce(r.force);
	return r;
}

/** \brief EnergyAndAllForce() on c; counts evaluations in which some species used Type 2 forces.
 *  \param p  Potential.
 *  \param c  Configuration to bind first.
 *  \return   Energy and forces. */
Result Fused(MultiStealthyPotential & p, const Configuration & c)
{
	Result r;
	p.SetConfiguration(c);
	r.energy = p.EnergyAndAllForce(r.force);
	for (const auto & sp : p.species_constraints)
		if (sp.CCPotential->GetLastForceBackendUsed() == ForceBackend::FinufftType2) {
			++type2Evaluations;
			break;
		}
	return r;
}

/** \brief Largest number of periodic images of one other particle that a particle is in contact with.
 *  \param c      Configuration.
 *  \param radii  Particle radii.
 *  \param rc     Neighbor cutoff of the potential.
 *  \return       The largest count over all ordered pairs (i, j != i). */
size_t MaxContactImages(const Configuration & c, const std::vector<double> & radii, double rc)
{
	size_t most = 0;
	for (size_t i = 0; i < c.NumParticle(); i++) {
		std::map<size_t, size_t> count;
		c.IterateThroughNeighbors(c.GetRelativeCoordinates(i), rc,
			[&](const GeometryVector & shift, const GeometryVector &, const signed long *, size_t j) {
			if (j != i && std::sqrt(shift.Modulus2()) < radii[i] + radii[j]) ++count[j];
		});
		for (const auto & kv : count) most = std::max(most, kv.second);
	}
	return most;
}

/** \brief Run every check for one case.
 *  \param tc  Test case. */
void Run(const Case & tc)
{
	std::ostringstream name;
	name << tc.dim << "D" << (tc.skew ? " skew" : "") << " S=" << tc.sizes.size() << " e0=" << tc.e0
	     << (tc.unequal ? " unequal" : "") << (tc.coincident ? " coincident" : "") << (tc.images ? " images" : "")
	     << (tc.fallback ? " fallback=" + std::to_string(tc.fallback) : std::string());
	const Configuration c = Positions(tc);
	Configuration moved(c);
	RandomGenerator gen(7 + c.NumParticle());
	for (size_t i = 0; i < c.NumParticle(); i++) {
		GeometryVector x = moved.GetRelativeCoordinates(i);
		for (DimensionType a = 0; a < tc.dim; a++) x.x[a] += 0.02 * gen.RandomDouble();
		moved.MoveParticle(i, x);
	}
	if (tc.images) {
		const auto p = Make(c, tc, Backend::Direct, 1);
		Check(MaxContactImages(c, Radii(tc, c.NumParticle()), p->GetRc()) >= 2, name.str() + ": no pair in contact through two images");
	}
	const Result ref = Separate(*Make(c, tc, Backend::Direct, 1), c);
	double scale = 0;
	for (const auto & f : ref.force)
		for (DimensionType d = 0; d < tc.dim; d++) scale = std::max(scale, std::abs(f.x[d]));

	for (Backend backend : {Backend::Direct, Backend::FinufftType1})
		for (int threads : {4, 8}) {
			const std::string label = name.str() + (backend == Backend::Direct ? " direct" : " finufft") + " T=" + std::to_string(threads);
			auto sep = Make(c, tc, backend, threads), fus = Make(c, tc, backend, threads);
			const Result a = Separate(*sep, c), b = Fused(*fus, c);
			Check(Same(a, b, tc.dim), label + ": cold fused differs from Energy + AllForce");
			const Result am = Separate(*sep, moved), bm = Fused(*fus, moved);
			Check(Same(am, bm, tc.dim), label + ": warm fused (moved particles) differs from Energy + AllForce");
			Result again;
			again.energy = fus->EnergyAndAllForce(again.force);  // same configuration, cached rho
			Check(Same(again, bm, tc.dim), label + ": repeated fused call without SetConfiguration differs");
			Check(Same(Separate(*fus, c), a, tc.dim), label + ": Energy + AllForce after fused calls changed");

			auto md = Make(c, tc, backend, threads);
			std::vector<GeometryVector> cold;
			md->SetConfiguration(c);
			md->AllForce(cold);
			const double e = md->Energy();
			Check(Same(cold, a.force, tc.dim), label + ": cold AllForce (MD order) differs");
			Check(std::memcmp(&e, &a.energy, sizeof(double)) == 0, label + ": energy after a cold AllForce differs");
			auto eo = Make(c, tc, backend, threads);
			eo->SetConfiguration(c);
			const double e1 = eo->Energy();
			Check(std::memcmp(&e1, &a.energy, sizeof(double)) == 0, label + ": energy-only call differs");

			double df = 0;
			for (size_t i = 0; i < ref.force.size(); i++)
				for (DimensionType d = 0; d < tc.dim; d++) df = std::max(df, std::abs(b.force[i].x[d] - ref.force[i].x[d]));
			Check(std::abs(b.energy - ref.energy) <= 1e-10 * std::abs(ref.energy),
			      label + ": energy rel diff " + std::to_string(std::abs(b.energy - ref.energy) / std::abs(ref.energy)));
			Check(df <= 1e-10 * std::max(1.0, scale), label + ": force diff " + std::to_string(df));
		}
}

/** \brief MultiStealthyPotential that answers EnergyAndAllForce with the Potential default (Energy, then AllForce). */
struct SeparatePotential : MultiStealthyPotential {
	using MultiStealthyPotential::MultiStealthyPotential;
	/** \brief The Potential default: Energy() followed by AllForce(forces).
	 *  \param forces  Output per-particle forces.
	 *  \return        Total potential energy. */
	double EnergyAndAllForce(std::vector<GeometryVector> & forces) override { return Potential::EnergyAndAllForce(forces); }
};

/** \brief MultiStealthyPotential counting Energy and EnergyAndAllForce calls.
 *  \tparam Base  MultiStealthyPotential or SeparatePotential. */
template <class Base>
struct Counting : Base {
	size_t energy = 0, fused = 0;  //< calls of Energy() and EnergyAndAllForce()
	using Base::Base;
	/** \brief Count and forward.
	 *  \return Total potential energy. */
	double Energy() override { ++energy; return Base::Energy(); }
	/** \brief Count and forward.
	 *  \param forces  Output per-particle forces.
	 *  \return        Total potential energy. */
	double EnergyAndAllForce(std::vector<GeometryVector> & forces) override { ++fused; return Base::EnergyAndAllForce(forces); }
};

/** \brief Fingerprint of all relative coordinates (FNV-1a over their bits).
 *  \param c  Configuration.
 *  \return   Hash value. */
std::uint64_t Hash(const Configuration & c)
{
	std::uint64_t h = 1469598103934665603ULL;
	for (size_t i = 0; i < c.NumParticle(); i++)
		for (DimensionType d = 0; d < c.GetDimension(); d++) {
			unsigned char bytes[sizeof(double)];
			std::memcpy(bytes, &c.GetRelativeCoordinates(i).x[d], sizeof(double));
			for (unsigned char byte : bytes) { h ^= byte; h *= 1099511628211ULL; }
		}
	return h;
}

/** \brief RelaxStructure_NLOPT with the fused callback ends where the separate calls end, after as many calls.
 *  \param tc       Test case.
 *  \param backend  Collective-density backend.
 *  \param threads  ParallelNumThread. */
void Optimize(const Case & tc, Backend backend, int threads)
{
	const std::string label = "optimization " + std::to_string(tc.dim) + "D S=" + std::to_string(tc.sizes.size())
		+ (backend == Backend::Direct ? " direct" : " finufft") + " T=" + std::to_string(threads);
	const Configuration start = Positions(tc);
	Configuration a(start), b(start);
	auto sep = Make<Counting<SeparatePotential>>(start, tc, backend, threads);
	auto fus = Make<Counting<MultiStealthyPotential>>(start, tc, backend, threads);
	RelaxStructure_NLOPT(a, *sep, 0.0, 0, 0.0, 400);
	RelaxStructure_NLOPT(b, *fus, 0.0, 0, 0.0, 400);
	Check(fus->fused > 0, label + ": the optimizer never called EnergyAndAllForce");
	// The default EnergyAndAllForce of sep calls Energy, which its counter sees too.
	Check(sep->fused == fus->fused && sep->energy == fus->energy + fus->fused,
	      label + ": call counts differ (" + std::to_string(sep->energy) + "+" + std::to_string(sep->fused) + " vs "
	      + std::to_string(fus->energy) + "+" + std::to_string(fus->fused) + ")");
	Check(Hash(a) == Hash(b), label + ": end points differ");
}

/** \brief The Potential default on a potential without an override equals Energy then AllForce bitwise. */
void Default()
{
	Case tc{2, true, {300}, 0.4, 0.0};
	const Configuration c = Positions(tc);
	ShiftedCCPotential p(2), q(2);
	const double K = 2 * pi * std::pow(2 * 2 * tc.chi / HyperSphere_Volume(2, 1.0), 0.5);
	for (const auto & k : GetKs(c, K, K, 1)) {
		p.AddConstraint(k, {1.0});
		q.AddConstraint(k, {1.0});
	}
	p.ParallelNumThread = q.ParallelNumThread = 4;
	p.SetConfiguration(c);
	q.SetConfiguration(c);
	Result a, b;
	a.energy = p.Energy();
	p.AllForce(a.force);
	b.energy = q.EnergyAndAllForce(b.force);
	Check(Same(a, b, 2), "Potential default: EnergyAndAllForce differs from Energy + AllForce");
}

}  // namespace

/** \brief Run the suite.
 *  \return 0 if every check passed, 1 otherwise. */
int main()
{
	if (!ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1)) {
		std::cout << "FINUFFT backend not compiled in; build with make USE_FINUFFT=1\n";
		return 1;
	}
	ProgramStart = std::time(nullptr);
	TimeLimit = std::numeric_limits<time_t>::max();
	std::vector<Case> cases = {
		{1, false, {400}, 0.3, 1.0},
		{1, false, {200, 150, 100}, 0.5, 1.0, 0.4, true, 3},
		{2, true, {300, 500}, 0.25, 1.0, 0.3, true},
		{2, false, std::vector<size_t>(5, 200), 0.4, 0.0},
		{3, false, std::vector<size_t>(5, 150), 0.2, 1.0, 0.3, true, 4},
		{3, true, {90, 100, 110, 120, 130, 140, 150, 160, 170, 180}, 0.2, 1.0},
		// species whose Type 2 forces the gate accepts (as in test_finufft_forces.cpp's combined cases)
		{2, true, std::vector<size_t>(6, 150), 0.2, 1.0, 0.2, true, 2},
		{1, false, {400, 300, 200}, 0.5, 1.0},
		{3, true, {500, 400}, 1.0, 1.0},
		// FINUFFT species falling back to Direct, all and mixed
		{2, false, {60, 80, 100, 120, 140, 160, 180, 200, 220, 240}, 0.3, 1.0, 0.2, false, 0, 0x3ff},
		{2, true, {300, 500}, 0.25, 1.0, 0.2, true, 1, 0x2},
		// tiny cells: contact through several images of a neighbor and of the particle itself
		{1, false, {3, 2}, 0.5, 1.0, 3.0, true, 1, 0, true},
		{2, true, {4, 3}, 0.5, 1.0, 3.0, true, 1, 0, true},
		{3, false, {5, 4}, 0.5, 1.0, 3.0, true, 1, 0, true},
		{3, true, {5, 4}, 0.5, 0.0, 3.0, true, 1},
	};
	for (size_t S = 1; S <= 10; S++)  // S = 1..10, alternating cell shape and contact
		cases.push_back({2, S % 2 == 0, std::vector<size_t>(S, 40 + 3 * S), 0.4, S % 3 ? 1.0 : 0.0, 0.3, true, S % 4});
	for (const auto & tc : cases)
		Run(tc);
	Check(type2Evaluations > 0, "no fused evaluation used Type 2 forces (gate constants changed?)");
	Default();
	for (Backend backend : {Backend::Direct, Backend::FinufftType1})
		Optimize({2, false, {150, 150, 100}, 0.4, 1.0}, backend, 4);
	Optimize({3, true, {120, 100}, 0.3, 1.0, 0.2, true}, Backend::FinufftType1, 8);
	std::cout << "fused evaluations with Type 2 forces: " << type2Evaluations << '\n';
	std::cout << (failures ? "FAILED: " + std::to_string(failures) + " of " : std::string("All ")) << checks << " fused contact checks" << (failures ? " failed\n" : " passed\n");
	return failures ? 1 : 0;
}
