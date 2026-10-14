/** \file test_finufft_multi.cpp
 * \brief Shared correctness oracle for FINUFFT acceleration of MultiStealthyPotential.
 *
 * Compares the FINUFFT backend (selected with MultiStealthyPotential::SetRhoBackend) against the
 * direct sums: energy to 1e-10 relative, every force component to 1e-10 times max(1, max|F|).
 * Cases cover dimensions 1-3, orthogonal and skew cells, unequal species sizes and mode sets,
 * k-dependent weights, contact on/off, and 4 and 8 threads. Each case also checks: a cold AllForce
 * (no prior Energy, the MD call order) equals the warm one bitwise; two fresh potentials give
 * bitwise-identical energies and forces (same-backend repeatability); results stay correct after
 * moving every particle and after a species reassignment round trip; and a lattice state with
 * rho(k) = 0 on all constraints has FINUFFT energy and forces within 1e-10 (absolute) of Direct.
 * Build and run: make USE_FINUFFT=1 test_finufft_multi.out && ./test_finufft_multi.out */

#include "MultiStealthyPotential.h"
#include <RandomGenerator.h>
#include <etc.h>

#include <cstring>
#include <sstream>
#include <stdexcept>

size_t Verbosity = 0;

namespace {

using Backend = ShiftedCCPotential::RhoBackend;

/** \brief One test case. */
struct Case {
	DimensionType dim;          //< spatial dimension
	bool skew;                  //< non-orthogonal cell
	std::vector<size_t> sizes;  //< particles per species
	double chi;                 //< base per-species chi; species s uses chi*(1+0.3 s)
	bool weighted;              //< V(k) = 1 + |k|^2 instead of flat
	double e0;                  //< contact strength (0 disables contact)
};

int failures = 0;

/** \brief Record a failed check without stopping the suite.
 *  \param ok    Check result.
 *  \param what  Description printed on failure. */
void Check(bool ok, const std::string & what)
{
	if (!ok) {
		++failures;
		std::cout << "FAIL: " << what << '\n';
	}
}

/** \brief Cell of volume N (unit number density), optionally skewed.
 *  \param d     Dimension.
 *  \param n     Particle number.
 *  \param skew  Use a non-orthogonal basis.
 *  \return      Empty configuration. */
Configuration Cell(DimensionType d, size_t n, bool skew)
{
	const double L = std::pow(double(n), 1.0 / d);
	GeometryVector basis[3] = {GeometryVector(d), GeometryVector(d), GeometryVector(d)};
	for (DimensionType i = 0; i < d; i++)
		basis[i].x[i] = L;
	if (skew && d >= 2) basis[1].x[0] = 0.31 * L;
	if (skew && d == 3) { basis[2].x[0] = 0.17 * L; basis[2].x[1] = -0.23 * L; }
	return Configuration(d, basis, 0.1 * L);
}

/** \brief Build a species-resolved potential on a configuration.
 *  \param c        Configuration (species blocks are contiguous in particle order).
 *  \param tc       Test case.
 *  \param backend  Collective-density backend.
 *  \param threads  ParallelNumThread.
 *  \param modes    Constraint wavevectors per species (empty: generate from chi).
 *  \return         The potential, bound to c. */
std::unique_ptr<MultiStealthyPotential> Make(const Configuration & c, const Case & tc, Backend backend, int threads,
                                             std::vector<std::vector<GeometryVector>> & modes)
{
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	auto p = std::make_unique<MultiStealthyPotential>(tc.dim, tc.e0);
	p->ParallelNumThread = threads;
	p->SetConfiguration(c);
	p->SetNumSpecies(tc.sizes.size());
	const double radius = std::pow(0.2 / HyperSphere_Volume(tc.dim, 1.0), 1.0 / tc.dim);
	p->SetSecureRadii(std::vector<double>(c.NumParticle(), radius));
	size_t first = 0;
	const bool generate = modes.empty();
	for (size_t s = 0; s < tc.sizes.size(); s++) {
		for (size_t i = 0; i < tc.sizes[s]; i++)
			p->SetSpeciesOfParticle(first + i, s);
		first += tc.sizes[s];
		if (generate) {
			const double rhoS = double(tc.sizes[s]) / c.PeriodicVolume(), chi = tc.chi * (1 + 0.3 * s);
			const double K = 2 * pi * std::pow(2 * tc.dim * rhoS * chi / HyperSphere_Volume(tc.dim, 1.0), 1.0 / tc.dim);
			modes.push_back(GetKs(c, K, K, 1));
		}
		for (const auto & k : modes[s])
			p->AddConstraint(s, k, {tc.weighted ? 1.0 + k.Modulus2() : 1.0});
	}
	p->SetRhoBackend(backend);
	p->SetConfiguration(c);
	std::cout.rdbuf(old);
	return p;
}

/** \brief Largest absolute force-component difference and largest reference component.
 *  \param a      Forces to test.
 *  \param b      Reference forces.
 *  \param dim    Dimension.
 *  \param[out] scale  max |b| over all components.
 *  \return       max |a - b|. */
double ForceDiff(const std::vector<GeometryVector> & a, const std::vector<GeometryVector> & b, DimensionType dim, double & scale)
{
	double diff = 0;
	scale = 0;
	for (size_t i = 0; i < b.size(); i++)
		for (DimensionType d = 0; d < dim; d++) {
			diff = std::max(diff, std::abs(a[i].x[d] - b[i].x[d]));
			scale = std::max(scale, std::abs(b[i].x[d]));
		}
	return diff;
}

/** \brief Bitwise equality of two force lists.
 *  \param a    First list.
 *  \param b    Second list.
 *  \param dim  Dimension.
 *  \return     true if every component has identical bits. */
bool Same(const std::vector<GeometryVector> & a, const std::vector<GeometryVector> & b, DimensionType dim)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); i++)
		if (std::memcmp(a[i].x, b[i].x, sizeof(double) * dim) != 0) return false;
	return true;
}

/** \brief Compare FINUFFT with Direct on configuration c.
 *  \param fast    FINUFFT potential (rebound to c here).
 *  \param direct  Direct potential (rebound to c here).
 *  \param c       Configuration.
 *  \param label   Case label for messages.
 *  \param absTol  If > 0, use this absolute tolerance for energy and forces (near-zero states). */
void Compare(MultiStealthyPotential & fast, MultiStealthyPotential & direct, const Configuration & c,
             const std::string & label, double absTol = 0.0)
{
	fast.SetConfiguration(c);
	direct.SetConfiguration(c);
	const double ef = fast.Energy(), ed = direct.Energy();
	std::vector<GeometryVector> ff, fd;
	fast.AllForce(ff);
	direct.AllForce(fd);
	double scale = 0;
	const double df = ForceDiff(ff, fd, c.GetDimension(), scale);
	if (absTol > 0) {
		Check(std::abs(ef - ed) <= absTol, label + ": near-zero energy differs " + std::to_string(ef - ed));
		Check(df <= absTol, label + ": near-zero force differs " + std::to_string(df));
	}
	else {
		Check(std::abs(ef - ed) <= 1e-10 * std::abs(ed), label + ": energy rel diff " + std::to_string(std::abs(ef - ed) / std::abs(ed)));
		Check(df <= 1e-10 * std::max(1.0, scale), label + ": force diff " + std::to_string(df) + " (scale " + std::to_string(scale) + ")");
	}
	for (auto & sp : fast.species_constraints)
		Check(sp.CCPotential->GetLastRhoBackendUsed() == Backend::FinufftType1, label + ": FINUFFT fell back to Direct");
}

/** \brief Run every check for one case at one thread count.
 *  \param tc       Test case.
 *  \param threads  ParallelNumThread. */
void Run(const Case & tc, int threads)
{
	std::ostringstream name;
	name << tc.dim << "D" << (tc.skew ? " skew" : "") << " S=" << tc.sizes.size() << " chi=" << tc.chi
	     << (tc.weighted ? " weighted" : "") << " e0=" << tc.e0 << " T=" << threads;
	const std::string label = name.str();
	size_t n = 0;
	for (size_t s : tc.sizes) n += s;
	Configuration c = Cell(tc.dim, n, tc.skew);
	RandomGenerator gen(1234 + n + tc.dim);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x(tc.dim);
		for (DimensionType d = 0; d < tc.dim; d++) x.x[d] = gen.RandomDouble();
		c.Insert("a", x);
	}
	std::vector<std::vector<GeometryVector>> modes;
	auto direct = Make(c, tc, Backend::Direct, threads, modes);
	auto fast = Make(c, tc, Backend::FinufftType1, threads, modes);
	auto again = Make(c, tc, Backend::FinufftType1, threads, modes);

	// Cold AllForce (MD order) must equal the warm result bitwise; fresh potentials must repeat bitwise.
	std::vector<GeometryVector> cold, warm, repeat;
	fast->AllForce(cold);
	const double e1 = fast->Energy();
	fast->AllForce(warm);
	const double e2 = again->Energy();
	again->AllForce(repeat);
	Check(Same(cold, warm, tc.dim), label + ": cold AllForce differs from warm AllForce");
	Check(std::memcmp(&e1, &e2, sizeof(double)) == 0 && Same(warm, repeat, tc.dim), label + ": fresh potentials not bitwise repeatable");

	Compare(*fast, *direct, c, label);

	// Move every particle: caches must be invalidated.
	Configuration moved(c);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x = moved.GetRelativeCoordinates(i);
		for (DimensionType d = 0; d < tc.dim; d++) x.x[d] += 0.05 * gen.RandomDouble();
		moved.MoveParticle(i, x);
	}
	Compare(*fast, *direct, moved, label + " after move");

	// Species reassignment round trip changes the internal particle order, not the physics.
	if (tc.sizes.size() > 1) {
		std::ostringstream sink;
		std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
		fast->SetSpeciesOfParticle(0, 1);
		fast->SetSpeciesOfParticle(0, 0);
		std::cout.rdbuf(old);
		Compare(*fast, *direct, moved, label + " after reassignment");
	}
}

/** \brief Lattice state with rho(k) = 0 on every constraint: FINUFFT must match Direct absolutely.
 *  \param dim      Dimension.
 *  \param threads  ParallelNumThread. */
void NearZero(DimensionType dim, int threads)
{
	const size_t side = dim == 3 ? 6 : (dim == 2 ? 14 : 150), species = 2;
	size_t ns = 1;
	for (DimensionType d = 0; d < dim; d++) ns *= side;
	Case tc{dim, false, std::vector<size_t>(species, ns), 0.0, true, 0.0};
	Configuration c = Cell(dim, species * ns, false);
	for (size_t s = 0; s < species; s++)
		for (size_t j = 0; j < ns; j++) {
			GeometryVector x(dim);
			size_t r = j;
			for (DimensionType d = 0; d < dim; d++) { x.x[d] = (double(r % side) + 0.37 * s + 0.11) / side; r /= side; }
			c.Insert("a", x);
		}
	// Each species is a simple lattice of `side` points per axis, so rho(k) = 0 for every reciprocal
	// vector of the cell with all integer indices |n_d| < side; constrain all with |k| < 2 pi (side-1)/L.
	const double L = std::pow(double(species * ns), 1.0 / dim), K = 2 * pi * (side - 1) / L;
	std::vector<std::vector<GeometryVector>> modes(species, GetKs(c, K, K, 1));
	auto direct = Make(c, tc, Backend::Direct, threads, modes);
	auto fast = Make(c, tc, Backend::FinufftType1, threads, modes);
	Compare(*fast, *direct, c, std::to_string(dim) + "D lattice rho=0 T=" + std::to_string(threads), 1e-10);
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
	const std::vector<Case> cases = {
		{1, false, {300, 500}, 0.3, false, 1.0},
		{2, false, {400, 400, 400, 400, 400}, 0.4, false, 1.0},
		{2, true, {150, 400, 250}, 0.25, true, 1.0},
		{2, false, {500}, 0.1, false, 0.0},
		{3, false, {300, 300}, 0.2, false, 1.0},
		{3, true, {200, 350, 120}, 0.35, true, 1.0},
	};
	for (int threads : {4, 8}) {
		for (const auto & tc : cases)
			Run(tc, threads);
		for (DimensionType d = 1; d <= 3; d++)
			NearZero(d, threads);
	}
	std::cout << (failures ? "FAILED: " + std::to_string(failures) + " check(s)\n" : std::string("All FINUFFT multi-species checks passed\n"));
	return failures ? 1 : 0;
}
