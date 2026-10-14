/** \file test_indexed_view.cpp
 * \brief Checks of the indexed species views that MultiStealthyPotential binds instead of sub-Configurations.
 *
 * For every species, rho(k), the species energy and every species force must be bitwise equal to those of a
 * ShiftedCCPotential bound the former way, to a sub-Configuration built as UpdateSubconfigs used to build it
 * (Insert in species idx order). With contact off, the total energy and AllForce must also be bitwise equal to
 * the sums assembled from those references. Cases: dimensions 1-3, orthogonal and skew cells, unequal species
 * sizes assigned in reverse and interleaved order, S = 1-4, Direct and FINUFFT, 4 threads; each is checked
 * cold (AllForce before Energy, the MD order) and after repeated SetConfiguration with moved particles and
 * after species reassignment. Also checks the binding policy: derived classes are rejected unless they opt in,
 * an out-of-range index is rejected without changing the object, and SetConfiguration unbinds the view.
 * Build and run: make test-indexed-view */

#include "MultiStealthyPotential.h"
#include <RandomGenerator.h>
#include <etc.h>

#include <algorithm>
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
	double chi;                 //< per-species chi
	double e0;                  //< contact strength (0 disables contact and enables the total checks)
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

/** \brief Bitwise equality of two doubles.
 *  \param a  First value.
 *  \param b  Second value.
 *  \return   true if the bits are identical. */
bool Same(double a, double b)
{
	return std::memcmp(&a, &b, sizeof(double)) == 0;
}

/** \brief Bitwise equality of two arrays of doubles.
 *  \param a  First array.
 *  \param b  Second array.
 *  \return   true if sizes and all bits are identical. */
bool Same(const std::vector<double> & a, const std::vector<double> & b)
{
	return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), sizeof(double) * a.size()) == 0);
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

/** \brief Random configuration in a cell of volume n (unit number density), optionally skewed.
 *  \param d     Dimension.
 *  \param n     Particle number.
 *  \param skew  Use a non-orthogonal basis.
 *  \param gen   Random generator.
 *  \return      The configuration. */
Configuration RandomConfig(DimensionType d, size_t n, bool skew, RandomGenerator & gen)
{
	const double L = std::pow(double(n), 1.0 / d);
	GeometryVector basis[3] = {GeometryVector(d), GeometryVector(d), GeometryVector(d)};
	for (DimensionType i = 0; i < d; i++)
		basis[i].x[i] = L;
	if (skew && d >= 2) basis[1].x[0] = 0.31 * L;
	if (skew && d == 3) { basis[2].x[0] = 0.17 * L; basis[2].x[1] = -0.23 * L; }
	Configuration c(d, basis, 0.1 * L);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x(d);
		for (DimensionType k = 0; k < d; k++) x.x[k] = gen.RandomDouble();
		c.Insert("a", x);
	}
	return c;
}

/** \brief The sub-configuration the former UpdateSubconfigs built: particles idx of c, inserted in that order.
 *  \param c    Full configuration.
 *  \param idx  Global indices of the species.
 *  \return     The sub-configuration. */
Configuration OldSubconfig(const Configuration & c, const std::vector<size_t> & idx)
{
	Configuration sub(c);
	sub.RemoveParticles();
	for (size_t i : idx)
		sub.Insert("a", c.GetRelativeCoordinates(i));
	return sub;
}

/** \brief Build a potential whose species are interleaved and assigned in reverse particle order.
 *  \param c        Configuration.
 *  \param tc       Test case.
 *  \param backend  Collective-density backend.
 *  \param gen      Random generator for the species labels.
 *  \return         The potential, bound to c. */
std::unique_ptr<MultiStealthyPotential> Make(const Configuration & c, const Case & tc, Backend backend, RandomGenerator & gen)
{
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	auto p = std::make_unique<MultiStealthyPotential>(tc.dim, tc.e0);
	p->ParallelNumThread = 4;
	p->SetConfiguration(c);
	p->SetNumSpecies(tc.sizes.size());
	const double radius = std::pow(0.2 / HyperSphere_Volume(tc.dim, 1.0), 1.0 / tc.dim);
	p->SetSecureRadii(std::vector<double>(c.NumParticle(), radius));
	std::vector<size_t> label;
	for (size_t s = 0; s < tc.sizes.size(); s++)
		label.insert(label.end(), tc.sizes[s], s);
	for (size_t i = label.size(); i > 1; i--)
		std::swap(label[i - 1], label[size_t(gen.RandomDouble() * i) % i]);
	for (size_t i = label.size(); i-- > 0;)
		p->SetSpeciesOfParticle(i, label[i]);
	for (size_t s = 0; s < tc.sizes.size(); s++) {
		const double rhoS = double(tc.sizes[s]) / c.PeriodicVolume();
		const double K = 2 * pi * std::pow(2 * tc.dim * rhoS * tc.chi / HyperSphere_Volume(tc.dim, 1.0), 1.0 / tc.dim);
		for (const auto & k : GetKs(c, K, K, 1))
			p->AddConstraint(s, k, {1.0 + 0.1 * k.Modulus2()});
	}
	p->SetRhoBackend(backend);
	p->SetConfiguration(c);
	std::cout.rdbuf(old);
	return p;
}

/** \brief Compare p on c with references bound the former way (sub-Configuration per species).
 *  \param p           Potential, bound to c (rebound here if rebind).
 *  \param c           Configuration.
 *  \param rebind      Call p.SetConfiguration(c) first; false checks the binding p already holds.
 *  \param backend     Expected backend of every species.
 *  \param forceFirst  Call AllForce before Energy (cold forces, the MD order).
 *  \param label       Message label. */
void CompareWithOldBinding(MultiStealthyPotential & p, const Configuration & c, bool rebind, Backend backend, bool forceFirst, const std::string & label)
{
	const DimensionType dim = c.GetDimension();
	const size_t n = c.NumParticle();
	if (rebind) p.SetConfiguration(c);
	std::vector<GeometryVector> total;
	double energy = 0.0;
	if (forceFirst) p.AllForce(total);
	energy = p.Energy();
	if (!forceFirst) p.AllForce(total);

	double cc = 0.0;
	std::vector<GeometryVector> ccForce(n, GeometryVector(dim));
	for (size_t s = 0; s < p.species_constraints.size(); s++) {
		const std::string tag = label + " species " + std::to_string(s);
		ShiftedCCPotential & view = *p.species_constraints[s].CCPotential;
		const std::vector<size_t> & idx = p.species_constraints[s].idx;
		const Configuration sub = OldSubconfig(c, idx);
		ShiftedCCPotential ref(view);	// same constraints, backend, thread count
		ref.SetConfiguration(sub);
		ref.GetRho();
		Check(Same(ref.RhoReal, view.RhoReal) && Same(ref.RhoImag, view.RhoImag), tag + ": rho differs from the sub-Configuration binding");
		Check(view.GetLastRhoBackendUsed() == backend, tag + ": unexpected backend used");
		const double es = ref.Energy();
		Check(Same(es, view.Energy()), tag + ": species energy differs");
		std::vector<GeometryVector> fr, fv;
		ref.AllForce(fr);
		view.AllForce(fv);
		Check(fr.size() == idx.size() && Same(fr, fv, dim), tag + ": species forces differ");
		// Assemble the totals exactly as MultiStealthyPotential does (valid when contact is off).
		cc += ((double)idx.size() / (double)n) * es;
		for (size_t j = 0; j < idx.size() && j < fr.size(); j++) {
			GeometryVector f = fr[j];
			f.MultiplyFrom(p.u0 * (double)idx.size() / (double)n);
			GeometryVector zero(dim);
			zero.MultiplyFrom(p.e0);
			ccForce[idx[j]] = zero + f;
		}
	}
	if (p.e0 == 0.0) {
		Check(Same(energy, p.u0 * cc + p.e0 * 0.0), label + ": total energy differs");
		Check(Same(total, ccForce, dim), label + ": total forces differ");
	}
}

/** \brief Run every check for one case and backend.
 *  \param tc       Test case.
 *  \param backend  Collective-density backend. */
void Run(const Case & tc, Backend backend)
{
	std::ostringstream name;
	name << tc.dim << "D" << (tc.skew ? " skew" : "") << " S=" << tc.sizes.size() << " e0=" << tc.e0
	     << (backend == Backend::Direct ? " Direct" : " FINUFFT");
	const std::string label = name.str();
	size_t n = 0;
	for (size_t s : tc.sizes) n += s;
	RandomGenerator gen(4321 + n + tc.dim);
	Configuration c = RandomConfig(tc.dim, n, tc.skew, gen);
	auto p = Make(c, tc, backend, gen);
	CompareWithOldBinding(*p, c, false, backend, true, label);

	// Repeated SetConfiguration with moved particles: the views must follow the new configuration.
	Configuration moved(c);
	for (int round = 0; round < 2; round++) {
		for (size_t i = 0; i < n; i++) {
			GeometryVector x = moved.GetRelativeCoordinates(i);
			for (DimensionType d = 0; d < tc.dim; d++) x.x[d] += 0.05 * gen.RandomDouble();
			moved.MoveParticle(i, x);
		}
		CompareWithOldBinding(*p, moved, true, backend, round == 1, label + " after move " + std::to_string(round));
	}

	// Species reassignment without SetConfiguration changes membership and local order (erase, then append).
	if (tc.sizes.size() > 1) {
		std::ostringstream sink;
		std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
		const size_t S = tc.sizes.size();
		for (size_t i : {size_t(0), n / 3, n / 2, n - 1})
			for (size_t s = 0; s < S; s++) {
				const auto & idx = p->species_constraints[s].idx;
				if (std::find(idx.begin(), idx.end(), i) != idx.end()) {
					p->SetSpeciesOfParticle(i, (s + 1) % S);
					break;
				}
			}
		std::cout.rdbuf(old);
		CompareWithOldBinding(*p, moved, false, backend, true, label + " after reassignment");
	}
}

/** \brief A derived class that opts in to the indexed view (it reads particles only through the base class). */
struct OptIn : ShiftedCCPotential
{
	/** \brief Construct an empty potential.
	 *  \param d  Dimension. */
	explicit OptIn(DimensionType d) : ShiftedCCPotential(d) {}

	/** \brief Opt in to SetIndexedConfiguration.
	 *  \return true. */
	bool SupportsIndexedView(void) const override { return true; }
};

/** \brief Binding policy: rejection of derived classes and bad indices, opt-in, and unbinding. */
void Policy()
{
	RandomGenerator gen(99);
	Configuration c = RandomConfig(2, 60, true, gen);
	std::vector<size_t> idx;
	for (size_t i = 59; i >= 5; i -= 3) idx.push_back(i);
	const Configuration sub = OldSubconfig(c, idx);
	std::vector<size_t> rest;
	for (size_t i = 0; i < 60; i += 4) rest.push_back(i);
	const Configuration sub2 = OldSubconfig(c, rest);

	ShiftedCCPotential base(2);
	for (const auto & k : GetKs(c, 4.0, 4.0, 1))
		base.AddConstraint(k, {1.0});
	ShiftedCCPotential old(base), view(base), fresh(base);
	old.SetConfiguration(sub);
	view.SetIndexedConfiguration(c, idx);
	const double eOld = old.Energy(), eView = view.Energy();
	Check(Same(eOld, eView), "policy: view energy differs from the sub-Configuration binding");

	bool threw = false;
	try { view.SetIndexedConfiguration(c, {0, 60}); } catch (const std::out_of_range &) { threw = true; }
	Check(threw, "policy: out-of-range index accepted");
	Check(view.RhoValid && Same(view.Energy(), eView), "policy: rejected binding changed the object");

	ShiftedCCPotential_S2 s2(2);
	threw = false;
	try { s2.SetIndexedConfiguration(c, idx); } catch (const std::logic_error &) { threw = true; }
	Check(threw, "policy: derived class without opt-in accepted the view");

	OptIn optIn(2);
	optIn.constraints = base.constraints;
	optIn.SetIndexedConfiguration(c, idx);
	Check(Same(optIn.Energy(), eOld), "policy: opted-in derived class differs");

	view.SetConfiguration(sub2);
	fresh.SetConfiguration(sub2);
	std::vector<GeometryVector> fv, ff;
	view.AllForce(fv);
	fresh.AllForce(ff);
	Check(Same(view.Energy(), fresh.Energy()) && Same(fv, ff, 2), "policy: SetConfiguration did not unbind the view");
}

}  // namespace

/** \brief Run the suite.
 *  \return 0 if every check passed, 1 otherwise. */
int main()
{
	const std::vector<Case> cases = {
		{1, false, {120, 300}, 0.3, 1.0},
		{1, false, {250}, 0.2, 0.0},
		{2, false, {200, 50, 400, 120}, 0.35, 0.0},
		{2, true, {150, 400, 250}, 0.25, 1.0},
		{3, false, {300, 90}, 0.2, 0.0},
		{3, true, {200, 350, 120}, 0.3, 1.0},
	};
	std::vector<Backend> backends = {Backend::Direct};
	if (ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1))
		backends.push_back(Backend::FinufftType1);
	else
		std::cout << "FINUFFT not compiled in; checking the Direct backend only\n";
	for (Backend backend : backends)
		for (const auto & tc : cases)
			Run(tc, backend);
	Policy();
	std::cout << (failures ? "FAILED: " + std::to_string(failures) + " check(s)\n" : std::string("All indexed-view checks passed\n"));
	return failures ? 1 : 0;
}
