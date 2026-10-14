/** \file test_finufft_forces.cpp
 * \brief Checks of the FINUFFT Type 2 collective-coordinate forces in ShiftedCCPotential::AllForce.
 *
 * Single ShiftedCCPotential, 1-3D, orthogonal and skew cells, +k and -k constraints, with and without
 * duplicated constraints (two weights on one mode), sized so the cost gate accepts Type 2: Type 2 forces
 * match Direct to 1e-10 max(1, max|F|); a cold AllForce (no prior Energy) equals the warm one and a fresh
 * object bitwise; a forced fallback (mode cap) reports Direct and equals Direct bitwise;
 * ShiftedCCPotential_S2 and LogCCPotential keep their own per-particle Force loop. Cost gate: a tiny 3D
 * system stays on the direct loop (bitwise equal to per-particle Force), and a mid-size one switches to
 * Type 2 with one thread but not with eight. MultiStealthyPotential with one heavy and two light species:
 * only the heavy one uses Type 2, and AllForce agrees with per-particle Force within the same tolerance;
 * with Direct it equals Force bitwise.
 * Build and run: make USE_FINUFFT=1 test_finufft_forces.out && ./test_finufft_forces.out */

#include "MultiStealthyPotential.h"
#include <RandomGenerator.h>
#include <etc.h>

#include <cstring>
#include <sstream>

size_t Verbosity = 0;

namespace {

using Backend = ShiftedCCPotential::RhoBackend;
using ForceBackend = ShiftedCCPotential::ForceBackend;

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

/** \brief Random configuration in a cell of volume n (unit number density), optionally skewed.
 *  \param d     Dimension.
 *  \param n     Particle number.
 *  \param skew  Use a non-orthogonal basis.
 *  \return      Configuration with n uniformly random particles. */
Configuration RandomCell(DimensionType d, size_t n, bool skew)
{
	const double L = std::pow(double(n), 1.0 / d);
	GeometryVector basis[3] = {GeometryVector(d), GeometryVector(d), GeometryVector(d)};
	for (DimensionType i = 0; i < d; i++)
		basis[i].x[i] = L;
	if (skew && d >= 2) basis[1].x[0] = 0.31 * L;
	if (skew && d == 3) { basis[2].x[0] = 0.17 * L; basis[2].x[1] = -0.23 * L; }
	Configuration c(d, basis, 0.1 * L);
	RandomGenerator gen(77 + n + d);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x(d);
		for (DimensionType a = 0; a < d; a++) x.x[a] = gen.RandomDouble();
		c.Insert("a", x);
	}
	return c;
}

/** \brief Reciprocal-lattice vectors with |k| < K(chi) at unit density, each with its negative.
 *  \param c    Configuration providing the cell.
 *  \param chi  Stealthiness parameter fixing K.
 *  \return     Wavevectors, +k block followed by the -k block. */
std::vector<GeometryVector> Modes(const Configuration & c, double chi)
{
	const DimensionType d = c.GetDimension();
	const double K = 2 * pi * std::pow(2 * d * chi / HyperSphere_Volume(d, 1.0), 1.0 / d);
	std::vector<GeometryVector> ks = GetKs(c, K, K, 1);
	const size_t n = ks.size();
	for (size_t i = 0; i < n; i++)
		ks.push_back(-1.0 * ks[i]);
	return ks;
}

/** \brief Fill a potential with weighted constraints.
 *  \param p          Potential to fill.
 *  \param ks         Wavevectors.
 *  \param duplicate  Add every wavevector a second time with weight 0.5 (same FINUFFT mode). */
void AddModes(ShiftedCCPotential & p, const std::vector<GeometryVector> & ks, bool duplicate)
{
	for (const auto & k : ks)
		p.AddConstraint(k, {1.0 + k.Modulus2()});
	if (duplicate)
		for (const auto & k : ks)
			p.AddConstraint(k, {0.5});
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
 *  \return     true if both have the same size and every component has identical bits. */
bool Same(const std::vector<GeometryVector> & a, const std::vector<GeometryVector> & b, DimensionType dim)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); i++)
		if (std::memcmp(a[i].x, b[i].x, sizeof(double) * dim) != 0) return false;
	return true;
}

/** \brief Per-particle Force of every particle, in particle order.
 *  \param p  Potential bound to a configuration.
 *  \param n  Particle number.
 *  \return   Forces from Force(f, i). */
template <class P>
std::vector<GeometryVector> PerParticle(P & p, size_t n)
{
	std::vector<GeometryVector> f(n);
	for (size_t i = 0; i < n; i++)
		p.Force(f[i], i);
	return f;
}

/** \brief Single-potential checks: Type 2 vs Direct, repeatability, fallback, derived classes.
 *  \param d          Dimension.
 *  \param skew       Non-orthogonal cell.
 *  \param duplicate  Duplicate every constraint.
 *  \param threads    ParallelNumThread. */
void Single(DimensionType d, bool skew, bool duplicate, int threads)
{
	std::ostringstream name;
	name << d << "D" << (skew ? " skew" : "") << (duplicate ? " duplicated" : "") << " T=" << threads;
	const std::string label = name.str();
	const size_t n = d == 1 ? 400 : (d == 2 ? 600 : 500);
	const Configuration c = RandomCell(d, n, skew);
	// chi large enough for a clear predicted Type 2 win at 8 threads (see ShiftedCCPotential::FinufftForcesPay)
	const std::vector<GeometryVector> ks = Modes(c, d == 1 ? 0.6 : (d == 2 ? 0.3 : 1.0));

	ShiftedCCPotential direct(d), fast(d), again(d);
	for (ShiftedCCPotential * p : {&direct, &fast, &again}) {
		AddModes(*p, ks, duplicate);
		p->ParallelNumThread = threads;
	}
	fast.SetRhoBackend(Backend::FinufftType1);
	again.SetRhoBackend(Backend::FinufftType1);
	ShiftedCCPotential capped(fast);
	capped.SetFinufftMaxModes(1);
	for (ShiftedCCPotential * p : {&direct, &fast, &again, &capped})
		p->SetConfiguration(c);

	std::vector<GeometryVector> fd, cold, warm, repeat, fb;
	Check(fast.UseFinufftForces(), label + ": cost gate declined a case chosen to use Type 2");
	fast.AllForce(cold);
	Check(fast.GetLastRhoBackendUsed() == Backend::FinufftType1 && fast.GetLastForceBackendUsed() == ForceBackend::FinufftType2,
	      label + ": cold AllForce did not use FINUFFT Type 1 rho and Type 2 forces");
	fast.Energy();
	fast.AllForce(warm);
	again.Energy();
	again.AllForce(repeat);
	Check(Same(cold, warm, d), label + ": cold AllForce differs from warm AllForce");
	Check(Same(warm, repeat, d), label + ": fresh potentials not bitwise repeatable");

	direct.AllForce(fd);
	Check(direct.GetLastForceBackendUsed() == ForceBackend::Direct, label + ": Direct backend used Type 2");
	double scale = 0;
	const double diff = ForceDiff(warm, fd, d, scale);
	Check(diff <= 1e-10 * std::max(1.0, scale), label + ": Type 2 force diff " + std::to_string(diff) + " (scale " + std::to_string(scale) + ")");

	capped.AllForce(fb);
	Check(capped.GetLastForceBackendUsed() == ForceBackend::Direct && capped.GetLastRhoBackendUsed() == Backend::Direct,
	      label + ": mode cap did not force a Direct fallback");
	Check(Same(fb, fd, d), label + ": fallback forces differ from Direct bitwise");

	// Derived classes override Force; AllForce must stay their per-particle loop (single thread inside).
	ShiftedCCPotential_S2 s2(d);
	for (const auto & k : ks)
		s2.AddConstraint(k, {1.0, 0.5});
	s2.SetRhoBackend(Backend::FinufftType1);
	s2.SetConfiguration(c);
	LogCCPotential logcc(fast);
	logcc.ParallelNumThread = 1;
	logcc.SetConfiguration(c);
	std::vector<GeometryVector> fs2, flog;
	s2.AllForce(fs2);
	logcc.AllForce(flog);
	Check(s2.GetLastForceBackendUsed() == ForceBackend::Direct && Same(fs2, PerParticle(s2, n), d),
	      label + ": ShiftedCCPotential_S2 AllForce is not its own Force loop");
	Check(logcc.GetLastForceBackendUsed() == ForceBackend::Direct && Same(flog, PerParticle(logcc, n), d),
	      label + ": LogCCPotential AllForce is not its own Force loop");
}

/** \brief MultiStealthyPotential: bulk species forces versus per-particle Force.
 *  \param d        Dimension.
 *  \param backend  Collective-density backend.
 *  \param threads  ParallelNumThread. */
void Multi(DimensionType d, Backend backend, int threads)
{
	const bool finufft = backend == Backend::FinufftType1;
	const std::string label = "multi " + std::to_string(d) + "D " + (finufft ? "FINUFFT" : "Direct") + " T=" + std::to_string(threads);
	const size_t n = 800, species = 3;
	const Configuration c = RandomCell(d, n, d >= 2);

	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	MultiStealthyPotential p(d, 1.0);
	p.ParallelNumThread = threads;
	p.SetConfiguration(c);
	p.SetNumSpecies(species);
	p.SetSecureRadii(std::vector<double>(n, std::pow(0.2 / HyperSphere_Volume(d, 1.0), 1.0 / d)));
	// Interleave species (sizes 200, 400, 200) so global and local particle indices differ.
	for (size_t i = 0; i < n; i++)
		p.SetSpeciesOfParticle(i, i % 4 == 0 ? 0 : (i % 4 == 1 ? 2 : 1));
	for (size_t s = 0; s < species; s++)
		for (const auto & k : Modes(c, s == 1 ? 0.5 : 0.02))  // species 1 heavy, 0 and 2 light
			p.AddConstraint(s, k, {1.0});
	p.SetRhoBackend(backend);
	p.SetConfiguration(c);
	std::cout.rdbuf(old);

	std::vector<GeometryVector> all;
	p.AllForce(all);  // cold: MD order, no prior Energy
	for (size_t s = 0; s < species; s++) {
		const ShiftedCCPotential & cc = *p.species_constraints[s].CCPotential;
		const bool type2 = finufft && s == 1;
		Check(cc.UseFinufftForces() == type2 && cc.GetLastForceBackendUsed() == (type2 ? ForceBackend::FinufftType2 : ForceBackend::Direct),
		      label + ": unexpected force backend of species " + std::to_string(s));
	}
	const std::vector<GeometryVector> each = PerParticle(p, n);
	if (finufft) {
		double scale = 0;
		const double diff = ForceDiff(all, each, d, scale);
		Check(diff <= 1e-10 * std::max(1.0, scale), label + ": AllForce vs Force diff " + std::to_string(diff) + " (scale " + std::to_string(scale) + ")");
	}
	else
		Check(Same(all, each, d), label + ": Direct AllForce differs from per-particle Force bitwise");
}

/** \brief Cost gate decisions: tiny system stays direct; the direct loop's thread count is weighed.
 *  \param threads  ParallelNumThread of the tiny case. */
void Gate(int threads)
{
	const std::string label = "gate T=" + std::to_string(threads);
	// Tiny 3D system: about 40 us of direct work against a Type 2 call of several hundred us.
	const Configuration tiny = RandomCell(3, 64, true);
	ShiftedCCPotential small(3);
	AddModes(small, Modes(tiny, 0.3), false);
	small.ParallelNumThread = threads;
	small.SetRhoBackend(Backend::FinufftType1);
	small.SetConfiguration(tiny);
	std::vector<GeometryVector> f;
	small.AllForce(f);
	Check(!small.UseFinufftForces() && small.GetLastForceBackendUsed() == ForceBackend::Direct, label + ": tiny case did not stay on the direct loop");
	Check(small.GetLastRhoBackendUsed() == Backend::FinufftType1, label + ": tiny case rho did not use FINUFFT");
	Check(Same(f, PerParticle(small, 64), 3), label + ": declined AllForce differs from per-particle Force bitwise");

	// Mid-size 3D system: Type 2 wins clearly against one direct thread, not against eight.
	const Configuration mid = RandomCell(3, 500, false);
	ShiftedCCPotential one(3);
	AddModes(one, Modes(mid, 0.3), false);
	one.SetRhoBackend(Backend::FinufftType1);
	ShiftedCCPotential eight(one);
	one.ParallelNumThread = 1;
	eight.ParallelNumThread = 8;
	one.SetConfiguration(mid);
	eight.SetConfiguration(mid);
	one.AllForce(f);
	Check(one.UseFinufftForces() && one.GetLastForceBackendUsed() == ForceBackend::FinufftType2, label + ": one-thread mid case did not use Type 2");
	eight.AllForce(f);
	Check(!eight.UseFinufftForces() && eight.GetLastForceBackendUsed() == ForceBackend::Direct, label + ": eight-thread mid case used Type 2");
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
	for (int threads : {4, 8}) {
		Gate(threads);
		for (DimensionType d = 1; d <= 3; d++) {
			for (bool skew : {false, true})
				if (!skew || d >= 2)
					for (bool duplicate : {false, true})
						Single(d, skew, duplicate, threads);
			Multi(d, Backend::FinufftType1, threads);
			Multi(d, Backend::Direct, threads);
		}
	}
	std::cout << (failures ? "FAILED: " + std::to_string(failures) + " check(s)\n" : std::string("All FINUFFT Type 2 force checks passed\n"));
	return failures ? 1 : 0;
}
