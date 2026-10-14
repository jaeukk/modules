/** \file test_finufft_forces.cpp
 * \brief Checks of the FINUFFT Type 2 collective-coordinate forces in ShiftedCCPotential::AllForce.
 *
 * Single ShiftedCCPotential, 1-3D, orthogonal and skew cells, +k and -k constraints, with and without
 * duplicated constraints (two weights on one mode), sized so the cost gate accepts Type 2: Type 2 forces
 * match Direct to 1e-10 max(1, max|F|); a cold AllForce (no prior Energy) equals the warm one and a fresh
 * object bitwise; a forced fallback (mode cap) reports Direct and equals Direct bitwise;
 * ShiftedCCPotential_S2 and LogCCPotential keep their own per-particle Force loop. Cost gate: a tiny 3D
 * system stays on the direct loop (bitwise equal to per-particle Force); a 3D system whose Type 2 work is too small to
 * split switches to Type 2 with one thread but not with eight; a mid-size one runs Type 2 as one block with one
 * thread and as several blocks, in parallel, with eight, with bitwise equal forces. MultiStealthyPotential with one heavy and two light species:
 * the heavy one uses Type 2, every species uses the backend of the joint gate (FinufftForceGate), and AllForce
 * agrees with per-particle Force within the same tolerance; with Direct it equals Force bitwise. Combined check of
 * the integrated design (see Combined): indexed views, persistent plans, species-concurrent Type 1 and Type 2, and
 * the joint gate match Direct, choose the expected paths, and repeat bitwise at 4 and 8 threads.
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
	// chi large enough for a clear predicted Type 2 win at 8 threads (see ShiftedCCPotential::FinufftForceGate)
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
	// Species run concurrently, so the gate decides them jointly; the heavy species alone must already pay.
	std::vector<const ShiftedCCPotential *> pots;
	for (const auto & sp : p.species_constraints)
		pots.push_back(sp.CCPotential.get());
	const std::vector<size_t> gate = ShiftedCCPotential::FinufftForceGate(pots, threads);
	Check(!finufft || (gate[1] && p.species_constraints[1].CCPotential->UseFinufftForces()), label + ": heavy species declined");
	for (size_t s = 0; s < species; s++) {
		const ShiftedCCPotential & cc = *p.species_constraints[s].CCPotential;
		const bool type2 = finufft && gate[s];
		Check(cc.GetLastForceBackendUsed() == (type2 ? ForceBackend::FinufftType2 : ForceBackend::Direct),
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

	// Small dense 3D system, Type 2 too small to split: it wins clearly against one direct thread, not against eight.
	const Configuration dense = RandomCell(3, 200, false);
	ShiftedCCPotential one(3);
	AddModes(one, Modes(dense, 1.0), false);
	one.SetRhoBackend(Backend::FinufftType1);
	ShiftedCCPotential eight(one);
	one.ParallelNumThread = 1;
	eight.ParallelNumThread = 8;
	one.SetConfiguration(dense);
	eight.SetConfiguration(dense);
	one.AllForce(f);
	Check(one.UseFinufftForces() && one.GetLastForceBackendUsed() == ForceBackend::FinufftType2, label + ": one-thread dense case did not use Type 2");
	Check(ShiftedCCPotential::FinufftForceGate({&eight}, 64)[0] <= 1, label + ": dense case splits its Type 2 work");
	eight.AllForce(f);
	Check(!eight.UseFinufftForces() && eight.GetLastForceBackendUsed() == ForceBackend::Direct, label + ": eight-thread dense case used Type 2");

	// Mid-size 3D system: one Type 2 block on one thread, several in parallel on eight; the forces agree bitwise.
	const Configuration mid = RandomCell(3, 500, false);
	ShiftedCCPotential serial(3);
	AddModes(serial, Modes(mid, 0.3), false);
	serial.SetRhoBackend(Backend::FinufftType1);
	ShiftedCCPotential split(serial);
	serial.ParallelNumThread = 1;
	split.ParallelNumThread = 8;
	serial.SetConfiguration(mid);
	split.SetConfiguration(mid);
	std::vector<GeometryVector> fs;
	serial.AllForce(f);
	split.AllForce(fs);
	Check(ShiftedCCPotential::FinufftForceGate({&serial}, 1)[0] == 1 && ShiftedCCPotential::FinufftForceGate({&split}, 8)[0] > 1,
	      label + ": mid case not one Type 2 block at T=1 and several at T=8");
	Check(serial.GetLastForceBackendUsed() == ForceBackend::FinufftType2 && split.GetLastForceBackendUsed() == ForceBackend::FinufftType2,
	      label + ": mid case did not use Type 2");
	Check(Same(f, fs, 3), label + ": Type 2 forces depend on the block count");
}

/** \brief One case of the combined check: species sizes and chi, cell, contact, and the expected gate decisions. */
struct Combo {
	DimensionType d;              //< dimension
	bool skew;                    //< non-orthogonal cell
	std::vector<size_t> sizes;    //< particles per species
	std::vector<double> chi;      //< per-species chi (at the species density)
	double e0;                    //< contact strength (0 disables contact)
	std::string expect[2];        //< per species at T=4 and T=8: '0' direct, '1' Type 2, '2' Type 2 only jointly
};

/** \brief Species-resolved potential on c with interleaved species (indexed views differ from contiguous blocks).
 *  \param c        Configuration.
 *  \param tc       Case.
 *  \param backend  Collective-density backend.
 *  \param threads  ParallelNumThread.
 *  \return         The potential, bound to c. */
std::unique_ptr<MultiStealthyPotential> MakeCombo(const Configuration & c, const Combo & tc, Backend backend, int threads)
{
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	auto p = std::make_unique<MultiStealthyPotential>(tc.d, tc.e0);
	p->ParallelNumThread = threads;
	p->SetConfiguration(c);
	p->SetNumSpecies(tc.sizes.size());
	p->SetSecureRadii(std::vector<double>(c.NumParticle(), std::pow(0.2 / HyperSphere_Volume(tc.d, 1.0), 1.0 / tc.d)));
	std::vector<size_t> label;
	for (size_t s = 0; s < tc.sizes.size(); s++)
		label.insert(label.end(), tc.sizes[s], s);
	RandomGenerator gen(99 + label.size());
	for (size_t i = label.size(); i > 1; i--)
		std::swap(label[i - 1], label[size_t(gen.RandomDouble() * i) % i]);
	for (size_t i = 0; i < label.size(); i++)
		p->SetSpeciesOfParticle(i, label[i]);
	for (size_t s = 0; s < tc.sizes.size(); s++) {
		const double rhoS = double(tc.sizes[s]) / c.PeriodicVolume();
		const double K = 2 * pi * std::pow(2 * tc.d * rhoS * tc.chi[s] / HyperSphere_Volume(tc.d, 1.0), 1.0 / tc.d);
		for (const auto & k : GetKs(c, K, K, 1))
			p->AddConstraint(s, k, {1.0 + 0.1 * k.Modulus2()});
	}
	p->SetRhoBackend(backend);
	p->SetConfiguration(c);
	std::cout.rdbuf(old);
	return p;
}

/** \brief Combined check of the integrated design on one case: exact-type species bound as indexed views, persistent
 *  Type 1 and Type 2 plans, species-concurrent rho and Type 2 forces, and the joint cost gate. At T=4 and T=8:
 *  every species' force path equals FinufftForceGate over the concurrent species and the expected decision
 *  (expectations follow the calibrated constants: recheck them after a refit); a cold AllForce (MD order) equals
 *  the warm one bitwise; energy and forces match Direct within 1e-10 (relative, and times max(1, max|F|)); repeated
 *  evaluations, a fresh potential, and a moved configuration repeat bitwise; and when the decisions at T=4 and T=8
 *  agree, so do the results, bitwise.
 *  \param tc  Case. */
void Combined(const Combo & tc)
{
	size_t n = 0;
	for (size_t ns : tc.sizes) n += ns;
	const Configuration c = RandomCell(tc.d, n, tc.skew);
	Configuration moved(c);
	RandomGenerator gen(5 + n);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x = moved.GetRelativeCoordinates(i);
		for (DimensionType a = 0; a < tc.d; a++) x.x[a] += 0.05 * gen.RandomDouble();
		moved.MoveParticle(i, x);
	}
	std::string decided[2];
	double energy[2] = {0, 0};
	std::vector<GeometryVector> force[2];
	for (int t = 0; t < 2; t++) {
		const int threads = t ? 8 : 4;
		std::ostringstream name;
		name << "combined " << tc.d << "D" << (tc.skew ? " skew" : "") << " S=" << tc.sizes.size() << " e0=" << tc.e0 << " T=" << threads;
		const std::string label = name.str();
		auto direct = MakeCombo(c, tc, Backend::Direct, threads);
		auto fast = MakeCombo(c, tc, Backend::FinufftType1, threads);
		auto again = MakeCombo(c, tc, Backend::FinufftType1, threads);

		std::vector<GeometryVector> cold, warm, repeat, fd;
		fast->AllForce(cold);  // MD order: no prior Energy
		const double e1 = fast->Energy();
		fast->AllForce(warm);
		Check(Same(cold, warm, tc.d), label + ": cold AllForce differs from warm AllForce");

		std::vector<const ShiftedCCPotential *> pots;
		for (const auto & sp : fast->species_constraints)
			pots.push_back(sp.CCPotential.get());
		const std::vector<size_t> gate = ShiftedCCPotential::FinufftForceGate(pots, threads);
		for (size_t s = 0; s < tc.sizes.size(); s++) {
			const ShiftedCCPotential & cc = *fast->species_constraints[s].CCPotential;
			const bool type2 = cc.GetLastForceBackendUsed() == ForceBackend::FinufftType2;
			decided[t] += type2 ? '1' : '0';
			const std::string tag = label + " species " + std::to_string(s);
			Check(cc.GetLastRhoBackendUsed() == Backend::FinufftType1, tag + ": rho fell back");
			Check(type2 == (gate[s] != 0), tag + ": force path differs from FinufftForceGate");
			Check(type2 == (tc.expect[t][s] != '0'), tag + ": gate decision differs from the expected one");
			if (tc.expect[t][s] == '2')
				Check(!cc.UseFinufftForces(), tag + ": accepted alone, expected only jointly");
		}

		const double ed = direct->Energy();
		direct->AllForce(fd);
		double scale = 0;
		const double df = ForceDiff(warm, fd, tc.d, scale);
		Check(std::abs(e1 - ed) <= 1e-10 * std::abs(ed), label + ": energy rel diff " + std::to_string(std::abs(e1 - ed) / std::abs(ed)));
		Check(df <= 1e-10 * std::max(1.0, scale), label + ": force diff " + std::to_string(df) + " (scale " + std::to_string(scale) + ")");

		const double e2 = again->Energy();
		again->AllForce(repeat);
		Check(std::memcmp(&e1, &e2, sizeof(double)) == 0 && Same(warm, repeat, tc.d), label + ": fresh potential not bitwise repeatable");
		for (int rep = 0; rep < 3; rep++) {
			fast->SetConfiguration(c);
			const double e = fast->Energy();
			fast->AllForce(repeat);
			Check(std::memcmp(&e, &e1, sizeof(double)) == 0 && Same(warm, repeat, tc.d), label + ": repeated evaluation not bitwise identical");
		}

		// Moved particles: plans reused, angles refilled; agreement and repeatability again (MD order on one of them).
		std::vector<GeometryVector> fm, fm2, fmd;
		fast->SetConfiguration(moved);
		again->SetConfiguration(moved);
		direct->SetConfiguration(moved);
		fast->AllForce(fm);
		const double em = fast->Energy(), em2 = again->Energy(), emd = direct->Energy();
		again->AllForce(fm2);
		direct->AllForce(fmd);
		Check(std::memcmp(&em, &em2, sizeof(double)) == 0 && Same(fm, fm2, tc.d), label + ": moved, not bitwise repeatable");
		Check(std::abs(em - emd) <= 1e-10 * std::abs(emd) && ForceDiff(fm, fmd, tc.d, scale) <= 1e-10 * std::max(1.0, scale), label + ": moved, differs from Direct");

		energy[t] = e1;
		force[t] = warm;
	}
	if (decided[0] == decided[1])
		Check(std::memcmp(&energy[0], &energy[1], sizeof(double)) == 0 && Same(force[0], force[1], tc.d),
		      "combined " + std::to_string(tc.d) + "D S=" + std::to_string(tc.sizes.size()) + ": same decisions at T=4 and T=8 but different bits");
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
	const std::vector<Combo> combos = {
		// uniform mid-size species: Type 2 at T=4 on their own; at T=8 only jointly (F = 6 concurrent species)
		{2, true, std::vector<size_t>(6, 150), std::vector<double>(6, 0.2), 1.0, {"111111", "222222"}},
		// one heavy and two light species: only the heavy one pays, alone or jointly
		{2, false, {700, 120, 120}, {0.4, 0.05, 0.05}, 0.0, {"100", "100"}},
		{1, false, {400, 300, 200}, {0.5, 0.5, 0.5}, 1.0, {"111", "111"}},
		// 3D: both alone at T=4; at T=8 the smaller species only jointly
		{3, true, {500, 400}, {1.0, 1.0}, 1.0, {"11", "12"}},
		// tiny 3D species: Type 2 never pays
		{3, false, std::vector<size_t>(4, 60), std::vector<double>(4, 0.2), 0.0, {"0000", "0000"}},
	};
	for (const auto & tc : combos)
		Combined(tc);
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
