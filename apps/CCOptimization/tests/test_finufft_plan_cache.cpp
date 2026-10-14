/** \file test_finufft_plan_cache.cpp
 * \brief Checks of the persistent FINUFFT Type 1 plan and mode map of ShiftedCCPotential (option B).
 *
 * The program links with -Wl,--wrap=finufft_makeplan,--wrap=finufft_destroy,--wrap=finufft_setpts, so every
 * plan that ShiftedCCPotential creates or destroys is observed. Checks:
 *  - rho is bitwise identical to the simple interface finufft{1,2,3}d1 with the same options, and within
 *    FINUFFT tolerance of the simple interface with default options (the previous code path);
 *  - repeated evaluations with moved particles reuse the plan (no makeplan);
 *  - in-place edits of the public constraints vector, a rescaled or different cell, a particle-count change,
 *    and changes of FinufftTolerance, FinufftMaxModes and FinufftNumThreads rebuild the map or plan as needed,
 *    and rho always matches a fresh Direct evaluation;
 *  - an ineligible (fallback) result and a failed makeplan are cached and re-checked when an input changes;
 *  - the Type 2 force plan: copies start without it, assignment, a declining gate and destruction release it,
 *    and a failed Type 2 setpts falls back to the direct forces for that call only;
 *  - copies (copy constructor, clone, derived converting constructors) start empty, assignment drops the
 *    target's plan, and every plan is destroyed exactly once whatever the destruction order (no leak, no
 *    double destroy); MultiStealthyPotential makes one Type 1 plan per species, one Type 2 plan per species whose
 *    forces use Type 2, and no more on later callbacks.
 * Build and run: make test-finufft-plan. AddressSanitizer run (from apps/CCOptimization, after make deps):
 *   g++ -std=c++17 -g -O1 -fopenmp -fsanitize=address -DCCPOTENTIAL_USE_FINUFFT -I../../core -I../../pairstat
 *     -I../../potential -I. -I$HOME/install/include tests/test_finufft_plan_cache.cpp MultiStealthyPotential.cpp
 *     ../../core/obj/?*.o ../../pairstat/obj/?*.o ../../potential/obj/Potential.o
 *     -Wl,--wrap=finufft_makeplan,--wrap=finufft_destroy,--wrap=finufft_setpts
 *     -L$HOME/install/lib -lfinufft -lgsl -lgslcblas -lnlopt -Wl,-rpath,$HOME/install/lib -o /tmp/plan_asan
 *   && /tmp/plan_asan */

#include "MultiStealthyPotential.h"

#include <cstring>
#include <memory>
#include <random>
#include <set>
#include <sstream>

size_t Verbosity = 0;

extern "C" {
int __real_finufft_makeplan(int type, int dim, const int64_t * n_modes, int iflag, int ntrans, double tol,
                            finufft_plan * plan, const finufft_opts * opts);
int __real_finufft_destroy(finufft_plan plan);
int __real_finufft_setpts(finufft_plan plan, int64_t M, const double * x, const double * y, const double * z,
                          int64_t N, const double * s, const double * t, const double * u);
}

namespace obs {
size_t makeplans = 0;            //< finufft_makeplan calls
size_t makeplans2 = 0;           //< finufft_makeplan calls of Type 2 (force) plans
size_t setpts = 0;               //< finufft_setpts calls
size_t badDestroys = 0;          //< finufft_destroy of a handle that is not live (double destroy)
int lastNthreads = -1;           //< opts->nthreads of the last makeplan
double lastUpsampfac = -1;       //< opts->upsampfac of the last makeplan
int lastSort = -1;               //< opts->spread_sort of the last makeplan
std::set<finufft_plan> live;     //< plans made and not yet destroyed
int injectMakeplan = 0;          //< nonzero: the next makeplan returns this code (after making a real plan)
int injectSetpts = 0;            //< nonzero: the next setpts returns this code without running
}  // namespace obs

/** \brief Observe finufft_makeplan, optionally reporting an injected warning code with a valid plan.
 *  \return The real status, or the injected code. */
extern "C" int __wrap_finufft_makeplan(int type, int dim, const int64_t * n_modes, int iflag, int ntrans, double tol,
                                       finufft_plan * plan, const finufft_opts * opts)
{
	obs::makeplans++;
	obs::makeplans2 += type == 2;
	obs::lastNthreads = opts->nthreads;
	obs::lastUpsampfac = opts->upsampfac;
	obs::lastSort = opts->spread_sort;
	int status = __real_finufft_makeplan(type, dim, n_modes, iflag, ntrans, tol, plan, opts);
	if (*plan != nullptr)
		obs::live.insert(*plan);
	if (obs::injectMakeplan != 0) {
		status = obs::injectMakeplan;
		obs::injectMakeplan = 0;
	}
	return status;
}

/** \brief Observe finufft_destroy; a handle that is not live is counted and not passed on.
 *  \return The real status, or 1 for a handle that is not live. */
extern "C" int __wrap_finufft_destroy(finufft_plan plan)
{
	if (obs::live.erase(plan) == 0) {
		obs::badDestroys++;
		return 1;
	}
	return __real_finufft_destroy(plan);
}

/** \brief Observe finufft_setpts, optionally failing once with an injected code.
 *  \return The real status, or the injected code. */
extern "C" int __wrap_finufft_setpts(finufft_plan plan, int64_t M, const double * x, const double * y, const double * z,
                                     int64_t N, const double * s, const double * t, const double * u)
{
	obs::setpts++;
	if (obs::injectSetpts != 0) {
		const int status = obs::injectSetpts;
		obs::injectSetpts = 0;
		return status;
	}
	return __real_finufft_setpts(plan, M, x, y, z, N, s, t, u);
}

namespace {

using Backend = ShiftedCCPotential::RhoBackend;
using Mode = std::vector<long>;
using Complex = std::complex<double>;

int failures = 0;
int checks = 0;

/** \brief Record one check.
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

/** \brief Configuration with uniformly random fractional coordinates.
 *  \param rows  Basis vectors, one row each.
 *  \param n     Particle count.
 *  \param rng   Random engine.
 *  \return      The configuration. */
Configuration MakeConfig(const std::vector<std::vector<double>> & rows, size_t n, std::mt19937_64 & rng)
{
	const DimensionType dim = rows.size();
	GeometryVector basis[3] = {GeometryVector(dim), GeometryVector(dim), GeometryVector(dim)};
	for (DimensionType d = 0; d < dim; d++)
		for (DimensionType c = 0; c < dim; c++)
			basis[d].x[c] = rows[d][c];
	Configuration config(dim, basis, 1.0);
	std::uniform_real_distribution<double> u(0.0, 1.0);
	for (size_t j = 0; j < n; j++) {
		GeometryVector s(dim);
		for (DimensionType d = 0; d < dim; d++) s.x[d] = u(rng);
		config.Insert("a", s);
	}
	return config;
}

/** \brief Reciprocal-lattice vector sum_d m_d b_d of a configuration.
 *  \param c  Configuration.
 *  \param m  Integer indices.
 *  \return   The wavevector. */
GeometryVector LatticeK(const Configuration & c, const Mode & m)
{
	GeometryVector k(c.GetDimension());
	for (DimensionType d = 0; d < c.GetDimension(); d++)
		k.AddFrom(double(m[d]) * c.GetReciprocalBasisVector(d));
	return k;
}

/** \brief Nonzero integer modes in [-K, K]^dim.
 *  \param dim  Dimension.
 *  \param K    Largest index.
 *  \return     The modes. */
std::vector<Mode> Modes(DimensionType dim, long K)
{
	std::vector<Mode> out;
	const long side = 2 * K + 1;
	long total = 1;
	for (DimensionType d = 0; d < dim; d++) total *= side;
	for (long t = 0; t < total; t++) {
		Mode m(dim);
		long r = t;
		bool zero = true;
		for (DimensionType d = 0; d < dim; d++) { m[d] = r % side - K; r /= side; zero = zero && m[d] == 0; }
		if (!zero) out.push_back(m);
	}
	return out;
}

/** \brief Potential with constraints at the given modes of c, bound to c, FINUFFT selected.
 *  \param c      Configuration.
 *  \param modes  Integer modes.
 *  \return       The potential. */
ShiftedCCPotential MakePot(const Configuration & c, const std::vector<Mode> & modes)
{
	ShiftedCCPotential pot(c.GetDimension());
	for (const Mode & m : modes)
		pot.constraints.emplace_back(LatticeK(c, m), 1.0 + 0.01 * double(pot.constraints.size() % 7));
	pot.SetRhoBackend(Backend::FinufftType1);
	pot.SetConfiguration(c);
	return pot;
}

/** \brief Simple-interface Type 1 rho at the given modes (the previous code path when explicit is false).
 *  \param c         Configuration.
 *  \param modes     Integer modes.
 *  \param tol       FINUFFT tolerance.
 *  \param explicitOpts  Also set upsampfac=2.0 and spread_sort=1, as the persistent plan does.
 *  \return          rho at each mode. */
std::vector<Complex> SimpleRho(const Configuration & c, const std::vector<Mode> & modes, double tol, bool explicitOpts)
{
	const DimensionType dim = c.GetDimension();
	const size_t np = c.NumParticle();
	std::int64_t K[3] = {0, 0, 0}, ms[3] = {1, 1, 1};
	for (const Mode & m : modes)
		for (DimensionType d = 0; d < dim; d++) K[d] = std::max<std::int64_t>(K[d], std::abs(m[d]));
	for (DimensionType d = 0; d < dim; d++) ms[d] = 2 * K[d] + 1;
	std::vector<double> a[3];
	for (DimensionType d = 0; d < dim; d++) a[d].resize(np);
	for (size_t j = 0; j < np; j++)
		for (DimensionType d = 0; d < dim; d++) {
			const double s = c.GetRelativeCoordinates(j).x[d];
			a[d][j] = 2.0 * ::pi * (s - std::floor(s));
		}
	std::vector<Complex> one(np, Complex(1.0, 0.0)), f(ms[0] * ms[1] * ms[2]);
	finufft_opts opts;
	finufft_default_opts(&opts);
	opts.modeord = 0;
	opts.nthreads = 1;
	if (explicitOpts) { opts.upsampfac = 2.0; opts.spread_sort = 1; }
	int status;
	if (dim == 1) status = finufft1d1(np, a[0].data(), one.data(), +1, tol, ms[0], f.data(), &opts);
	else if (dim == 2) status = finufft2d1(np, a[0].data(), a[1].data(), one.data(), +1, tol, ms[0], ms[1], f.data(), &opts);
	else status = finufft3d1(np, a[0].data(), a[1].data(), a[2].data(), one.data(), +1, tol, ms[0], ms[1], ms[2], f.data(), &opts);
	Check(status == 0, "simple interface status");
	std::vector<Complex> rho;
	for (const Mode & m : modes) {
		std::int64_t idx = 0;
		for (int d = int(dim) - 1; d >= 0; d--) idx = idx * ms[d] + (m[d] + K[d]);
		rho.push_back(f[idx]);
	}
	return rho;
}

/** \brief Bitwise equality of a potential's cached rho with a reference.
 *  \param pot  Evaluated potential.
 *  \param ref  Reference rho.
 *  \return     true if every real and imaginary part has identical bits. */
bool SameBits(const ShiftedCCPotential & pot, const std::vector<Complex> & ref)
{
	if (pot.RhoReal.size() != ref.size()) return false;
	for (size_t i = 0; i < ref.size(); i++) {
		const double re = ref[i].real(), im = ref[i].imag();
		if (std::memcmp(&pot.RhoReal[i], &re, sizeof(double)) != 0 || std::memcmp(&pot.RhoImag[i], &im, sizeof(double)) != 0)
			return false;
	}
	return true;
}

/** \brief Largest |rho| difference between a potential's cached rho and a reference.
 *  \param pot  Evaluated potential.
 *  \param ref  Reference rho.
 *  \return     max_i |rho_i - ref_i| (infinity on a size mismatch). */
double MaxDiff(const ShiftedCCPotential & pot, const std::vector<Complex> & ref)
{
	if (pot.RhoReal.size() != ref.size()) return INFINITY;
	double diff = 0;
	for (size_t i = 0; i < ref.size(); i++)
		diff = std::max(diff, std::abs(Complex(pot.RhoReal[i], pot.RhoImag[i]) - ref[i]));
	return diff;
}

/** \brief Re-evaluate pot (bound to c) and compare it with a fresh Direct copy.
 *  \param pot       Potential; its cached rho is invalidated first.
 *  \param c         Configuration pot is bound to.
 *  \param expected  Backend pot must report.
 *  \param label     Case name. */
void CheckDirect(ShiftedCCPotential & pot, const Configuration & c, Backend expected, const std::string & label)
{
	pot.RhoValid = false;
	pot.GetRho();
	Check(pot.GetLastRhoBackendUsed() == expected, label + ": expected backend");
	if (expected == Backend::FinufftType1)
		Check(pot.GetLastFinufftStatus() == 0, label + ": status 0");
	ShiftedCCPotential ref(pot);
	ref.SetRhoBackend(Backend::Direct);
	ref.SetConfiguration(c);
	ref.GetRho();
	std::vector<Complex> r;
	for (size_t i = 0; i < ref.RhoReal.size(); i++) r.emplace_back(ref.RhoReal[i], ref.RhoImag[i]);
	const double diff = MaxDiff(pot, r);
	if (expected == Backend::FinufftType1)
		Check(diff <= 1e-10 * double(c.NumParticle()), label + ": rho differs from Direct by " + std::to_string(diff));
	else
		Check(diff == 0.0, label + ": fallback rho differs from Direct");
}

/** \brief Bitwise agreement with the simple interface and plan reuse across moves, for one cell.
 *  \param rows   Basis rows.
 *  \param n      Particle count.
 *  \param K      Largest mode index.
 *  \param label  Case name. */
void TestAgreement(const std::vector<std::vector<double>> & rows, size_t n, long K, const std::string & label)
{
	std::mt19937_64 rng(11 + n);
	Configuration c = MakeConfig(rows, n, rng);
	std::vector<Mode> modes = Modes(c.GetDimension(), K);
	std::shuffle(modes.begin(), modes.end(), rng);
	ShiftedCCPotential pot = MakePot(c, modes);
	const size_t made = obs::makeplans;
	pot.GetRho();
	Check(obs::makeplans == made + 1, label + ": one makeplan on first use");
	Check(obs::lastNthreads == 1 && obs::lastUpsampfac == 2.0 && obs::lastSort == 1, label + ": plan options");
	Check(pot.GetLastRhoBackendUsed() == Backend::FinufftType1, label + ": FINUFFT used");
	Check(SameBits(pot, SimpleRho(c, modes, 1e-12, true)), label + ": not bitwise equal to simple interface with the same options");
	const std::vector<Complex> old = SimpleRho(c, modes, 1e-12, false);
	Check(MaxDiff(pot, old) <= 1e-10 * double(n), label + ": differs from the default-option simple interface");
	std::cout << "  " << label << ": bitwise equal to default-option simple interface (previous path): "
	          << (SameBits(pot, old) ? "yes" : "no") << '\n';

	// moved particles, rebound copies of the configuration: same plan, map not rebuilt
	for (int step = 0; step < 3; step++) {
		Configuration moved(c);
		std::uniform_real_distribution<double> u(-0.05, 0.05);
		for (size_t j = 0; j < n; j++) {
			GeometryVector s = moved.GetRelativeCoordinates(j);
			for (DimensionType d = 0; d < c.GetDimension(); d++) s.x[d] += u(rng);
			moved.MoveParticle(j, s);
		}
		pot.SetConfiguration(moved);
		pot.GetRho();
		Check(SameBits(pot, SimpleRho(moved, modes, 1e-12, true)), label + ": moved, not bitwise equal to simple interface");
		c = moved;
		pot.SetConfiguration(c);
	}
	Check(obs::makeplans == made + 1, label + ": moves made a new plan");
}

/** \brief Invalidation: constraints edited in place, cell changes, particle count and option changes. */
void TestInvalidation()
{
	std::mt19937_64 rng(77);
	Configuration c = MakeConfig({{9.0, 0.0}, {2.7, 8.1}}, 250, rng);
	std::vector<Mode> modes = Modes(2, 5);
	ShiftedCCPotential pot = MakePot(c, modes);
	CheckDirect(pot, c, Backend::FinufftType1, "baseline");
	const std::vector<double> baseReal = pot.RhoReal;
	size_t made = obs::makeplans;

	// in-place edits of the public vector (RhoValid is reset by CheckDirect, as a caller must do):
	// swap two k (same grid), change one k inside the grid
	std::swap(pot.constraints[0].k, pot.constraints[17].k);
	pot.constraints[3].k = LatticeK(c, {-4, 2});
	CheckDirect(pot, c, Backend::FinufftType1, "edited k, same grid");
	Check(obs::makeplans == made, "edited k, same grid: plan reused");
	// a k outside the grid: larger plan
	pot.constraints[5].k = LatticeK(c, {7, -6});
	CheckDirect(pot, c, Backend::FinufftType1, "edited k, larger grid");
	Check(obs::makeplans == ++made, "edited k, larger grid: new plan");
	// a non-lattice k: cached fallback, no FINUFFT call, then repaired
	const GeometryVector keep = pot.constraints[8].k;
	pot.constraints[8].k = 1.37 * c.GetReciprocalBasisVector(0);
	const size_t live = obs::live.size(), sets = obs::setpts;
	CheckDirect(pot, c, Backend::Direct, "non-lattice k");
	Check(pot.GetLastFinufftStatus() == ShiftedCCPotential::FinufftNotCalled, "non-lattice k: FINUFFT not called");
	Check(obs::live.size() == live - 1, "non-lattice k: plan released");
	CheckDirect(pot, c, Backend::Direct, "non-lattice k, cached");
	Check(obs::setpts == sets && obs::makeplans == made, "non-lattice k: cached ineligibility calls no FINUFFT");
	pot.constraints[8].k = keep;
	CheckDirect(pot, c, Backend::FinufftType1, "repaired k");
	Check(obs::makeplans == ++made, "repaired k: new plan");
	// appended and removed constraints
	pot.constraints.emplace_back(LatticeK(c, {1, 1}), 0.5);
	CheckDirect(pot, c, Backend::FinufftType1, "appended k");
	pot.constraints.pop_back();
	pot.constraints.pop_back();
	CheckDirect(pot, c, Backend::FinufftType1, "removed k");

	// rescaled cell with k rescaled to the same integer modes: map rebuilt, plan reused
	std::vector<Mode> current;
	for (const auto & con : pot.constraints) {
		Mode m(2);
		for (DimensionType d = 0; d < 2; d++) m[d] = std::lround(con.k.Dot(c.GetBasisVector(d)) / (2 * ::pi));
		current.push_back(m);
	}
	Configuration scaled(c);
	scaled.Resize(1.21 * c.PeriodicVolume());	// every length times 1.1
	for (size_t i = 0; i < current.size(); i++) pot.constraints[i].k = LatticeK(scaled, current[i]);
	pot.SetConfiguration(scaled);
	made = obs::makeplans;
	CheckDirect(pot, scaled, Backend::FinufftType1, "rescaled cell, rescaled k");
	Check(obs::makeplans == made, "rescaled cell, same modes: plan reused");
	// the same k on a twice larger cell: integer modes double, larger grid
	Configuration doubled(scaled);
	doubled.Resize(4.0 * scaled.PeriodicVolume());
	pot.SetConfiguration(doubled);
	CheckDirect(pot, doubled, Backend::FinufftType1, "doubled cell");
	Check(obs::makeplans == ++made, "doubled cell: new plan");
	// a cell where the k are no longer lattice vectors: fallback
	Configuration other = MakeConfig({{9.3, 0.0}, {1.1, 8.0}}, 250, rng);
	pot.SetConfiguration(other);
	CheckDirect(pot, other, Backend::Direct, "unrelated cell");
	Check(pot.GetLastFinufftStatus() == ShiftedCCPotential::FinufftNotCalled, "unrelated cell: FINUFFT not called");
	pot.SetConfiguration(doubled);
	CheckDirect(pot, doubled, Backend::FinufftType1, "back to lattice cell");

	// particle count change: map re-checked, plan reused
	Configuration more(doubled);
	more.Insert("a", GeometryVector(0.3, 0.6));
	pot.SetConfiguration(more);
	made = obs::makeplans;
	CheckDirect(pot, more, Backend::FinufftType1, "particle added");
	Check(obs::makeplans == made, "particle added: plan reused");

	// options: each change rebuilds; returning to the defaults reproduces rho bitwise
	pot.SetConfiguration(c);
	pot.constraints.resize(modes.size(), pot.constraints[0]);
	for (size_t i = 0; i < modes.size(); i++) pot.constraints[i].k = LatticeK(c, modes[i]);
	CheckDirect(pot, c, Backend::FinufftType1, "restored");
	Check(pot.RhoReal == baseReal, "restored: rho not bitwise repeatable");
	made = obs::makeplans;
	pot.SetFinufftTolerance(1e-9);
	pot.GetRho();
	Check(obs::makeplans == ++made && pot.GetLastRhoBackendUsed() == Backend::FinufftType1, "tolerance: new plan");
	pot.SetFinufftNumThreads(2);
	pot.GetRho();
	Check(obs::makeplans == ++made && obs::lastNthreads == 2, "threads=2: new plan with nthreads=2");
	pot.SetFinufftNumThreads(1);
	pot.SetFinufftTolerance(1e-12);
	pot.GetRho();
	Check(obs::makeplans == ++made && obs::lastNthreads == 1, "defaults: new plan with nthreads=1");
	Check(pot.RhoReal == baseReal, "defaults: rho not bitwise repeatable");
	pot.SetFinufftMaxModes(120);	// grid is 11 x 11 = 121
	CheckDirect(pot, c, Backend::Direct, "max modes 120");
	pot.SetFinufftMaxModes(121);
	CheckDirect(pot, c, Backend::FinufftType1, "max modes 121");
	Check(obs::makeplans == ++made, "max modes 121: new plan");

	// a makeplan failure (warning code with a valid plan) is cached with its code and the plan released
	const size_t liveBefore = obs::live.size();
	obs::injectMakeplan = 1;
	pot.SetFinufftTolerance(1e-11);
	CheckDirect(pot, c, Backend::Direct, "makeplan warning");
	Check(pot.GetLastFinufftStatus() == 1, "makeplan warning: code reported");
	Check(obs::live.size() == liveBefore - 1, "makeplan warning: plans released");
	made = obs::makeplans;
	CheckDirect(pot, c, Backend::Direct, "makeplan warning, cached");
	Check(pot.GetLastFinufftStatus() == 1 && obs::makeplans == made, "makeplan warning: cached with its code");
	pot.SetFinufftTolerance(1e-12);
	CheckDirect(pot, c, Backend::FinufftType1, "after makeplan warning");
	// a setpts failure falls back for that call only and keeps the plan
	made = obs::makeplans;
	obs::injectSetpts = 7;
	CheckDirect(pot, c, Backend::Direct, "setpts failure");
	Check(pot.GetLastFinufftStatus() == 7, "setpts failure: code reported");
	CheckDirect(pot, c, Backend::FinufftType1, "after setpts failure");
	Check(obs::makeplans == made && pot.RhoReal == baseReal, "after setpts failure: plan reused, rho repeatable");
}

/** \brief Copies, clone, converting constructors and assignment never share a plan; any destruction order. */
void TestCopies()
{
	std::mt19937_64 rng(99);
	Configuration c = MakeConfig({{6.0, 0.0, 0.0}, {1.0, 6.5, 0.0}, {-0.7, 0.9, 5.8}}, 300, rng);
	const std::vector<Mode> modes = Modes(3, 3);
	const size_t live0 = obs::live.size();
	{
		ShiftedCCPotential a = MakePot(c, modes);
		a.GetRho();
		const std::vector<double> ref = a.RhoReal;
		auto b = std::make_unique<ShiftedCCPotential>(a);	// copy of an object with a live plan
		Check(!b->RhoValid || b->RhoReal == ref, "copy: rho copied");
		b->RhoValid = false;
		b->GetRho();
		Check(obs::live.size() == live0 + 2 && b->RhoReal == ref, "copy: own plan, same rho");
		std::unique_ptr<ShiftedCCPotential> cl(b->clone());
		cl->RhoValid = false;
		cl->GetRho();
		Check(obs::live.size() == live0 + 3 && cl->RhoReal == ref, "clone: own plan, same rho");
		b.reset();	// destroy the source before its clone
		cl->RhoValid = false;
		cl->GetRho();
		Check(obs::live.size() == live0 + 2 && cl->RhoReal == ref, "clone outlives source");
		ShiftedCCPotential_ForMC mc(a);
		LogCCPotential lg(a);
		mc.RhoValid = false;
		lg.RhoValid = false;
		mc.GetRho();
		lg.GetRho();
		Check(obs::live.size() == live0 + 4 && mc.RhoReal == ref && lg.RhoReal == ref, "derived copies: own plans");
		// assignment drops the target's plan and never shares the source's; self-assignment is harmless
		ShiftedCCPotential & alias = mc;
		alias = mc;
		Check(obs::live.size() == live0 + 3, "self-assignment releases only its own plan");
		ShiftedCCPotential d = MakePot(c, {{1, 0, 0}, {0, 2, -1}});
		d.GetRho();
		d = a;
		Check(obs::live.size() == live0 + 3, "assignment released the target plan");
		d.SetConfiguration(c);
		d.GetRho();
		mc.RhoValid = false;
		mc.GetRho();
		Check(obs::live.size() == live0 + 5 && d.RhoReal == ref && mc.RhoReal == ref, "assigned objects rebuild their own plans");
	}
	Check(obs::live.size() == live0, "copies: every plan destroyed at scope exit");
}

/** \brief The Type 2 plan: copies start without it, assignment, a declining gate and destruction release it, and a
 *  failed Type 2 setpts falls back to the direct forces for that call only. */
void TestType2Plan()
{
	std::mt19937_64 rng(123);
	Configuration c = MakeConfig({{14.0, 0.0}, {3.0, 13.0}}, 500, rng);
	const size_t live0 = obs::live.size();
	{
		ShiftedCCPotential a = MakePot(c, Modes(2, 8));
		ShiftedCCPotential direct(a);
		direct.SetRhoBackend(Backend::Direct);
		direct.SetConfiguration(c);
		std::vector<GeometryVector> fa, fb, fd;
		direct.AllForce(fd);
		a.PrepareFinufft();
		Check(a.PrepareFinufftForces(true) && obs::live.size() == live0 + 2, "type 2: both plans live");
		Check(a.TryFinufftAllForce(fa) && a.GetLastForceBackendUsed() == ShiftedCCPotential::ForceBackend::FinufftType2, "type 2: forces");
		double diff = 0, scale = 0;
		for (size_t j = 0; j < fd.size(); j++)
			for (DimensionType d = 0; d < 2; d++) {
				diff = std::max(diff, std::abs(fa[j].x[d] - fd[j].x[d]));
				scale = std::max(scale, std::abs(fd[j].x[d]));
			}
		Check(fa.size() == fd.size() && diff <= 1e-10 * std::max(1.0, scale), "type 2: forces differ from Direct by " + std::to_string(diff));
		ShiftedCCPotential b(a);	// copy of an object holding both plans
		Check(obs::live.size() == live0 + 2, "type 2 copy: shares no plan");
		b.RhoValid = false;
		b.PrepareFinufft();
		b.PrepareFinufftForces(true);
		Check(b.TryFinufftAllForce(fb) && obs::live.size() == live0 + 4, "type 2 copy: own plans");
		bool same = fb.size() == fa.size();
		for (size_t j = 0; same && j < fa.size(); j++)
			same = std::memcmp(fa[j].x, fb[j].x, 2 * sizeof(double)) == 0;
		Check(same, "type 2 copy: forces not bitwise repeatable");
		// a failed Type 2 setpts (rho is valid, so the next setpts is the Type 2 one) falls back for that call only
		const size_t made = obs::makeplans;
		obs::injectSetpts = 5;
		Check(!b.TryFinufftAllForce(fb) && b.GetLastForceBackendUsed() == ShiftedCCPotential::ForceBackend::Direct, "type 2 setpts failure: direct fallback");
		Check(b.TryFinufftAllForce(fb) && obs::makeplans == made, "type 2 setpts failure: plan kept and reused");
		b = a;	// assignment releases both of b's plans
		Check(obs::live.size() == live0 + 2, "type 2 assignment released both plans");
		a.PrepareFinufftForces(false);	// the gate declined: the Type 2 plan goes, the Type 1 plan stays
		Check(obs::live.size() == live0 + 1 && !a.TryFinufftAllForce(fb), "declined gate releases only the Type 2 plan");
	}
	Check(obs::live.size() == live0, "type 2: every plan destroyed at scope exit");
}

/** \brief MultiStealthyPotential: one plan per species, none on later optimizer-style callbacks. */
void TestMulti()
{
	std::mt19937_64 rng(5);
	const size_t n = 600, S = 3;
	Configuration c = MakeConfig({{24.0, 0.0}, {7.0, 25.0}}, n, rng);
	std::ostringstream sink;
	std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
	MultiStealthyPotential fast(2, 1.0), direct(2, 1.0);
	for (MultiStealthyPotential * p : {&fast, &direct}) {
		p->ParallelNumThread = 4;
		p->SetConfiguration(c);
		p->SetNumSpecies(S);
		p->SetSecureRadii(std::vector<double>(n, 0.3));
		for (size_t i = 0; i < n; i++) p->SetSpeciesOfParticle(i, i % S);
		for (size_t s = 0; s < S; s++)
			for (const Mode & m : Modes(2, 3 + long(s)))
				p->AddConstraint(s, LatticeK(c, m), {1.0});
		p->SetConfiguration(c);
	}
	std::cout.rdbuf(old);
	fast.SetRhoBackend(Backend::FinufftType1);
	const size_t made = obs::makeplans, made2 = obs::makeplans2;
	for (int step = 0; step < 5; step++) {
		Configuration moved(c);
		for (size_t j = 0; j < n; j++) {
			GeometryVector s = moved.GetRelativeCoordinates(j);
			s.x[0] += 0.01 * step;
			s.x[1] -= 0.007 * step;
			moved.MoveParticle(j, s);
		}
		fast.SetConfiguration(moved);
		direct.SetConfiguration(moved);
		const double ef = fast.Energy(), ed = direct.Energy();
		std::vector<GeometryVector> ff, fd;
		fast.AllForce(ff);
		direct.AllForce(fd);
		double df = 0, scale = 0;
		for (size_t i = 0; i < n; i++)
			for (DimensionType d = 0; d < 2; d++) {
				df = std::max(df, std::abs(ff[i].x[d] - fd[i].x[d]));
				scale = std::max(scale, std::abs(fd[i].x[d]));
			}
		Check(std::abs(ef - ed) <= 1e-10 * std::abs(ed) && df <= 1e-10 * std::max(1.0, scale), "multi: matches Direct at step " + std::to_string(step));
	}
	// AllForce adds one persistent Type 2 plan per species whose forces the gate sent to Type 2.
	size_t type2 = 0;
	for (const auto & sp : fast.species_constraints)
		type2 += sp.CCPotential->GetLastForceBackendUsed() == ShiftedCCPotential::ForceBackend::FinufftType2;
	Check(obs::makeplans - made - (obs::makeplans2 - made2) == S, "multi: one Type 1 plan per species over five callbacks (made " + std::to_string(obs::makeplans - made) + ")");
	Check(obs::makeplans2 - made2 == type2, "multi: one Type 2 plan per Type 2 species over five callbacks");
}

}  // namespace

/** \brief Run every check.
 *  \return 0 if every check passed, 1 otherwise. */
int main()
{
	if (!ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1)) {
		std::cout << "FINUFFT backend not compiled in; build with make USE_FINUFFT=1\n";
		return 1;
	}
	TestAgreement({{7.3}}, 200, 40, "1D");
	TestAgreement({{12.0, 0.0}, {0.0, 12.0}}, 400, 10, "2D");
	TestAgreement({{13.0, 0.0}, {4.0, 11.5}}, 400, 12, "2D skew");
	TestAgreement({{6.0, 0.0, 0.0}, {1.0, 6.5, 0.0}, {-0.7, 0.9, 5.8}}, 300, 4, "3D skew");
	TestInvalidation();
	TestCopies();
	TestType2Plan();
	TestMulti();
	Check(obs::live.empty(), "plans still live at exit: " + std::to_string(obs::live.size()));
	Check(obs::badDestroys == 0, "double or unknown destroys: " + std::to_string(obs::badDestroys));
	std::cout << (failures ? "FAILED: " + std::to_string(failures) + " of " + std::to_string(checks) + " check(s)\n"
	                       : "All " + std::to_string(checks) + " plan-cache checks passed\n");
	return failures ? 1 : 0;
}
