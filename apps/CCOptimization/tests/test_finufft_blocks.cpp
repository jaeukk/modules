/** \file test_finufft_blocks.cpp
 * \brief Checks of the particle-block split of the FINUFFT transforms (ShiftedCCPotential, MultiStealthyPotential).
 *
 * The program links with -Wl,--wrap=finufft_makeplan,--wrap=finufft_destroy,--wrap=finufft_setpts,--wrap=finufft_execute,
 * so every plan and execution is observed. Checks:
 *  - single potentials (1D, 2D skew, 3D skew) large enough for several Type 1 blocks: rho is bitwise the sum, in block
 *    order, of FINUFFT simple-interface Type 1 transforms of each block's particles (contiguous ranges b N / B); it is
 *    within FINUFFT tolerance of the unsplit simple-interface transform and of Direct; it is bitwise identical for
 *    ParallelNumThread 1, 4 and 8 and on repeated calls; the Type 2 forces, one block at T=1 and several at T=8, are
 *    bitwise those of the unsplit simple-interface Type 2 transform and agree with Direct;
 *  - MultiStealthyPotential with mixed species (split, unsplit, FINUFFT falling back, Direct) at T = 1, 4, 8: every
 *    species uses the expected backend, every Type 1 block executes once per evaluation, energy and forces agree with
 *    Direct, and the results are bitwise identical across T, on repetition, for a cold AllForce (MD order) and for
 *    EnergyAndAllForce against Energy then AllForce;
 *  - plan lifecycle: the Type 1 plans follow the block count when the particle count changes (trimmed or completed,
 *    never remade), a tolerance change remakes them all, copies start without plans and get their own, assignment
 *    releases them; Type 2 plans follow the requested block count; a makeplan failure in a later block leaves no plan
 *    and a cached Direct fallback; a failed setpts in one Type 1 or Type 2 block falls back for that call only;
 *    no plan is made inside a parallel region, none leaks and none is destroyed twice.
 * Build and run: make test-finufft-blocks */

#include "MultiStealthyPotential.h"

#include <cstring>
#include <map>
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
int __real_finufft_execute(finufft_plan plan, std::complex<double> * weights, std::complex<double> * result);
}

namespace obs {
std::map<finufft_plan, int> live;  //< plans made and not yet destroyed, with their type
size_t made[3] = {0, 0, 0};        //< finufft_makeplan calls by type (index 1, 2)
size_t executed[3] = {0, 0, 0};    //< finufft_execute calls by plan type
size_t parallelMakeplans = 0;      //< makeplan calls inside an OpenMP parallel region
size_t badDestroys = 0;            //< destroys of a handle that is not live
int failMakeplanAfter = -1;        //< >= 0: count down successful makeplans, then fail one with code 1
int injectSetpts = 0;              //< nonzero: the next setpts returns this code without running
}  // namespace obs

/** \brief Observe finufft_makeplan, optionally failing one call after a number of successes.
 *  \return The real status, or 1 for the injected failure (its plan stays live until the caller destroys it). */
extern "C" int __wrap_finufft_makeplan(int type, int dim, const int64_t * n_modes, int iflag, int ntrans, double tol,
                                       finufft_plan * plan, const finufft_opts * opts)
{
	obs::made[type]++;
	obs::parallelMakeplans += omp_in_parallel() != 0;
	int status = __real_finufft_makeplan(type, dim, n_modes, iflag, ntrans, tol, plan, opts);
	if (*plan != nullptr)
		obs::live[*plan] = type;
	if (obs::failMakeplanAfter == 0)
		status = 1;
	if (obs::failMakeplanAfter >= 0)
		obs::failMakeplanAfter--;
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
	int inject = 0;
#pragma omp critical(observeSetpts)
	{
		inject = obs::injectSetpts;
		obs::injectSetpts = 0;
	}
	return inject != 0 ? inject : __real_finufft_setpts(plan, M, x, y, z, N, s, t, u);
}

/** \brief Count finufft_execute calls by plan type.
 *  \return The real status. */
extern "C" int __wrap_finufft_execute(finufft_plan plan, std::complex<double> * weights, std::complex<double> * result)
{
	int type = 0;
#pragma omp critical(observeExecute)
	{
		const auto it = obs::live.find(plan);
		type = it == obs::live.end() ? 0 : it->second;
		obs::executed[type]++;
	}
	return __real_finufft_execute(plan, weights, result);
}

namespace {

using Backend = ShiftedCCPotential::RhoBackend;
using ForceBackend = ShiftedCCPotential::ForceBackend;
using Mode = std::vector<long>;
using Complex = std::complex<double>;

int failures = 0, checks = 0;

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

/** \brief Configuration of unit number density with uniformly random fractional coordinates.
 *  \param dim   Dimension.
 *  \param n     Particle count.
 *  \param skew  Non-orthogonal cell (dim >= 2).
 *  \param seed  Random seed.
 *  \return      The configuration. */
Configuration MakeConfig(DimensionType dim, size_t n, bool skew, unsigned seed)
{
	const double L = std::pow(double(n), 1.0 / dim);
	GeometryVector basis[3] = {GeometryVector(dim), GeometryVector(dim), GeometryVector(dim)};
	for (DimensionType d = 0; d < dim; d++)
		basis[d].x[d] = L;
	if (skew && dim >= 2) basis[1].x[0] = 0.29 * L;
	if (skew && dim == 3) { basis[2].x[0] = -0.17 * L; basis[2].x[1] = 0.23 * L; }
	Configuration config(dim, basis, 0.1 * L);
	std::mt19937_64 rng(seed);
	std::uniform_real_distribution<double> u(0.0, 1.0);
	for (size_t j = 0; j < n; j++) {
		GeometryVector s(dim);
		for (DimensionType d = 0; d < dim; d++) s.x[d] = u(rng);
		config.Insert("a", s);
	}
	return config;
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

/** \brief Bitwise equality of two force lists.
 *  \param a    First list.
 *  \param b    Second list.
 *  \param dim  Dimension.
 *  \return     true if both have the same size and identical bits. */
bool Same(const std::vector<GeometryVector> & a, const std::vector<GeometryVector> & b, DimensionType dim)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); i++)
		if (std::memcmp(a[i].x, b[i].x, sizeof(double) * dim) != 0) return false;
	return true;
}

/** \brief Largest force-component difference and the largest reference component.
 *  \param a          Forces.
 *  \param ref        Reference forces.
 *  \param dim        Dimension.
 *  \param[out] scale max |ref|.
 *  \return           max |a - ref|. */
double ForceDiff(const std::vector<GeometryVector> & a, const std::vector<GeometryVector> & ref, DimensionType dim, double & scale)
{
	double diff = a.size() == ref.size() ? 0.0 : INFINITY;
	scale = 0;
	for (size_t i = 0; i < std::min(a.size(), ref.size()); i++)
		for (DimensionType d = 0; d < dim; d++) {
			diff = std::max(diff, std::abs(a[i].x[d] - ref[i].x[d]));
			scale = std::max(scale, std::abs(ref[i].x[d]));
		}
	return diff;
}

/** \brief FINUFFT options of the persistent plans (MakeFinufftPlans).
 *  \return The options. */
finufft_opts PlanOptions()
{
	finufft_opts opts;
	finufft_default_opts(&opts);
	opts.modeord = 0;
	opts.nthreads = 1;
	opts.upsampfac = 2.0;
	opts.spread_sort = 1;
	return opts;
}

/** \brief A potential on integer modes, with everything the simple-interface references need. */
struct Setup {
	Configuration config;            //< bound configuration
	std::vector<Mode> modes;         //< integer mode of each constraint
	std::int64_t K[3] = {0, 0, 0};   //< largest |mode| per dimension
	std::int64_t ms[3] = {1, 1, 1};  //< grid sizes 2K+1
	std::vector<double> angles[3];   //< 2 pi s of every particle
	ShiftedCCPotential pot;          //< the potential (FINUFFT backend)

	/** \brief Bind constraints at the given modes, with distinct weights, to config.
	 *  \param c  Configuration.
	 *  \param m  Modes. */
	Setup(const Configuration & c, const std::vector<Mode> & m) : config(c), modes(m), pot(c.GetDimension())
	{
		const DimensionType dim = c.GetDimension();
		for (const Mode & mode : modes) {
			pot.constraints.emplace_back(LatticeK(config, mode), 1.0 + 0.01 * double(pot.constraints.size() % 7));
			for (DimensionType d = 0; d < dim; d++) K[d] = std::max<std::int64_t>(K[d], std::abs(mode[d]));
		}
		for (DimensionType d = 0; d < dim; d++) {
			ms[d] = 2 * K[d] + 1;
			angles[d].resize(c.NumParticle());
			for (size_t j = 0; j < c.NumParticle(); j++) {
				const double s = c.GetRelativeCoordinates(j).x[d];
				angles[d][j] = 2.0 * ::pi * (s - std::floor(s));
			}
		}
		pot.SetRhoBackend(Backend::FinufftType1);
		pot.SetConfiguration(config);
	}

	/** \brief Position of a constraint's mode in one modeord=0 grid.
	 *  \param i  Constraint index.
	 *  \return   The slot. */
	size_t Slot(size_t i) const
	{
		std::int64_t idx = 0;
		for (int d = int(config.GetDimension()) - 1; d >= 0; d--) idx = idx * ms[d] + (modes[i][d] + K[d]);
		return size_t(idx);
	}

	/** \brief Simple-interface Type 1 modes of the particles [lo, hi).
	 *  \param lo  First particle.
	 *  \param hi  One past the last particle.
	 *  \return    The mode grid. */
	std::vector<Complex> Type1(size_t lo, size_t hi) const
	{
		const DimensionType dim = config.GetDimension();
		const finufft_opts opts = PlanOptions();
		std::vector<Complex> one(hi - lo, Complex(1.0, 0.0)), f(ms[0] * ms[1] * ms[2]);
		const double tol = pot.GetFinufftTolerance();
		int status;
		if (dim == 1) status = finufft1d1(hi - lo, &angles[0][lo], one.data(), +1, tol, ms[0], f.data(), &opts);
		else if (dim == 2) status = finufft2d1(hi - lo, &angles[0][lo], &angles[1][lo], one.data(), +1, tol, ms[0], ms[1], f.data(), &opts);
		else status = finufft3d1(hi - lo, &angles[0][lo], &angles[1][lo], &angles[2][lo], one.data(), +1, tol, ms[0], ms[1], ms[2], f.data(), &opts);
		Check(status == 0, "simple-interface Type 1 status");
		return f;
	}

	/** \brief rho at every constraint from Type 1 transforms of `blocks` contiguous particle blocks, summed in block order.
	 *  \param blocks  Block count B (1: the unsplit transform).
	 *  \return        rho per constraint. */
	std::vector<Complex> BlockRho(size_t blocks) const
	{
		const size_t n = config.NumParticle();
		std::vector<Complex> rho(modes.size());
		for (size_t b = 0; b < blocks; b++) {
			const std::vector<Complex> f = Type1(b * n / blocks, (b + 1) * n / blocks);
			for (size_t i = 0; i < modes.size(); i++)
				rho[i] = b == 0 ? f[Slot(i)] : rho[i] + f[Slot(i)];
		}
		return rho;
	}

	/** \brief Unsplit simple-interface Type 2 forces from a given rho, as RunForceBlock defines them.
	 *  \param rho  rho per constraint.
	 *  \return     Forces. */
	std::vector<GeometryVector> Type2(const std::vector<Complex> & rho) const
	{
		const DimensionType dim = config.GetDimension();
		const size_t n = config.NumParticle(), total = ms[0] * ms[1] * ms[2];
		std::vector<Complex> coef(dim * total, Complex(0.0, 0.0)), out(dim * n);
		for (size_t i = 0; i < modes.size(); i++) {
			const Complex w = pot.constraints[i].V * std::conj(rho[i]);
			for (DimensionType a = 0; a < dim; a++)
				coef[a * total + Slot(i)] += pot.constraints[i].k.x[a] * w;
		}
		const finufft_opts opts = PlanOptions();
		const double tol = pot.GetFinufftTolerance();
		int status;
		if (dim == 1) status = finufft1d2many(1, n, angles[0].data(), out.data(), +1, tol, ms[0], coef.data(), &opts);
		else if (dim == 2) status = finufft2d2many(2, n, angles[0].data(), angles[1].data(), out.data(), +1, tol, ms[0], ms[1], coef.data(), &opts);
		else status = finufft3d2many(3, n, angles[0].data(), angles[1].data(), angles[2].data(), out.data(), +1, tol, ms[0], ms[1], ms[2], coef.data(), &opts);
		Check(status == 0, "simple-interface Type 2 status");
		const double scale = 2.0 / config.PeriodicVolume();
		std::vector<GeometryVector> f(n, GeometryVector(dim));
		for (size_t j = 0; j < n; j++)
			for (DimensionType a = 0; a < dim; a++) f[j].x[a] = scale * out[a * n + j].imag();
		return f;
	}
};

/** \brief Bitwise equality of a potential's cached rho with a reference.
 *  \param pot  Evaluated potential.
 *  \param ref  Reference rho.
 *  \return     true if every real and imaginary part has identical bits. */
bool SameRho(const ShiftedCCPotential & pot, const std::vector<Complex> & ref)
{
	if (pot.RhoReal.size() != ref.size()) return false;
	for (size_t i = 0; i < ref.size(); i++) {
		const double re = ref[i].real(), im = ref[i].imag();
		if (std::memcmp(&pot.RhoReal[i], &re, sizeof(double)) != 0 || std::memcmp(&pot.RhoImag[i], &im, sizeof(double)) != 0)
			return false;
	}
	return true;
}

/** \brief Largest |rho - ref| of a potential's cached rho.
 *  \param pot  Evaluated potential.
 *  \param ref  Reference rho.
 *  \return     The difference (infinity on a size mismatch). */
double RhoDiff(const ShiftedCCPotential & pot, const std::vector<Complex> & ref)
{
	if (pot.RhoReal.size() != ref.size()) return INFINITY;
	double diff = 0;
	for (size_t i = 0; i < ref.size(); i++)
		diff = std::max(diff, std::abs(Complex(pot.RhoReal[i], pot.RhoImag[i]) - ref[i]));
	return diff;
}

/** \brief Single potential: block structure of rho, independence of the team, Type 2 against the unsplit transform.
 *  \param dim    Dimension.
 *  \param n      Particle count.
 *  \param K      Largest mode index.
 *  \param label  Case name. */
void Single(DimensionType dim, size_t n, long K, const std::string & label)
{
	Setup s(MakeConfig(dim, n, dim >= 2, 31 + n), Modes(dim, K));
	ShiftedCCPotential & pot = s.pot;
	const size_t made1 = obs::made[1];
	pot.ParallelNumThread = 1;
	pot.GetRho();
	const size_t blocks = obs::made[1] - made1;
	Check(blocks > 1, label + ": expected several Type 1 blocks, got " + std::to_string(blocks));
	Check(pot.GetLastRhoBackendUsed() == Backend::FinufftType1 && pot.GetLastFinufftStatus() == 0, label + ": FINUFFT not used");
	Check(SameRho(pot, s.BlockRho(blocks)), label + ": rho is not the block-order sum of the block transforms");
	const std::vector<Complex> unsplit = s.BlockRho(1);
	Check(RhoDiff(pot, unsplit) <= 1e-10 * double(n), label + ": rho differs from the unsplit transform by " + std::to_string(RhoDiff(pot, unsplit)));
	std::cout << "  " << label << ": " << blocks << " Type 1 blocks; bitwise equal to the unsplit transform: "
	          << (SameRho(pot, unsplit) ? "yes" : "no") << '\n';
	ShiftedCCPotential direct(pot);
	direct.SetRhoBackend(Backend::Direct);
	direct.SetConfiguration(s.config);
	direct.GetRho();
	std::vector<Complex> rd;
	for (size_t i = 0; i < direct.RhoReal.size(); i++) rd.emplace_back(direct.RhoReal[i], direct.RhoImag[i]);
	Check(RhoDiff(pot, rd) <= 1e-10 * double(n), label + ": rho differs from Direct");

	const std::vector<double> re = pot.RhoReal, im = pot.RhoImag;
	for (int threads : {4, 8, 1}) {
		pot.ParallelNumThread = threads;
		for (int rep = 0; rep < 3; rep++) {
			pot.RhoValid = false;
			pot.GetRho();
			Check(pot.RhoReal == re && pot.RhoImag == im, label + ": rho not bitwise repeatable at T=" + std::to_string(threads));
		}
	}

	// Type 2: one block at T=1, several at T=8; both bitwise the unsplit simple-interface transform of the same rho.
	std::vector<Complex> rho;
	for (size_t i = 0; i < re.size(); i++) rho.emplace_back(re[i], im[i]);
	const std::vector<GeometryVector> ref = s.Type2(rho);
	std::vector<GeometryVector> fd, f;
	direct.AllForce(fd);
	for (int threads : {1, 8}) {
		pot.ParallelNumThread = threads;
		const size_t forceBlocks = ShiftedCCPotential::FinufftForceGate({&pot}, threads)[0];
		Check(threads == 1 ? forceBlocks == 1 : forceBlocks > 1, label + ": unexpected Type 2 block count " + std::to_string(forceBlocks) + " at T=" + std::to_string(threads));
		const size_t exec2 = obs::executed[2];
		pot.AllForce(f);
		Check(pot.GetLastForceBackendUsed() == ForceBackend::FinufftType2 && obs::executed[2] - exec2 == forceBlocks, label + ": Type 2 blocks not run (" + std::to_string(obs::executed[2] - exec2) + " executions, " + std::to_string(obs::executed[0]) + " unknown)");
		Check(Same(f, ref, dim), label + ": Type 2 forces in " + std::to_string(forceBlocks) + " blocks differ from the unsplit transform");
		double scale = 0;
		const double diff = ForceDiff(f, fd, dim, scale);
		Check(diff <= 1e-10 * std::max(1.0, scale), label + ": Type 2 forces differ from Direct by " + std::to_string(diff));
	}
}

/** \brief Mixed species in MultiStealthyPotential: split, unsplit, FINUFFT fallback and Direct, at T = 1, 4, 8. */
void Mixed()
{
	const DimensionType dim = 3;
	const std::vector<size_t> sizes = {2000, 300, 400, 300};	// split, unsplit, fallback (mode cap), Direct
	const std::vector<double> chi = {0.05, 0.4, 0.3, 0.3};
	size_t n = 0;
	for (size_t ns : sizes) n += ns;
	const Configuration c = MakeConfig(dim, n, true, 7);
	Configuration moved(c);
	std::mt19937_64 rng(3);
	std::uniform_real_distribution<double> u(-0.02, 0.02);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x = moved.GetRelativeCoordinates(i);
		for (DimensionType d = 0; d < dim; d++) x.x[d] += u(rng);
		moved.MoveParticle(i, x);
	}
	auto make = [&](bool finufft, int threads) {
		std::ostringstream sink;
		std::streambuf * old = std::cout.rdbuf(sink.rdbuf());
		auto p = std::make_unique<MultiStealthyPotential>(dim, 1.0);
		p->ParallelNumThread = threads;
		p->SetConfiguration(c);
		p->SetNumSpecies(sizes.size());
		p->SetSecureRadii(std::vector<double>(n, std::pow(0.2 / HyperSphere_Volume(dim, 1.0), 1.0 / dim)));
		std::vector<size_t> label;	// interleaved species: local and global indices differ
		for (size_t s = 0; s < sizes.size(); s++) label.insert(label.end(), sizes[s], s);
		std::shuffle(label.begin(), label.end(), std::mt19937_64(11));
		for (size_t i = 0; i < n; i++) p->SetSpeciesOfParticle(i, label[i]);
		for (size_t s = 0; s < sizes.size(); s++) {
			const double K = 2 * pi * std::pow(2 * dim * double(sizes[s]) / c.PeriodicVolume() * chi[s] / HyperSphere_Volume(dim, 1.0), 1.0 / dim);
			for (const auto & k : GetKs(c, K, K, 1))
				p->AddConstraint(s, k, {1.0 + 0.1 * k.Modulus2()});
		}
		if (finufft) {
			p->SetRhoBackend(Backend::FinufftType1);
			p->species_constraints[2].CCPotential->SetFinufftMaxModes(1);
			p->species_constraints[3].CCPotential->SetRhoBackend(Backend::Direct);
		}
		p->SetConfiguration(c);
		std::cout.rdbuf(old);
		return p;
	};
	auto direct = make(false, 1);
	const double ed = direct->Energy();
	std::vector<GeometryVector> fd;
	direct->AllForce(fd);
	direct->SetConfiguration(moved);
	const double edm = direct->Energy();
	std::vector<GeometryVector> fdm;
	direct->AllForce(fdm);

	double e0 = 0, em0 = 0;
	std::vector<GeometryVector> f0, fm0;
	std::string decided0;
	for (int threads : {1, 4, 8}) {
		const std::string label = "mixed T=" + std::to_string(threads);
		auto p = make(true, threads);
		std::vector<GeometryVector> cold, warm, fused, again;
		const size_t exec1 = obs::executed[1];
		p->AllForce(cold);	// MD order
		const size_t blocks1 = obs::executed[1] - exec1;
		const double e = p->Energy();
		Check(obs::executed[1] - exec1 == blocks1, label + ": Energy after AllForce ran Type 1 again");
		p->AllForce(warm);
		Check(Same(cold, warm, dim), label + ": cold AllForce differs from warm AllForce");
		const auto & sp = p->species_constraints;
		Check(sp[0].CCPotential->GetLastRhoBackendUsed() == Backend::FinufftType1 && sp[1].CCPotential->GetLastRhoBackendUsed() == Backend::FinufftType1
		      && sp[2].CCPotential->GetLastRhoBackendUsed() == Backend::Direct && sp[3].CCPotential->GetLastRhoBackendUsed() == Backend::Direct,
		      label + ": unexpected rho backends");
		Check(blocks1 > 2, label + ": the large species did not split its Type 1 work (" + std::to_string(blocks1) + " Type 1 executions)");
		Check(sp[2].CCPotential->GetLastForceBackendUsed() == ForceBackend::Direct && sp[3].CCPotential->GetLastForceBackendUsed() == ForceBackend::Direct,
		      label + ": a Direct species used Type 2");
		Check(sp[0].CCPotential->GetLastForceBackendUsed() == ForceBackend::FinufftType2, label + ": the large species did not use Type 2");
		std::string decided;
		for (const auto & species : sp)
			decided += species.CCPotential->GetLastForceBackendUsed() == ForceBackend::FinufftType2 ? '1' : '0';
		p->SetConfiguration(c);
		const size_t exec1b = obs::executed[1], exec2b = obs::executed[2];
		const double ef = p->EnergyAndAllForce(fused);
		const size_t blocks2 = obs::executed[2] - exec2b;
		Check(obs::executed[1] - exec1b == blocks1, label + ": fused evaluation did not run every Type 1 block once");
		Check(threads == 1 ? blocks2 == 2 : blocks2 > 2, label + ": unexpected Type 2 executions " + std::to_string(blocks2));
		Check(std::memcmp(&e, &ef, sizeof(double)) == 0 && Same(warm, fused, dim), label + ": EnergyAndAllForce differs from Energy + AllForce");
		for (int rep = 0; rep < 2; rep++) {
			p->SetConfiguration(c);
			const double er = p->Energy();
			p->AllForce(again);
			Check(std::memcmp(&e, &er, sizeof(double)) == 0 && Same(warm, again, dim), label + ": repeated evaluation not bitwise identical");
		}
		double scale = 0;
		Check(std::abs(e - ed) <= 1e-10 * std::abs(ed), label + ": energy differs from Direct");
		Check(ForceDiff(warm, fd, dim, scale) <= 1e-10 * std::max(1.0, scale), label + ": forces differ from Direct");
		p->SetConfiguration(moved);
		std::vector<GeometryVector> fm;
		const double em = p->EnergyAndAllForce(fm);
		Check(std::abs(em - edm) <= 1e-10 * std::abs(edm) && ForceDiff(fm, fdm, dim, scale) <= 1e-10 * std::max(1.0, scale), label + ": moved, differs from Direct");
		std::cout << "  " << label << ": " << blocks1 << " Type 1 and " << blocks2 << " Type 2 executions per evaluation, Type 2 species " << decided << '\n';
		if (threads == 1) { e0 = e; f0 = warm; em0 = em; fm0 = fm; decided0 = decided; }
		else if (decided == decided0)	// the same gate decisions: Type 1 blocks are fixed by the inputs, Type 2 blocks change no bit
			Check(std::memcmp(&e, &e0, sizeof(double)) == 0 && Same(warm, f0, dim) && std::memcmp(&em, &em0, sizeof(double)) == 0 && Same(fm, fm0, dim),
			      label + ": results differ from T=1");
	}
}

/** \brief Plan lifecycle with several blocks: particle-count and tolerance changes, copies, assignment, Type 2 counts,
 *  injected makeplan and setpts failures. */
void Lifecycle()
{
	const size_t live0 = obs::live.size();
	{
		const Configuration c = MakeConfig(3, 2000, true, 19);
		Setup s(c, Modes(3, 4));
		ShiftedCCPotential & pot = s.pot;
		pot.ParallelNumThread = 4;
		const size_t made1 = obs::made[1];
		pot.GetRho();
		const size_t blocks = obs::made[1] - made1;
		Check(blocks > 2 && obs::live.size() == live0 + blocks, "lifecycle: several Type 1 plans");
		const std::vector<double> ref = pot.RhoReal;
		// moved particles: plans reused
		Configuration moved(c);
		for (size_t j = 0; j < c.NumParticle(); j++) {
			GeometryVector x = moved.GetRelativeCoordinates(j);
			x.x[1] += 0.01;
			moved.MoveParticle(j, x);
		}
		pot.SetConfiguration(moved);
		pot.GetRho();
		pot.SetConfiguration(c);
		pot.GetRho();
		Check(obs::made[1] == made1 + blocks && pot.RhoReal == ref, "lifecycle: moves made plans or changed rho");
		// fewer particles: fewer blocks, plans trimmed, none made; then back: completed, rho bitwise as before
		Configuration fewer(c);
		while (fewer.NumParticle() > 600) fewer.DeleteParticle(fewer.NumParticle() - 1);
		pot.SetConfiguration(fewer);
		pot.GetRho();
		const size_t trimmed = obs::live.size() - live0;
		Check(trimmed >= 1 && trimmed < blocks && obs::made[1] == made1 + blocks, "lifecycle: fewer particles did not trim the Type 1 plans");
		ShiftedCCPotential check(pot);
		check.SetRhoBackend(Backend::Direct);
		check.SetConfiguration(fewer);
		check.GetRho();
		double diff = 0;
		for (size_t i = 0; i < check.RhoReal.size(); i++)
			diff = std::max(diff, std::abs(Complex(pot.RhoReal[i], pot.RhoImag[i]) - Complex(check.RhoReal[i], check.RhoImag[i])));
		Check(diff <= 1e-10 * 600.0, "lifecycle: fewer particles, rho differs from Direct");
		pot.SetConfiguration(c);
		pot.GetRho();
		Check(obs::live.size() == live0 + blocks && obs::made[1] == made1 + blocks + (blocks - trimmed) && pot.RhoReal == ref,
		      "lifecycle: restored particles did not complete the plans or changed rho");
		// tolerance change: every plan remade
		const size_t made1b = obs::made[1];
		pot.SetFinufftTolerance(1e-11);
		pot.GetRho();
		Check(obs::made[1] - made1b >= 1 && obs::live.size() == live0 + (obs::made[1] - made1b), "lifecycle: tolerance change did not remake the plans");
		pot.SetFinufftTolerance(1e-12);
		pot.GetRho();
		Check(pot.RhoReal == ref, "lifecycle: default tolerance did not restore rho bitwise");
		// copies start empty and make their own plans; assignment releases the target's
		{
			ShiftedCCPotential copy(pot);
			Check(obs::live.size() == live0 + blocks, "lifecycle: copy shares plans");
			copy.RhoValid = false;
			copy.GetRho();
			Check(obs::live.size() == live0 + 2 * blocks && copy.RhoReal == ref, "lifecycle: copy did not make its own plans");
			copy = pot;
			Check(obs::live.size() == live0 + blocks, "lifecycle: assignment did not release the plans");
		}
		Check(obs::live.size() == live0 + blocks, "lifecycle: copy destruction released the source's plans");
		// Type 2 plans follow the requested count; forces equal whatever the count
		std::vector<GeometryVector> f1, f;
		Check(pot.PrepareFinufftForces(1) == 1 && pot.TryFinufftAllForce(f1), "lifecycle: one Type 2 block");
		for (size_t count : {5, 3, 7}) {
			Check(pot.PrepareFinufftForces(count) == count && obs::live.size() == live0 + blocks + count, "lifecycle: Type 2 plans not " + std::to_string(count));
			Check(pot.TryFinufftAllForce(f) && Same(f, f1, 3), "lifecycle: Type 2 forces depend on the block count");
		}
		// a failed setpts in one Type 2 block falls back for that call only
		obs::injectSetpts = 9;
		Check(!pot.TryFinufftAllForce(f) && pot.GetLastForceBackendUsed() == ForceBackend::Direct, "lifecycle: failed Type 2 block did not fall back");
		Check(pot.TryFinufftAllForce(f) && Same(f, f1, 3), "lifecycle: Type 2 not back after a failed block");
		Check(pot.PrepareFinufftForces(0) == 0 && obs::live.size() == live0 + blocks, "lifecycle: declined Type 2 kept plans");
		// a failed setpts in one Type 1 block: Direct rho for that call only, with the injected status
		obs::injectSetpts = 7;
		pot.RhoValid = false;
		pot.GetRho();
		Check(pot.GetLastRhoBackendUsed() == Backend::Direct && pot.GetLastFinufftStatus() == 7, "lifecycle: failed Type 1 block did not fall back");
		pot.RhoValid = false;
		pot.GetRho();
		Check(pot.GetLastRhoBackendUsed() == Backend::FinufftType1 && pot.RhoReal == ref, "lifecycle: Type 1 not back after a failed block");
		// a makeplan failure in the second block: no plan stays, the fallback is cached until an input changes
		obs::failMakeplanAfter = 1;
		pot.SetFinufftTolerance(1e-11);
		pot.GetRho();
		Check(pot.GetLastRhoBackendUsed() == Backend::Direct && pot.GetLastFinufftStatus() == 1 && obs::live.size() == live0,
		      "lifecycle: failed block makeplan left plans or did not fall back");
		const size_t madeAll = obs::made[1];
		pot.RhoValid = false;
		pot.GetRho();
		Check(obs::made[1] == madeAll && pot.GetLastFinufftStatus() == 1, "lifecycle: failed makeplan not cached");
		obs::failMakeplanAfter = -1;
		pot.SetFinufftTolerance(1e-12);
		pot.GetRho();
		Check(pot.GetLastRhoBackendUsed() == Backend::FinufftType1 && pot.RhoReal == ref && obs::live.size() == live0 + blocks, "lifecycle: no recovery after a failed makeplan");
	}
	Check(obs::live.size() == live0, "lifecycle: plans left at scope exit");
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
	Single(1, 20000, 400, "1D");
	Single(2, 6000, 12, "2D skew");
	Single(3, 2000, 5, "3D skew");
	Mixed();
	Lifecycle();
	Check(obs::live.empty(), "plans still live at exit: " + std::to_string(obs::live.size()));
	Check(obs::badDestroys == 0, "double or unknown destroys: " + std::to_string(obs::badDestroys));
	Check(obs::parallelMakeplans == 0, "plans made inside a parallel region: " + std::to_string(obs::parallelMakeplans));
	std::cout << (failures ? "FAILED: " + std::to_string(failures) + " of " + std::to_string(checks) + " check(s)\n"
	                       : "All " + std::to_string(checks) + " block checks passed\n");
	return failures ? 1 : 0;
}
