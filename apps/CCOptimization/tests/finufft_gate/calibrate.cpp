/** \file calibrate.cpp
 * \brief Timings for the cost model of ShiftedCCPotential: the Type 2 force gate (FinufftForceGate, FinufftForceCosts)
 *  and the Type 1 particle blocks (FinufftBlocks, FinufftExecutionCosts).
 *
 * For each case (a cell of unit number density, orthogonal or skew, with N uniformly random particles and a
 * constraint at every reciprocal-lattice vector with |k| < K(chi)) it times, at the default eps and FINUFFT
 * nthreads 1, the minimum over repetitions of:
 *  - one Type 2 force execution with the persistent plan (coefficients, setpts, execute, output copy) on one
 *    thread, whatever the gate would decide;
 *  - the direct per-particle Force() loop on one thread, divided by the number of particle-constraint pairs;
 *  - the direct Force() loop over teams of 4 and 8 threads (wall time, used only by the margin check);
 *  - one Type 1 particle block (angles, setpts, execute: RunRhoBlock) of n1 = N / B1 points, B1 the case's Type 1
 *    block count, on one thread.
 * It prints one line per case on stdout (redirect it to the data area, never into the repository):
 *   skew chi dim N M G type2_us pair_us direct4_us direct8_us n1 type1_us
 * with M the number of constraints and G the FINUFFT mode-grid size. fit_gate.py turns the lines into the
 * constants of CollectiveCoordinatePotential.h. Results are machine-specific (see README.md).
 * Usage (from apps/CCOptimization, after make finufft-gate):
 *   ./calibrate_gate.out                     the default list of 142 cases (several minutes)
 *   ./calibrate_gate.out DIM N CHI [SKEW]    one case */

#ifndef CCPOTENTIAL_USE_FINUFFT
#error "calibrate.cpp needs the FINUFFT backend: build it with make finufft-gate"
#endif

#include <CollectiveCoordinatePotential.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <etc.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

size_t Verbosity = 0;

namespace {

using Clock = std::chrono::steady_clock;
constexpr int Teams[2] = {4, 8};	//< team sizes of the direct-loop wall times

/** \brief ShiftedCCPotential that exposes its FINUFFT mode-grid size. */
struct Probe : ShiftedCCPotential
{
	/** \brief Construct an empty potential.
	 *  \param d  Dimension. */
	explicit Probe(DimensionType d) : ShiftedCCPotential(d) {}

	/** \brief Mode-grid size of the bound configuration.
	 *  \return G, or 0 if FINUFFT is not eligible. */
	long long GridModes() const { return this->FinufftGridModes(); }
};

/** \brief Smallest wall time of repeated calls.
 *  \param call    Work to time.
 *  \param reps    Minimum number of calls.
 *  \param budget  Keep calling, up to 20 times reps, while the summed time is below this [s].
 *  \return        Smallest time of one call [s]. */
template <class F>
double MinTime(F call, int reps, double budget)
{
	double best = 1e300, spent = 0;
	for (int r = 0; r < reps || (spent < budget && r < 20 * reps); r++) {
		const auto t = Clock::now();
		call();
		const double s = std::chrono::duration<double>(Clock::now() - t).count();
		best = std::min(best, s);
		spent += s;
	}
	return best;
}

/** \brief Time one case and print its line (or a comment line if FINUFFT is not eligible).
 *  \param d     Dimension (1-3).
 *  \param n     Particle number.
 *  \param chi   Constraint density: every reciprocal-lattice k with |k| < K(chi).
 *  \param skew  Use a non-orthogonal cell (d >= 2). */
void Case(int d, size_t n, double chi, bool skew)
{
	const double L = std::pow(double(n), 1.0 / d);
	GeometryVector basis[3] = {GeometryVector(d), GeometryVector(d), GeometryVector(d)};
	for (int i = 0; i < d; i++)
		basis[i].x[i] = L;
	if (skew && d >= 2) basis[1].x[0] = 0.31 * L;
	if (skew && d == 3) { basis[2].x[0] = 0.17 * L; basis[2].x[1] = -0.23 * L; }
	Configuration c(d, basis, 0.1 * L);
	RandomGenerator gen(7 + n);
	for (size_t i = 0; i < n; i++) {
		GeometryVector x(d);
		for (int k = 0; k < d; k++) x.x[k] = gen.RandomDouble();
		c.Insert("a", x);
	}
	const double K = 2 * pi * std::pow(2 * d * chi / HyperSphere_Volume(d, 1.0), 1.0 / d);
	Probe p(d);
	for (const auto & k : GetKs(c, K, K, 1))
		p.AddConstraint(k, {1.0});
	p.SetRhoBackend(ShiftedCCPotential::RhoBackend::FinufftType1);
	p.SetConfiguration(c);
	p.ParallelNumThread = 1;
	p.GetRho();
	// one block makes the unsplit Type 2 plan whatever the gate decides; with rho valid, TryFinufftAllForce runs only Type 2
	if (p.GetLastRhoBackendUsed() != ShiftedCCPotential::RhoBackend::FinufftType1 || !p.PrepareFinufft() || !p.PrepareFinufftForces(true)) {
		std::printf("# skipped (FINUFFT not eligible): skew %d chi %g dim %d N %zu\n", int(skew), chi, d, n);
		return;
	}
	std::vector<GeometryVector> f;
	const double type2 = MinTime([&] { p.TryFinufftAllForce(f); }, 15, 0.3);
	p.RhoValid = false;
	const size_t blocks1 = p.StartRho(false);
	const double type1 = MinTime([&] { p.RunRhoBlock(0); }, 15, 0.3);
	p.FinishRho();
	const size_t m = p.constraints.size();
	const size_t sample = std::min<size_t>(n, std::max<size_t>(16, 200000 / std::max<size_t>(m, 1)));
	std::vector<GeometryVector> g(sample);
	const double pair = MinTime([&] { for (size_t i = 0; i < sample; i++) p.Force(g[i], i); }, 5, 0.2) / double(sample * m);
	double wall[2];
	for (int t = 0; t < 2; t++) {
		const int team = Teams[t];
		std::vector<GeometryVector> fa(n);
		wall[t] = MinTime([&] {
#pragma omp parallel for num_threads(team)
			for (long i = 0; i < long(n); i++)
				p.Force(fa[i], i);
		}, 5, 0.3);
	}
	std::printf("%d %g %d %zu %zu %lld %.3f %.6f %.3f %.3f %zu %.3f\n", int(skew), chi, d, n, m, p.GridModes(),
	            1e6 * type2, 1e6 * pair, 1e6 * wall[0], 1e6 * wall[1], n / std::max<size_t>(blocks1, 1), 1e6 * type1);
	std::fflush(stdout);
}

}  // namespace

/** \brief Run the default case list or one case from the command line.
 *  \param argc  Argument count.
 *  \param argv  Nothing, or DIM N CHI [SKEW].
 *  \return      0 on success, 1 on bad arguments. */
int main(int argc, char ** argv)
{
	const int dim = argc >= 4 ? std::atoi(argv[1]) : 0;
	if ((argc != 1 && argc != 4 && argc != 5) || (argc > 1 && (dim < 1 || dim > 3))) {
		std::fprintf(stderr, "usage: calibrate_gate.out [DIM N CHI [SKEW]]\n");
		return 1;
	}
	std::printf("# skew chi dim N M G type2_us pair_us direct%d_us direct%d_us n1 type1_us\n", Teams[0], Teams[1]);
	if (argc > 1) {
		Case(dim, std::stoul(argv[2]), std::atof(argv[3]), argc > 4 && std::atoi(argv[4]) != 0);
		return 0;
	}
	for (int d = 1; d <= 3; d++)
		for (size_t n : {100, 300, 1000, 3000, 10000})
			for (double chi : {0.02, 0.05, 0.1, 0.2, 0.4, 0.8})
				for (bool skew : {false, true})
					if ((!skew || d >= 2) && !(d == 3 && n == 10000 && chi > 0.05))	// 3D N=1e4 dense direct loops take minutes
						Case(d, n, chi, skew);
	return 0;
}
