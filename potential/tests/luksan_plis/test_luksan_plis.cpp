/** \file test_luksan_plis.cpp
 * \brief LuksanPLIS against the NLopt library and against itself (see README.md).
 * Strict mode (fast = 0) must give NLopt's NLOPT_LD_LBFGS bit for bit: every point passed to the objective, the
 * number of evaluations, the return value, the minimizer and the minimum. Fast mode must reach the same minimum,
 * repeat bitwise and not depend on the number of threads. */
#include "LuksanPLIS.h"

#include <nlopt.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

namespace {
/** \brief Record of one optimization: a hash of every point evaluated, and the evaluation count. */
struct Trace {
	std::uint64_t hash = 1469598103934665603ull;	//< FNV-1a over the bits of every x passed to the objective
	int calls = 0;	//< objective evaluations
	int kind = 0;	//< 0: chained Rosenbrock (nonconvex valley), 1: stiff quadratic plus quartic
};

/** \brief Test objective with gradient; mixes every evaluated point into the trace.
 *  \param n     Number of variables.
 *  \param x     Variables.
 *  \param grad  Receives the gradient.
 *  \param data  Trace.
 *  \return Objective value. */
double Objective(unsigned n, const double * x, double * grad, void * data)
{
	Trace & trace = *static_cast<Trace *>(data);
	for (unsigned i = 0; i < n; ++i) {
		std::uint64_t bits;
		std::memcpy(&bits, &x[i], sizeof(bits));
		trace.hash = (trace.hash ^ bits) * 1099511628211ull;
	}
	++trace.calls;
	double f = 0.0;
	for (unsigned i = 0; i < n; ++i)
		grad[i] = 0.0;
	if (trace.kind == 0) {
		for (unsigned i = 0; i + 1 < n; ++i) {
			const double a = x[i + 1] - x[i] * x[i], b = 1.0 - x[i];
			f += 100.0 * a * a + b * b;
			grad[i] += -400.0 * a * x[i] - 2.0 * b;
			grad[i + 1] += 200.0 * a;
		}
	}
	else {
		for (unsigned i = 0; i < n; ++i) {
			const double w = 1.0 + 999.0 * double(i) / double(n), d = x[i] - std::cos(0.37 * i);
			f += 0.5 * w * d * d + 0.25 * d * d * d * d;
			grad[i] = w * d + d * d * d;
		}
	}
	return f;
}

/** \brief One NLopt run with the settings RelaxStructure_NLOPT uses.
 *  \param[in,out] x  Start point; the minimizer on return.
 *  \param trace      Receives the evaluations.
 *  \param history    Number of correction pairs.
 *  \param stopval    Stop value.
 *  \param maxeval    Evaluation limit.
 *  \param[out] minf  Minimum found.
 *  \return NLopt's result code. */
int RunNLopt(std::vector<double> & x, Trace & trace, int history, double stopval, int maxeval, double & minf)
{
	nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, static_cast<unsigned>(x.size()));
	nlopt_set_vector_storage(opt, static_cast<unsigned>(history));
	nlopt_set_min_objective(opt, Objective, &trace);
	nlopt_set_stopval(opt, stopval);
	nlopt_set_maxeval(opt, maxeval);
	const int status = nlopt_optimize(opt, x.data(), &minf);
	nlopt_destroy(opt);
	return status;
}
}

/** \brief Run every case; print one line per case.
 *  \return Number of failed cases. */
int main()
{
	int failures = 0;
	const double stops[2] = {-HUGE_VAL, 1e-6};
	for (int kind = 0; kind < 2; ++kind)
		for (int n : {2, 37, 1000, 30000})
			for (int history : {1, 3, 40})
				for (int maxeval : {25, 4000})
					for (double stopval : stops) {
						std::vector<double> start(n);
						for (int i = 0; i < n; ++i)
							start[i] = kind == 0 ? -1.2 + 0.01 * (i % 7) : 2.0 * std::sin(1.3 * i);
						Trace nlopt, strict, fast1, fast1again, fast4;
						nlopt.kind = strict.kind = fast1.kind = fast1again.kind = fast4.kind = kind;
						std::vector<double> xn(start), xs(start), xf(start), xr(start), xt(start);
						double fn = 0, fs = 0, ff = 0, fr = 0, ft = 0;
						int evals = -1;
						const int sn = RunNLopt(xn, nlopt, history, stopval, maxeval, fn);
						const int ss = LuksanPLIS(n, Objective, &strict, xs.data(), &fs, stopval, maxeval, history, &evals, 0);
						const bool same = sn == ss && nlopt.calls == strict.calls && evals == strict.calls && nlopt.hash == strict.hash
							&& std::memcmp(xn.data(), xs.data(), sizeof(double) * n) == 0 && std::memcmp(&fn, &fs, sizeof(double)) == 0;
						const int sf = LuksanPLIS(n, Objective, &fast1, xf.data(), &ff, stopval, maxeval, history, nullptr, 1);
						LuksanPLIS(n, Objective, &fast1again, xr.data(), &fr, stopval, maxeval, history, nullptr, 1);
						LuksanPLIS(n, Objective, &fast4, xt.data(), &ft, stopval, maxeval, history, nullptr, 4);
						const bool repeat = fast1.hash == fast1again.hash && fast1.hash == fast4.hash && fast1.calls == fast4.calls
							&& std::memcmp(xf.data(), xt.data(), sizeof(double) * n) == 0 && std::memcmp(&ff, &ft, sizeof(double)) == 0;
						// Converged runs (success or stop value) of both modes must agree on the minimum; capped runs need not.
						const bool converged = (sn == 1 || sn == 2) && (sf == 1 || sf == 2);
						const bool minimum = !converged || std::abs(ff - fn) <= 1e-6 * std::max(1.0, std::abs(fn)) + (stopval > 0 ? stopval : 0.0);
						const bool ok = same && repeat && minimum;
						failures += !ok;
						std::cout << (ok ? "ok   " : "FAIL ") << "kind " << kind << " n " << n << " history " << history << " maxeval " << maxeval
							<< " stopval " << stopval << ": NLopt " << sn << " in " << nlopt.calls << " evaluations, f " << fn << "; strict "
							<< (same ? "identical" : "DIFFERENT") << "; fast " << sf << " in " << fast1.calls << ", f " << ff
							<< (repeat ? ", repeats for 1 and 4 threads" : ", NOT repeatable") << (minimum ? "" : ", minimum differs") << '\n';
					}
	// Invalid arguments are refused.
	{
		Trace trace;
		double x = 0.5, f = 0.0;
		const bool ok = LuksanPLIS(0, Objective, &trace, &x, &f, 0.0, 10, 5, nullptr, 0) == -2
			&& LuksanPLIS(1, Objective, &trace, &x, &f, 0.0, 10, 0, nullptr, 0) == -2 && trace.calls == 0;
		failures += !ok;
		std::cout << (ok ? "ok   " : "FAIL ") << "invalid arguments\n";
	}
	std::cout << (failures == 0 ? "ALL PASSED\n" : "FAILURES\n");
	return failures;
}
