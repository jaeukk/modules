/** \file LuksanPLIS.h
 * \brief Limited-memory BFGS (Luksan's PLIS, the algorithm of NLopt's NLOPT_LD_LBFGS) with a ring-buffer history. */

#ifndef LUKSANPLIS_INCLUDED
#define LUKSANPLIS_INCLUDED

/** \brief Objective with the signature of NLopt's nlopt_func.
 *  \param n         Number of variables.
 *  \param x         Variables.
 *  \param gradient  Receives the gradient (never null here).
 *  \param data      The caller's pointer. */
typedef double (*LuksanObjective)(unsigned n, const double *x, double *gradient, void *data);

/** \brief Minimize f by Luksan's PLIS, as NLopt's NLOPT_LD_LBFGS does for unbounded variables with only a stop
 *  value and an evaluation limit set (NLopt's defaults otherwise: relative x tolerance 1e-16 and relative f
 *  tolerance 1e-14, each on two consecutive iterations, and gradient maximum norm 1e-8). The history of mf
 *  correction pairs is a ring buffer instead of NLopt's two matrices that are shifted by one column in every
 *  iteration. With fast = 0 every arithmetic operation keeps NLopt's order, so the iterates are bitwise NLopt's.
 *  With fast != 0 the dot products of the two Strang recurrences, the bulk of the optimizer's own work, are summed
 *  block by block in eight interleaved partial sums, and with fast > 1 a team of that many threads runs them:
 *  deterministic and independent of the number of threads, but the iterates then differ from NLopt's at rounding
 *  level.
 *  \param n            Number of variables (>= 1).
 *  \param f            Objective and gradient.
 *  \param f_data       Passed to f.
 *  \param[in,out] x    Initial guess; the last iterate on return.
 *  \param[out] minf    Objective at x.
 *  \param stopval      Stop when the objective is at or below this value.
 *  \param maxeval      Evaluation limit (<= 0: none).
 *  \param mf           Number of correction pairs (>= 1).
 *  \param[out] evaluations  Number of objective evaluations, if not null.
 *  \param fast         0 for NLopt's arithmetic; otherwise the reordered dot products on this many threads (1: serial).
 *  \return NLopt's nlopt_result value: 1 success (gradient tolerance), 2 stop value, 3 f tolerance, 4 x tolerance,
 *          5 evaluation limit, -1 failure, -2 invalid arguments, -3 out of memory. */
int LuksanPLIS(int n, LuksanObjective f, void *f_data, double *x, double *minf, double stopval, int maxeval, int mf,
	int *evaluations, int fast);

#endif
