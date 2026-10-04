#ifndef __GN_H__
#define __GN_H__

/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: February 2024 */

/** \file gn.h
 * \brief The gn class: n-body correlation function g_n (n = 2, 3, 4) of 1D point
 *        patterns, computed directly from a ConfigurationPack.
 *
 * Estimator. For every particle x_1 of a configuration, the relative positions
 * x_{j1} = x_j - x_1 of all neighbors (periodic images included) are histogrammed
 * into bins of width dx centered at i*dx, |i| <= num_half = floor(Rmax/dx); each
 * bin is fully sampled (neighbors are searched up to (num_half + 1/2)*dx).  With
 * N_i the occupancy of bin i around x_1 (x_1 itself sits in bin 0), g_n at the bin
 * tuple (x_21, ..., x_n1) counts ordered tuples of *distinct* particles:
 * \f[ g_n \approx \frac{1}{N\,\rho^{n-1} dx^{n-1}} \sum_{x_1}
 *     \prod_{k=2}^{n} \big(N_{i_k} - \#\{j<k : i_j = i_k\}\big), \f]
 * with i_1 = 0 and rho = N/L.  This is the bin average of g_n over the box of side
 * dx.  The ensemble mean and its standard error over configurations are accumulated.
 * g_2 is symmetrized over x_21 -> -x_21.
 *
 *
 * Periodic images. Neighbors are all periodic images within (num_half + 1/2)*dx.  If
 * that exceeds L/2, images of the same particle (including x_1's own) count as distinct
 * particles, as in the legacy cell-list search.
 *
 * Storage. Only x_21 >= 0 is kept, so the flat array has (num_half+1) * len_row^{n-2}
 * entries, len_row = 2*num_half+1, laid out as
 * idx = i2*len_row^{n-2} + (i3+num_half)*len_row^{n-3} + ... (row-major).  For n = 4
 * the counts are symmetric under x_31 <-> x_41, so only x_31 <= x_41 is accumulated and
 * the other half is read through the mirrored index.
 * Memory grows like (Rmax/dx)^{n-1}; time like Nc * N * (2*rho*Rmax)^{n-1}.
 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
//-------Core part--------
#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>

/** \brief Ensemble-averaged n-body correlation function g_n of 1D point patterns,
 *  with checkpointing (binary file) so an interrupted run can be resumed. */
struct gn{
	size_t n = 0;				///< order n of g_n (2, 3, or 4)
	double Rmax = 0.;			///< largest bin center |x_j1| [length]
	double dx = 0.;				///< bin width [length]
	int num_half = 0;			///< number of bins on each side of x_1 = floor(Rmax/dx)
	size_t len_row = 0;			///< 2*num_half+1, number of bins along x_31, x_41
	size_t num_threads = 1;		///< OpenMP threads
	size_t curr_c_idx = 0;		///< number of configurations already accumulated
	size_t targ_c_idx = 0;		///< number of configurations to accumulate in total
	double Rmin = 5.0;			///< lower cutoff of the "long-range" tau_n (Compute_tau_n_excluded) [length]
	bool print_gn = true;		///< write the g_n text files at the end of Run
	std::string savename;		///< prefix of the output text files
	std::string bin_name;		///< checkpoint file name (<prefix>.g<n>)

	/** \brief Default constructor (empty, unusable until parameters are set). */
	gn() {};

	/** \brief Interactive constructor: reads n, the checkpoint prefix, Rmax, dx, the
	 *  number of threads, the input ConfigPack, the number of configurations, and the
	 *  output prefix from ifile (prompts go to ofile).  If the checkpoint
	 *  <prefix>.g<n> exists, Rmax/dx are taken from it and the two input tokens are ignored.
	 *  \param ofile  Output stream for prompts.
	 *  \param ifile  Input stream with the parameters.
	 *  \throw std::invalid_argument if n is not 2, 3, or 4, dR <= 0, or Rmax < dR. */
	gn(std::ostream & ofile, std::istream & ifile);

	/** \brief Write the checkpoint (progress, grid, and per-bin ensemble means of g_n and
	 *  g_n^2) in binary; the layout is unchanged from the legacy code.
	 *  \param ofile  Binary output stream. */
	void WriteBinary(std::ostream & ofile ) const;

	/** \brief Read a checkpoint written by WriteBinary and rebuild the derived grid members.
	 *  \param ifile  Binary input stream. */
	void ReadBinary(std::istream & ifile );

	/** \brief Accumulate g_n of one more configuration (multi-threaded; the result does not
	 *  depend on the number of threads).
	 *  \param c  A 1D configuration.
	 *  \throw std::invalid_argument if c is not one-dimensional, is empty, or has a
	 *         non-positive box length. */
	void AddConfig(const Configuration & c);

	/** \brief Write g_n and its standard error as text.
	 *  n = 2, 3: one file <prefix>_g<n>.txt; n = 4: one file per x_21 bin,
	 *  <prefix>_g4_x21_<x21>.txt.
	 *  \param filename  Output prefix. */
	void Write_gn(const std::string & filename);

	/** \brief Set the number of OpenMP threads.
	 *  \param num_threads  Number of threads. */
	void SetNumThreads(size_t num_threads){
		this->num_threads = num_threads;
		omp_set_num_threads(this->num_threads);
	}

	/** \brief Do not write the g_n text files at the end of Run (tau_n is still written). */
	void DisablePrint_gn(){
		this->print_gn = false;
		std::cout << "Disabled printing gn data\n";
	}

	/** \brief Accumulate configurations until done or the time limit is near, save the
	 *  checkpoint, write g_n (if finished), and write tau_n and the long-range tau_n.
	 *  \param ofile      Output stream for progress messages.
	 *  \param TimeLimit  Wall-clock deadline (absolute time_t) [s]. */
	void Run(std::ostream & ofile, time_t TimeLimit);

	/** \brief tau_n(L) = integral of [g_n - 1]^2 over the box max|x_j1| < L, with its error.
	 *  \param[out] result  result[i] = (L = (i+1)*dx, tau_n(L), error), i < num_half. */
	void Compute_tau_n(std::vector<GeometryVector> & result);

	/** \brief Variant of tau_n restricted to tuples whose pair separations all exceed Rmin
	 *  (excludes short-range correlations), as a function of L.
	 *  \param[out] result  (L, tau_n', error) for L >= Rmin. */
	void Compute_tau_n_excluded(std::vector<GeometryVector> & result);

private:
	std::vector<double> sum1;			///< per bin: sum over configurations of g_n
	std::vector<double> sum2;			///< per bin: sum over configurations of g_n^2
	std::vector<std::vector<double>> slab;	///< per-thread scratch for one x_21 slab (n >= 3)
	ConfigurationPack CP;				///< input ensemble

	/** \brief Index of the bin that holds the data of bin i: for n = 4, (i2, i3, i4) with
	 *  i3 > i4 maps to its mirror (i2, i4, i3); otherwise i itself.
	 *  \param i  Flat bin index.
	 *  \return Flat index of the accumulated bin. */
	size_t Canonical(size_t i) const;

	/** \brief Ensemble mean of g_n in bin i and its standard error.
	 *  \param i          Flat bin index.
	 *  \param[out] mean  Mean over the curr_c_idx configurations (0 if none).
	 *  \param[out] err   Standard error of the mean (0 if fewer than two configurations). */
	void Moments(size_t i, double & mean, double & err) const;
};

#endif // __GN_H__
