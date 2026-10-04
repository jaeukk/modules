/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: March 2022 */

/**
 * \file Computation.h
 * \brief Header file for a Computation class.
 * It simplified the input and output for some calculations involving an ensemble of point configurations.
 * Computation PairStatisticsCLI code written by Ge.
 */


#ifndef COMPUTATION_H__
#define COMPUTATION_H__

#include <set>
#include <cmath>
#include <list>
#include <vector>
#include <omp.h>

/* header files in $(cores) */
#include <etc.h>
#include <PeriodicCellList.h>
#include <GeometryVector.h>

#include "WeightedStructureFactor.h"
#include "PairCorrelation.h"
#include "LocalVariance.h"

typedef PeriodicCellList<std::complex<double>> MarkedPointPattern;



/**
 * \brief Abstract base for a statistic computed over an ensemble of marked point configurations.
 * Each derived class implements one pair statistic; the CLI driver creates the objects, calls
 * Compute(), then Write() (and optionally Plot()).
 */
class Computation
{
public:
	int num_threads = 1;	///< Number of OpenMP threads used by Compute().

	/**
	 * \brief Compute the statistic over the ensemble.
	 * \param GetConfigsFunction	A lambda returning the i-th marked configuration.
	 * \param NumConfig			The number of configurations to process.
	 */
	virtual void Compute(std::function<const MarkedPointPattern (size_t i)> GetConfigsFunction, size_t NumConfig) =0;

	/** \brief Write the computed result to files whose names begin with \p OutputPrefix. */
	virtual void Write(const std::string OutputPrefix) =0;

	/** \brief Write a Grace (.agr) plot of the result with the given \p Title. */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title) =0;

	/**
	 * \brief Handle a CLI command not understood by PairStatisticsCLI().
	 * Any such command is forwarded here so a computation can expose extra options.
	 * \param option	The unrecognized command token.
	 * \param input	The CLI input stream (for reading option arguments).
	 * \param output	The CLI output stream (for messages).
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output)
	{
		//by default, accept no additional option.
		output << "Unrecognized command!\n";
	}

	virtual ~Computation()
	{
	}

	/** \brief Set the number of OpenMP threads used by Compute(). */
	void SetNumThreads(int n){
		this->num_threads = n;
	}
};

/** \brief Computes the weighted pair correlation function g2 as a function of radius r. */
class g2Computation : public Computation
{
public:
	std::vector<GeometryVector> result;	///< (r, Re g2, Im g2, err) rows produced by Compute().
	double g2rmax, dr;					///< Largest radius and bin width read from the CLI.
	double resolution;
	bool single_sided = false;			///< If true, compute the single-sided variant instead of g2.

	/** \brief Read g2 R_max and dr from the CLI streams. */
	g2Computation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"g2 R_max=";
		ifile>>g2rmax;
		ofile<<"dr = ";
		ifile>>dr;
//		resolution=1.0;
	}

	/** \brief Compute g2 (or the single-sided variant) over the ensemble. */
	virtual void Compute(std::function<const MarkedPointPattern(size_t i)> GetConfigsFunction, size_t NumConfig);

	/** \brief Write the result to "<prefix>_s2" (or "<prefix>_single-sided-s2"). */
	virtual void Write(const std::string OutputPrefix);

	/** \brief Write a Grace plot of g2(r). */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/** \brief Handle the ResolutionPreference, ManualBins, and SingleSided options. */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	virtual ~g2Computation()
	{
	}
};


/** \brief Computes the weighted structure factor as a function of the wavenumber |k|. */
class SkComputation : public Computation
{
public:
	std::vector<GeometryVector> result;					///< (|k|, S(|k|), err.|k|, err.S) rows produced by Compute().
	double CircularKMax, LinearKMax, KPrecision, CircularKMin = 0.0;	///< |k| range and bin width read from the CLI.
	double SampleProbability;							///< Fraction of wavevectors sampled above CircularKMax.
	size_t average_option = 0;							///< 0 = over counts, 1 = over reciprocal-space volume, 2 = Bragg peaks.

	/** \brief Read circular K_max, linear K_max, and the binning width from the CLI streams. */
	SkComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"circular K_max=";
		ifile>>CircularKMax;
		ofile<<"linear K_max=";
		ifile>>LinearKMax;
		ofile<<"K precision (binning width)=";
		ifile>>KPrecision;

		SampleProbability = 1.0;
	}

	/** \brief Compute the angularly-averaged weighted structure factor over the ensemble. */
	virtual void Compute(std::function<const MarkedPointPattern(size_t i)> GetConfigsFunction, size_t NumConfig);

	/** \brief Write the result to "<prefix>_chif". */
	virtual void Write(const std::string OutputPrefix);

	/** \brief Write a Grace plot of S(k). */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/** \brief Handle the SampleProbability, AverageOver*, BraggPeaks, and CircularKMin options. */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


/** \brief Computes the local variance of the particle weights as a function of the window radius R. */
class LocalVariance : public Computation
{
public:
	std::vector<GeometryVector> result;	///< (R, variance, SE, 0) rows produced by Compute().
	std::vector<GeometryVector> mean_NR;	///< Mean weighted count per window: (R-tagged) Re/Im mean and their errors.
	double Rmax, dR, Rmin = -1.0;		///< Largest window radius, radial increment, and inner exclusion radius.
	size_t num_samp_centers, seed;		///< Number of random window centers and the RNG seed (from the CLI).
	RandomGenerator rng;
	bool use_prt_centers = false;		///< If true, use the particle positions as window centers instead of random ones.

	/** \brief Read Rmax, dR, the number of random sampling centers, and the random seed from the CLI streams. */
	LocalVariance(std::istream & ifile, std::ostream & ofile)
	{
		ofile << "Rmax = ";
		ifile >> Rmax;
		ofile << "dR = ";
		ifile >> dR;
		ofile << "number of random sampling centers = ";
		ifile >> num_samp_centers;
		ofile << "random seed = ";
		ifile >> seed;
		rng.seed(seed);
	}

	/** \brief Compute the local variance of the weights over the ensemble. */
	virtual void Compute(std::function<const MarkedPointPattern(size_t i)> GetConfigsFunction, size_t NumConfig);

	/** \brief Write the variance to "<prefix>_varf[_self]" and the means to "<prefix>_means[_self].txt". */
	virtual void Write(const std::string OutputPrefix);

	/** \brief Write a Grace plot of the local variance vs R. */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/** \brief Handle the UseParticleCenters and Rmin options. */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


#endif // COMPUTATION_H__