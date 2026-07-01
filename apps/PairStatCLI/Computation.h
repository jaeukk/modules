/**
 *	Author	: Ge Zhang
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

#include "PairCorrelation.h"
#include "StructureFactor.h"
#include "NumberVariance.h"
#include "ContactStatistics.h"

/**
 * @brief Base Computation class to interface computating codes to PairStatCLI. 
 */
class Computation
{
public:
	int num_threads = 1;

	/**
	 * \brief Run the computation over an ensemble of configurations.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig) =0;

	/**
	 * \brief Write the computed result to a file with the given output prefix.
	 * \param OutputPrefix  Prefix of the output file name.
	 */
	virtual void Write(const std::string OutputPrefix) =0;

	/**
	 * \brief Plot the computed result in Grace format.
	 * \param OutputPrefix  Prefix of the output file name.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title) =0;

	/**
	 * \brief Process an additional, computation-specific command. By default no additional option is accepted.
	 * Some computations allow additional options, 
	 * any command not understood by PairStatisticsCLI() is passed to here.
	 * \param option  The command keyword.
	 * \param input   Stream from which the command's arguments are read.
	 * \param output  Stream to which messages are written.
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output)
	{
		//by default, accept no additional option.
		output << "Unrecognized command!\n";
	}

	/** \brief Destructor. */
	virtual ~Computation()
	{
	}

	/**
	 * \brief Set the number of OpenMP threads used by the computation.
	 * \param n  The number of threads.
	 */
	void SetNumThreads(int n){
		this->num_threads = n;
	}
};

/* class to compute pair correlation functions as a function of radius $r$ */
class g2Computation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double g2rmax;
	double resolution;
	HistogramGenerator * pHGen;

	/**
	 * \brief Construct a g2 computation class, reading the maximal radius rmax from the input stream.
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
	g2Computation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"g2 R_max=";
		ifile>>g2rmax;
		resolution=1.0;
		pHGen = nullptr;
	}

	/**
	 * \brief Compute the isotropic g2(r) over the ensemble of configurations.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	
	/**
	 * \brief Write the g2(r) result to "<OutputPrefix>_gofr".
	 * \param OutputPrefix  Prefix of the output file name.
	 */
	virtual void Write(const std::string OutputPrefix);

	/**
	 * \brief Plot the g2(r) result in Grace format.
	 * \param OutputPrefix  Prefix of the output file name.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/**
	 * \brief Process g2-specific options (ResolutionPreference, ManualBins).
	 * \param option  The command keyword.
	 * \param input   Stream from which the command's arguments are read.
	 * \param output  Stream to which messages are written.
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	/** \brief Destructor; releases the optional histogram generator.*/
	virtual ~g2Computation()
	{
		if (pHGen != nullptr)
			delete pHGen;
	}
};

/* class to compute pair correlation functions as a function of a displacement vector */
class Directionalg2Computation : public Computation
{
public:
	double g2rmax;
	size_t NumThetaBins;
	std::vector< std::vector<GeometryVector> > results;

	/**
	 * \brief Construct a directional g2 computation, reading g2rmax and the number of angular bins.
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
	Directionalg2Computation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"g2 R_max=";
		ifile>>g2rmax;
		ofile<<"Num Theta Bins=";
		ifile>>NumThetaBins;
	}

	/**
	 * \brief Compute the directional g2(r, theta) over the ensemble of configurations.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	
	/**
	 * \brief Write the directional g2 result to "<OutputPrefix>_x/_y/_g2" text files.
	 * \param OutputPrefix  Prefix of the output file names.
	 */
	virtual void Write(const std::string OutputPrefix);

	/**
	 * \brief Plotting is not implemented for the directional g2 computation.
	 * \param OutputPrefix  Prefix of the output file name (unused).
	 * \param Title         Title of the plot (unused).
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title)
	{
		std::cerr<<"Directionalg2Computation::Plot is not implemented!\n";
	}
};

/* a class to compute the structure factor as a function of a wavenumber */
class SkComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double CircularKMax, LinearKMax, KPrecision, CircularKMin = 0.0;
	double SampleProbability;
	size_t average_option = 0;
	
	/**
	 * \brief Construct an S(k) computation, reading the circular and linear K maxima and the binning width.
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
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

	/**
	 * \brief Compute the isotropic structure factor S(k) over the ensemble of configurations.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	
	/**
	 * \brief Write the S(k) result to "<OutputPrefix>_Sofk".
	 * \param OutputPrefix  Prefix of the output file name.
	 */
	virtual void Write(const std::string OutputPrefix);
	
	/**
	 * \brief Plot the S(k) result in Grace format.
	 * \param OutputPrefix  Prefix of the output file name.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/**
	 * \brief Process S(k)-specific options (SampleProbability, AverageOver*, BraggPeaks, CircularKMin).
	 * \param option  The command keyword.
	 * \param input   Stream from which the command's arguments are read.
	 * \param output  Stream to which messages are written.
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};



/* a class to compute the nearest-neighbor distribution as a function of a distance */
class NearestNeighbors : public Computation
{
public:
	std::vector<GeometryVector> distribution, min_distribution;
	std::vector<double> minPairs;
	size_t SampleDistanceSize;
	double resolution;
	bool logscale_bin = false;
	bool Output_AllMinPairDistances = false;
	
	/**
	 * \brief Construct a nearest-neighbor statistics computation, reading the sample distance size used to build bins.
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
	NearestNeighbors(std::istream & ifile, std::ostream & ofile)
	{
		ofile << "Sample distance sizes (to generate bins) ";
		ifile >> SampleDistanceSize;
		resolution=1.0;
	}

	/**
	 * \brief Compute the nearest-neighbor distribution and the distribution of minimal nearest-neighbor distances.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);

	/**
	 * \brief Write the nearest-neighbor results (and optionally all minimal pair distances) to file.
	 * \param OutputPrefix  Prefix of the output file names.
	 */
	virtual void Write(const std::string OutputPrefix);

	/**
	 * \brief Plot the nearest-neighbor distributions in Grace format.
	 * \param OutputPrefix  Prefix of the output file names.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/**
	 * \brief Process nearest-neighbor-specific options (ResolutionPreference, LogScale, PrintAll).
	 * \param option  The command keyword.
	 * \param input   Stream from which the command's arguments are read.
	 * \param output  Stream to which messages are written.
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


/* class to compute contact numbers */
class CoordinationNumberComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double rmin, rmax;
	int bins;
	
	/**
	 * \brief Construct a coordination-number computation, reading rmin, rmax, and the number of bins.
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
	CoordinationNumberComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"Coordination numbers between rmin and rmax = ";
		ifile>>rmin;
		ifile>>rmax;
		ofile<<"Number of bins = ";
		ifile>>bins;
	}

	/**
	 * \brief Compute the distribution of coordination numbers between rmin and rmax.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);

	/**
	 * \brief Write the coordination-number distribution to "<OutputPrefix>_Z".
	 * \param OutputPrefix  Prefix of the output file name.
	 */
	virtual void Write(const std::string OutputPrefix);

	/**
	 * \brief Plot the coordination-number distribution in Grace format.
	 * \param OutputPrefix  Prefix of the output file name.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/**
	 * \brief Process additional options; this computation accepts none.
	 * \param option  The command keyword.
	 * \param input   Stream from which the command's arguments are read.
	 * \param output  Stream to which messages are written.
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	/** \brief Destructor.*/
	virtual ~CoordinationNumberComputation()
	{}
};

/** Class to compute the cumulative coordination number as a function of a separation */
class CumulativeCoordinationNumberComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double rmin, rmax, dr;
	bool logscale = false;
	
	/**
	 * \brief Construct a cumulative coordination-number computation, reading rmin, rmax, and dr.
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
	CumulativeCoordinationNumberComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"Cumulative Coordination numbers between rmin and rmax = ";
		ifile>>rmin;
		ifile>>rmax;
		ofile<<"dr = ";
		ifile>>dr;
	}

	/**
	 * \brief Compute the cumulative coordination number Z(r) between rmin and rmax.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);

	/**
	 * \brief Write the cumulative coordination number Z(r) to "<OutputPrefix>_Zofr".
	 * \param OutputPrefix  Prefix of the output file name.
	 */
	virtual void Write(const std::string OutputPrefix);

	/**
	 * \brief Plot the cumulative coordination number Z(r) in Grace format.
	 * \param OutputPrefix  Prefix of the output file name.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/**
	 * \brief Process cumulative coordination-number options (logscale).
	 * \param option  The command keyword.
	 * \param input   Stream from which the command's arguments are read.
	 * \param output  Stream to which messages are written.
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	/** \brief Destructor.*/
	virtual ~CumulativeCoordinationNumberComputation()
	{}
};


/* a class to compute the local number variance as a function of a wavenumber */
class LocalNumberVariance : public Computation
{
public:
	std::vector<GeometryVector> result;
	double Rmax, dR;
	size_t num_samp_centers, seed;
	RandomGenerator rng;
	bool use_prt_centers = false;

	/**
	 * \brief Construct a local number variance computation, reading Rmax, dR, the number of random sampling centers, and the random seed.
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
	LocalNumberVariance(std::istream & ifile, std::ostream & ofile)
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

	/**
	 * \brief Compute the local number variance via Monte Carlo sampling over the ensemble.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);

	/**
	 * \brief Write the local number variance to "<OutputPrefix>_nv" (or "_nv_self" when particle centers are used).
	 * \param OutputPrefix  Prefix of the output file name.
	 */
	virtual void Write(const std::string OutputPrefix);

	/**
	 * \brief Plot the local number variance in Grace format.
	 * \param OutputPrefix  Prefix of the output file name.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);

	/**
	 * \brief Process local number variance options (UseParticleCenters).
	 * \param option  The command keyword.
	 * \param input   Stream from which the command's arguments are read.
	 * \param output  Stream to which messages are written.
	 */
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


/* Class to calculate average cluster size for a bunch of radii, until a particle is connected to its periodic image */
class AverageClusterSizeComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	size_t SampleSize;
	int seed;
	std::string prefix;
	double PercolationRadius;//a particle has been found to connect to its periodic image at this radius, therefore calculation for higher radii is futile

	/**
	 * \brief Construct an average cluster size computation, reading the min/increment/max radii, sample size, and random seed. 
	 * \param ifile  Input stream supplying parameters.
	 * \param ofile  Output stream for prompts.
	 */
	AverageClusterSizeComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile << "Enter min, increment, and radius:";
		double min, inc, max;
		ifile >> min;
		ifile >> inc;
		ifile >> max;
		ofile << "Enter sample size per configuration per radius:";
		ifile >> SampleSize;
		ofile << "Enter Random Seed:";
		ifile >> seed;
		for (double x = min; x<max; x += inc)
		{
			result.push_back(GeometryVector(x, 0.0));
		}
		std::stringstream ss;
		ss << "_AverageClusterSize";
		prefix = ss.str();
		PercolationRadius = ::MaxDistance;
	}

	/**
	 * \brief Compute the average cluster size at each radius until a particle connects to its periodic image.
	 * \param GetConfigsFunction  A lambda function returning the configuration of index i.
	 * \param NumConfig           The number of configurations to process.
	 */
	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);

	/**
	 * \brief Write the average cluster size result to "<OutputPrefix>_AverageClusterSize".
	 * \param OutputPrefix  Prefix of the output file name.
	 */
	virtual void Write(const std::string OutputPrefix);

	/**
	 * \brief Plot the average cluster size result in Grace format.
	 * \param OutputPrefix  Prefix of the output file name.
	 * \param Title         Title of the plot.
	 */
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
};

#endif // COMPUTATION_H__