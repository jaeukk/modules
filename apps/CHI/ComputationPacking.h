/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: March 2022 */

/** \file ComputationPacking.h 
 * \brief Header file for a ComputationPacking class.
 * It simplified the input and output for some calculations involving an ensemble of packing configurations.
 * Computation PairStatisticsCLI code written by Ge.*/


#ifndef COMPUTATIONPACKING_H__
#define COMPUTATIONPACKING_H__

#include <set>
#include <cmath>
#include <list>
#include <vector>
#include <omp.h>
#include <string.h>

// ------ core ------
#include <etc.h>
#include <PeriodicCellList.h>
#include <GeometryVector.h>

// ------ in the directory ----
//#include <StructureFactor.h>

//#include "NumberVariance.h"
#include <LinealPathFunction.h>
#include <SpectralDensity.h>

//TODO: implement functions for the dipsersion class
//TODO: implement various times of calculations (S2, local volume fraction variance, and so on.)

class ComputationPacking
{
public:
	int numThreads = 1;

	/** \brief Run the computation over an ensemble of configurations. (CLAUDE)
	 *  \param GetConfigsFunction  Functor returning the i-th SpherePacking configuration.
	 *  \param NumConfig           Number of configurations to process. */
	virtual void Compute(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfig) =0;

	/** \brief Write the computed result to output files. (CLAUDE)
	 *  \param OutputPrefix  Prefix used to build the output file name(s). */
	virtual void Write(const std::string OutputPrefix) =0;
	//virtual void Plot(const std::string OutputPrefix, const std::string & Title) =0;

	//some computations allow additional options
	//any command not understood by PairStatisticsCLI() is passed to here
	/** \brief Handle an extra command not recognized by PairStatisticsCLI(). (CLAUDE)
	 *  \param option  The unrecognized command token.
	 *  \param input   Stream to read any further parameters from.
	 *  \param output  Stream to report messages to. */
	virtual void ProcessAdditionalOption(const std::string &, std::istream &, std::ostream & output)
	{
		//by default, accept no additional option.
		output << "Unrecognized command!\n";
	}

	/** \brief Virtual destructor. (CLAUDE) */
	virtual ~ComputationPacking()
	{
	}

	/** \brief Set the number of OpenMP threads to use. (CLAUDE)
	 *  \param n  Number of threads. */
	void SetNumThreads(int n){
		this->numThreads = n;
	}
};

class chivComputation : public ComputationPacking{
public:
	std::vector<GeometryVector> result;
	double CircularKMax, LinearKMax, KPrecision;

	/** \brief Construct by reading the spectral-density parameters. (CLAUDE)
	 *  \param ifile  Input stream supplying CircularKMax, LinearKMax and KPrecision.
	 *  \param ofile  Output stream echoing the parsed parameters. */
	chivComputation(std::istream & ifile, std::ostream & ofile) {
		ofile << "circular K_max=";
		ifile >> CircularKMax; ofile << CircularKMax << "\n";
		ofile << "linear K_max=";
		ifile >> LinearKMax; ofile << LinearKMax << "\n";
		ofile << "K precision (binning width)=";
		ifile >> KPrecision;	ofile << KPrecision << "\n";
		//SampleProbability = 1.0;
	}

	/** \brief Compute the isotropic spectral density over the ensemble. (CLAUDE)
	 *  \param GetConfig   Functor returning the i-th SpherePacking configuration.
	 *  \param numConfigs  Number of configurations to process. */
	virtual void Compute(std::function<const SpherePacking(size_t i)> GetConfig, size_t numConfigs) {
		result.clear();
		omp_set_num_threads(this->numThreads);
		IsotropicSpectralDensity(GetConfig, numConfigs, CircularKMax, LinearKMax, result, KPrecision);
	}

	/** \brief Write the spectral density to "<OutputPrefix>_chiv". (CLAUDE)
	 *  \param OutputPrefix  Prefix used to build the output file name. */
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_chiv");
		WriteFunction(result, name);
	}
};

//Compute the loca volume-fraction variance of particle phase
// class LocalVolumeFractionVariance : public ComputationPacking {
// public:
// 	int numSamplingPts;
// 	double Rmax, dR;
// 	std::vector<GeometryVector> result;

// 	LocalVolumeFractionVariance(std::istream & ifile, std::ostream & ofile) {
// 		ofile << "Rmax = ";
// 		ifile >> Rmax;	ofile << "\t" << Rmax << "\n";
// 		ofile << "dR = ";
// 		ifile >> dR;	ofile << "\t" << dR << "\n";
// 		ofile << "the number of sampling windows = ";
// 		ifile >> numSamplingPts;	ofile << "\t" << numSamplingPts << "\n";
// 	}
// 	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
// 		omp_set_num_threads(this->numThreads);
// 		MC_vfv(GetConfigs, numConfigs, result, Rmax, dR, numSamplingPts);
// 	}
// 	virtual void Write(const std::string OutputPrefix) {
// 		std::string name = OutputPrefix + std::string("_localVolVar");
// 		WriteFunction(result, name);
// 	}
// };


//Compute the Hole Probability Function for a given phase.
class HoleProbability : public ComputationPacking {
public:
	int numSamplingPts;
	int phase=0;
	double Rmax, dR;
	std::vector<GeometryVector> result;

	/** \brief Construct by reading the hole-probability parameters. (CLAUDE)
	 *  \param ifile  Input stream supplying Rmax, dR, numSamplingPts and phase.
	 *  \param ofile  Output stream echoing the parsed parameters. */
	HoleProbability(std::istream & ifile, std::ostream & ofile) {
		ofile << "Rmax = ";
		ifile >> Rmax;	ofile << "\t" << Rmax << "\n";
		ofile << "dR = ";
		ifile >> dR;	ofile << "\t" << dR << "\n";
		ofile << "the number of sampling windows = ";
		ifile >> numSamplingPts;	ofile << "\t" << numSamplingPts << "\n";
		ofile << "Choice of phase (0=matrix/1=particle) = ";
		ifile >> phase;
		ofile << "\t" << phase << "\n";
	}

	/** \brief Compute the hole probability function over the ensemble (not yet implemented). (CLAUDE)
	 *  \param GetConfigs  Functor returning the i-th SpherePacking configuration.
	 *  \param numConfigs  Number of configurations to process. */
	virtual void Compute(std::function<const SpherePacking(size_t)>, size_t) {
		omp_set_num_threads(this->numThreads);
		// ToDo

	}

	/** \brief Write the hole probability function to a phase-specific file. (CLAUDE)
	 *  \param OutputPrefix  Prefix used to build the output file name. */
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_Ev_");
		if (this->phase == 0)
			name += std::string("matrix");
		else if (this->phase == 1)
			name += std::string("particle");
		WriteFunction(result, name);
	}
};


//Compute the lineal path Function for a given phase.
class LinealPath : public ComputationPacking {
public:
	int numSamplingPts;
	int rng_seed=0;
	double Rmax, dR;
	std::vector<GeometryVector> result0;
	std::vector<GeometryVector> result1;

	/** \brief Construct by reading the lineal-path parameters. (CLAUDE)
	 *  \param ifile  Input stream supplying Rmax, dR, numSamplingPts and rng_seed.
	 *  \param ofile  Output stream echoing the parsed parameters. */
	LinealPath(std::istream & ifile, std::ostream & ofile) {
		ofile << "Rmax = ";
		ifile >> Rmax;	ofile << "\t" << Rmax << "\n";
		ofile << "dR = ";
		ifile >> dR;	ofile << "\t" << dR << "\n";
		ofile << "the number of sampling windows = ";
		ifile >> numSamplingPts;	ofile << "\t" << numSamplingPts << "\n";
		ofile << "RNG seed = ";
		ifile >> rng_seed;
		ofile << "\t" << rng_seed << "\n";
	}

	/** \brief Compute the lineal-path function for both phases over the ensemble. (CLAUDE)
	 *  \param GetConfigs  Functor returning the i-th SpherePacking configuration.
	 *  \param numConfigs  Number of configurations to process. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		omp_set_num_threads(this->numThreads);
		LinealPathFunction(GetConfigs, numConfigs, result0, result1, this->Rmax, this->dR, this->numSamplingPts, this->rng_seed);
	}

	/** \brief Write the matrix- and particle-phase lineal-path functions to files. (CLAUDE)
	 *  \param OutputPrefix  Prefix used to build the output file names. */
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_Lz_");
		std::vector<std::string> comment1 {"Lineal-path function via MC method\n", 
			"Line segments inside the matrix phase.\n", 
			"number of sampling lines = " + std::to_string(this->numSamplingPts) + "\n", 
			"z\tL(z)\terror"};
		std::vector<std::string> comment2 {"Lineal-path function via MC method\n", 
			"Line segments inside the particle phase.\n", 
			"number of sampling lines = " + std::to_string(this->numSamplingPts)+"\n", 
			"z\tL(z)\terror"};
		WriteFunction(result0, (name + std::string("matrix")), comment1);
		WriteFunction(result1, (name + std::string("particle")), comment2);

	}
};


//Compute the chord-length distribution.
class ChordLength : public ComputationPacking {
public:
	int numSamplingPts;
	int rng_seed=0;
	double Rmax, dR;
	std::vector<GeometryVector> result0;
	std::vector<GeometryVector> result1;

	/** \brief Construct by reading the chord-length parameters. (CLAUDE)
	 *  \param ifile  Input stream supplying Rmax, dR, numSamplingPts and rng_seed.
	 *  \param ofile  Output stream echoing the parsed parameters. */
	ChordLength(std::istream & ifile, std::ostream & ofile) {
		ofile << "Rmax = ";
		ifile >> Rmax;	ofile << "\t" << Rmax << "\n";
		ofile << "dR = ";
		ifile >> dR;	ofile << "\t" << dR << "\n";
		ofile << "the number of sampling windows = ";
		ifile >> numSamplingPts;	ofile << "\t" << numSamplingPts << "\n";
		ofile << "RNG seed = ";
		ifile >> rng_seed;
		ofile << "\t" << rng_seed << "\n";
	}

	/** \brief Compute the chord-length distribution for both phases over the ensemble. (CLAUDE)
	 *  \param GetConfigs  Functor returning the i-th SpherePacking configuration.
	 *  \param numConfigs  Number of configurations to process. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		omp_set_num_threads(this->numThreads);
		ChordLengthDistribution(GetConfigs, numConfigs, result0, result1, this->Rmax, this->dR, this->numSamplingPts, this->rng_seed);
	}

	/** \brief Write the matrix- and particle-phase chord-length distributions to files. (CLAUDE)
	 *  \param OutputPrefix  Prefix used to build the output file names. */
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_pz_");
		std::vector<std::string> comment1 {"Chord-length distribution function via MC method\n", 
			"Chords inside the matrix phase.\n", 
			"number of sampling lines = " + std::to_string(this->numSamplingPts) +"\n", 
			"z\tp(z)\terror"};
		std::vector<std::string> comment2 {"Chord-length distribution function via MC method\n", 
			"Chords inside the particle phase.\n", 
			"number of sampling lines = " + std::to_string(this->numSamplingPts)+"\n", 
			"z\tp(z)\terror"};

		WriteFunction(result0, (name + std::string("matrix")), comment1);
		WriteFunction(result1, (name + std::string("particle")), comment2);
	}
};



#endif // COMPUTATIONPACKING_H__
