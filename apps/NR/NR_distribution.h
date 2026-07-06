/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	August, 2020 */

 /** \file NR_distribution.h
  *	\brief Header file to compute the distribution of local number from Configuration objects. */


#ifndef __NR_DISTRIBUTION_H__
#define __NR_DISTRIBUTION_H__

#include <vector>
#include <omp.h>
#include <deque>		// For LocalNumber class
#include <memory>		// For uniuqe_ptr
#include <functional>
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>	//For Gaussian distribution.

#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <etc.h>

/** \brief Write a table of GeometryVector data to an output stream, with header fields. (CLAUDE)
 *  \param result  The table of GeometryVector rows to write.
 *  \param ofile   The output stream to write into.
 *  \param fields  Header/comment lines prepended to the output. */
void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile, const std::vector<std::string> & fields);

/** \brief Write a table of GeometryVector data to a file named by a prefix, with header fields. (CLAUDE)
 *  \param result         The table of GeometryVector rows to write.
 *  \param OutFilePrefix  Prefix of the output file name.
 *  \param fields         Header/comment lines prepended to the output. */
void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::vector<std::string> & fields);


/** \brief A class to store the empirical distribution of an integer random varianble.
 *	This class is designed to efficiently store raw data of the empirical distribution of integral random numbers,
 extract some estimators and empirical distribution. */
class LocalNumber {
protected:
	// Variables 
	int min,		//< The smallest value of N(R).
		max;		//< The largest value of N(R).
	size_t	totalCount;	//<	Total number of counts.
	std::deque<size_t> counts;	//< \brief Container for the frequency of N(R).
								//< counts[i] =  the number of observations that N(R) = i + min.
								//< deque<T> cases is used for efficiently modifying containers.

	mutable std::vector<double> CDF;	//< cumulative distribution function obtained from empirical distribution.
	mutable bool CDF_uptodate;			//< True = this->CDF is up to date.



	/** \brief Give an index of random number according to the empirical distribution.
	 *	\param rng	A random number generator.
	 *	\return	An index for the random number. Random number = index + min. */
	size_t GetRandomIndex(RandomGenerator & rng) const;
public:
	/** \brief A default constructor. */
	LocalNumber() : min(0), max(0), totalCount(0), CDF_uptodate(false) {};

	/** \brief An empty constructor for a given range.
	 * \param min	The lower bound on data.
	 * \param max	The upper bound on data. */
	LocalNumber(int min, int max) {
		this->min = min;	this->max = max;
		totalCount = 0;
		counts.resize(max - min + 1, 0);
		CDF_uptodate = false;
	}

	/** \brief A copy constructor. */
	LocalNumber(const LocalNumber & source) :
		min(source.min), max(source.max), totalCount(source.totalCount), counts(source.counts),
		CDF(source.CDF), CDF_uptodate(source.CDF_uptodate)
	{}

	/** \brief A constructor by loading data written by WriteBinary().
	 *	\param ifile	An input stream positioned at data written by WriteBinary(). */
	LocalNumber(std::istream & ifile);

	/** \brief A member function to compare all counts.
	 *	\return True, if totalCount == sum of all counts. */
	bool CheckSum();

	/** \brief A member function to add a datum.
	 * This is not thread-safe. Please use "critical" clause before using it.
	 * \param num	A random number to add. */
	void report(int num);

	/** \brief A member function to store data in a binary file.
	 *	Data will be stored in the following order.
		min
		max
		totalCount
		counts[0]
		counts[1]
		....
		counts[..]
	 *	\param ofile	The binary output stream to write into. */
	void WriteBinary(std::ostream & ofile) const;

	/** \return The smallest element. */
	int GetMin() { return min; }
	int GetMin() const { return min; }

	/** \return The largest element. */
	int GetMax() { return max; }
	int GetMax() const { return max; }

	//====== Member functions to obtain statistics ======
	/** \return Sample size.*/
	size_t GetTotalCount() const { return totalCount; }

	/** \brief Get an estimated mean of this distribution. */
	double GetMean() const;

	/** \brief Get a standard error in the estimated mean.
	 * \param[out]	mean	An estimated mean value.
	 * \return a standard error in mean. */
	double GetSE_Mean(double & mean) const;

	/** \brief Get an estimated variance. */
	double GetVariance() const;

	/** \brief Get a standard error in the estimated variance.
	 * Bootstrap method is used.
	 * \param Bootstrap_Size	The number of bootstrap replica.
	 * \param seed	A seed for the random number generator used to draw bootstrap samples. (CLAUDE)
	 * \return Standard error of variance via the Bootstrap algorithm. */
	double GetSE_Variance(int Bootstrap_Size = 50, int seed = 0) const;

	/**	\brief Draw a random number according to the empirical distribution. (CLAUDE)
	 *	\param	rng	A random number generator.
	 *	\return A random number between this->min to this->max
	according to the empirical distribution. */
	size_t GetRandomNumber(RandomGenerator & rng) const {
		return this->max + GetRandomIndex(rng);
	}

	/** \brief Compute CDF again if CDF_uptodate = false. */
	void UpdateCDF();
	void UpdateCDF() const;

	/** \brief Generate a Bootstrap sample of this empirical distribution.
	 *	Since its return object is a smart pointer, users don't need to deallocate memory.
	 *	\param[in] seed	A seed for a random number generator.
	 *	\param[in] sample_size	The size of a Bootstrap sample. ("0" means the BootstrapSample has the same size as the original sample.)
	 *	\return	A unique pointer of a Bootstrap sample.  */
	std::unique_ptr<LocalNumber> GetBootstrapSample(int seed = 0, size_t sample_size = 0) {
		std::unique_ptr<LocalNumber> BootstrapSample = std::make_unique<LocalNumber>(this->min, this->max);
		RandomGenerator rng(seed);
		if(sample_size <= this->totalCount){
			size_t SampleSize = sample_size;
			if(SampleSize == 0)
				SampleSize = this->totalCount;

			for (size_t i = 0; i < SampleSize; i++) {
				BootstrapSample->counts[this->GetRandomIndex(rng)]++;
			}
			BootstrapSample->totalCount = SampleSize;
		}
		else
		{
			std::cerr<<"sample_size must be smaller than the totalCount\n";
		}
		
		return BootstrapSample;
	}

	/** \brief Generate a Bootstrap sample of this empirical distribution (const overload). (CLAUDE)
	 *	Since its return object is a smart pointer, users don't need to deallocate memory.
	 *	\param[in] seed	A seed for a random number generator.
	 *	\param[in] sample_size	The size of a Bootstrap sample. ("0" means the BootstrapSample has the same size as the original sample.)
	 *	\return	A unique pointer of a Bootstrap sample. */
	std::unique_ptr<LocalNumber> GetBootstrapSample(int seed = 0, size_t sample_size = 0) const {
		std::unique_ptr<LocalNumber> BootstrapSample = std::make_unique<LocalNumber>(this->min, this->max);
		RandomGenerator rng(seed);
		if(sample_size <= this->totalCount){
			size_t SampleSize = sample_size;
			if(SampleSize == 0)
				SampleSize = this->totalCount;

			for (size_t i = 0; i < SampleSize; i++) {
				BootstrapSample->counts[this->GetRandomIndex(rng)]++;
			}
			BootstrapSample->totalCount = SampleSize;
		}
		else
		{
			std::cerr<<"sample_size must be smaller than the totalCount\n";
		}

		return BootstrapSample;
	}

	/** \brief Compute an average of f = <f(i+min)> .
	 *	\param f	A Lambda function to compute an average.
	 *	\return A value of the estimator computed from this LocalNumber object. */
	inline double ComputeAverageOf(const std::function<double(double)> & f) const {
		double val = 0;
		for (size_t i = 0; i < counts.size(); i++) {
			val += f(i) * (double)counts[i] / totalCount;
		}
		return val;
	}

	/** \brief Compute an average of f = <f(i+min)> and its standard error.
	 *	\param f	A Lambda function to compute an average.
	 *	\param[out] mean	The bootstrap-averaged value of the estimator.
	 *	\param[out] SE	The standard error of the estimator from bootstrap samples.
	 *	\param Bootstrap_Size	The number of bootstrap samples used. */
	void ComputeEstimator_Bootstrap(const std::function<double(double)> & f, double & mean, double & SE, size_t Bootstrap_Size = 50) const {
		mean = 0.0; SE = 0.0;
#pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < Bootstrap_Size; i++) {
			double val = GetBootstrapSample(i)->ComputeAverageOf(f);
#pragma omp atomic
			mean += val;
#pragma omp atomic
			SE += val * val;
		}
		mean /= Bootstrap_Size;
		SE = SE / Bootstrap_Size - mean * mean;
		SE = (SE > 0) ? sqrt(SE * (double)Bootstrap_Size / (Bootstrap_Size - 1)) : 0.0;
	}

	/** \brief Give a histogram.
	 * \param result	A vector of GeometryVectors that describe the table of histogram.
						result[i].x[0] = x value
						result[i].x[1] = probability density
						result[i].x[2] = error in x
						result[i].x[3] = error in y (via multinomial distribution)
	 * \param width	Bin width in which counts are averaged.	 */
	void GetHistogram(std::vector<GeometryVector> & result, size_t width) const;

	/** \brief Give a CDF.
	 * \param result	A vector of GeometryVectors that describe the table of CDF.
		result[i].x[0] = x value
		result[i].x[1] = value of CDF */
	void GetCDF(std::vector<GeometryVector> & result) {
		result.clear();
		if (totalCount > 0) {
			UpdateCDF();
			result.reserve(this->totalCount);
			for (size_t i = 0; i < this->CDF.size(); i++)
				result.emplace_back(this->min + i, this->CDF[i]);

		}
		else {
			std::cerr << "Empirical distribution is empty.\n";
			return;
		}
	}

	/** \brief Give a CDF.
	 * \param result	A vector of GeometryVectors that describe the table of CDF.
		result[i].x[0] = x value
		result[i].x[1] = value of CDF */
	void GetCDF(std::vector<GeometryVector> & result) const {
		result.clear();
		if (totalCount > 0 && CDF_uptodate) {
			result.reserve(this->totalCount);
			for (size_t i = 0; i < this->CDF.size(); i++)
				result.emplace_back(this->min + i, this->CDF[i]);

		}
		else {
			std::cerr << "Empirical distribution is not ready.\n";
			return;
		}
	}
};

/** \brief A member function to add a datum.*/
inline void LocalNumber::report(int num) {
	//Add to existing empirical distribution.
	if (totalCount > 0) {
		// change the size of counts if necessary.
		// new datum is larger than the upper bound.
		if (num > max) {
			int num_add_end = num - max;
			for (int i = 0; i < num_add_end; i++)
				counts.push_back(0);

			max = num;
		}
		// new datum is smaller than the lower bound.
		else if (num < min) {
			int num_add_front = min - num;
			for (int i = 0; i < num_add_front; i++)
				counts.push_front(0);

			min = num;
		}
	}
	//The first datum.
	else {
		min = num;	max = num;
		counts.push_back(0);
	}

	totalCount++;
	counts[num - min]++;
	CDF_uptodate = false;
}



/** \brief Class to store information about some basis computations on LocalNumber.
 *  An object of this class can be used as an input argument for MC_nv() functions to perform prescribed Bootstrap computations. */
struct BootstrapSetup {
	size_t BootstrapSize;	//< The number of Bootstrap samples to generate
	bool PrintCDF,			//<	True = save CDF
		PrintPDF,			//< True = save PDF
		ComputeMean,		//< True = compute mean
		ComputeVariance,	//< True = compute variance 
		ComputeSkewness,	//< True = compute third moment.
		ComputeKurtosis,	//< True = compute excess kurtosis (= fourth moment - 3). 
		SaveData,			//< True = preserve data...
		ComputeL2norm,		//< True = compute L2 norm
		ComputeKL;			//< True = compute Kullback-Leibler distance. 
	std::string prefix;		//< Prefix of files in which data will be saved.
	std::string output_prefix;	//< Prefix of output (post-processed) files.
							
							/** \brief A default constructor*/
	BootstrapSetup() {
		BootstrapSize = 50;
		PrintCDF = false;	PrintPDF = false;	
		ComputeMean = true; ComputeVariance = true;
		ComputeSkewness = false;  ComputeKurtosis = false;
		prefix = "";	SaveData = false;
	}

	/** \brief A constructor by prefix of files
	 *	\param file_prefix	Prefix used for both the data and output file names. */
	BootstrapSetup(const std::string & file_prefix) : BootstrapSetup() {
		prefix = std::string(file_prefix);
		output_prefix = std::string(file_prefix);
	}

	/** \brief Generate LocalNumber object using the prefix and radius.
	 *	\param[out]	 data	If a designated file exists, a LocalNumber object corresponding is returned.
	 *						If the file doesn't exist, an empty object is returned.
	 *	\param[in]	R	Window radius. If R < 0, this argument will be ignored.*/
	void Load(LocalNumber & data, double R = -1.0) const;

	/** \brief Record data in a LocalNumber object in a binary file.
	 *	\param LN		A LocalNumber object to store.
	 *	\param R		Radius of windows that will be written in the binary file. If R < 0, this argument will be ignored.*/
	void Save(const LocalNumber & LN, double R = -1.0) const;

	/** \brief Reset output prefix for the post-processed data.
	 *  \param[in] prefix	The new prefix for post-processed output files. */
	void SetOutputPrefix(const std::string & prefix) {
		output_prefix = std::string(prefix);
	}

	/** \brief Compute requested statistics.
	 *	Parallelized with OpenMP.
	 *	\param[in] data		A LocalNumber object that stores the empirical distribution to compute statistics.
	 *	\param[out] mean	If ComputeVariance = true, mean = {mean, SE}.
	 *	\param[out] variance	If ComputeVariance = true, variance = {variance, SE}.
	 *	\param[out]	skewness	If ComputeSkewness = true, skewness = {skewness, SE}.
	 *	\param[out]	kurtosis	If ComputeKurtosis = true, kurtosis = {kurtosis, SE}.
	 *	\param[in] R		Window radius.
	 *	\param[in] SampleSize	If SampleSize = 0, the size of a Bootstrap sample is
	 							identical to the original size. */
	void Compute(const LocalNumber & data, GeometryVector & mean, GeometryVector & variance,
		GeometryVector & skewness, GeometryVector & kurtosis, double R = -1.0, size_t SampleSize = 0) const;

	/** \brief Compute distance metrices from the normal (or Gaussian) distribution.
	 *	Parallelized with OpenMP.
	 *  \param[in] data		A LocalNumber object that stores the empirical distribution to compute statistics.
	 *	\param[in] metrics	metric[i].x[0] = metric_i, metric[i].x[1] = SE_metric_i.
	 *	\param[in] R		Window radius.
	 *	\param[in] SampleSize	If SampleSize = 0, the size of a Bootstrap sample is
	 							identical to the original size. */
	void ComputeMetrics(const LocalNumber & data, std::vector<GeometryVector> & metrics,
		double R = -1.0, size_t SampleSize = 0) const;

};

/** \brief Compute Local number variance and other statistics by taking advantages of Bootstrap method.
 *	This function is parallelized in the for loops of particles.
 *	\param[in] GetConfig	A lambda function to construct Configuration objects.
 *	\param[in] NumConfigs	The number of configurations.
 *	\param[out] Result		A table of local number variance. \see MC_nv()
 *	\param[in] Rmax			Largest window radius.
 *	\param[in] dR			Bin width of window radius
 *	\param[in] WHAT2DO		A BootstrapSetup class that stores information about what to compute and aids certain bootstrap computations.
 *	\param[in] WindowCenters	A list of sampling centers in relative coordinates.
								Users can prescribe the centers of sampling windows through this argumement.  */
void MC_nv_distribution(const std::function<Configuration(size_t i)> & GetConfig, size_t NumConfigs,
	std::vector<GeometryVector> & Result, double Rmax, double dR, const BootstrapSetup & WHAT2DO, const std::vector<GeometryVector> & WindowCenters);


/** \brief Compute Local number variance and other statistics by taking advantages of Bootstrap method.
 *	This function is parallelized in the for loops of particles.
 *	This function implements a post-sampling method to reduce the finite-size effects \
 *	(currently disabled because it does not reduce the finite-size effects.)
 *	in the void MC_nv_distribution(...) function.
 *	\param[in] GetConfig	A lambda function to construct Configuration objects.
 *	\param[in] NumConfigs	The number of configurations.
 *	\param[out] Result		A table of local number variance. \see MC_nv()
 *	\param[in] Rmax			Largest window radius.
 *	\param[in] dR			Bin width of window radius
 *	\param[in] WHAT2DO		A BootstrapSetup class that stores information about what to compute and aids certain bootstrap computations.
 *	\param[in] WindowCenters	A list of sampling centers in relative coordinates.
								Users can prescribe the centers of sampling windows through this argumement.  */
void MC_nv_distribution_PostSampling(const std::function<Configuration(size_t i)> & GetConfig, size_t NumConfigs,
	std::vector<GeometryVector> & Result, double Rmax, double dR, const BootstrapSetup & WHAT2DO, const std::vector<GeometryVector> & WindowCenters);


/** \brief Compute the L2-norm of a test distrbituion from a reference distribution.
 *  The size of two input parameters must be the same.
 *  \param[in] dist_test	A test distribution function.
 *  \param[in] dist_ref		A reference distribution function.
 *  \param[in] dx			A bin size.
 *  \return L2-norm.		*/
double Compute_L2_norm(const std::vector<GeometryVector> & dist_test, const std::vector<GeometryVector> & dist_ref, double dx = 1.0);

/** \brief Compute the Kullback-Leibler distance (~relative entropy) of a test distrbituion from a reference distribution.
 *  The size of two input parameters must be the same.
 *  \param[in] dist_test	A test distribution function.
 *  \param[in] dist_ref		A reference distribution function.
 *  \param[in] dx 			A bin size.
 *  \return Kullback-Leibler distance.		*/
double Compute_KullbackLeibler_distance(const std::vector<GeometryVector> & dist_test,
										const std::vector<GeometryVector> & dist_ref, double dx = 1.0);



/** \brief	A function to give uncorrelated sampling points.
 *	\param dim			Spatial dimension of the sampling points.
 *	\param NumPoints	The number of sampling points to generate.
 *	\param seed			A seed for the random number generator.
 *	\return A vector of NumPoints uncorrelated points in relative coordinates. */
std::vector<GeometryVector> GetSamplingPoints(DimensionType dim, size_t NumPoints, int seed);



#endif