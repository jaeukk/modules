/** Author	: Ge Zhang
 *  Email	:
 *	Date	:	*/

/** \file PairCorrelation.h
 *	\brief Header file to compute g2 function for point configurations, and S2 function for sphere packings. */

#ifndef PAIRCORRELATION_INCLUDED
#define PAIRCORRELATION_INCLUDED

#include "PeriodicCellList.h"
#include <vector>

/** \brief Class to generate a histogram from a list of values. */
class HistogramGenerator
{
public:

	/** A structure to store an individual bin for Histograms. */
	struct bin
	{
		size_t Count;	//<	The frequency in this bin.
		double Start;	//<	smallest value of the range of this bin.

		/** A constructor for given data.
		 * \param[in] (Count, Start)	Generate a bin whose smallest value in the range is "Start" and frequency is "Count".*/
		bin(size_t Count, double Start) : Count(Count), Start(Start){}

		/** A default constructor for an empty bin. */
		bin() : Count(0), Start(0){}
	};
	std::vector<bin> bins;	//<  A histrogram consists of a series of bins.

private:

	/** \brief Recursively subdivide a sorted range of values into adaptive histogram bins. (CLAUDE)
	 *  \param MinDist      Smallest value of the current range.
	 *  \param MaxDist      Largest value of the current range.
	 *  \param begin        Iterator to the first value of the current range.
	 *  \param end          Iterator past the last value of the current range.
	 *  \param MoreData     Factor by which more data are expected to come, used to estimate the relative error.
	 *  \param MinBinWidth  Minimum allowed bin width; ranges narrower than this are not subdivided.
	 *  \param TotalRange   Total span of all data, used to normalize the x-direction error. */
	void Iteration(double MinDist, double MaxDist, std::vector<double>::iterator begin, std::vector<double>::iterator end,
		double MoreData, double MinBinWidth, double TotalRange)
	{
		double count=static_cast<double>(end-begin);
		double yError = 1.0/std::sqrt(count*MoreData);//if treat these distances as a point in g_2, this is the estimated relative error in y direction
		//double Center = (MaxDist+MinDist)/2.0;
		double xError = (MaxDist-MinDist)/TotalRange/2.0;
		//if(yError>xError || count<3 || MaxDist==MinDist || MaxDist-MinDist<MinBinWidth)

		if(yError>xError || count<3 || MaxDist-MinDist<MinBinWidth)
		{
			//treat this chunk of distances as a point
			if(bins.empty())
			{
				//create a new bin, we also want to insert empty bins before the first non-empty bin
				size_t EmptyBinCount = 100;

				if(EmptyBinCount>0)
				{
					double EmptyBinWidth=0.2*TotalRange/EmptyBinCount;
					for(size_t i=0; i<EmptyBinCount; i++)
					{
						bin temp;
						temp.Count=0;
						temp.Start=MinDist-0.2*TotalRange+i*EmptyBinWidth;
						bins.push_back(temp);
					}
				}
				else if(MinDist>0)
				{
					//always have a bin at 0
					bin temp;
					temp.Count=0;
					temp.Start=0.0;
					bins.push_back(temp);
				}
				bin temp;
				temp.Count=end-begin;
				temp.Start=MinDist;
				bins.push_back(temp);
			}
			else if(MinDist==bins.back().Start /*|| count<5*/)
			{
				//merge this into previous bin because either:
				//the bin is identical to the previous bin, or
				//the count is too small, the bin cannot be precise.
				bins.back().Count+=end-begin;
			}
			else
			{
				//debug temp
				//if (MaxDist - MinDist>10000) std::cout << MaxDist << " \t" << MinDist << " \t" << count << '\n';
				//create a new bin
				bin temp;
				temp.Count=end-begin;
				temp.Start=MinDist;
				bins.push_back(temp);
			}

		}
		else
		{
			//divide this chunk
			std::vector<double>::iterator mid=begin+(end-begin)/2;
			std::nth_element(begin, mid, end);
			Iteration(MinDist, *mid, begin, mid, MoreData, MinBinWidth, TotalRange);
			Iteration(*mid, MaxDist, mid, end, MoreData, MinBinWidth, TotalRange);
		}
	}
public:
	/** \brief Count all data, generate bins using adaptive method, set the Count of each bin to their occurence in data.
	 *  The vector data is changed upon return.
	 *  \param data                 The values to histogram; reordered in place upon return.
	 *  \param MoreData             MoreData=x indicates that x times more data are expected to come.
	 *  \param ResolutionPreference >1 to get better x resolution, <1 to get better y resolution.
	 *  \param MinBinWidth          Minimum allowed bin width. */
	void GenerateBins(std::vector<double> & data, double MoreData, double ResolutionPreference, double MinBinWidth=0)
	{
		if(data.empty())
		{
			std::cerr<<"Error in HistogramGenerator::GenerateBins : no data!\n";
			return;
		}
		//std::nth_element(data.begin(), data.begin(), data.end());
		//std::nth_element(data.begin()+1, data.end()-1, data.end());

		//debug temp
		//for(auto iter=data.begin(); iter!=data.end(); iter++)
		//	std::cout<<*iter<<" \t";

		double Max=*std::max_element(data.begin(), data.end());
		double Min=*std::min_element(data.begin(), data.end());
		Iteration(Min, Max, data.begin(), data.end(), MoreData*ResolutionPreference, MinBinWidth, Max-Min); 

		//The bins generated by the previous iteration process is often too small. Merge them here.
		std::vector<bin> newBins;
		for(auto iter=this->bins.begin();;)
		{
			if(iter==bins.end())
				break;
			else if(iter==bins.end()-1)
			{
				newBins.push_back(*iter);
				break;
			}
			else
			{
				bin temp;
				temp.Count=iter->Count+(iter+1)->Count;
				temp.Start=iter->Start;
				newBins.push_back(temp);
				iter+=2;
			}
		}
		double dist = 0.2*(Max - Min) / 50.0;
		for(int i=0; i<50; i++)
		{
			bin temp;
			temp.Start=Max+i*dist;
			temp.Count=0;
			newBins.push_back(temp);
		}
		std::swap(newBins, this->bins);
	}

	/** \brief Report a new data, increase the corresponding bin count by 1.
	 *  Call it after calling GenerateBins. This member function is thread safe.
	 *  \param Data  The value whose bin count is incremented. */
	void Report(const double & Data)
	{
		//if user call Report() before calling GenerateBins(), assertion will fail.
		assert(this->bins.empty()==false);

		if(Data<bins[0].Start)
		{
			std::cerr<<"Warning in HistogramGenerator::Report : Data is out of range!\n";
			return;
		}
		bin temp;
		temp.Start=Data;
		std::vector<bin>::iterator loc=std::upper_bound(bins.begin(), bins.end(), temp, [](const bin & left, const bin & right)->bool{return left.Start<right.Start;});
		assert(loc!=bins.begin());
		loc--;// decrease it by one because std::upper_bound returns the first iterator which is GREATER. However, our bin::Start is the start location of the bin

#pragma omp atomic
		loc->Count++;
	}
};


//deprecated. Use IsotropicTwoPairCorrelation with large SampleDistanceSize instead
//void IsotropicTwoPairCorrelation_FullAdaptive(std::vector<const Configuration *> Configs, double MaxDistance, std::vector<GeometryVector> & Result, double ResolutionPreference=1.0);

//clear Result, then fill it with the pair correlation data calculated from Config
//Result is filled with 4-dimensional GeometryVectors, the elements are:
//(r, g_2(r), \delta r, \delta g_2(r) )
//gets SampleDistanceSize number of sample pair distances, use it to generate bins, then count all pair distances
//ResolutionPreference : >1 to get better r resolution, <1 to get better g_2 resolution
/** \brief Compute an isotropic g2 function from a number of point configurations.
 * \param[in] (GetConfigsFunction, NumConfigs)	A Lambda function to generate "NumConfigs" point configurations.
 * \param[in] MaxDistance	Up to this distance, the g2 function is computed.
 * \param[out] Result		A resulting g2 function.
							First clear this vector and store data as follows:
							Result[i] = (radius, g2(radius), error in radius, error in g2).
 * \param[in] SampleDistanceSize	The number of pairs to determine a histogram.
									The larger this argument, the finer bins the histogram has. This sampling is performed in a single thread.
 * \param[in] ResolutionPreference	>1 to get better r resolution, <1 to get better g_2 resolution. */
void IsotropicTwoPairCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs,
	double MaxDistance, std::vector<GeometryVector> & Result, size_t SampleDistanceSize=5000000, double ResolutionPreference=1.0);

/** \brief Compute an isotropic g2 function from a number of point configurations.
 * \param[in] (GetConfigsFunction, NumConfigs)	A Lambda function to generate "NumConfigs" point configurations.
 * \param[in] MaxDistance	Up to this distance, the g2 function is computed.
 * \param[out] Result		A resulting g2 function.
							First clear this vector and store data as follows.
							Result[i] = (radius, g2(radius), error in radius, error in g2).
 * \param[in] HGen			A customized histogram to compute the g2 function. */
void IsotropicTwoPairCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, 
	double MaxDistance, std::vector<GeometryVector> & Result, HistogramGenerator & HGen);


//clear Result, then fill it with the pair correlation data calculated from Config
//Result is filled with 3-dimensional GeometryVectors, the elements are:
//(r, theta, g_2(r) )
//r bins are adaptively generated
//gets SampleDistanceSize number of sample pair distances, use it to generate bins, then count all pair distances
//ResolutionPreference : >1 to get better r resolution, <1 to get better g_2 resolution
//theta bins are fixed. there are NumDirectionalBins of them.
//if pError is not nullptr, it is filled with 3-dimensional GeometryVectors, the elements are:
//(delta r, delta theta = pi/(2*NumDirectionalBins), delta g_2(r) )
/** \brief Compute the directional (2D anisotropic) pair correlation function g2(r, theta) from a number of point configurations. (CLAUDE)
 *  \param GetConfigsFunction  A lambda function to generate "NumConfigs" point configurations.
 *  \param NumConfigs          The number of point configurations.
 *  \param MaxDistance         Up to this distance, the g2 function is computed.
 *  \param Result              Cleared, then filled with (r, theta, g2(r,theta)) per directional bin.
 *  \param NumDirectionalBins  Number of fixed angular (theta) bins.
 *  \param pError              If non-null, filled with (delta r, delta theta = pi/(2*NumDirectionalBins), delta g2(r)).
 *  \param SampleDistanceSize  Number of sample pair distances used to generate the adaptive r bins.
 *  \param ResolutionPreference  >1 to get better r resolution, <1 to get better g2 resolution. */
void TwoPairCorrelation_2DAnisotropic(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs,
	double MaxDistance, std::vector<std::vector<GeometryVector> > & Result, size_t NumDirectionalBins,
	std::vector<std::vector<GeometryVector> > * pError=nullptr, size_t SampleDistanceSize=2000000, double ResolutionPreference=1.0);

//TODO: interface this function to Computation class.
//Probability density function of the nearest neighbor distance
//clear Result, then fill it with the PDF calculated from Config
//Result is filled with 4-dimensional GeometryVectors, the elements are:
//(r, p(r), \delta r, \delta g_2(r) )
//gets SampleDistanceSize number of sample pair distances, use it to generate bins, then count all pair distances
//ResolutionPreference : >1 to get better r resolution, <1 to get better g_2 resolution
/** \brief Compute the probability density function H_p(r) of the nearest-neighbor distance from a number of point configurations. (CLAUDE)
 *  \param GetConfigsFunction  A lambda function to generate "NumConfigs" point configurations.
 *  \param NumConfigs          The number of point configurations.
 *  \param Result              Cleared, then filled with (r, p(r), delta r, delta p(r)).
 *  \param SampleDistanceSize  Number of sample distances used to generate the adaptive bins.
 *  \param ResolutionPreference  >1 to get better r resolution, <1 to get better p(r) resolution. */
void NearestNeighborDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleDistanceSize=500000, double ResolutionPreference=1.0);

//calculate Hv
//Each configuration will be sampled OverSampling*NumParticle() times
/** \brief Calculate the void nearest-neighbor distribution function H_v(r); each configuration is sampled OverSampling*NumParticle() times. (CLAUDE)
 *  \param GetConfigsFunction  A lambda function to generate "NumConfigs" point configurations.
 *  \param NumConfigs          The number of point configurations.
 *  \param Result              Cleared, then filled with (r, p(r), delta r, delta p(r)).
 *  \param SampleDistanceSize  Number of sample distances used to generate the adaptive bins.
 *  \param ResolutionPreference  >1 to get better r resolution, <1 to get better p(r) resolution.
 *  \param OverSampling        Number of random sampling points per particle. */
void HvDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleDistanceSize=500000, double ResolutionPreference=1.0, double OverSampling=1.0);


//return the minimum pair distance in the configuration
//if pi and pj are not nullptr, *pi and *pj will be set to the particles that are closest to each other
/** \brief Return the minimum pair distance in the configuration. If pi and pj are not nullptr, *pi and *pj will be set to the particles that are closest to each other.
 *  \param a   The point configuration.
 *  \param pi  If non-null, set to the index of one particle of the closest pair.
 *  \param pj  If non-null, set to the index of the other particle of the closest pair.
 *  \return    The minimum pair distance. */
double MinDistance(const Configuration & a, size_t * pi=nullptr, size_t * pj=nullptr);


/** \brief Return the mean nearest-neighbor distance over all particles in the configuration. (CLAUDE)
 *  \param a  The point configuration.
 *  \return   The mean nearest-neighbor distance. */
double MeanNearestNeighborDistance(const Configuration & a);


//convert histogram to probability density function or cumulative distribution function
//return total count the HistogramGenerator has
/** \brief Convert histogram to probability density function.
 *  \param HGen    The histogram to convert.
 *  \param result  Cleared, then filled with the probability density function (r, p(r), delta r, delta p(r)).
 *  \return        Total count the HistogramGenerator has. */
size_t HistogramToPDF(const HistogramGenerator & HGen, std::vector<GeometryVector> & result);
/** \brief Convert histogram to cumulative distribution function.
 *  \param HGen    The histogram to convert.
 *  \param result  Cleared, then filled with the cumulative distribution function (r, P(r), delta r, delta P(r)).
 *  \return        Total count the HistogramGenerator has. */
size_t HistogramToCDF(const HistogramGenerator & HGen, std::vector<GeometryVector> & result);



#endif