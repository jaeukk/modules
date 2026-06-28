#ifndef ANALYZEDATA_INCLUDED
#define ANALYZEDATA_INCLUDED

#include <vector>

const double SampleAutocorrelationLimit=0.01;//when calculating the error, we need to integrate autocorrelation, when abs(autocorrelation) is smaller than this value, cutoff the integral

/** \brief Compute the mean, statistical error, and integrated autocorrelation time of a data series. (CLAUDE)
 *  \param data      Input sample series (e.g. a time series of a measured observable).
 *  \param average   Output: sample mean of the data.
 *  \param error     Output: statistical error of the mean, corrected for autocorrelation.
 *  \param corrtime  Output: integrated autocorrelation time of the series. */
inline void AnalyzeData(const std::vector<double> & data, double & average, double & error, double & corrtime)
{
	//calculate <data> and <data^2>
	double a2=0;
	average=0;
	for(auto iter=data.begin(); iter!=data.end(); iter++)
	{
		average+=*iter;
		a2+=(*iter) * (*iter);
	}
	average/=data.size();
	a2/=data.size();

	//calculate autocorrelation function
	corrtime=0.5;
	for(size_t i=0; i<data.size(); i++)
	{
		double autocor=0.0;
		for(auto iter=data.begin()+i; iter!=data.end(); iter++)
			autocor+=(*(iter-i)-average)*(*(iter)-average);
		autocor/=(data.size()-i);
		autocor=(autocor)/(a2-average*average);

		corrtime+=autocor;

		if(std::abs(autocor)<SampleAutocorrelationLimit) break;
		if(autocor!=autocor) break;//this is true if and only if autocor=NaN
	}

	error=std::sqrt((a2-average*average)*(1+2*corrtime)/data.size());
}

#endif