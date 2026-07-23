#ifndef ANALYZEDATA_INCLUDED
#define ANALYZEDATA_INCLUDED

#include <vector>
#include <cmath>

const double SampleAutocorrelationLimit=0.01;//when calculating the error, we need to integrate autocorrelation, when abs(autocorrelation) is smaller than this value, cutoff the integral

/** \brief Compute the sample mean, this routine's autocorrelation-adjusted error estimate, and its truncated autocorrelation sum.
 *  \param data Input sample series.
 *  \param[out] average Sample mean.
 *  \param[out] error Square root of variance times (1+2*corrtime)/sample count.
 *  \param[out] corrtime Value initialized to 0.5 and augmented by normalized autocorrelations starting at lag one until the magnitude falls below SampleAutocorrelationLimit or becomes NaN. */
inline void AnalyzeData(const std::vector<double> & data, double & average, double & error, double & corrtime)
{
	//empty input has no defined mean/error: return zeros instead of dividing by zero
	if(data.empty())
	{
		average=0.0;
		error=0.0;
		corrtime=0.0;
		return;
	}

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

	double variance=a2-average*average;
	//constant data has zero variance: the normalized autocorrelation is undefined,
	//so report zero error and the bare tau_int prefactor without summing autocorrelations
	if(variance<=0.0)
	{
		corrtime=0.5;
		error=0.0;
		return;
	}

	//integrated autocorrelation time tau_int = 0.5 + sum_{lag>=1} normalized autocorrelation;
	//the loop starts at lag 1 because the lag-0 term is already folded into the 0.5 offset
	corrtime=0.5;
	for(size_t i=1; i<data.size(); i++)
	{
		double autocor=0.0;
		for(auto iter=data.begin()+i; iter!=data.end(); iter++)
			autocor+=(*(iter-i)-average)*(*(iter)-average);
		autocor/=(data.size()-i);
		autocor=(autocor)/variance;

		corrtime+=autocor;

		if(std::abs(autocor)<SampleAutocorrelationLimit) break;
		if(autocor!=autocor) break;//this is true if and only if autocor=NaN
	}

	error=std::sqrt(variance*(1+2*corrtime)/data.size());
}

#endif