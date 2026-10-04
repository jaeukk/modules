/** Author	: Jaeuk Kim
 *  Email	: phy000.kim@gmail.com
 *	Date	: 2024. 10.	*/

/**
 * \file PairCorrelation.h
 * \brief Header file to compute g2 function for point configurations, and S2 function for sphere packings.
 */

#ifndef PAIRCORRELATION_INCLUDED
#define PAIRCORRELATION_INCLUDED

#include <PeriodicCellList.h>
#include <vector>
#include <complex>

//void Weighted IsotropicTwoPairCorrelation_FullAdaptive(std::vector<const Configuration *> Configs, double MaxDistance, std::vector<GeometryVector> & Result, double ResolutionPreference=1.0);

/**
 * \brief Compute the weighted (two-sided) isotropic pair correlation function over an ensemble.
 * For each pair the contribution is the product of conjugated weights
 * \f$ \overline{w_i}\, w_j \f$, so the result reduces to the ordinary \f$ g_2(r) \f$ when all weights are 1.
 * \param[in] GetConfigsFunction	A lambda returning the i-th marked configuration.
 * \param[in] NumConfigs	The number of configurations to average over.
 * \param[in] MaxDistance	Up to this distance, the g2 function is computed.
 * \param[out] Result		Cleared, then filled with 4-dimensional GeometryVectors:
 * Result[i] = (radius, Re g2(radius), Im g2(radius), error in g2).
 * \param[in] dR			Bin width in r (controls the histogram resolution).
 */
void WeightedIsotropicTwoPairCorrelation(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs,
	double MaxDistance, std::vector<GeometryVector> & Result, double dR = 0.1);


/**
 * \brief Compute the weighted single-sided isotropic pair correlation function over an ensemble.
 * Unlike WeightedIsotropicTwoPairCorrelation, each neighbor contributes only the weight of the
 * source particle (not the conjugated product), giving the one-sided weighted neighbor density.
 * \param[in] GetConfigsFunction	A lambda returning the i-th marked configuration.
 * \param[in] NumConfigs	The number of configurations to average over.
 * \param[in] MaxDistance	Up to this distance, the function is computed.
 * \param[out] Result		Cleared, then filled with 4-dimensional GeometryVectors:
 * Result[i] = (radius, Re value(radius), Im value(radius), error).
 * \param[in] dR			Bin width in r.
 */
void WeightedIsotropicSingleSidedTwoPairCorrelation(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs,
	double MaxDistance, std::vector<GeometryVector> & Result, double dR = 0.1);



#endif