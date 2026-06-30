/**
 *	Author	: Ge Zhang
 *	Email	: 
 *	Date	:  */



#ifndef WeightedSTRUCTUREFACTOR_INCLUDED
#define WeightedSTRUCTUREFACTOR_INCLUDED

#include <functional>
#include <vector>

/* header fileds in $(cores)*/
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>
/**
 * \file WeightedStructureFactor.h
 * \brief Header file for computing the structure factor of point configurations with weights
 */


/**
 * \brief Compute the weighted structure factor of a single marked configuration at one wavevector.
 * Here \c Config.GetCharacteristics(i) is the (complex) weight on the i-th particle, so the returned
 * quantity is \f$ |\sum_i w_i e^{i k\cdot r_i}|^2 / V \f$ rather than the unweighted structure factor.
 * \param[in] Config	A marked point pattern; the characteristic of each particle is its weight.
 * \param[in] k		The wavevector at which the spectral function is evaluated.
 * \return The weighted spectral function at \p k.
 */
double WeightedStructureFactor(const PeriodicCellList<std::complex<double>> & Config, const GeometryVector & k);

/**
 * \brief Compute the angularly-averaged weighted structure factor over an ensemble of configurations.
 * \param[in] GetConfigsFunction	A lambda returning the i-th marked configuration.
 * \param[in] NumConfigs			The number of configurations to average over.
 * \param[in] CircularKMax		Below this |k|, every reciprocal-lattice wavevector is enumerated exactly.
 * \param[in] LinearKMax			The largest |k| at which the spectral function is evaluated.
 * \param[out] Results			Cleared, then filled with GeometryVector(|k|, S(|k|), err.|k|, err.S).
 * \param[in] KPrecision			Bin width in |k| (0.0 auto-selects the reciprocal basis spacing).
 * \param[in] SampleProbability	Fraction of wavevectors sampled above CircularKMax (1.0 = all).
 * \param[in] option				Averaging mode: 0 = over wavevector counts, 1 = over reciprocal-space
 * 								volume, 2 = individual Bragg peaks (no binning).
 * \param[in] CircularKMin		Exclude wavevectors with |k| below this value.
 */
void IsotropicStructureFactor_weighted(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision = 0.01, double SampleProbability = 1.0, size_t option = 0, double CircularKMin = 0.0);


#endif