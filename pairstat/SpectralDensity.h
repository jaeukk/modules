/**
 *	Author	: Ge Zhang
 *	Email	:
 *	Date	:  */

/**
 * \file SpectralDensity.h
 * \brief Spectral density chi_V(k) of sphere packings.
 *
 * Split out of StructureFactor.{h,cpp}: the spectral density is the
 * packing-volume analog of the structure factor, built from the sphere form
 * factor m~(k;R) (Mtilde, which remains in StructureFactor.h). Spherical
 * packings only — the nonspherical/USE_SHAPES variants are not included.
 */


#ifndef SPECTRALDENSITY_INCLUDED
#define SPECTRALDENSITY_INCLUDED

#include <functional>
#include <vector>

/* header files in $(cores) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>


/**
 * \brief Spectral density chi_V(k) of a single sphere packing at one wavevector.
 * chi_V(k) = |sum_i m~(k;R_i) e^{i k.x_i}|^2 / V.
 * (Ported from the full top-level StructureFactor; used by the tiling app.)
 * \param Config  A sphere packing.
 * \param k       The wavevector at which chi_V(k) is evaluated.
 * \return        The spectral density chi_V(k).
 */
double SpectralDensity(const SpherePacking & Config, const GeometryVector & k);

/**
 * \brief Low-storage variant of SpectralDensity: recomputes the form factor per
 * particle instead of caching it in a map (for very large, highly polydisperse
 * packings where the cache would be costly).
 * \param Config  A sphere packing.
 * \param k       The wavevector at which chi_V(k) is evaluated.
 * \return        The spectral density chi_V(k).
 */
double SpectralDensity_Lstorage(const SpherePacking & Config, const GeometryVector & k);

/**
 * \brief Ensemble-averaged isotropic spectral density chi_V(|k|) of sphere packings.
 * Results[i] = ( |k|, chi_V, error in |k|, error in chi_V ), binned in |k|.
 * \param GetConfigsFunction  Returns the sphere packing of index i when called.
 * \param NumConfigs          Number of configurations to average over.
 * \param CircularKMax        Upper bound on |k| for the densely sampled region.
 * \param LinearKMax          Upper bound on |k| up to which the multiples are taken.
 * \param Results             Output table (see above).
 * \param KPrecision          Binning width in |k|.
 */
void IsotropicSpectralDensity(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision);

/**
 * \brief Low-storage variant of IsotropicSpectralDensity (uses SpectralDensity_Lstorage).
 * \param GetConfigsFunction  Returns the sphere packing of index i when called.
 * \param NumConfigs          Number of configurations to average over.
 * \param CircularKMax        Upper bound on |k| for the densely sampled region.
 * \param LinearKMax          Upper bound on |k| up to which the multiples are taken.
 * \param Results             Output table ( |k|, chi_V, error in |k|, error in chi_V ).
 * \param KPrecision          Binning width in |k|.
 */
void IsotropicSpectralDensity_Lstorage(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision);


#endif
