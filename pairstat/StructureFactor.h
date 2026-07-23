/**
 *	Author	: Ge Zhang
 *	Editor	: Jaeuk Kim (phy000.kim@gmail.com)
 *	Email	: 
 *	Date	:  */



#ifndef STRUCTUREFACTOR_INCLUDED
#define STRUCTUREFACTOR_INCLUDED

#include <functional>
#include <vector>

/* header fileds in $(cores)*/
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>
/**
 * \file StructureFactor.h
 * \brief Header file for computing the structure factor of point configurations,
 * spectral density of sphere packings, and spectral density for packings of nonspherical particles (optional).
 */





/**
 * \brief Fourier transform of a d-dimensional sphere indicator of radius R
 * (the "form factor" m~(k;R)) in dimension d = 1..4.
 * (Ported from the full top-level StructureFactor; GetKs already lives in core's
 * PeriodicCellList, so only Mtilde is added here. Used by FFH / CCOptimization.)
 */
double Mtilde(double k, double R, DimensionType d);

// The spectral density chi_V(k) of sphere packings — SpectralDensity[_Lstorage]
// and IsotropicSpectralDensity[_Lstorage] — was split into SpectralDensity.h.
// (It is built from Mtilde, which stays here as the shared sphere form factor.)

/**
 * \brief Calculate the structure factor at a certain k.
 * \param Config  A point configuration.
 * \param k       The wavevector at which S(k) is evaluated.
 * \return        The structure factor S(k).
 */
double StructureFactor(const Configuration & Config, const GeometryVector & k);

/**
 * \brief Calculate the structure factor with weights; here config.getCharacteristics(i) = the weight on the i-th particle.
 * \param Config  A point configuration whose per-particle characteristics hold the complex weights.
 * \param k       The wavevector at which S(k) is evaluated.
 * \return        The weighted structure factor S(k).
 */
double WeightedStructureFactor(const PeriodicCellList<std::complex<double>> & Config, const GeometryVector & k);

/**
 * \brief Calculate structure factor of all ks in the range CircularKMax (i.e. for all ks that abs(k)<CircularKMax) and all of their multiples in the range LinearKMax.
 * Results has the form ( abs(k), S(k), KPrecision, \delta S(k) ). If SampleProbability is <1, then for any k point
 * satisfying the above condition there is SampleProbability probability that it will be used to calculate S(k)
 * and (1-SampleProbability) probability that it will not be used.
 * \param GetConfigsFunction  Returns the configuration of index i when called.
 * \param NumConfigs          Number of configurations to average over.
 * \param CircularKMax        Upper bound on |k| for the densely sampled (circular) region.
 * \param LinearKMax          Upper bound on |k| up to which the multiples of the circular ks are taken.
 * \param Results             Output table; each entry is ( |k|, S(k), KPrecision, error in S(k) ).
 * \param KPrecision          Binning width in |k|.
 * \param SampleProbability   Probability that each eligible k point is included in the average.
 * \param option              Averaging mode: 0 average over reciprocal lattice points, 1 average over reciprocal-space volume, 2 structure factor at Bragg peaks.
 * \param CircularKMin        Lower bound on |k| below which wavevectors are excluded. 0 by default.
 */
void IsotropicStructureFactor(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision=0.01, double SampleProbability=1.0, size_t option = 0, double CircularKMin = 0.0);

// DirectionalStructureFactor(...) was declared here but never defined (its
// implementation was lost when the top-level StructureFactor was split); the
// declaration is removed to avoid link failures for callers.  The 2D directional
// S(k) is available via visualization/Plots' PlotDirectionalStructureFactor2D.

#endif