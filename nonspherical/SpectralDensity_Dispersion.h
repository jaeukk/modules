/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com	
*	Date	:	Jan. 2021 */

#ifndef SPECTRALDENSITY_DISPERSION_INCLUDED
#define SPECTRALDENSITY_DISPERSION_INCLUDED

#include <functional>
#include <vector>

#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include "Dispersions.h"
#include "ParticleShapes.h"

/** \file SpectralDensity_Dispersion.h
 * \brief Header file for computing the spectral density of a Dispersion class.
 */



/** A function to compute spectral density of a dispersion at a wavevector k.
* It assumes that particles are not overlapping.
* It follows the convention (Fourier transform  = \int d{\vect{r}} e^{i k}    )
* @param config	A Dispersion object that describes a packing of (presumably) nonoverlapping nonspherical particles.
* @param k	A wavevector at which we compute the spectral density.
* @return A value of spectral density at k. */
double SpectralDensity_Dispersion(const Dispersion & config, const GeometryVector &k);


/** A function to compute an orientationally averaged spectral density.
 * We assume that all dispersions are nonoverlapping. @see SpectralDensity_nonspherical().
 * It computes values of the spectral density at all wavevectors searched by two parameters (CircularKmax, LinearKmax), and
	then computes the averages of spectral density at wavevectors whose magnitude is between [n*KPrecision, (n+1)*KPrecision].
 * @param[in] GetConfiguration	A Lambda function to generate Dispersion objects that we have interest in.
 * @param[in] NumConfigs	The number of configurations used in computation.
 * @param[in] (CircularKmax, LinearKmax)	Searching radii for wavevectors.
				Search all wavevectors (except for their parity pairs) up to k=CircularKmax,
				and search all wavevectors that are scalar multiplications of the aforementioned wavevectors up to k=LinearKmax.
 * @param[out] result		A table of spectral density.
				result[i].x[0] = wavenumber.
				result[i].x[1] = spectral density.
				result[i].x[2] = standard error of k
				result[i].x[3] = standard error of spectral density.
 * @param[in] KPrecision	Bin witdth of wavenumbers.*/
void IsotropicSpectralDensity_Dispersion(std::function<const Dispersion(size_t)> GetConfiguration, size_t NumConfigs,
	double CircularKmax, double LinearKMax, std::vector<GeometryVector> & result, double KPrecision = 0.01);


/** A function to compute spectral density to draw a contour plot.
 * @param[in] GetConfiguration	A Lambda function to generate Dispersion objects that we have interest in.
 * @param[in] NumConfigs	The number of configurations used in computation.
 * @param[in] limit	The number of wavevectors to consider. Wavevectors are k = i*k_1 + j*k_2 + ..., where i, j, ... are integers between -limit and limit.
 * @param[out] (ks, chik)	ks[i] is the ith wavevector, and chik[i] is the spectral density at this wavevector.
 * @param[in] dk	The number of pixels over which the spectral density is averaged.*/
void DirectionalSpectralDensity_Dispersion(std::function<const Dispersion(size_t)> GetConfigsFunction, size_t NumConfigs,
	long limit, std::vector<GeometryVector> & ks, std::vector<GeometryVector> & chik, double dk = 1);



#endif