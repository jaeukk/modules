/**
 * \file SpectralDensity.h
 * \author Ge Zhang. [Modified by Jaeuk Kim (phy000.kim@gmail.com)]
 * \brief  Header file for computing the spectral density of configurations of sphere packing.
 * \version 0.1
 * \date 2024-08-27
 *
 *
 */

#ifndef __SPECTRALDENSITY_H__
#define __SPECTRALDENSITY_H__

#include <functional>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
/* ----- core source files */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>

/** \brief Fourier transform of the indicator function of a hypersphere of radius R
 *  \param k Wavenumber
 *  \param R Sphere radius
 *  \param d Space dimensions
 *  \return real-valued Fourier transform */
double Mtilde(double k, double R, DimensionType d);


/**
 * \brief Compute a spectral density of a sphere packing at a wavevector. Since this code stores a list of values of Mtilde, it is efficient for a sphere packing whose particle radii are similar or narrowly distributed.
 * \param Config A sphere packing
 * \param k Wavevector
 * \return Spectral density	*/

double SpectralDensity(const SpherePacking & Config, const GeometryVector & k);


/** \brief Compute spectral densities of sphere packings for all wavevectors based on CircularKMax and LinearKMax.
 *  \param GetConfigsFunction A function to load (or generate) sphere packings
 *  \param NumConfigs Number of configurations
 *  \param CircularKMax All reciprocal lattice vectors are considered if their magnitudes are smaller than this value.
 *  \param LinearKMax All multiples of the wavevectors whose magnitudes are smaller than CircularKMax are considered up to this value.
 *  \param Results Non-negative real-valued spectral densities.
 *  \param KPrecision Bin size of wavenumbers. */
void IsotropicSpectralDensity(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision);


/** \brief A low-memory version of SpectralDensity(...). It saves memory space significantly for extremely large sphere packings with a high polydispersity in size.
 *  \param Config A sphere packing
 *  \param k Wavevector
 *  \return Spectral density */
inline double SpectralDensity_Lstorage(const SpherePacking & Config, const GeometryVector & k) {
	std::complex<double> rho;
	size_t NumParticle = Config.NumParticle();
	double kl = std::sqrt(k.Modulus2());
	DimensionType d = Config.GetDimension();
	double mt;
	for (size_t i = 0; i<NumParticle; i++)
	{
		std::complex<double> a = std::exp(std::complex<double>(0, k.Dot(Config.GetCartesianCoordinates(i))));
		double R = Config.GetCharacteristics(i);
		mt = Mtilde(kl, R, d);
		rho += a * mt;
	}
	return (rho.real()*rho.real() + rho.imag()*rho.imag()) / Config.PeriodicVolume();
}



/** \brief Low-memory version of IsotropicSpectralDensity(...). It computes spectral densities of sphere packings for all wavevectors based on CircularKMax and LinearKMax.
 *  \param GetConfigsFunction A function to load (or generate) sphere packings
 *  \param NumConfigs Number of configurations
 *  \param CircularKMax All reciprocal lattice vectors are considered if their magnitudes are smaller than this value.
 *  \param LinearKMax All multiples of the wavevectors whose magnitudes are smaller than CircularKMax are considered up to this value.
 *  \param Results Non-negative real-valued spectral densities.
 *  \param KPrecision Bin size of wavenumbers. */
void IsotropicSpectralDensity_Lstorage(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision);


#endif // __SPECTRALDENSITY_H__