/** Author	: Ge Zhang
 *  Email	:
 *	Date	:	*/

/**
 * \file SpectralDensity.cpp
 * \brief Implementation of the spectral density chi_V(k) of sphere packings
 * (split out of StructureFactor.cpp).
 */

#include "SpectralDensity.h"
#include "StructureFactor.h"	// Mtilde (sphere form factor)

#include <complex>
#include <cmath>
#include <map>
#include <etc.h>		// Verbosity, progress_display, MaxDimension


/**
 * \struct SkBin
 * \brief Per-bin accumulator for isotropic averaging in |k|.
 * (Mirrors the file-local SkBin in StructureFactor.cpp; kept local here so the
 * two translation units stay independent.)
 */
struct SkBin
{
	double Sum1,	///< count
		SumK1,		///<	sum of the magnitudes of wavevectors (|k|)
		SumK2,		///< sum of squared magnitudes of wavevectors (|k|^2)
		SumS,		///< sum of spectral functions
		SumS2;		///< sum of squared spectral functions
	/** \brief Construct an empty bin with all accumulators set to zero. */
	SkBin() : Sum1(0), SumK1(0), SumK2(0), SumS(0), SumS2(0)
	{}
};


// ---- Spectral density chi_V(k) of sphere packings -------------------------
//      Ported verbatim from the full top-level StructureFactor (spherical only;
//      no nonspherical/USE_SHAPES baggage). Needed by the tiling app.
double SpectralDensity(const SpherePacking & Config, const GeometryVector & k)
{
	std::complex<double> rho;
	std::map<double, double> mtildes;
	size_t NumParticle = Config.NumParticle();
	double kl = std::sqrt(k.Modulus2());
	DimensionType d = Config.GetDimension();
	for (size_t i = 0; i<NumParticle; i++)
	{
		std::complex<double> a = std::exp(std::complex<double>(0, k.Dot(Config.GetCartesianCoordinates(i))));
		double mt;
		double R = Config.GetCharacteristics(i);
		auto iter = mtildes.find(R);
		if (iter != mtildes.end())
			mt = iter->second;
		else
		{
			mt = Mtilde(kl, R, d);
			mtildes.insert(std::make_pair(R, mt));
		}
		rho += a*mt;
	}
	return (rho.real()*rho.real() + rho.imag()*rho.imag()) / Config.PeriodicVolume();
}

double SpectralDensity_Lstorage(const SpherePacking & Config, const GeometryVector & k)
{
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

void IsotropicSpectralDensity(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision)
{
	Results.clear();
	if (!(CircularKMax>0))
		return;
	if (!(LinearKMax>0))
		return;

	if (KPrecision == 0.0)
	{
		std::cerr << "Warning in IsotropicSpectralDensity : This version does not support KPrecision==0.0. Auto choosing this quantity!\n";
		SpherePacking c = GetConfigsFunction(0);
		KPrecision = std::sqrt(c.GetReciprocalBasisVector(0).Modulus2());
	}
	size_t NumBin = std::floor(LinearKMax / KPrecision) + 1;
	std::vector<SkBin> vSkBin(NumBin, SkBin());

	GeometryVector prevBasis[::MaxDimension];
	std::vector<GeometryVector> ks;
	if (Verbosity>1)
		std::cout << "Computing chi(k)";

	{
		SpherePacking CurrentConfig = GetConfigsFunction(0);
		ks = GetKs(CurrentConfig, CircularKMax, LinearKMax);
		if (Verbosity > 0)
			std::cout << "the number of considered k-points = " << ks.size() << std::endl;
	}

	progress_display pd(NumConfigs);
	for (size_t j = 0; j<NumConfigs; j++)
	{
		SpherePacking CurrentConfig = GetConfigsFunction(j);
		if (CurrentConfig.GetDimension() == 0)
			continue;
		if (j != 0)
		{
			bool SameBasis = true;
			for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
				if (!(prevBasis[i] == CurrentConfig.GetBasisVector(i)))
					SameBasis = false;

			if (SameBasis == false)
				ks=GetKs(CurrentConfig, CircularKMax, LinearKMax);
		}
		signed long end = ks.size();
#pragma omp parallel for schedule(guided)
		for (signed long i = 0; i<end; i++)
		{
			double s = SpectralDensity(CurrentConfig, ks[i]);
			double k2 = ks[i].Modulus2();
			size_t Bin = std::floor(std::sqrt(k2) / KPrecision);
#pragma omp atomic
			vSkBin[Bin].Sum1 += 1.0;
#pragma omp atomic
			vSkBin[Bin].SumK1 += sqrt(k2);
#pragma omp atomic
			vSkBin[Bin].SumK2 += k2;
#pragma omp atomic
			vSkBin[Bin].SumS += s;
#pragma omp atomic
			vSkBin[Bin].SumS2 += s*s;
		}
		for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
			prevBasis[i] = CurrentConfig.GetBasisVector(i);
		pd++;
	}
	for (auto iter = vSkBin.begin(); iter != vSkBin.end(); iter++)
	{
		if (iter->Sum1 != 0.0)
		{
			GeometryVector temp(4);
			temp.x[0] = iter->SumK1 / iter->Sum1;
			temp.x[1] = iter->SumS / iter->Sum1;
			double var = (iter->SumK2 / (iter->Sum1) - temp.x[0] * temp.x[0]) / (iter->Sum1);
			temp.x[2] = (var > 0) ? std::sqrt(var) : 0;
			temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1));
			Results.push_back(temp);
		}
	}

	if (Verbosity>2)
		std::cout << "done!\n";
}

void IsotropicSpectralDensity_Lstorage(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision)
{
	Results.clear();
	if (!(CircularKMax>0))
		return;
	if (!(LinearKMax>0))
		return;

	if (KPrecision == 0.0)
	{
		std::cerr << "Warning in IsotropicSpectralDensity : This version does not support KPrecision==0.0. Auto choosing this quantity!\n";
		SpherePacking c = GetConfigsFunction(0);
		KPrecision = std::sqrt(c.GetReciprocalBasisVector(0).Modulus2());
	}
	size_t NumBin = std::floor(LinearKMax / KPrecision) + 1;
	std::vector<SkBin> vSkBin(NumBin, SkBin());

	GeometryVector prevBasis[::MaxDimension];
	std::vector<GeometryVector> ks;
	if (Verbosity>1)
		std::cout << "Computing chi(k)";

	{
		SpherePacking CurrentConfig = GetConfigsFunction(0);
		ks = GetKs(CurrentConfig, CircularKMax, LinearKMax);
		if (Verbosity > 0)
			std::cout << "the number of considered k-points = " << ks.size() << std::endl;
	}

	progress_display pd(NumConfigs);
	for (size_t j = 0; j<NumConfigs; j++)
	{
		SpherePacking CurrentConfig = GetConfigsFunction(j);
		if (CurrentConfig.GetDimension() == 0)
			continue;
		if (j != 0)
		{
			bool SameBasis = true;
			for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
				if (!(prevBasis[i] == CurrentConfig.GetBasisVector(i)))
					SameBasis = false;

			if (SameBasis == false)
				ks = GetKs(CurrentConfig, CircularKMax, LinearKMax);
		}

		signed long end = ks.size();
#pragma omp parallel for schedule(guided)
		for (signed long i = 0; i<end; i++)
		{
			double s = SpectralDensity_Lstorage(CurrentConfig, ks[i]);
			double k2 = ks[i].Modulus2();
			size_t Bin = std::floor(std::sqrt(k2) / KPrecision);
#pragma omp atomic
			vSkBin[Bin].Sum1 += 1.0;
#pragma omp atomic
			vSkBin[Bin].SumK1 += sqrt(k2);
#pragma omp atomic
			vSkBin[Bin].SumK2 += k2;
#pragma omp atomic
			vSkBin[Bin].SumS += s;
#pragma omp atomic
			vSkBin[Bin].SumS2 += s * s;
		}
		for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
			prevBasis[i] = CurrentConfig.GetBasisVector(i);
		pd++;
	}
	for (auto iter = vSkBin.begin(); iter != vSkBin.end(); iter++)
	{
		if (iter->Sum1 != 0.0)
		{
			GeometryVector temp(4);
			temp.x[0] = iter->SumK1 / iter->Sum1;
			temp.x[1] = iter->SumS / iter->Sum1;
			double var = (iter->SumK2 / (iter->Sum1) - temp.x[0] * temp.x[0]) / (iter->Sum1);
			temp.x[2] = (var > 0) ? std::sqrt(var) : 0;
			temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1));
			Results.push_back(temp);
		}
	}

	if (Verbosity>2)
		std::cout << "done!\n";
}
