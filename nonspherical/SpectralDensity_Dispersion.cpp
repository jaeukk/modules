/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com	
*	Date	:	Jan. 2021 */

#include "SpectralDensity_Dispersion.h"
#include <etc.h>
#include <ParameterInput.h>
#include <complex>
#include <cmath>
#include <algorithm>

/** \file SpectralDensity_Dispersion.cpp
 *  \brief Function implementations to compute spectral density of Dispersion objects. */


/** \struct A class to store data in each bin. */
struct SkBin
{
	double Sum1,	//< count
		SumK1,		//<	sum of the magnitudes of wavevectors (|k|)
		SumK2,		//< sum of squared magnitudes of wavevectors (|k|^2)
		SumS,		//< sum of spectral functions
		SumS2;		//< sum of squared spectral functions
	SkBin() : Sum1(0), SumK1(0), SumK2(0), SumS(0), SumS2(0)
	{}
};


/** Compute a spectral density for a given configuration. */
double SpectralDensity_Dispersion(const Dispersion & config, const GeometryVector &k) {
	//Check Dimensions
	if (config.GetDimension() == k.Dimension) {
		std::complex <double> tildeI(0.0, 0.0);
		/*Note that the phase factor is exp(-i k.x)
		 for the consistency with the Fourier transform employed
		 in the form factor. */
		for (size_t i = 0; i < config.NumParticle(); i++) {
			tildeI += config.GetFormFactor(i, k) * std::exp(std::complex<double>(0.0, -k.Dot(config.GetCartesianCoordinates(i))));
		}

		// |tildeI|^2 / V_F
		return std::norm(tildeI) / config.PeriodicVolume();	
	}
	else {
		std::cerr << "SpectralDensity_nonspherical() have incorrect arguements." << std::endl;
	}
	return 0.0;
}

void IsotropicSpectralDensity_Dispersion(std::function<const Dispersion(size_t)> GetConfiguration, size_t NumConfigs,
	double CircularKmax, double LinearKMax, std::vector<GeometryVector> & Result, double KPrecision){

	double dK = KPrecision, SampleProbability = 1.0;
	if (KPrecision == 0) {
					std::cerr << "Warning in IsotropicStructureFactor : This version does not support KPrecision==0.0. Auto choosing this quantity!\n";
		dK = std::sqrt(GetConfiguration(0).GetReciprocalBasisVector(0).Modulus2());
	}

	size_t NumBin = std::floor(LinearKMax / dK) + 1;
	std::vector<SkBin> vSkBin(NumBin, SkBin());		//	data in each bin
	std::vector<GeometryVector> ks;	//	Wavevectors
	GeometryVector prevBasis[::MaxDimension];	//	An array of basis vectors of the fundamental cell.
	DimensionType d = 0;

	if (Verbosity > 1)
		std::cout << "Computing tildeChik(k)";
	progress_display pd(NumConfigs);

	for (size_t i = 0; i < NumConfigs; i++) {
		//Get a Dispersion object
		const Dispersion CurrConfig = GetConfiguration(i);
		
		//For the first configuration,
		//obtain wavevectors.
		//If dK = 0, assign a new value.
		if (i == 0) {


			d = CurrConfig.GetDimension();
			ks = GetKs(CurrConfig, CircularKmax, LinearKMax, SampleProbability);	
		}
		else {
			//GetKs() function can be very time-consuming.
			//So save time, we recycle wavevectors used in the previous configuration, 
			//if the current one has the same basis vectors.
			bool SameBasis = true;
			if (d == CurrConfig.GetDimension()) {
				for (size_t idx = 0; idx < d && SameBasis; idx++)
					SameBasis = prevBasis[idx] == CurrConfig.GetBasisVector(idx);
			}
			else {
				std::cerr << "Dimensions of configurations are inconsistent." << std::endl;
				return; 
			}
			
			//if the basis vectors are different, obtain a new set of basis vectors.
			if (SameBasis == false)
				ks = GetKs(CurrConfig, CircularKmax, LinearKMax, SampleProbability);
		}

		//Parallelize computations with respect to wavevectors.
		signed long end = ks.size();
#pragma omp parallel for schedule(guided)
		for (signed long j = 0; j < end; j++) {

			double s = SpectralDensity_Dispersion(CurrConfig, ks[j]);
			double k2 = ks[j].Modulus2();
			size_t Bin = std::min(static_cast<size_t>(std::floor(std::sqrt(k2) / dK)), NumBin - 1);

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
		//save the basis vectors of the current configuration.
		for (size_t idx = 0; idx < d; idx++)
			prevBasis[idx] = CurrConfig.GetBasisVector(idx);
		pd++;
	}

	//Summarize data in vSkBin
	for (auto iter = vSkBin.begin(); iter != vSkBin.end(); iter++)
	{
		if (iter->Sum1 != 0.0)
		{
			GeometryVector temp(4);
			temp.x[0] = iter->SumK1 / iter->Sum1;
			temp.x[1] = iter->SumS / iter->Sum1;
			double var = (iter->SumK2 / (iter->Sum1) - temp.x[0] * temp.x[0]) / (iter->Sum1);
			temp.x[2] = (var > 0) ? std::sqrt(var) : 0;
			temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1)); // I modified it
			Result.push_back(temp);
		}
	}

	if (Verbosity>2)
		std::cout << "done!\n";	
}


void DirectionalSpectralDensity_Dispersion(std::function<const Dispersion(size_t)> GetConfigsFunction, size_t NumConfigs,
	long limit, std::vector<GeometryVector> & ks, std::vector<GeometryVector> & chik, double dk){

	ks.clear();
	chik.clear();

	GeometryVector prevBasis[::MaxDimension];
	DimensionType dim;
	std::vector<SkBin> Bins;

	if (NumConfigs == 0) {
		std::cerr << "DirectionalSpectralDensity_nonspherical:: no configurations \n";
		return;
	}
	else {

		for (size_t i = 0; i < NumConfigs; i++) {
			const Dispersion CurrConfig = GetConfigsFunction(i);

			//Check wavevectors;
			//For the first configuration, obtain ks.
			if (i == 0) {
				dim = CurrConfig.GetDimension();
				if (dim == 0) {
					std::cerr << "Realizations should have positive dimensions.\n";
					return;
				}

				for (size_t j = 0; j < dim; j++)
					prevBasis[j] = CurrConfig.GetReciprocalBasisVector(j);

				//Obtain nontrivial wavevectors (parity pairs are excluded).
				{
					int num = 1;	/* integer (2*limit+1)^dim, avoids pow() double round-off/overflow */
					for (size_t j = 0; j < dim; j++)
						num *= (2 * limit + 1);
					n_naryNumber rel_coord(dim, 2 * limit + 1);

					for (int idx = 0; idx < num; idx++) {
						GeometryVector k;	k.SetDimension(dim);
						bool nontrivial = true;
						int component;

						for (int j = 0; j < dim && nontrivial; j++) {
							component = rel_coord.GetNthDigit(j) - limit;
							if (component > 0)
								nontrivial = true;
							else if (component < 0)
								nontrivial = false;

							k = k + component * prevBasis[j];
						}
						if (nontrivial)
							ks.push_back(k);

						rel_coord++;
					}
				}
				//data
				Bins.assign(ks.size(), SkBin());
			}
			else {
				//Check whether this system has identical basis vectors as the previous ones.

				if (dim != CurrConfig.GetDimension()) {
					std::cerr << "Realizations should have identical space dimensions.\n";
					return;
				}

				bool sameBasis = true;
				for (int j = 0; j < dim; j++)
					sameBasis &= (prevBasis[j] == CurrConfig.GetReciprocalBasisVector(j));

				if (!sameBasis) {
					std::cerr << "Realizations should have the same basis vectors\n";
					return;
				}

			}

			if (dk == 1.0) {

				int end = std::min(ks.size(), Bins.size());
#pragma omp parallel for schedule (guided)
				for (int n = 0; n < end; n++) {
					double s = SpectralDensity_Dispersion(CurrConfig, ks[n]);

#pragma omp atomic
					Bins[n].Sum1++;
#pragma omp atomic
					Bins[n].SumS += s;
#pragma omp atomic
					Bins[n].SumS2 += s * s;
				}
			}
			else {
				std::cerr << "It isn't implemented\n";
				return;
			}



		}

		//Summarize data
		chik.resize(Bins.size(), GeometryVector());
		for (size_t i = 0; i < Bins.size(); i++) {
			chik[i].x[0] = Bins[i].SumS / Bins[i].Sum1;
			chik[i].x[1] = Bins[i].SumS2 / Bins[i].Sum1 - chik[i].x[0] * chik[i].x[0];
			chik[i].x[1] = std::sqrt(std::abs(chik[i].x[1]) / Bins[i].Sum1);
		}
	}
}
