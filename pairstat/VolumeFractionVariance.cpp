/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file VolumeFractionVariance.cpp
 *	\brief Implementation of the local volume-fraction variance sigma_V^2(R) of
 *	sphere packings (split out of NumberVariance.cpp). */

#include "VolumeFractionVariance.h"
#include <etc.h>		// v2(), HyperSphere_Volume() — used by the local-volume helpers
#include <RandomGenerator.h>	// random sampling-point generation
#include <ctime>


// ==========================================================================
//   Local volume-fraction variance of sphere packings
//   (ported from the full top-level NumberVariance; I_1/I_config are the
//    per-window local-volume helpers, kept file-local here.)
// ==========================================================================

// local particle-phase volume inside a spherical window of radius dR*(i+1)
// centered at x0_relative; I[i] is the cumulative local volume up to that radius.
static void I_1(const SpherePacking & packing, std::vector<double> &I, double Rmax, double dR, const GeometryVector & x0_relative) {
	I.clear();
	I.resize((size_t)floor(round(Rmax / dR)));

	double r_M = packing.GetMaxRadius();
	DimensionType d = packing.GetDimension();

	packing.IterateThroughNeighbors(x0_relative, Rmax + r_M,
		[&I, dR, &packing, d](const GeometryVector & x_i, const GeometryVector & LatticeShift, const signed long * PereodicShift, const size_t i) ->void {
		double abs_x_i = sqrt(x_i.Modulus2());
		double r_prt = packing.GetCharacteristics(i);
		size_t rU = std::min(I.size() - 1, (size_t)floor((abs_x_i + r_prt) / dR));
		size_t rL = (size_t)std::max(0.0, floor((abs_x_i - r_prt) / dR));

		double v_int_prev = 0;
		double v_int_curr = 0;
		double r_win = 0;
		for (size_t j = rL; j <= rU; j++) {
			r_win = dR*(j + 1.0);
			v_int_curr = v2(d, r_win, r_prt, abs_x_i);
			I[j] += (v_int_curr - v_int_prev);
			v_int_prev = v_int_curr;
		}
	});

	for (size_t i = 1; i < I.size(); i++)
		I[i] += I[i - 1];
}

// phi_inter[i].x = { <phi(R_i)>, <phi(R_i)^2>, sigma_V^2(R_i), <phi_theor> }
static void I_config(const SpherePacking & packing, std::vector<GeometryVector> &phi_inter, double Rmax, double dR, const std::vector<GeometryVector> & x0_relative) {
	double phi_theor = packing.PackingFraction();
	double N = x0_relative.size();
	DimensionType d = packing.GetDimension();
	phi_inter.clear();
	{
		GeometryVector t;	t.Dimension = 4;
		phi_inter.resize((size_t)ceil(Rmax / dR), t);
		for (size_t i = 0; i < phi_inter.size(); i++)
			phi_inter[i].x[3] = phi_theor*HyperSphere_Volume(d, dR*(i + 1.0));
	}

	for (size_t n = 0; n < x0_relative.size(); n++) {
		std::vector<double> I_n;
		I_1(packing, I_n, Rmax, dR, x0_relative[n]);

		for (size_t i = 0; i < I_n.size(); i++) {
			double dV = I_n[i];
			phi_inter[i].x[0] = (1.0 - 1.0 / (n + 1.0))* phi_inter[i].x[0] + dV / (n + 1.0);
			phi_inter[i].x[1] = (1.0 - 1.0 / (n + 1.0))* phi_inter[i].x[1] + dV*dV / (n + 1.0);
		}
	}

	for (size_t i = 0; i < phi_inter.size(); i++) {
		double R = dR*(i + 1.0);
		double Vr = HyperSphere_Volume(d, R);
		phi_inter[i].x[0] /= Vr;
		phi_inter[i].x[1] /= Vr*Vr;
		phi_inter[i].x[2] = N / (N - 1.0)*(phi_inter[i].x[1] - phi_inter[i].x[0] * phi_inter[i].x[0]);
		phi_inter[i].x[3] /= Vr;
	}
}

//result[i].x = { R, sigma_V^2(R), error, running average }
void MC_vfv(const std::function<SpherePacking(size_t i)> &GetConfig, size_t numConfig, std::vector<GeometryVector> & result, double Rmax, double dR, std::vector<GeometryVector> &samplingPoints, size_t num_threads) {
	int max = (int)ceil(Rmax / dR);
	size_t numSampling = samplingPoints.size();
	std::vector<GeometryVector> sum;
	GeometryVector t;	t.Dimension = 4;
	//---Initialization
	{
		result.clear();
		sum.reserve(max);
		result.reserve(max);
		for (int i = 0; i < max; i++) {
			result.push_back(t);
			sum.push_back(t);
			result[i].x[0] = dR*(i + 1);
		}
	}
	double count = 0;

	if (num_threads > 1) {
#pragma omp parallel for schedule (guided) shared(count) num_threads(num_threads)
		for (int i = 0; i < (int)numConfig; i++) {
			std::vector<GeometryVector> phi1;
			SpherePacking conf(GetConfig(i));
			I_config(conf, phi1, Rmax, dR, samplingPoints);
#pragma omp critical
			{
				count++;
				for (size_t j = 0; j < phi1.size(); j++) {
					sum[j].x[0] = (1.0 - 1.0 / count)*sum[j].x[0] + phi1[j].x[0] / count;
					sum[j].x[1] = (1.0 - 1.0 / count)*sum[j].x[1] + phi1[j].x[1] / count;
					if (numConfig > 1)
						sum[j].x[2] = (1.0 - 1.0 / count)*sum[j].x[2] + phi1[j].x[2] * phi1[j].x[2] / count;
				}
			}
		}
	}
	else {
		for (size_t i = 0; i < numConfig; i++) {
			std::vector<GeometryVector> phi1;
			SpherePacking conf(GetConfig(i));
			I_config(conf, phi1, Rmax, dR, samplingPoints);
			{
				count++;
				for (size_t j = 0; j < phi1.size(); j++) {
					sum[j].x[0] = (1.0 - 1.0 / count)*sum[j].x[0] + phi1[j].x[0] / count;
					sum[j].x[1] = (1.0 - 1.0 / count)*sum[j].x[1] + phi1[j].x[1] / count;
					if (numConfig > 1)
						sum[j].x[2] = (1.0 - 1.0 / count)*sum[j].x[2] + phi1[j].x[2] * phi1[j].x[2] / count;
				}
			}
		}
	}

	//analyze data
	for (size_t i = 0; i < result.size(); i++) {
		double av = sum[i].x[0];
		result[i].x[1] = (sum[i].x[1] - av*av) * (1.0 + 1.0 / (numConfig*numSampling - 1.0));

		if (numConfig > 1)
			result[i].x[2] = sqrt(std::abs(sum[i].x[2] - result[i].x[1] * result[i].x[1]) / (numConfig - 1.0));
		//running average
		if (i == 0)
			result[i].x[3] = result[i].x[1];
		else
			result[i].x[3] = (1.0 - 1.0 / (i + 1.0))*result[i - 1].x[3] + result[i - 1].x[1] / (i + 1.0);
	}
}

//with random sampling points
void MC_vfv(const std::function<SpherePacking(size_t i)> &GetConfig, size_t numConfig, std::vector<GeometryVector> & result, double Rmax, double dR, size_t numSampling, size_t num_threads) {
	std::cout << "Computing Sigma_V^2(R)\n";
	std::vector<GeometryVector> samplingPoints;
	RandomGenerator g(std::time(NULL));
	DimensionType d = GetConfig(0).GetDimension();
	GeometryVector x(static_cast<int> (d));

	for (size_t i = 0; i < numSampling; i++) {
		for (int j = 0; j < d; j++)
			x.x[j] = g.RandomDouble();
		samplingPoints.push_back(x);
	}
	MC_vfv(GetConfig, numConfig, result, Rmax, dR, samplingPoints, num_threads);
}
