#include "FFH.h"
#include <GeometryVector.h>
#include <vector>

/** \brief Build a random checkerboard two-phase dispersion. (CLAUDE)
 *  \param[out] dispersion  The resulting two-phase dispersion.
 *  \param num_domains      Number of checkerboard domains per side (must be even).
 *  \param num_pixels       Number of grid pixels per side.
 *  \param rng              Random number generator. */
void RandomCheckerBoard(TwoPhaseDispersion &, int num_domains, int num_pixels, RandomGenerator & rng);

/** \brief Build a periodic checkerboard two-phase dispersion. (CLAUDE)
 *  \param[out] dispersion  The resulting two-phase dispersion.
 *  \param num_domains      Number of checkerboard domains per side (must be even).
 *  \param num_pixels       Number of grid pixels per side. */
void PeriodicCheckerBoard(TwoPhaseDispersion &, int num_domains, int num_pixels);

/** \brief Build a periodic alternating-strip two-phase dispersion. (CLAUDE)
 *  \param[out] dispersion  The resulting two-phase dispersion.
 *  \param num_domains      Number of strip domains per side (must be even).
 *  \param num_pixels       Number of grid pixels per side. */
void PeriodicStrips(TwoPhaseDispersion &, int num_domains, int num_pixels);

size_t FFH_verbosity = 2;
#include <StructureFactor.h>
#include <vector>
#include <functional>

#include "FFH_padded.h"

/** \brief A function to measure the four-point correlation function for 1D point patterns at a given set of points,
 * i.e., x_{21} = x0.x[0] + x.x[0] * r, x_{31} = x0.x[1] + x.x[1] * r, x_{41} = x0.x[2] + x.x[2] * r,
 * \param[in]	getc	A function to generate configurations.
 * \param[in]	num_config	The number of configurations.
 * \param[in]	x		direction of increase.
 * \param[in]	dx		increment in r.
 * \param[in]	Rmax	largest value of r.
 * \param[in,out]	result	result[i].x[0] = r, result[i].x[1] = g4, result[i].x[2] = standard error in g4
 * \param[in]	x0		The origin. */
void Measure_g4(std::function<const Configuration(size_t i)> getc, size_t num_config,
	const GeometryVector & x, double dx, double Rmax, std::vector<GeometryVector> & result, const GeometryVector & x0 = GeometryVector(0, 0, 0)) {

	result.clear();
	result.resize((size_t)ceil(Rmax / dx), GeometryVector(0, 0, 0, 0));
#pragma omp parallel for
	for (int i = 0; i < num_config; i++) {
		Configuration c = getc(i);
		std::vector<GeometryVector> sum;
		sum.resize((size_t)ceil(Rmax / dx), GeometryVector(0, 0, 0, 0));

		for (int j = 1; j <= sum.size(); j++) {
			double prob = 0.0L;
			GeometryVector x21(x.x[0] * dx *(double)j + x0.x[0]),
				x31(x.x[1] * dx*(double)j + x0.x[1]), x41(x.x[2] * dx*(double)j + x0.x[2]);
			x21 = c.CartesianCoord2RelativeCoord(x21);
			x31 = c.CartesianCoord2RelativeCoord(x31);
			x41 = c.CartesianCoord2RelativeCoord(x41);
			for (int prt = 0; prt < c.NumParticle(); prt++) {
				bool find_prt = true;
				GeometryVector x1 = c.GetRelativeCoordinates(prt);
				GeometryVector x2 = x1 + x21,
					x3 = x1 + x31, x4 = x1 + x41;
				if (!((x2 == x3) || (x3 == x4) || (x4 == x2))) {

					find_prt = (c.NearestParticleDistance(x2) < 0.5*dx)
						&& (c.NearestParticleDistance(x3) < 0.5*dx)
						&& (c.NearestParticleDistance(x4) < 0.5*dx);
				}
				else
					find_prt = false;

				if (find_prt)
					prob++;
			}

			sum[j - 1].x[0] = prob / (double)c.NumParticle() / (dx*dx*dx);
		}

		for (int j = 0; j < result.size(); j++) {
#pragma omp atomic
			result[j].x[1] += sum[j].x[0];
#pragma omp atomic
			result[j].x[2] += sum[j].x[0] * sum[j].x[0];
		}
	}

	for (int j = 0; j < result.size(); j++) {
		result[j].x[0] = dx * (double)(j + 1);
		result[j].x[1] /= (double)(num_config);
		result[j].x[2] /= (double)(num_config);
		result[j].x[2] -= result[j].x[1] * result[j].x[1];
		result[j].x[2] = sqrt(result[j].x[2] / (double)(num_config - 1));
	}
}


/** \brief A function to measure the three-point correlation function for 1D point patterns at a given set of points,
 * i.e., x_{21} and x_{31} in [-R_max, R_max]^2,
 * \param[in]	getc	A function to generate configurations.
 * \param[in]	num_config	The number of configurations.
 * \param[in]	dx		increment in r and bin size.
 * \param[in]	Rmax	largest value of r.
 * \param[in,out]	g3	result[i].x[0] = x_{21}, result[i].x[1] = x_{31}, result[i].x[2] = g3, result[i].x[3] = standard error in g3. */
void Measure_g3(std::function<const Configuration(size_t i)> getc, size_t num_config,
	double dx, double Rmax, std::vector<GeometryVecotr> & g3) {

	int num_half = floor(Rmax / dx);
	int num_points = 2*num_half+1;
	g3.resize((size_t)num_points*num_points, GeometryVector(0,0,0,0));
	
#pragma omp parallel for
	for (int i = 0; i < num_config; i++) {
		Configuration c = getc(i);
		/*a point index for a configuration. */
		int idx = 0;
		/*int idx2, idx3;*/
		for (int j = -num_half; j <= num_half; j++) {
			double x2 = dx * (double)j;
			/* idx2 = j;// (num_half + 1) + j;
			if (idx2 < 0)
				idx2 += num_points; */
			for (int k = -num_half; k <= num_half; k++) {
				double x3 = dx * (double)k;
			/* 	idx3 = k;// (num_half + 1) + k;
				if (idx3 < 0)
					idx3 += num_points; */

				GeometryVector X21(static_cast<double>(dx*(double)j));
				GeometryVector X31(static_cast<double>(dx*(double)k));
				GeometryVector X2(X21);
				GeometryVector X3(X31);

				X21 = c.CartesianCoord2RelativeCoord(X21);
				X31 = c.CartesianCoord2RelativeCoord(X31);
				double prob = 0.0L;
				for (int prt = 0; prt < c.NumParticle(); prt++) {

					GeometryVector X1 = c.GetRelativeCoordinates(prt);
					bool find_prt = false;

					X2 = X21 + X1;
					X3 = X31 + X1;

					if (j != k && j!=0 && k!=0) {
						find_prt = (c.NearestParticleDistance(X2) < 0.5*dx)
							&& (c.NearestParticleDistance(X3) < 0.5*dx);
					}
					

					if (find_prt)
						prob++;

				}

				prob /= ((double)c.NumParticle() * (dx*dx));
				if(i==0){
#pragma omp atomic
					g3[idx].x[0] = X21.x[0];
#pragma omp atomic					
					g3[idx].x[1] = X31.x[0];	
				}

#pragma omp atomic
				g3[idx].x[2] += prob;
#pragma omp atomic
				g3[idx].x[3] += (prob*prob);
			
				idx++;
			}
		}

	}

#pragma omp parallel for
	for (int idx = 0; idx < g3.num_total_pixels; idx++) {
	
		g3[idx].x[2] /= (double)num_config;
		g3[idx].x[3] /= (double)num_config;
		g3[idx].x[3] = sqrt(abs(g3[idx].x[3] - g3[idx].x[2] * g3[idx].x[2]) / ((double)num_config-1.0) );
		
	}

}


/** \brief Measure the three-point correlation function and store it in a TensorField. (CLAUDE)
 *  \param[in] getc        A function to generate configurations.
 *  \param[in] num_config  The number of configurations.
 *  \param[in] dx          Increment in r and bin size.
 *  \param[in] Rmax        Largest value of r.
 *  \param[in,out] g3      Scalar tensor field receiving g3 and its standard error. */
void Measure_g3(std::function<const Configuration(size_t i)> getc, size_t num_config,
	double dx, double Rmax, TensorField & g3) {

	int num_half = floor(Rmax / dx);
	int num_points = 2*num_half+1;
	g3 = TensorField(2, 0, std::vector<size_t>({ (size_t)num_points, (size_t)num_points }));

#pragma omp parallel for
	for (int i = 0; i < num_config; i++) {
		Configuration c = getc(i);
		int idx2, idx3;
		for (int j = -num_half; j <= num_half; j++) {
			double x2 = dx * (double)j;
			idx2 = j;// (num_half + 1) + j;
			if (idx2 < 0)
				idx2 += num_points;
			for (int k = -num_half; k <= num_half; k++) {
				double x3 = dx * (double)k;
				idx3 = k;// (num_half + 1) + k;
				if (idx3 < 0)
					idx3 += num_points;

				GeometryVector X21(static_cast<double>(dx*(double)j));
				GeometryVector X31(static_cast<double>(dx*(double)k));
				GeometryVector X2(X21);
				GeometryVector X3(X31);

				X21 = c.CartesianCoord2RelativeCoord(X21);
				X31 = c.CartesianCoord2RelativeCoord(X31);
				double prob = 0.0L;
				for (int prt = 0; prt < c.NumParticle(); prt++) {

					GeometryVector X1 = c.GetRelativeCoordinates(prt);
					bool find_prt = false;

					X2 = X21 + X1;
					X3 = X31 + X1;

					if (j != k && j!=0 && k!=0) {
						find_prt = (c.NearestParticleDistance(X2) < 0.5*dx)
							&& (c.NearestParticleDistance(X3) < 0.5*dx);
					}
					

					if (find_prt)
						prob++;

				}

				prob /= ((double)c.NumParticle() * (dx*dx));
#pragma omp atomic
				g3[0][idx3 + num_points * idx2][0] += prob;
#pragma omp atomic
				g3[0][idx3 + num_points * idx2][1] += (prob*prob);
			}
		}


	}

#pragma omp parallel for
	for (int idx = 0; idx < g3.num_total_pixels; idx++) {
	
		g3[0][idx][0] /= (double)num_config;
		g3[0][idx][1] /= (double)num_config;
		g3[0][idx][1] = sqrt(abs(g3[0][idx][1] - g3[0][idx][0]* g3[0][idx][0]) / ((double)num_config-1.0) );
		
	}

}

/** \brief compute F(Q) for the Zd
 *  \param d               Spatial dimension.
 *  \param packingfraction Packing fraction of the lattice.
 *  \param k_max           Largest wavenumber Q.
 *  \param k_min           Smallest wavenumber Q (and step).
 *  \param[out] FQ         Resulting (Q, F(Q)) pairs. */
void ComputeAttenuationFn_Zd(int d, double packingfraction, double k_max, double k_min, std::vector<GeometryVector> & FQ) {

	double a = pow(packingfraction / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);
	double phi2 = packingfraction;

		double k_limit;
		if (d == 2)
			k_limit = 1000;
		else if (d == 3)
			k_limit = 200;

		k_limit *= (2.0*pi);

		
		Configuration c = GetUnitCubicBox(d, 1);
		std::vector<GeometryVector> qs = GetKs(c, k_limit, k_limit);

		FQ.clear();
		FQ.reserve((int)floor(k_max / k_min) + 1);

		double factor = -pow(2.0, 0.5*d)*tgamma(0.5*d) / pi;
		for (double Q = k_min; Q <= k_max; Q += k_min) {
			double val_FQ = 0;
#pragma omp parallel for
			for (int i = 0; i < qs.size(); i++) {
				double q = sqrt(qs[i].Modulus2()),
					mtilde = Mtilde(q, a, d),
					val = (mtilde*mtilde*2.0*q*q / (q*q*q*q - 4.0*Q*Q*qs[i].x[0]* qs[i].x[0]));

#pragma omp atomic
				val_FQ += val;
			}
			val_FQ *= (factor*Q*Q);
			FQ.emplace_back(Q, val_FQ, 0.0L);
		}
}

/** \brief Compute the local attenuation function F(Q) for the Zd lattice. (CLAUDE)
 *  \param d               Spatial dimension.
 *  \param packingfraction Packing fraction of the lattice.
 *  \param k_max           Largest wavenumber Q.
 *  \param k_min           Smallest wavenumber Q (and step).
 *  \param[out] FQ         Resulting (Q, F(Q)) pairs. */
void ComputeLocalAttenuationFn_Zd(int d, double packingfraction, double k_max, double k_min, std::vector<GeometryVector> & FQ) {

	double a = pow(packingfraction / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);
	double phi2 = packingfraction;

	double k_limit;
	if (d == 2)
		k_limit = 1000;
	else if (d == 3)
		k_limit = 200;

	k_limit *= (2.0*pi);


	Configuration c = GetUnitCubicBox(d, 1);
	std::vector<GeometryVector> qs = GetKs(c, k_limit, k_limit);

	FQ.clear();
	FQ.reserve((int)floor(k_max / k_min) + 1);

	double factor = -pow(2.0, 0.5*d)*tgamma(0.5*d) / pi;
	for (double Q = k_min; Q <= k_max; Q += k_min) {
		double val_FQ = 0;
#pragma omp parallel for
		for (int i = 0; i < qs.size(); i++) {
			double q = sqrt(qs[i].Modulus2()),
				mtilde = Mtilde(q, a, d),
				val = (mtilde*mtilde*2.0 / (q*q - Q*Q));

#pragma omp atomic
			val_FQ += val;
		}
		val_FQ *= (factor*Q*Q);
		FQ.emplace_back(Q, val_FQ, 0.0L);
	}
}



/** \brief Compute the effective dielectric constant of a cubic lattice via the strong-contrast formula. (CLAUDE)
 *  \param d               Spatial dimension.
 *  \param packingfraction Packing fraction of the lattice.
 *  \param k_max           Largest wavenumber q.
 *  \param k_min           Smallest wavenumber q (and step).
 *  \param epsilon2        Permittivity of the inclusion phase.
 *  \param[out] result     Resulting (q, epsilon_e) pairs. */
void ComputeDielectricCubicLattice(int d, double packingfraction, double k_max, double k_min, double epsilon2, std::vector<GeometryVector> & result) {

	double a = pow(packingfraction / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);
	double phi2 = packingfraction;
	double beta = (epsilon2 - 1.0) / (epsilon2 + d-1.0);
	std::vector<double> Z;
	std::vector<double> k;
	
	{
		//double k_limit = 2.0*pi*pow(5000000, 1.0 / (double)d);
		double k_limit = 2.0*pi*pow(100, 1.0 / (double)d);

		Configuration c = GetUnitCubicBox(d, 1);
		std::vector<GeometryVector> qs = GetKs(c, k_limit, k_limit);
		std::vector<double> k2;
		k2.resize(qs.size());
#pragma omp parallel for
		for (int i = 0; i < qs.size(); i++) {
			k2[i] = qs[i].Modulus2();
		}
		qs.clear();
		std::sort(k2.begin(), k2.end());
		
		int index = 0, count=0;
		for (auto iter = k2.begin(); iter != k2.end()-1; iter++) {
			if (iter == k2.begin()) {
				k.push_back(*iter);
				Z.push_back(0);
				count++;
			}
			else {
				count++;
				if (*iter != *(iter + 1)) {
					k.push_back(*(iter + 1));
					
					Z[Z.size() - 1] = 2*count;
					Z.push_back(0);
					index++;
					count = 0;
				}
			}
		}
		k2.clear();
	
		for (int i = 0; i < k.size(); i++) {
			k[i] = sqrt(k[i]);
		}
	
	}
	result.clear();
	for (double q = k_min; q <= k_max; q += k_min) {
		double A2 = 0, epsilon_e = 0;
#pragma omp parallel for
		for (int i = 0; i < k.size(); i++) {
			double val = Mtilde(k[i], a, d);
			val = phi2 * val*val / (k[i] * k[i] - q*q);
#pragma omp atomic
			A2 += val;
		}
		A2 *= phi2 * q*q / (2.0*pi*pi);
		
		epsilon_e = 1.0 - (double)d*beta*phi2*phi2 / (A2*beta - phi2*(1.0-phi2*beta));
		
		result.emplace_back(q, epsilon_e);
		
	}



}


/** \brief Estimate the pore-size distribution moments of a sphere packing via Monte Carlo sampling. (CLAUDE)
 *  \param src           The sphere packing.
 *  \param[out] Moments  The first two moments of the scaled pore radius.
 *  \param num_sampling  Number of random sample points. */
void MeasurePoreSize(SpherePacking & src, GeometryVector & Moments, int num_sampling = 10000){
	
	double dx = 0.5*src.NearestParticleDistance(0);
	double tolerance = dx * sqrt(src.GetDimension());
	double a = sqrt(src.GetBasisVector(0).Modulus2());
	RandomGenerator rng;
	Moments = GeometryVector(0, 0);

	int count = 0;
#pragma omp parallel
	{
		SpherePacking dig(src);
#pragma omp critical
		{
			int idx = omp_get_thread_num();
			std::cout << "idx for threads = " << idx << "\n";
			RandomGenerator rng(idx);
		}
#pragma omp for
		for (int i = 0; i < num_sampling; i++) {
			GeometryVector x0(static_cast<int>(dig.GetDimension()));
			for (int j = 0; j < x0.Dimension; j++)
				x0.x[j] = rng.RandomDouble();

			double delta = dig.NearestParticleDistance(x0);
			if (delta > tolerance) {
#pragma omp atomic
				count++;
				double scaled_delta = ((delta - tolerance) / a);
#pragma omp atomic
				Moments.x[0] += scaled_delta;
#pragma omp atomic
				Moments.x[1] += scaled_delta * scaled_delta;
			}
		}
	}
	Moments.x[0] /= (double)count;
	Moments.x[1] /= (double)count;
	std::cout << "count = " << count<<"\n";
	std::cout << "End Measurement\n";
}

/** \brief Compute a normalized measure of the curl of the electric field in Fourier space. (CLAUDE)
 *  \param src  A static dielectric problem with computed fields.
 *  \return     The normalized RMS magnitude of the transverse (curl) component. */
double curl(const StaticDielectricProblem & src) {
	TensorField E, kE;
	src.ComputeField_Fourier(E);
	InnerProduct(kE, src.ks, E, 1);
	//kE.Print_shifted_origin("test", 0, 0);
	//kE.Print_shifted_origin("test", 0, 1);
	//src.ks.Print("test_ks", 0, 1);
	double value = 0.0L, mean = 0.0L;
#pragma omp parallel for
	for (int x = 0; x < kE.num_total_pixels; x++) {
		double k2 = 0.0L, E2 = 0.0L, kE2 = 0.0L, temp = 0.0L;
		for (int j = 0; j < src.d; j++) {
			k2 += Modulus2(src.ks[j][x]);
		}
		for (int j = 0; j < src.d; j++) {
			E2 += Modulus2(E[j][x]);
		}
		kE2 = Modulus2(kE[0][x]);
		
		if(k2 > 0)
			temp = (E2 - kE2 / k2);
		if ((log10(E2) - log10(abs(temp))) > 15.0) {
			temp = 0;
		}
		
		k2 = k2*temp;
		



		kE[0][x][0] = k2;
		kE[0][x][1] = 0;
		//if (k2 < 0){
		//	ofile << "something is wrong\n";
		//else {
#pragma omp atomic
		value += k2;
#pragma omp atomic
		mean += E2;
		//}
	}
	value /= mean;
//	kE.Print_shifted_origin("test", 0, 0);
//	kE.Print_shifted_origin("test", 0, 1);
	return sqrt(value);
}

/** \brief Test/scratch driver exercising the various homogenization solvers. (CLAUDE)
 *  \return 0 on success. */
int main() {
	char tempstring[1000] = {};
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	char dir[] = "./results/n1/";
	omp_set_num_threads(4);
	RandomGenerator rng(0);
	fftw_init_threads();
	std::string dump;

	int mode;
	ofile << "mode (1:dielectric/2:elastic)= ";
	ifile >> mode;

	{       /* Periodic Checkerboard array */
		int d = 2;
		FFH_verbosity = 2;
		double epsilon1 = 1.0L, phi = 0.5;
		double epsilon2 = 10.0L;
		int num_pixels = 11;
		std::vector<GeometryVector> epsilon_eff, errors;


		

		for (num_pixels = 2; num_pixels < 30; num_pixels++) {
			//ifile >> epsilon2;
			SpherePacking c(GetUnitCubicBox(d, 1), 0.0);
			TwoPhaseDispersion packing(c, num_pixels);      //It gives a discretized version of a sphere packing.
			PeriodicCheckerBoard(packing, 2, num_pixels);

			epsilon2 = 10.0L;
			StaticDielectricProblem test(packing, epsilon1, epsilon2);
			test.SetMethod(true);


			// define the problem.
			ConstantTensor E0(d, 1);
			E0[0] = 1.0L;
			test.DefineProperties(epsilon1, epsilon2);
			test.Prepare(E0);

			// Iteratively compute the resulting electric field
			test.Run(1e-4, 1e-10);
			test.Print("periodic??");

			std::complex<double> epsilone = test.ComputeEffectiveEpsilon();
			ofile << "epsilon2/epsilon1 = " << epsilon2 << "\n";
			ofile << "Effective dielectric constant = " << epsilone << "\n";
			ofile << "next...\n";
			epsilon_eff.emplace_back(num_pixels, sqrt(epsilon2), epsilone.real(), epsilone.imag());
//			errors.emplace_back(epsilon2, test.Tol_Field, test.Tol_div, test.m);

			sprintf(name, "strong_%.2f", epsilon2);
			test.Print(name);
			test.Print_Fourier(name);


		}
		sprintf(name, "PeriodicCheckerBoard_epsilone_strong");
		WriteFunction(epsilon_eff, name);
//		sprintf(name, "PeriodicCheckerBoard_errors_strong");
//		WriteFunction(errors, name);
	}
	return 0;
	{       /* Periodic Checkerboard array */
		int d = 2;
		FFH_verbosity = 2;
		double epsilon1 = 1.0L, phi = 0.5;
		double epsilon2 = 10.0L;
		int num_pixels = 11;
		std::vector<GeometryVector> epsilon_eff, errors;


		SpherePacking c(GetUnitCubicBox(d, 1), 0.0);
		TwoPhaseDispersion packing(c, num_pixels);      //It gives a discretized version of a sphere packing.
		PeriodicCheckerBoard(packing, 2, num_pixels);


		StaticDielectricProblem test(packing, epsilon1, epsilon2);
		test.SetMethod(true);
		for (epsilon2 = 10.0L; epsilon2 <= 100L; epsilon2 += 10) {
			//ifile >> epsilon2;
			ofile << "epsilon2 = " << epsilon2 << "\n";


			// define the problem.
			ConstantTensor E0(d, 1);
			E0[0] = 1.0L;
			test.DefineProperties(epsilon1, epsilon2);
			test.Prepare(E0);

			// Iteratively compute the resulting electric field
			test.Run(1e-4, 1e-10);
			test.Print("periodic??");

			std::complex<double> epsilone = test.ComputeEffectiveEpsilon();
			ofile << "epsilon2/epsilon1 = " << epsilon2 << "\n";
			ofile << "Effective dielectric constant = " << epsilone << "\n";
			ofile << "next...\n";
			epsilon_eff.emplace_back(epsilon2, epsilone.real(), epsilone.imag());
			errors.emplace_back(epsilon2, test.Tol_Field, test.Tol_div, test.m);

			sprintf(name, "strong_%.2f", epsilon2);
			test.Print(name);
			test.Print_Fourier(name);
		}
		sprintf(name, "PeriodicCheckerBoard_epsilone_strong");
		WriteFunction(epsilon_eff, name);
		sprintf(name, "PeriodicCheckerBoard_errors_strong");
		WriteFunction(errors, name);
	}
	return 0;


	if (mode == 1) {
		ofile << "\nDynamic Dielectric\n";
		/* periodic array of spheres (dynamic dielectric). */
		FFH_verbosity = 2;
		int d = 3;
		std::string method("strong"), filename ("Z3");
		int pixels = 201, num_threads = 4;
		double epsilon1 = 1.0L, phi = 0.25;
		double epsilon2 = 10.01L; double k = 1.0;
		ofile << "compute the effective dielectric constant of 3D two-phase media (export data via ReadConfiguration)\n";
		ofile << "planar external electric field incident along the x-axis\n";
		ofile << "Green's function is modified\n";
		ofile << "eletric field is in the y-direction\n";
		//ofile << "file name = "; ifile >> filename; ofile << filename << "\n";
		ofile << "num pixles = "; ifile >> pixels;	ofile << pixels << "\n";
		ofile << "epsilon2 = "; ifile >> epsilon2;	ofile << epsilon2 << "\n";
		ofile << "phi2 = "; ifile >> phi;		ofile << phi << "\n";
		ofile << "num threads = "; ifile >> num_threads; ofile << num_threads << "\n"; omp_set_num_threads(num_threads);
		double k_min = 0.01, k_max = 10, dk = 0.1;

		ofile << "k_min = "; ifile >> k_min;	ofile << k_min;
		ofile << "k_max = "; ifile >> k_max;	ofile << k_max;
		ofile << "dk = "; ifile >> dk; 		ofile << dk;
		double a = pow(phi / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);

		SpherePacking c2;
		{
			Configuration c = GetUnitCubicBox(d);
			c.Insert("a", GeometryVector(0.5, 0.5, 0.5));
			c2 = SpherePacking(c, a);
			//ReadConfiguration(c2, filename);
		}

		/* c4 is zero-padding...*/
		TwoPhaseDispersion packing(c2, pixels);
		packing.ComputeIndicator();
		DynamicDielectricProblem test(packing, epsilon1, epsilon2);
		test.UseOptimalReference(false);
		std::vector<GeometryVector> result;
		std::vector<GeometryVector> error;

		test.DefineProperties(epsilon1, epsilon2, GeometryVector(static_cast<int>(d)));
		sprintf(tempstring, "Dimension = %d, Packing fraction = %0.3f, Particle number = %d, \n# resolutions = %d^%d", d, phi, packing.packing.NumParticle(), pixels, d);
		std::string system(tempstring);
		sprintf(tempstring, "epsilon1 = %0.3f, epsilon2 = %0.3f, epsilon0 = %0.3f", epsilon1, epsilon2, test.epsilon0);
		std::string phase_properties(tempstring);
		sprintf(tempstring, "method = %s + regularized", method.c_str());
		std::string solver(tempstring);

		std::vector<std::string> label1 = { system, phase_properties, tempstring,
			"\n#k1", "Re[epsilon_e]", "Im[epsilon_e]" };
		std::vector<std::string> label2 = { system, phase_properties, tempstring,
			"\n#k1", "m", "Tol_field", "Tol_div" };

		//system = std::string("Z3_sphere-packing_");

		for (double k = k_min; k < k_max; k += dk) {
			ofile << "k = " << k << "\n";

			std::complex<double> epsilon_eff;
			GeometryVector k1;
			{
				GeometryVector e(static_cast<int>(d));
				e.x[0] = 1.0;
				k1 = k * e;
			}

			ofile << "epsilon2 = " << epsilon2 << "\n";
			ofile << "k1 a = " << k1 << "\n";

			ConstantTensor E0(d, 1);
			E0[1] = 1.0L;
			test.SetMethod(method.compare("strong") == 0);
			test.DefineProperties(epsilon1, epsilon2, k1);
			test.Prepare(E0);
			//test.ComputeG_Fourier_test3();
			test.Run(1e-2, 1e-9);
			/*
			epsilon_eff = test.ComputeEffectiveEpsilon();
			ofile << "epsilon_eff = " << epsilon_eff << "\n";
			if (test.m <= num_limit && test.IsCovergeningSolution()) {
				method2.emplace_back(k, epsilon_eff.real(), epsilon_eff.imag());
			}*/

			epsilon_eff = test.ComputeEffectiveEpsilon();

			if (test.m <= num_limit && test.IsCovergeningSolution() && test.Tol_Field < 1e-4) {
				result.emplace_back(k, epsilon_eff.real(), epsilon_eff.imag());
			}
			error.emplace_back(k, test.m, test.Tol_Field, test.Tol_div);

			sprintf(tempstring, "%s_%0.2f_epsilon2__%0.2f__%s_%d", filename.c_str(), phi, epsilon2, method.c_str(), pixels);
			WriteFunction(result, tempstring, label1);
			//sprintf(tempstring, "%s_%0.2f_epsilon2__%0.2f__%s_%d_error", filename.c_str(), phi, epsilon2, method.c_str(), pixels);
			//WriteFunction(error, tempstring, label2);

			ofile << "k_1 = " << k1 << ",\tepsilon_eff = " << epsilon_eff << "\n";



		}
	}
	else if (mode == 2) {
		ofile << "\nDynamic Elastic\n";
		/* periodic array of spheres. */
		FFH_verbosity = 2;
		int d = 3;
		std::string method("strong");
		int pixels = 151, num_threads = 4;

		double phi;
		double K1 = 1.0, G1 = 3.0 / 8.0;
		double K2 = 2.0, G2 = K2 * G1;


		ofile << "compute the effective dielectric constant of 3D cubic sphere packing\n";
		ofile << "planar external elastic field incident along the x-axis\n";
		ofile << "num pixles = "; ifile >> pixels;      ofile << pixels << "\n";
		ofile << "K2 = "; ifile >> K2;      ofile << K2 << "\n";
		ofile << "K1 = "; ifile >> K1;      ofile << K1 << "\n";
		ofile << "G2 = "; ifile >> G2;      ofile << G2 << "\n";
		ofile << "G1 = "; ifile >> G1;      ofile << G1 << "\n";
		ofile << "phi2 = "; ifile >> phi;               ofile << phi << "\n";
		ofile << "num threads = "; ifile >> num_threads; ofile << num_threads << "\n"; omp_set_num_threads(num_threads);
		double k_min = 0.01, k_max = 10, dk = 0.1;

		ofile << "k_min = "; ifile >> k_min;    ofile << k_min;
		ofile << "k_max = "; ifile >> k_max;    ofile << k_max;
		ofile << "dk = "; ifile >> dk;          ofile << dk;
		double a = pow(phi / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);

		SpherePacking c(GetUnitCubicBox(d, 1), 0.0), c2, c4;
		{
			GeometryVector x(static_cast<int> (d));
			for (int i = 0; i < d; i++)
				x.x[i] = 0.5;
			c.Insert(a, x);
		}

		TwoPhaseDispersion packing(c, pixels);  //It gives a discretized version of a sphere packing.
		packing.ComputeIndicator();
		DynamicElasticProblem test(packing, K1, G1, K2, G2);
		test.SetMethod(method.compare("strong") == 0);
		test.UseOptimalReference(false);
		std::vector<GeometryVector> method1;
		std::vector<GeometryVector> method2;

		sprintf(tempstring, "Dimension = %d, Packing fraction = %0.3f, Particle number = %d, \n# resolutions = %d^%d", d, phi, packing.packing.NumParticle(), pixels, d);
		std::string system(tempstring);
		sprintf(tempstring, "K1 = %0.3f, G1 = %0.3f, K2 = %0.3f, G2 = %0.3f", K1, G1, K2, G2);
		std::string phase_properties(tempstring);
		sprintf(tempstring, "method = %s", method.c_str());
		std::string solver(tempstring);

		//k = 0.09 * 2 * pi; k < 1.0 * 2 * pi; k += 0.025 * 2 * pi
		std::vector<std::string> label1 = { system, phase_properties, tempstring,
				"a = " + std::to_string(a),
				"\n#k1", "Re[K_e]", "Im[K_e]" };
		std::vector<std::string> label2 = { system, phase_properties, tempstring,
				"a = " + std::to_string(a),
				"\n#k1", "Re[G_e]", "Im[G_e]" };
		sprintf(&system[0], "Z%d_packing_%0.2f_%0.2f_%0.2f", d, phi, K2/K1, G2/G1);


		for (double k = k_min; k < k_max; k += dk)
		{
			ofile << "k = ?";// ifile >> dump;
			ofile << k;

			GeometryVector k1;
			{
				GeometryVector e(static_cast<int>(d));
				e.x[0] = 1.0;
				k1 = k * e;
			}
			/* longitudinal strain. */
			ConstantTensor E0(d, 2);
			E0[0] = 1.0L;

			test.DefineProperties(K1, G1, K2, G2);
			test.SetWavevector(k1);
			test.Prepare(E0);
			test.Run(1e-8, 1e-5);
			std::vector<std::complex<double>> result;
			test.ComputeEffectiveModuli(result);

			/* longitudinal strain -> bulk modulus. */
			if (test.m <= num_limit && test.IsCovergeningSolution()) {
				method1.emplace_back(k, result[0].real(), result[0].imag());
			}


			/* transverse strain. */
			E0[0] = 0.0L;
			E0[1 + 0 * d] = 1.0L;   E0[0 + d * 1] = 1.0L;

			test.DefineProperties(K1, G1, K2, G2);
			test.SetWavevector(k1);
			test.Prepare(E0);
			test.Run(1e-6, 1e-5);
			test.ComputeEffectiveModuli(result);

			/* transverse strain -> shear modulus. */
			if (test.m <= num_limit && test.IsCovergeningSolution()) {
				method2.emplace_back(k, result[1].real(), result[1].imag());
			}

			sprintf(tempstring, "%s_bulk_%s", system.c_str(), method.c_str());
			WriteFunction(method1, tempstring, label1);
			sprintf(tempstring, "%s_shear_%s", system.c_str(), method.c_str());
			WriteFunction(method2, tempstring, label2);
		}
	}
	return 0;



	{
		double G1 = 1.0L, K1 = 1.0L;
		double G2, K2;
		double phi;
		ofile << "G2 = "; ifile >> G2; ofile << "\n";
		ofile << "K2 = "; ifile >> K2; ofile << "\n";
		ofile << "phi = "; ifile >> phi; ofile << "\n";
		std::string method("strong");
		bool use_sc = (method.compare("strong") == 0);
		int d = 3;

		int num_pixels, num_threads;
		ofile << "num_pixels = "; ifile >> num_pixels; ofile << "\n";
		ofile << "num_threads = "; ifile >> num_threads; ofile << "\n"; omp_set_num_threads(num_threads);
		std::vector<GeometryVector> summary, error_hydro, error_shear;
		FFH_verbosity = 2;

		std::vector<std::complex<double>> moduli;
		double a = pow(phi / (4 * 4 * 4 * HyperSphere_Volume(d, 1.0)), 1.0 / (double)d);
		SpherePacking c(GetUnitCubicBox(d, 1), 0.0);
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					c.Insert(a, GeometryVector(0.25*i, 0.25*j, 0.25*k));
				}
			}
		}
		TwoPhaseDispersion packing(c, num_pixels);
		packing.ComputeIndicator();

		ofile << "size max = " << SIZE_MAX << "\n";
		ofile << pow(num_pixels, 3)*size_fftw << "\n";

		//packing.ComputeIndicator2();

		ofile << "K2/K1 = " << K2 << ", G2/G1 = " << G2 << "\n";
		std::complex<double> Ke, Ge;
		Homo_ElasticProblem_static test(packing, K1, G1, K2, G2);
		//test.SetStrongContrast(use_sc);
		test.RedefineProblem(packing, K1, G1, K2, G2);
		test.SetStrongContrast(use_sc);

		// compression
		test.Prepare(HydrostaticStrain(d));
		test.Run(5e-7, 1e-10);
		test.ComputeEffectiveProperties(moduli);
		ofile << "Under hydrostatic compress\n";
		ofile << "Ke = " << moduli[0] << "\n";
		ofile << "Ge = " << moduli[1] << "\n";
		ofile << "m = " << test.m << ", error = " << test.Tol_Field << ", div = " << test.Tol_div << "\n";
		// shear
		test.Prepare(DeviatoricStrain(d));
		test.Run();
		test.ComputeEffectiveProperties(moduli);
		ofile << "Under deviatoric strain\n";
		ofile << "Ke = " << moduli[0] << "\n";
		ofile << "Ge = " << moduli[1] << "\n";
		ofile << "m = " << test.m << ", error = " << test.Tol_Field << ", div = " << test.Tol_div << "\n";
	}
	return 0;

	{
		std::vector<GeometryVector> result;
		ComputeLocalAttenuationFn_Zd(3, 0.05, 50, 0.01, result);
		WriteFunction(result, "Z3_packing_0.05__F");
	}
	return 0;

	


	


	{
		double phi = 0.2, a = pow(phi / (HyperSphere_Volume(3, 1.0)), 1.0 / 3.0);
		for (int i = 0; i < 3; i++) {
			Configuration c;
			sprintf(name, "configs/3D_Stealthy_Ka-0.0-1.5__1.0-0.92__L_10.0_Success_%d", i);
			ReadConfiguration(c, name);

			SpherePacking sp (c, a);
			sprintf(name, "configs/3D_Stealthy_%d", i);
			WriteConfiguration(sp, name);

		}
		
	}
	return 0;
	


	{       /* Simple cubic array of sphere packing (static dielectric) */
		int d = 3;
		FFH_verbosity = 2;
		double epsilon1 = 1.0L, phi = 0.50;
		double epsilon2 = 10.0L;
		int num_pixels = 101;// std::atoi(argv[1]);

		std::vector<GeometryVector> epsilon_eff, errors;

		sprintf(tempstring, "Z3-packing-%0.2f-%.2f", phi, epsilon2);
		std::string system(tempstring), method("strong");
		double a = pow(phi / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);
		SpherePacking c(GetUnitCubicBox(d, 1), 0.0);
		{
			GeometryVector pos(static_cast<int>(d));
			for (int i = 0; i < d; i++) {
				pos.x[i] = 0.5;
			}
			c.Insert(a, pos);
		}

		TwoPhaseDispersion packing(c, num_pixels);      //It gives a discretized version of a sphere packing.
		packing.ComputeIndicator();


		StaticDielectricProblem test(packing, epsilon1, epsilon2);
		test.SetMethod(true);
		ConstantTensor E0(d, 1);
		E0[0] = 1.0L;
		//	test.SetMethod(false);
		test.DefineProperties(epsilon1, epsilon2);
		test.Prepare(E0);

		sprintf(tempstring, "Dimension = %d, Packing fraction = %0.3f, Particle number = %d, \n# resolutions = %d^%d", d, phi, packing.packing.NumParticle(), num_pixels, d);
		std::string sys_descrip(tempstring);
		sprintf(tempstring, "epsilon1 = %0.3f, epsilon2 = %0.3f", epsilon1, epsilon2);
		std::string phase_properties(tempstring);
		sprintf(tempstring, "method = %s", method.c_str());
		std::string solver(tempstring);

		std::vector<std::string> label1 = { sys_descrip, phase_properties, tempstring, "The effective dielectric constant according to the tolerance in the divergence-free condition.\n",
				"tolerance", "Re[epsilon_e]", "Im[epsilon_e]" };

		std::vector<std::string> label2 = { sys_descrip, phase_properties, tempstring, "The effective dielectric constant according to the tolerance in the divergence-free condition.\n",
"tolerance", "m", "tol_field", "tol_curl" };
		int m_prev = 0;

		test.Run(1e-5, 1e-9);
		std::complex<double> epsilone = test.ComputeEffectiveEpsilon();
		ofile << "epsilon2/epsilon1 = " << epsilon2 << "\n";
		ofile << "Effective dielectric constant = " << epsilone << "\n";
		ofile << "done\n";

	}
	return 0;


	

	{
		FFH_verbosity = 2;
		int d = 3;
		std::string method("strong"), filename;
		int pixels = 201, num_threads = 4, L = 1;
		double epsilon1 = 1.0L, phi = 0.1;
		double epsilon2 = 10.0L;
		ofile << "compute the effective dielectric constant of 3D two-phase media (export data via ReadConfiguration)\n";
		ofile << "planar external electric field incident along the x-axis\n";
		ofile << "Green's function is modified\n";
		ofile << "electric field is in the y-direction\n";
		ofile << "file name = "; ifile >> filename; ofile << filename << "\n";
		ofile << "num pixles = "; ifile >> pixels;	ofile << pixels << "\n";
		ofile << "epsilon2 = "; ifile >> epsilon2;	ofile << epsilon2 << "\n";
		ofile << "phi2 = "; ifile >> phi;		ofile << phi << "\n";
		ofile << "num threads = "; ifile >> num_threads; ofile << num_threads << "\n"; omp_set_num_threads(num_threads);
		double a = pow(phi / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);
		SpherePacking c(GetUnitCubicBox(d, 1), 0.0), c2, c4;
		{
			GeometryVector x(static_cast<int> (d));
			for (int i = 0; i < d; i++)
				x.x[i] = 0.5;
			c.Insert(a, x);

			Configuration cc(c, "a");
			c2 = SpherePacking(MultiplicateStructure(cc, L), a);
		}
		ConstantTensor E0(d, 1);
		E0[1] = 1.0L;

		TwoPhaseDispersion packing(c2, pixels);
		packing.ComputeIndicator();
		StaticDielectricProblem test(packing, epsilon1, epsilon2);
		test.UseOptimalReference(false);
		test.SetMethod(method.compare("strong") == 0);
		test.DefineProperties(epsilon1, epsilon2);
		test.Prepare(E0);


		std::vector<GeometryVector> result;
		std::vector<GeometryVector> error;

		sprintf(tempstring, "Dimension = %d, Packing fraction = %0.3f, Particle number = %d, \n# resolutions = %d^%d", d, phi, packing.packing.NumParticle(), pixels, d);
		std::string system(tempstring);
		sprintf(tempstring, "epsilon1 = %0.3f, epsilon2 = %0.3f, epsilon0 = %0.3f", epsilon1, epsilon2, test.epsilon0);
		std::string phase_properties(tempstring);
		sprintf(tempstring, "method = %s + regularized", method.c_str());
		std::string solver(tempstring);

		std::vector<std::string> label1 = { system, phase_properties, tempstring,
			"\n#Re[epsilon_e]", "Im[epsilon_e]" };
		std::vector<std::string> label2 = { system, phase_properties, tempstring,
			"\nm", "Tol_field", "Tol_div" };
		system = std::string("Z3_sphere-packing_");

		test.Run(1e-5, 1e-9);


		std::complex<double> epsilon_eff;
		epsilon_eff = test.ComputeEffectiveEpsilon();

		if (test.m <= num_limit && test.IsCovergeningSolution() && test.Tol_Field < 1e-4) {
			result.emplace_back(epsilon_eff.real(), epsilon_eff.imag());
		}
		error.emplace_back(test.m, test.Tol_Field, test.Tol_div);

		sprintf(tempstring, "%s_%0.2f_epsilon2__%0.2f__%s_%d_corr", filename.c_str(), phi, epsilon2, method.c_str(), pixels);
		WriteFunction(result, tempstring, label1);
		sprintf(tempstring, "%s_%0.2f_epsilon2__%0.2f__%s_%d_corr_error", filename.c_str(), phi, epsilon2, method.c_str(), pixels);
		WriteFunction(error, tempstring, label2);
	}
	return 0;









	{
		/* periodic array of spheres (dynamic dielectric). */
		FFH_verbosity = 2;
		int d = 2;
		std::string method("strong");
		int pixels = 51, num_threads = 4, L = 1;
		double epsilon1 = 1.0L, phi = 0.05;
		double epsilon2 = 10.01L; double k = 1.0;
		ofile << "compute the effective dielectric constant of 3D cubic sphere packing\n";
		ofile << "planar external electric field incident along the x-axis\n";
		ofile << "eletric field is in the y-direction\n";
		ofile << "num pixles = "; ifile >> pixels;	ofile << pixels << "\n";
		ofile << "epsilon2 = "; ifile >> epsilon2;	ofile << epsilon2 << "\n";
		ofile << "phi2 = "; ifile >> phi;		ofile << phi << "\n";
		ofile << "num threads = "; ifile >> num_threads; ofile << num_threads << "\n"; omp_set_num_threads(num_threads);
		double k_min = 0.01, k_max = 10, dk = 0.1;

		ofile << "k_min = "; ifile >> k_min;	ofile << k_min;
		ofile << "k_max = "; ifile >> k_max;	ofile << k_max;
		ofile << "dk = "; ifile >> dk; 		ofile << dk;
		double a = pow(phi / HyperSphere_Volume(d, 1.0), 1.0 / (double)d);
		SpherePacking c(GetUnitCubicBox(d, 1), 0.0), c2, c4;
		{
			GeometryVector x(static_cast<int> (d));
			for (int i = 0; i < d; i++)
				x.x[i] = 0.5;
			c.Insert(a, x);

			Configuration cc(c, "a");
			c2 = SpherePacking(MultiplicateStructure(cc, L), a);
		}


		/* c4 is zero-padding...*/
		TwoPhaseDispersion packing(c2, pixels);
		packing.ComputeIndicator();
		DynamicDielectricProblem test(packing, epsilon1, epsilon2);
		test.UseOptimalReference(false);

		std::vector<GeometryVector> result, error;

		test.DefineProperties(epsilon1, epsilon2, GeometryVector(static_cast<int>(d)));
		sprintf(tempstring, "Dimension = %d, Packing fraction = %0.3f, Particle number = %d, \n# resolutions = %d^%d", d, phi, packing.packing.NumParticle(), pixels, d);
		std::string system(tempstring);
		sprintf(tempstring, "epsilon1 = %0.3f, epsilon2 = %0.3f, epsilon0 = %0.3f", epsilon1, epsilon2, test.epsilon0);
		std::string phase_properties(tempstring);
		sprintf(tempstring, "method = %s", method.c_str());
		std::string solver(tempstring);

		std::vector<std::string> label1 = { system, phase_properties, tempstring,
			"\n#k1", "Re[epsilon_e]", "Im[epsilon_e]" };
		std::vector<std::string> label2 = { system, phase_properties, tempstring,
			"\n#k1", "m", "Tol_field", "Tol_div" };
		system = std::string("Z3_sphere-packing_");

		for (double k = k_min; k < k_max; k += dk) {
			ofile << "k = " << k << "\n";

			std::complex<double> epsilon_eff;
			GeometryVector k1;
			{
				GeometryVector e(static_cast<int>(d));
				e.x[0] = 1.0;
				k1 = k * e;
			}

			ofile << "epsilon2 = " << epsilon2 << "\n";
			ofile << "k1 a = " << k1 << "\n";

			ConstantTensor E0(d, 1);
			E0[1] = 1.0L;
			test.SetMethod(method.compare("strong") == 0);
			test.DefineProperties(epsilon1, epsilon2, k1);
			test.Prepare(E0);
			//test.ComputeG_Fourier_test3();
			test.Run(0.1, 1e-10);

			epsilon_eff = test.ComputeEffectiveEpsilon();
			ofile << "epsilon_eff = " << epsilon_eff << "\n";
			error.emplace_back(k, test.m, test.Tol_Field, test.Tol_div);

			if (test.m <= num_limit && test.IsCovergeningSolution() && test.Tol_Field < 1e-4) {
				result.emplace_back(k, epsilon_eff.real(), epsilon_eff.imag());
			}


			sprintf(tempstring, "%s_%0.2f_epsilon2_corr_%0.2f__%s_%d", system.c_str(), phi, epsilon2, method.c_str(), pixels);
			WriteFunction(result, tempstring, label1);
			sprintf(tempstring, "%s_%0.2f_epsilon2_corr_%0.2f_error_%s_%d", system.c_str(), phi, epsilon2, method.c_str(), pixels);
			WriteFunction(error, tempstring, label2);

		}
	}
	return 0;

	

	{
		/* Load configurations ... */


		FFH_verbosity = 2;
		std::string method("weak");

		int d = 3, pixels, num_threads = 1;
		double epsilon1 = 1.0L, phi = 0.4;// 0.4*0.4*pi;
		double epsilon2 = 1.3L; double k = 1.0;
		double target_tol_field = 1e-6, target_tol_div = 1e-10;
		double k_min, k_max, dk;
		ofile << "Verbosity = "; ifile >> FFH_verbosity;
		ofile << "target tolerance in fields = "; ifile >> target_tol_field;
		ofile << "target tolerance in divergence = "; ifile >> target_tol_div;
		ofile << "method = ";   ifile >> method;
		ofile << "number of threads = ";        ifile >> num_threads;
		ofile << "k_min = "; ifile >> k_min;
		ofile << "k_max = "; ifile >> k_max;
		ofile << "dk = "; ifile >> dk;
		for (;;) {
			std::string load_prefix;
			ofile << "input the prefix of the sphere packing = "; ifile >> load_prefix;
			ofile << "input the number of pixels in each side = "; ifile >> pixels;
			if (load_prefix.compare("EXIT") == 0)
				return 0;

			SpherePacking c;
			ReadConfiguration(c, load_prefix);

			d = c.GetDimension();	phi = c.PackingFraction();
			int N = c.NumParticle(), num_chr;
			num_chr = sprintf(tempstring, "Dimension = %d, Packing fraction = %0.3f, Particle number = %d, \n# resolutions = %d^%d", d, phi, N, pixels, d);
			std::string system(tempstring, num_chr);
			ofile << system << "\n";

			num_chr = sprintf(tempstring, "epsilon1 = %0.3f, epsilon2 = %0.3f", epsilon1, epsilon2);
			std::string phase_properties(tempstring, num_chr);
			num_chr = sprintf(tempstring, "method = %s", method.c_str());
			std::string solver(tempstring, num_chr);

			TwoPhaseDispersion packing(c, pixels);	//It gives a discretized version of a sphere packing.
			packing.ComputeIndicator();
			DynamicDielectricProblem test(packing, epsilon1, epsilon2);

			double a = c.GetCharacteristics(0);
			std::vector<GeometryVector> dispersion_relation;

			for (k = k_min; k < k_max; k += dk) {
			std::complex<double> epsilon_eff;
			GeometryVector k1;
			{ GeometryVector e(static_cast<int>(d)); e.x[0] = 1.0L;
			k1 = k * e; }
			ofile << "k1 a = " << k1 << "\n";

			ConstantTensor E0(d, 1);
			E0[1] = 1.0L;
			test.DefineProperties(epsilon1, epsilon2, k1);
			//test.DefineProperties(epsilon1, epsilon2);
			test.SetMethod(method.compare("strong") == 0);
			test.Prepare(E0);
			test.ComputeG_Fourier_test();
			test.Print_GreenFn_2D("truncated_GreenFn");
			//test.ComputeG_Fourier_test4(1e-3/pow(c.NumParticle(), 1.0/(double)d));
			test.Run(1e-5, 1e-8);
			if (test.m <= num_limit && test.IsCovergeningSolution()) {
				epsilon_eff = test.ComputeEffectiveEpsilon();
				ofile << "epsilon_eff = " << epsilon_eff << "\n";
				dispersion_relation.emplace_back(k, epsilon_eff.real(), epsilon_eff.imag(), std::max(test.Tol_div, test.Tol_Field));
			}

			std::vector<std::string> label1 = { system, phase_properties, tempstring,
				"\n#k1", "Re[epsilon_e]", "Im[epsilon_e]", "Error" };
			sprintf(tempstring, "%s_epsilon_eff_%s", load_prefix.c_str(), method.c_str());
			WriteFunction(dispersion_relation, tempstring, label1);
			}
		}

	}
	return 0;
	

	{
		/* random checkerboard */

		double G1 = 1.0L, G2 = 1.3L;
		double K1 = 10000.0L, K2 = 10000.0L;
		double phi = 0.05;	//packing fraction
							//ifile >> K2;
		std::string prefix("8x8R_checkerboard3");
		std::string method("strong");
		bool use_sc = (method.compare("strong") == 0);


		int num_domains = 8, num_pixels = 400;
		int num_configs = 200;
		std::vector<GeometryVector> summary, error_hydro, error_shear;
		FFH_verbosity = 0;

		std::vector<std::complex<double>> moduli;
		SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
		c.Insert(pow(phi / HyperSphere_Volume(2, 1.0), 0.5), GeometryVector(0.5, 0.5));
		TwoPhaseDispersion packing(c, 201);	//It gives a discretized version of a sphere packing.
											//packing.ComputeIndicator();
											//PeriodicCheckerBoard(packing, num_domains, num_pixels);


		for (double value = 0.0; value <= 10.0; value += 0.5) {
			G2 = (value < 1e-5) ? 0.1 : value;
			ofile << "K2/K1 = " << K2 << ", G2/G1 = " << G2 << "\n";
			std::complex<double> Ke, Ge;
			std::complex<double> Ke2, Ge2;

			for (int i = 0; i < num_configs; i++) {
				RandomCheckerBoard(packing, num_domains, num_pixels, rng);
				Homo_ElasticProblem_static test(packing, K1, G1, K2, G2);
				test.SetStrongContrast(use_sc);

				// define the problem.
				test.RedefineProblem(packing, K1, G1, K2, G2);
				test.SetStrongContrast(use_sc);
				/* compression */
				test.Prepare(HydrostaticStrain(2));
				test.Run(5e-7, 1e-10);
				test.ComputeEffectiveProperties(moduli);
				ofile << "Under hydrostatic compression\n";
				ofile << "Ke = " << moduli[0] << "\n";
				ofile << "Ge = " << moduli[1] << "\n";
				Ke += moduli[0];
				Ke2 += moduli[0] * moduli[0];

				ofile << "m = " << test.m << ", error = " << test.Tol_Field << ", div = " << test.Tol_div << "\n";

				/*if (((int)round(value*20.0)) % 10 == 0) {
					sprintf(name, "%s_%0.2f__%.2f_%0.2f_hydro_%s", prefix.c_str(), phi, K2 / K1, G2 / G1, method.c_str());
					test.Print(name);
				}*/

				/* shear moduli */
				test.Prepare(DeviatoricStrain(2));
				test.Run();
				test.ComputeEffectiveProperties(moduli);
				ofile << "Under deviatoric strain\n";
				ofile << "Ke = " << moduli[0] << "\n";
				ofile << "Ge = " << moduli[1] << "\n";
				Ge += moduli[1];
				Ge2 += moduli[1] * moduli[1];

				//error_shear.emplace_back(K2, G2, (double)test.m + test.Tol_Field, test.Tol_div);

				if (((int)round(value*20.0)) % 10 == 0) {
					sprintf(name, "%s_%0.2f__%.2f_%0.2f_shear_%s", prefix.c_str(), phi, K2 / K1, G2 / G1, method.c_str());
					test.Print(name);
				}

			}
			Ke /= (double)num_configs;
			Ge /= (double)num_configs;
			Ke2 = abs(Ke2) / (double)num_configs - Ke.real()*Ke.real();
			Ge2 = abs(Ge2) / (double)num_configs - Ge.real()*Ge.real();
			Ke2 = sqrt(Ke2 / ((double)num_configs - 1.0));
			Ge2 = sqrt(Ge2 / ((double)num_configs - 1.0));

			ofile << "====================\n";
			ofile << " Ke = " << Ke << " +- " << Ke2 << "\n";
			ofile << " Ge = " << Ge << " +- " << Ge2 << "\n";
			ofile << "====================\n";
			summary.emplace_back(K2, G2, Ge.real(), Ge2.real());
			//error_hydro.emplace_back(K2, G2, (double)test.m + test.Tol_Field, test.Tol_div);

		}
		//std::vector<std::string> label({ "K1 = " + std::to_string(K1), "G1 = " + std::to_string(G1) + "\n",
		//	"K2", "G2", "Effective bulk modulus", "Effective shear modulus" });
		std::vector<std::string> label({ "K1 = " + std::to_string(K1), "G1 = " + std::to_string(G1) + "\n",
			"K2", "G2", "Effective shear modulus", "Standard error" });
		/*std::vector<std::string> label2({ "K1 = " + std::to_string(K1), "G1 = " + std::to_string(G1) + "\n",
			"K2", "G2", "iterations + error in field", "error in divergence of stress" });*/
		ofile << prefix << "\t" << name;
		sprintf(name, "%s_ShearModuli_%s", prefix.c_str(), method.c_str());
		WriteFunction(summary, name, label);
		/*sprintf(name, "%s_%0.2f_errors_hydro_%s", prefix.c_str(), phi, method.c_str());
		WriteFunction(error_hydro, name, label2);
		sprintf(name, "%s_%0.2f_errors_shear_%s", prefix.c_str(), phi, method.c_str());
		WriteFunction(error_shear, name, label2);*/

	}
	return 0;

	{

		double G1 = 1.0L, G2 = 1.3L;
		double K1 = 1.0L, K2 = 1.3L;
		double phi = 0.005;	//packing fraction
							//ifile >> K2;
		std::string prefix("square_array1");
		std::vector<GeometryVector> summary_hydrostatic, error_hydro;
		std::vector<GeometryVector> summary_shear, error_shear;
		FFH_verbosity = 0;
		bool use_sc = true;

		std::vector<std::complex<double>> moduli;
		SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
		c.Insert(pow(phi / HyperSphere_Volume(2, 1.0), 0.5), GeometryVector(0.5, 0.5));
		TwoPhaseDispersion packing(c, 401);	//It gives a discretized version of a sphere packing.
		packing.ComputeIndicator();
		Homo_ElasticProblem_static test(packing, K1, G1, K2, G2);
		test.SetStrongContrast(use_sc);

		for (double value = 0.0; value <= 10.0; value += 0.2) {
			G2 = (value < 1e-5) ? 0.01 : value;
			ofile << "K2/K1 = " << K2 << ", G2/G1 = " << G2 << "\n";

			// define the problem.
			test.RedefineProblem(packing, K1, G1, K2, G2);
			test.SetStrongContrast(use_sc);
			/* compression */
			test.Prepare(HydrostaticStrain(2));
			test.Run(5e-7, 1e-10);
			test.ComputeEffectiveProperties(moduli);
			ofile << "Under hydrostatic compression\n";
			ofile << "Ke = " << moduli[0] << "\n";
			ofile << "Ge = " << moduli[1] << "\n";
			summary_hydrostatic.emplace_back(K2, G2, moduli[0].real(), moduli[1].real());
			error_hydro.emplace_back(K2, G2, (double)test.m + test.Tol_Field, test.Tol_div);
			ofile << "m = " << test.m << ", error = " << test.Tol_Field << ", div = " << test.Tol_div << "\n";

			if (((int)round(value*20.0)) % 10 == 0) {
				sprintf(name, "%s_%0.2f__%.2f_%0.2f_hydro", prefix.c_str(), phi, K2 / K1, G2 / G1);
				test.Print(name);
			}

			/* shear moduli */
			test.Prepare(DeviatoricStrain(2));
			test.Run();
			test.ComputeEffectiveProperties(moduli);
			ofile << "Under deviatoric strain\n";
			ofile << "Ke = " << moduli[0] << "\n";
			ofile << "Ge = " << moduli[1] << "\n";
			summary_shear.emplace_back(K2, G2, moduli[0].real(), moduli[1].real());
			error_shear.emplace_back(K2, G2, (double)test.m + test.Tol_Field, test.Tol_div);

			if (((int)round(value*20.0)) % 10 == 0) {
				sprintf(name, "%s_%0.2f__%.2f_%0.2f_shear", prefix.c_str(), phi, K2 / K1, G2 / G1);
				test.Print(name);
			}
		}
		std::vector<std::string> label({ "K1 = " + std::to_string(K1), "G1 = " + std::to_string(G1) + "\n",
			"K2", "G2", "Effective bulk modulus", "Effective shear modulus" });
		std::vector<std::string> label2({ "K1 = " + std::to_string(K1), "G1 = " + std::to_string(G1) + "\n",
			"K2", "G2", "iterations + error in field", "error in divergence of stress" });
		sprintf(name, "%s_%0.2f_Moduli_hydro_strong", prefix.c_str(), phi);
		WriteFunction(summary_hydrostatic, name, label);
		sprintf(name, "%s_%0.2f_errors_hydro_strong", prefix.c_str(), phi);
		WriteFunction(error_hydro, name, label2);

		sprintf(name, "%s_%0.2f_Moduli_shear_strong", prefix.c_str(), phi);
		WriteFunction(summary_shear, name, label);
		sprintf(name, "%s_%0.2f_errors_shear_strong", prefix.c_str(), phi);
		WriteFunction(error_shear, name, label2);

	}
	return 0;


	




	{
		/* usage of the class TensorField*/
		int dimension = 3;
		int pixel = 1, rank1 = 3, rank2 = 2;
		TensorField A(dimension, rank1, std::vector<size_t>(dimension, pixel));
		ofile << "A: " << A.field[0] << "\n";

		TensorField B(dimension, rank2, std::vector<size_t>(dimension, pixel));
		ofile << "B: " << B.field[0] << "\n";

		TensorField C(A);
		ofile << "C: " << C.field[0] << "\n";

		C = A;
		ofile << "C: " << C.field[0] << "\n";

		C = TensorField(A);
		ofile << "C: " << C.field[0] << "\n";

	}
	return 0;

	{
		int dimension = 3;
		int pixel = 1, rank1 = 3, rank2 = 2;
		TensorField A(dimension, rank1, std::vector<size_t>(dimension, pixel));
		TensorField B(dimension, rank2, std::vector<size_t>(dimension, pixel));

		for (size_t i = 0; i < A.num_total_pixels; i++) {
			for (size_t j = 0; j < A.total_comp; j++) {
				A.field[j][i][0] = rng.RandomDouble();
				A.field[j][i][1] = rng.RandomDouble();
			}
		}

		for (size_t i = 0; i < B.num_total_pixels; i++) {
			for (size_t j = 0; j < B.total_comp; j++) {
				B.field[j][i][0] = rng.RandomDouble();
				B.field[j][i][1] = rng.RandomDouble();
			}
		}
		//TensorField C_(dimension, rank1 + rank2 - 2, std::vector<size_t>(dimension, pixel));
		TensorField C; InnerProduct(C, A, B, 1);


		std::vector<size_t> A_i(rank1, 0), B_i(rank2, 0), C_i(rank1 + rank2 - 2, 0);

		/*		for (size_t x = 0; x < C_.num_total_pixels; x++) {
		for (size_t i= 0; i < dimension; i++) {
		//		for (size_t j= 0; j < 3; j++) {
		//for (int m = 0; m < 3; m++) {
		for (size_t n = 0; n < dimension; n++) {
		ofile << "\n\n";
		for (size_t k = 0; k < dimension; k++) {
		C_i = std::vector<size_t>({ i,n });
		A_i = std::vector<size_t>({ i,k });
		B_i = std::vector<size_t>({ k,n });
		ofile << i + dimension * k << "," << k + dimension * n << "\n";
		//								C_i = std::vector<size_t>({ i,j,n });
		//								A_i = std::vector<size_t>({ i,j,k });
		//								B_i = std::vector<size_t>({k,n});


		//								ofile << B_i[0] << ", " << B_i[1] << ": "
		//									<< B.at(B_i)[x][0] << "\t" << B.at(B_i)[x][1] << "\n";

		//								ofile << C_i[0] <<", "<<C_i[1] <<", "<<C_i[2]<<": "
		//									<< C.at(C_i)[x][0] << "\t" << C.at(C_i)[x][1] << "\n";
		C_.at(C_i)[x][0] += (A.at(A_i)[x][0] * B.at(B_i)[x][0] - A.at(A_i)[x][1] * B.at(B_i)[x][1]);
		C_.at(C_i)[x][1] += (A.at(A_i)[x][0] * B.at(B_i)[x][1] + A.at(A_i)[x][1] * B.at(B_i)[x][0]);

		}

		}
		//}
		}
		//	}
		}*/
		/*for (size_t x = 0; x < C_.num_total_pixels; x++) {
		for (size_t i = 0; i < dimension; i++) {
		for (size_t j= 0; j < dimension; j++) {
		//for (int m = 0; m < 3; m++) {
		for (size_t n = 0; n < dimension; n++) {
		ofile << "\n\n";
		for (size_t k = 0; k < dimension; k++) {
		C_i = std::vector<size_t>({ i,j,n });
		A_i = std::vector<size_t>({ i,j,k });
		B_i = std::vector<size_t>({k,n});


		//								ofile << B_i[0] << ", " << B_i[1] << ": "
		//									<< B.at(B_i)[x][0] << "\t" << B.at(B_i)[x][1] << "\n";

		//								ofile << C_i[0] <<", "<<C_i[1] <<", "<<C_i[2]<<": "
		//									<< C.at(C_i)[x][0] << "\t" << C.at(C_i)[x][1] << "\n";
		C_.at(C_i)[x][0] += (A.at(A_i)[x][0] * B.at(B_i)[x][0] - A.at(A_i)[x][1] * B.at(B_i)[x][1]);
		C_.at(C_i)[x][1] += (A.at(A_i)[x][0] * B.at(B_i)[x][1] + A.at(A_i)[x][1] * B.at(B_i)[x][0]);

		}

		}
		//}
		}
		}
		}*/

		ofile << "check !!\n";
		for (int x = 0; x < C.num_total_pixels; x++) {
			ofile << "x = " << x << "\n";
			ofile << "Real[diff.]\tImag[diff.]\n";
			for (int i = 0; i < C.total_comp; i++) {
				//ofile << C.field[i][x][0] - C_.field[i][x][0] << "\t" << C.field[i][x][1] - C_.field[i][x][1] << "\n";
			}
		}

		ofile << "A = " << A.field[0] << "\n";
		ofile << "B = " << B.field[0] << "\n";
		ofile << "C = " << C.field[0] << "\n";
		//ofile << "C_ = " << C_.field[0] << "\n";

		for (int i = 0; i < 10; i++) {
			TensorField c(C);
		}
		ofile << "end?\n";
		//C.~TensorField();
	}
	return 0;
	{

		TensorField ff(3, 3, std::vector<size_t>(3, 1));
		for (int i = 0; i < 27; i++) {
			ff.field[i][0][0] = 1.0;
			ff.field[i][0][1] = i;
		}
		TensorField g = ff;

		std::vector<size_t> dd;
		for (size_t i = 0; i<3; i++)
			for (size_t j = 0; j<3; j++)
				for (size_t k = 0; k < 3; k++) {
					dd = std::vector<size_t>({ i,j,k });
					ofile << "ff_(" << i << ", " << j << ", " << k << "): " << ff.comp(dd)[0][0] << " +i" << (ff.comp(dd))[0][1] << "\n";
					ofile << "g_(" << i << ", " << j << ", " << k << "): " << g.comp(dd)[0][0] << " +i" << (g.comp(dd))[0][1] << "\n";

				}

	}
	return 0;
	{
		std::vector<int *> c;	c.reserve(10);
		int idx = 0;
		for (int i = 0; i < 10; i++) {
			int * temp; temp = (int *)malloc(5 * sizeof(int));
			for (int j = 0; j < 5; j++) {
				temp[j] = (idx++);
			}
			c.push_back(temp);
		}

		ofile << "show !!\n";
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 5; j++) {
				ofile << c[i][j] << "\t";
			}
			ofile << "\n";
		}

		ofile << "new values !!\n";
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 5; j++) {
				c[i][j] = 123;
			}
		}

		ofile << "show !!\n";
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 5; j++) {
				ofile << c[i][j] << "\t";
			}
			ofile << "\n";
		}

	}

	return 0;
	{
		std::vector<GeometryVector> b;
		b.emplace_back(1, 0);
		b.emplace_back(0, 1);

		SpherePacking c(2, &b[0], 1);
		for (int i = 0; i < 10; i++) {
			GeometryVector x(rng.RandomDouble(), rng.RandomDouble());
			c.Insert(0.1, x);
		}
		int num_pix = 10;
		WriteConfiguration(c, "test_config");
		TwoPhaseDispersion test(c, num_pix);
		test.ComputeIndicator();

		{
			std::vector<GeometryVector> x;
			for (int i = 0; i < test.prt_indicator.size(); i++) {
				std::vector<size_t> comp; test.Index2Components(test.prt_indicator[i], comp);
				x.emplace_back((comp[0] + 0.5) / (double)num_pix, (comp[1] + 0.5) / (double)num_pix);
			}
			WriteFunction(x, "discrete1");
		}

		num_pix = 50;
		{
			test.SetResolutions(num_pix);
			test.ComputeIndicator();
			std::vector<GeometryVector> x;
			for (int i = 0; i < test.prt_indicator.size(); i++) {
				std::vector<size_t> comp; test.Index2Components(test.prt_indicator[i], comp);
				x.emplace_back((comp[0] + 0.5) / (double)num_pix, (comp[1] + 0.5) / (double)num_pix);
			}
			WriteFunction(x, "discrete2");
		}
	}
	return 0;
}


void RandomCheckerBoard(TwoPhaseDispersion & out, int num_domains, int num_pixels, RandomGenerator & rng)
{
	SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
	out = TwoPhaseDispersion(c, num_pixels);
	//generate a random checkboard;
	if ((num_domains % 2) == 0)
	{
		std::vector<std::vector<size_t>> domains(num_domains*num_domains);
		int index = 0;
		for (int i = 0; i < num_domains; i++) {
			for (int j = 0; j < num_domains; j++) {
				domains[index].push_back(i);
				domains[index].push_back(j);
				index++;
			}
		}

		for (int i = 0; i < num_domains*num_domains / 2; i++) {
			domains.erase(domains.begin() + int(floor(domains.size()*rng.RandomDouble())));
		}


		std::vector<GeometryVector> basis;
		basis.emplace_back(out.packing.GetBasisVector(0));
		basis.emplace_back(out.packing.GetBasisVector(1));
		basis[0] = 1.0 / (double)num_pixels * basis[0];
		basis[1] = 1.0 / (double)num_pixels * basis[1];

		Configuration cell(2, &basis[0], 1);
		cell.Insert("a", GeometryVector(0.5, 0.5));

		Configuration grid = MultiplicateStructure(cell, std::vector<size_t>({ (size_t)num_pixels, (size_t)num_pixels }));
		out.prt_indicator.resize(grid.NumParticle(), 0);

		index = 0;
		for (int i = 0; i < grid.NumParticle(); i++) {
			GeometryVector pos = grid.GetRelativeCoordinates(i);
			int x = (int)floor(num_domains*pos.x[0]);
			int y = (int)floor(num_domains*pos.x[1]);
			bool inside = false;
			for (int j = 0; j < domains.size(); j++) {
				inside = inside || (domains[j][0] == x && domains[j][1] == y);
			}
			if (inside) {
				std::vector<size_t> comp({ (size_t)floor(pos.x[0] * num_pixels), (size_t)floor(pos.x[1] * num_pixels) });
				out.prt_indicator[index] = out.Components2Index(comp);
				index++;
			}
		}
		out.prt_indicator.resize(index);
		out.update_indicator = true;
	}

	else
		std::cerr << "num_domains must be an even number\n";
}

void PeriodicCheckerBoard(TwoPhaseDispersion & out, int num_domains, int num_pixels)
{
	SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
	out = TwoPhaseDispersion(c, num_pixels);
	if((num_domains % 2)==0)
	{
		std::vector<std::vector<size_t>> domains(num_domains*num_domains);
		int index = 0;
		for (int i = 0; i < num_domains; i++) {
			for (int j = 0; j < num_domains; j++) {
				domains[index].push_back(i);
				domains[index].push_back(j);
				index++;
			}
		}

		for (int i = 0; i < domains.size(); i++) {
			if ((domains[i][0] + domains[i][1]) % 2 == 0) {
				domains.erase(domains.begin() + i);
				i--;
			}
		}


		std::vector<GeometryVector> basis;
		basis.emplace_back(out.packing.GetBasisVector(0));
		basis.emplace_back(out.packing.GetBasisVector(1));
		basis[0] = 1.0 / (double)num_pixels * basis[0];
		basis[1] = 1.0 / (double)num_pixels * basis[1];

		Configuration cell(2, &basis[0], 1);
		cell.Insert("a", GeometryVector(0.5, 0.5));

		Configuration grid = MultiplicateStructure(cell, std::vector<size_t>({ (size_t)num_pixels, (size_t)num_pixels }));
		out.prt_indicator.resize(grid.NumParticle(), 0);

		index = 0;
		for (int i = 0; i < grid.NumParticle(); i++) {
			GeometryVector pos = grid.GetRelativeCoordinates(i);
			int x = (int)floor(num_domains*pos.x[0]);
			int y = (int)floor(num_domains*pos.x[1]);
			bool inside = false;
			for (int j = 0; j < domains.size(); j++) {
				inside = inside || (domains[j][0] == x && domains[j][1] == y);
			}
			if (inside) {
				std::vector<size_t> comp({ (size_t)floor(pos.x[0] * num_pixels), (size_t)floor(pos.x[1] * num_pixels) });
				out.prt_indicator[index] = out.Components2Index(comp);
				index++;
			}
		}
		out.prt_indicator.resize(index);
		out.update_indicator = true;
	}

	else {
		std::cout << "num_domains should be an even number\n";
	}
}


void PeriodicStrips(TwoPhaseDispersion & out, int num_domains, int num_pixels)
{
	SpherePacking c(GetUnitCubicBox(2, 1), 0.0);
	out = TwoPhaseDispersion(c, num_pixels);
	if ((num_domains % 2) == 0)
	{
		std::vector<std::vector<size_t>> domains(num_domains);
		int index = 0;
		for (int i = 0; i < num_domains; i+=2) {
				domains[index].push_back(i);
				index++;
		}

		std::vector<GeometryVector> basis;
		basis.emplace_back(out.packing.GetBasisVector(0));
		basis.emplace_back(out.packing.GetBasisVector(1));
		basis[0] = 1.0 / (double)num_pixels * basis[0];
		basis[1] = 1.0 / (double)num_pixels * basis[1];

		Configuration cell(2, &basis[0], 1);
		cell.Insert("a", GeometryVector(0.5, 0.5));

		Configuration grid = MultiplicateStructure(cell, std::vector<size_t>({ (size_t)num_pixels, (size_t)num_pixels }));
		out.prt_indicator.resize(grid.NumParticle(), 0);

		index = 0;
		for (int i = 0; i < grid.NumParticle(); i++) {
			GeometryVector pos = grid.GetRelativeCoordinates(i);
			int x = (int)floor(num_domains*pos.x[0]);
//			int y = (int)floor(num_domains*pos.x[1]);
			bool inside = ((x%2)==0);
//			for (int j = 0; j < domains.size(); j++) {
//				inside = inside || (domains[j][0] == x && domains[j][1] == y);
//			}
			if (inside) {
				std::vector<size_t> comp({ (size_t)floor(pos.x[0] * num_pixels), (size_t)floor(pos.x[1] * num_pixels) });
				out.prt_indicator[index] = out.Components2Index(comp);
				index++;
			}
		}
		out.prt_indicator.resize(index);
		out.update_indicator = true;
	}

	else {
		std::cout << "num_domains should be an even number\n";
	}
}