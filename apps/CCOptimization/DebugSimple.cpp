#include <fstream>
#include <sstream>
#include <iostream>
#include <complex>
#include <list>
#include <set>
#include <boost/math/special_functions.hpp>
#include <iomanip> // setprecision

#include "PeriodicCellList.h"
#include "Plots.h"
//#include "MC_System.h"
#include "MD_System.h"
#include "StructureOptimization.h"
//#include "PhononFrequency.h"
//#include "Voronoi.h"
#include "NumberVariance.h"
#include "PairCorrelation.h"
#include "StructureFactor.h"

#include "Potential.h"
#include "CollectiveCoordinatePotential.h"
//#include "LJPotential.h"

#include "RandomSequentialAddition.h"
//#include "TorquatoJiao.h"
#include "RandomOrganization.h"
//#include "Optimize.h"
//#include "Solvers.h"
//#include "ThreeBody.h"
//#include "Degeneracy.h"

#include "Discrepancy.h"
#include "etc1.h"

#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <limits>

#include "Fitting.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
//#include <windows.h>
#include <sstream>
#include <string>

#ifdef USE_CGAL
#include "Voronoi.h"
#endif
//#include "Discrepancy.h"


#include "CCInversionProblem.h"
//These functions are defined in PairStatisticsCLI.cpp, but other modules could also find them useful
//void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile);
//void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix);
void ReadFunction(std::vector<GeometryVector> & result, std::istream & ifile, size_t NumColumns=4);
void ReadFunction(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t NumColumns=4);
//void ReadFunction__(std::vector<GeometryVector> & result, std::istream & ifile, size_t DataColumns, size_t NumColumns);
//void ReadFunction__(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t DataColumns, size_t NumColumns);


//namespace
//{
//	Configuration getConfigsFunction(size_t i)
//	{
//		//read a configuration from a grace file
//			std::vector<std::vector<GeometryVector> > loc;
//		std::string prefix = "pts_config_VdC_N_1024";
//		ReadGraceData(loc, prefix);
//		Configuration c = GetUnitCubicBox(2);
//		for (int j = 0; j < i; j++) {
//			for (auto iter = loc[j].begin(); iter != loc[j].end(); iter++)
//				c.Insert("A", *iter);
//		}
//
//		return c;
//	}
//};




//Generate N_c configurations with the same collective coordinates at the first N_kPoints k points. Compare it with the target
//if inequivalent up to permutation, return false
//bool FrankInvCCUniqueTest(const Configuration & target, size_t Nc, size_t N_kPoints, double Tolerance = 1e-4, size_t NKPointsSkip = 0, std::fstream & outFile = logfile)
//{
//	DimensionType d = target.GetDimension();
//	size_t NumParticle = target.NumParticle();
//	RandomGenerator gen(1234567);
//
//	//generate k point list
//	double KLimit = 2 * pi*std::pow(1.2*N_kPoints / target.PeriodicVolume() / HyperSphere_Volume(d, 1.0), 1.0 / d);
//	std::vector<GeometryVector> kPoints;
//	while (kPoints.size() < N_kPoints + NKPointsSkip)
//	{
//		kPoints = GetKs(target, KLimit, KLimit);
//		KLimit *= 1.5;
//	}
//	auto compareFunc = [](const GeometryVector & left, const GeometryVector & right) -> bool
//	{ return left.Modulus2() < right.Modulus2(); };
//	std::sort(kPoints.begin(), kPoints.end(), compareFunc);
//	FrankInvertCCPotential pot(d);	pot.ParallelNumThread = 6;
//	for (size_t i = 0; i < N_kPoints; i++)
//	{
//		GeometryVector k = kPoints[i + NKPointsSkip];
//		std::complex<double> rho;
//		for (size_t i = 0; i<NumParticle; i++)
//			rho += std::exp(std::complex<double>(0, k.Dot(target.GetCartesianCoordinates(i))));
//		pot.constraints.push_back(FrankInvertCCPotential::KPointConstraint(k, 1.0, rho.real(), rho.imag()));
//	}
//
//	double maxDr = 0.0;
//	//generate configurations and test
//	for (size_t i = 0; i < Nc; i++)
//	{
//		double E = 100.0;
//		Configuration temp = target;
//		while (E > 1e-20)
//		{
//			temp.RemoveParticles();
//			for (size_t j = 0; j < NumParticle; j++)
//				temp.Insert("A", gen);
//			RelaxStructure_NLOPT(temp, pot, 0.0, 0, 0.0);
//			RelaxStructure_MINOP(temp, pot, 0.0, 0, 0.0);
//			pot.SetConfiguration(temp);
//			E = pot.Energy();
//		}
//
//		for (size_t j = 0; j < NumParticle; j++)
//		{
//			GeometryVector r = target.GetRelativeCoordinates(j);
//			double dr = temp.NearestParticleDistance(r);
//			maxDr = std::max(dr, maxDr);
//
//			//debug temp
//			//std::cout << "E=" << E << std::endl;
//			//if (maxDr > Tolerance)
//			//            std::cout << "nonuniqueness found\n";
//		}
//	}
//	outFile << "max dr=" << maxDr;
//	if (maxDr > Tolerance)
//	{
//		return false;
//	}
//	return true;
//}

//int Debug()
//{
//
//	{
//		const size_t Nc = 100;
//		size_t NKSkip;
//		std::cout << "NKSkip=";
//		std::cin >> NKSkip;
//		size_t dim;
//		std::cout << "dim=";
//		std::cin >> dim;
//		if (NKSkip>0)
//			std::cout << "Skip the central " << NKSkip << " k points.\n";
//
//
//		std::stringstream ss;
//		ss << "CC_Inv_Log_d_" << dim << "_NKSkip_" << NKSkip;
//		std::fstream outfile(ss.str(), std::fstream::out);
//
//		for (size_t N = 4; N < 30; N++)
//		{
//			Configuration c = GeneratePoisson(dim, N, RandomGenerator(1234));
//
//			outfile << "testing Nk for N=" << N;
//			outfile << std::endl;
//			for (size_t Nk = 1;; Nk++)
//			{
//				outfile << "Nk=" << Nk << ", testing...";
//				outfile.flush();
//
//				bool result = FrankInvCCUniqueTest(c, Nc, Nk, 1e-4, NKSkip, outfile);
//				if (result)
//				{
//					outfile << ", Configuration Unique.\n";
//					std::cout << N << " " << Nk << std::endl;
//
//					break;
//				}
//				else
//					outfile << ", Configuration not Unique. Increasing Nk.\n";
//			}
//		}
//
//		return 0;
//	}
//}


void ReadCoordinationNumber(std::ifstream &ifile, std::vector<double> & r2, std::vector<double> &Z);
void WriteToContourf(const std::vector<GeometryVector> & vec_k, const std::vector<GeometryVector> &s_k, const std::string &prefix, int length);
//int PoissonRNG(const boost::random::mt19937 &g, int mean);


void SetTrialMovement(SimpleHarmonic &pot, Configuration & conf, double &dx, double T, double lowAccept, double upperAccept, int randomSeed, int numSteps);

void SetTrialMovement(std::vector<SimpleHarmonic> &pot, std::vector<Configuration> & confList, std::vector<RandomGenerator> & g, double &dx, double T, double lowAccept, double upperAccept, int numSteps);

void SetTrialMovement(std::vector<SimpleHarmonic> &pot, std::vector<RandomGenerator> & g, double &dx, double T, double lowAccept, double upperAccept, int numSteps);

/** \brief Try displacement du for particle i and accept it by the Metropolis rule; keep going if acceptable.
 *  \param pot          Harmonic potential whose energy change is evaluated.
 *  \param conf         Configuration being evolved.
 *  \param g            Random number generator for the acceptance test.
 *  \param temperature  Temperature of the Metropolis acceptance.
 *  \param du           Trial displacement for particle i.
 *  \param i            Index of the particle to move.
 *  \return             1 if the move is accepted, 0 otherwise. */
double MCevolve(SimpleHarmonic &pot, Configuration & conf, RandomGenerator &g, const double temperature, GeometryVector & du, size_t i) {
	double result = 0;
	double delta_E = pot.Delta_energy(i, du);
	double Pr = g.RandomDouble();

	if (delta_E < 0 || exp(-delta_E / temperature) > Pr) {
		pot.updateUn(i, du);
		result = 1.0;
	}

	return result;
}

/** \brief Try displacement du for particle i using the potential's own configuration and accept it by the Metropolis rule. (CLAUDE)
 *  \param pot          Harmonic potential whose energy change is evaluated.
 *  \param g            Random number generator for the acceptance test.
 *  \param temperature  Temperature of the Metropolis acceptance.
 *  \param du           Trial displacement for particle i.
 *  \param i            Index of the particle to move.
 *  \return             1 if the move is accepted, 0 otherwise. */
double MCevolve(SimpleHarmonic &pot, RandomGenerator &g, const double temperature, GeometryVector & du, size_t i) {
	double result = 0;
	double delta_E = pot.Delta_energy(i, du);
	double Pr = g.RandomDouble();

	if (delta_E < 0 || exp(-delta_E / temperature) > Pr) {
		pot.updateUn(i, du);
		result = 1.0;
	}

	return result;
}

/** \brief Make one trial move for every particle in every configuration and return the overall acceptance rate. (CLAUDE)
 *  \param pot          List of harmonic potentials, one per configuration.
 *  \param g            List of random number generators, one per configuration.
 *  \param temperature  Temperature of the Metropolis acceptance.
 *  \param dx           Maximum magnitude of a trial displacement component.
 *  \return             Fraction of accepted moves over all particles and configurations. */
double MCevolve(std::vector<SimpleHarmonic> &pot, std::vector<RandomGenerator> &g, const double temperature, const double dx) {//make trial moves for all particles in all configurations
	if (pot.size() == g.size()) {
		double rate = 0;
		int d = pot[0].Dimension;
		int N = pot[0].confi.NumParticle();
		std::vector<double> temp_rate(pot.size(), 0);
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < pot.size(); i++) {
			temp_rate[i] = 0;
			GeometryVector du(static_cast<int> (d)); //trial displacement
			for (int j = 0; j < N; j++) {//evolve one particle in each configuration
				for (int k = 0; k < d; k++) { du.x[k] = dx*(g[i].RandomDouble()-0.5); } //random displacements
				temp_rate[i] += MCevolve(pot[i], g[i], temperature, du, j); 
			}
		}
		for (int i = 0; i < pot.size(); i++)
			rate += temp_rate[i];

		rate /= (double)(pot.size()*N); //
		return rate;
	}
	else {
		std::cout << "MCevolve:the number of pot and g are different \n";
		return 0;
	}
}



#include "MonteCarlo.h"

//MonteCarlo<SpherePacking, Potential>

#include "TilingProblem.h"

/** \brief Convert a monodisperse packing into a polydisperse one by sizing particles to their Voronoi-cell volumes. (CLAUDE)
 *  \param config       Input packing whose Voronoi cells are measured.
 *  \param resultConfig Output packing with radii set from Voronoi-cell volumes.
 *  \param targetPhi    Target packing fraction of the result.
 *  \param Volumes      Output Voronoi-cell volume of each particle.
 *  \param numThreads   Number of OpenMP threads to use.
 *  \param numVertices  Optional output count of vertices of each Voronoi cell (may be nullptr).
 *  \return             Ratio of the largest-sphere volume to the maximum Voronoi-cell volume. */
double MakeMono2Poly_parallel_test(const SpherePacking & config, SpherePacking & resultConfig, double targetPhi, std::vector<double> &Volumes, int numThreads, std::vector<int> * numVertices) {
	//make resultConfig empty
	resultConfig = SpherePacking(config);
	resultConfig.RemoveParticles();

	if (numVertices != nullptr) {
		(*numVertices).resize(config.NumParticle(), 0);
	}

	DimensionType d = config.GetDimension();
	int numPrts = config.NumParticle();
	double v_1 = HyperSphere_Volume(d, config.GetMaxRadius());
	Volumes.erase(Volumes.begin(), Volumes.end());
	Volumes.resize(numPrts, 0.0);

	double max_cellVol = 0;
	double r = 0;

	double factor = pow(targetPhi / HyperSphere_Volume(d, 1.0), 1.0 / d);
	//------------- measure volumes of Voronoi cells

#if defined (USE_VORO)
	int NMAX;
	switch (d) {
	case 2:
		NMAX = 100000; /*500 */  break;
	case 3:
		NMAX = 10000000; break;
	default:
		NMAX = 1;
		std::cout << "dimension is too high" << std::endl;
		break;
	}
		
	double fac = 0.5;
	int nx, ny, nz = 0;	bool trans;
	int mem = 30;
	std::vector<GeometryVector> basis;
	std::vector<container_leak> con_list;
	

	if (d == 2 || d == 3) {
		/*configuration is small enough*/
		if (numPrts < NMAX) {
			GetConversionParam(config, basis, trans, nx, ny, nz, mem);
			con_list.emplace_back(container_leak(config, basis, trans, nx, ny, nz, mem));
		}
		/*a configuration is too large => split the configuration*/
		else {
			double margin = 6.0;
			double rho = config.NumParticle() / config.PeriodicVolume();
			mem = 10;
			int numSplit = (int)ceil(pow(config.NumParticle(), 1.0 / d) / (pow(NMAX*fac, 1.0 / d) - 2.0*margin));
			std::vector<GeometryVector> basis_;//Basis vectors for the domain of our interst in each split
			std::vector<GeometryVector> basis_subdomain;
			std::vector<GeometryVector> centerDomain; //In Cartesian coordinates.

			for (int i = 0; i < d; i++)
				basis_.push_back((1.0 / (double)numSplit) * config.GetBasisVector(i));
			{/*Center of domains*/
				GeometryVector center;
				for (int i = 0; i < numSplit; i++) {
					for (int j = 0; j < numSplit; j++) {
						if (d == 2)
							center = (0.5 + i)*basis_[0] + (0.5 + j)*basis_[1];
						else if (d == 3) {
							for (int k = 0; k < numSplit; k++)
								center = (0.5 + i)*basis_[0] + (0.5 + j)*basis_[1] + (0.5 + k)*basis_[2];
						}
						centerDomain.push_back(center);
					}
				}
			}
			for (int i = 0; i < d; i++) {//Extend the basis vectors of subdomains.
				basis_subdomain.push_back( (1.0 + 2.0*margin / pow(rho, 1.0 / d) / sqrt(basis_[i].Modulus2()))*basis_[i]);
			}
			
			/*new basis vectors when basis_subdomain is not suitable for container_leak!*/
			std::vector<GeometryVector> basis_subdomain_standard;
			std::vector<GeometryVector> unitbasisVectors;
			for (int i = 0; i < d; i++)
				unitbasisVectors.push_back((1.0/sqrt(basis_[i].Modulus2()))*basis_[i]);

			Configuration subdomain_template(d, &basis_subdomain[0], 1.0);
			GetConversionParam(subdomain_template, (int)rho*subdomain_template.PeriodicVolume(),basis_subdomain_standard, trans, nx, ny, nz, mem);
			for (int i = 0; i < centerDomain.size(); i++) {
				container_leak con(subdomain_template, basis_subdomain_standard, trans, nx, ny, nz, mem);
				con_list.emplace_back(con);
				//con_list.emplace_back(container_leak(subdomain_template, basis_subdomain_standard, trans, nx, ny, nz, mem));
			}
			double r_max = 0.5 * subdomain_template.MaxDiagonalLength();
			PeriodicCellList<Empty> c(config);	c.PrepareIterateThroughNeighbors(r_max);
			
			//half lengths of basis vectors
			std::vector<double> region0;
			std::vector<double> region1;
			for (int i = 0; i < d; i++) {
				region0.push_back(0.5*sqrt(basis_[i].Modulus2()));
				region1.push_back(0.5*sqrt(basis_subdomain[i].Modulus2()));
			}
			/* return 0 = region of interest
			return 1 = marginal region
			return 2 = don't consider
			*/
			auto regionDet = [&unitbasisVectors, &region0, &region1](const GeometryVector & cartesianCoord) ->int {
				int region = 0;
				int region_idx;
				double projLength;
				for (int i = 0; i < unitbasisVectors.size(); i++) {
					projLength = fabs(unitbasisVectors[i].Dot(cartesianCoord));
					if (projLength <= region0[i])
						region_idx = 0;
					else if (projLength < region1[i])
						region_idx = 1;
					else
						region_idx = 2;

					region = std::max(region, region_idx);
				}
				
				return region;
			};
			GeometryVector ver2cen(basis_subdomain[0]);
			for (int i = 1; i < d; i++)
				ver2cen = ver2cen + basis_subdomain[i];
			ver2cen = 0.5*ver2cen;
			/*Fill particle positions!*/
#pragma omp parallel for num_threads(numThreads) schedule(guided)
			for (int i = 0; i < con_list.size(); i++) {
				GeometryVector x0 = c.CartesianCoord2RelativeCoord(centerDomain[i]);
				c.IterateThroughNeighbors(x0, r_max, 
					[i, &con_list, &regionDet, &ver2cen](const GeometryVector & shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, size_t Sourceparticle)->void {
					int region = regionDet(shift);
					if (region <= 1) {
						GeometryVector r = shift + ver2cen;
						if (con_list[i].transform) {
							r = con_list[i].c_.RelativeCoord2CartesianCoord(con_list[i].c0.CartesianCoord2RelativeCoord(r));
						}
						double x = r.x[0], y = r.x[1], z=0.5;
						if (shift.Dimension == 3)
							z = r.x[2];
						int prt_idx = (region == 0) ? Sourceparticle : -1;
						int ijk = con_list[i].put_ijk(prt_idx, x, y, z);
						IDs temp;	temp.blockId = ijk;	temp.order = (con_list[i].co[ijk]) - 1;
						con_list[i].id_leak.push_back(temp); 
					}
				}
				);
			}
		}
	}
	else {
		std::cout << "MakeMono2Poly_parallel::" << "Dimension is too high\n";
	}

	/*Compute volume of Voronoi cells*/
	/*Parameters

	*/
	int count_fail = 0;
#pragma omp parallel for num_threads(numThreads) schedule(guided)
	for (int i = 0; i < con_list.size(); i++) {
		voro::voronoicell vc;
		voro::c_loop_all_periodic routine(con_list[i]);
		int idx=0, count = 0;
		double cellVolume;
		std::vector<double> vertices;
		if (routine.start()) {
			do {
				idx = routine.pid();
				if (idx > -1) {/*Particles are in the inner region*/
					if (con_list[i].compute_cell(vc, routine)) {
						Volumes[idx] = vc.volume();
						if (numVertices != nullptr) {
							vc.vertices(vertices);
							int numV = vertices.size() / 3;
							if (d == 2)
								(*numVertices)[idx] = numV / 2;
							else if (d == 3)
								(*numVertices)[idx] = numV;
						}
					}
					else
						count++;
				}
			} while (routine.inc());
		}
#pragma omp atomic
		count_fail += count;
	}
	std::cout << "# of failure = " << count_fail << std::endl;
#else
#pragma omp parallel num_threads(numThreads)
	{
		SpherePacking c_copy(config);
		std::vector<GeometryVector> vertices;
		double temp = 0;
#pragma omp for
		for (int i = 0; i < numPrts; i++) {
			GetVoronoiCell(c_copy, i, vertices);
			if (numVertices != nullptr)
				numVertices->at(i) = (int)vertices.size();
			temp = PolygonVolume(vertices, d);
#pragma omp critical
			{
				max_cellVol = (max_cellVol > temp) ? max_cellVol : temp;
			}
			Volumes.at(i) = temp;
		}
	}
#endif

	for (int i = 0; i < numPrts; i++) {
		double temp = Volumes[i];
		if (d == 2)
			r = sqrt(temp);
		else if (d == 3)
			r = cbrt(temp);
		else { r = pow(temp, 1.0 / d); }

		r *= factor;
		resultConfig.Insert(r, config.GetRelativeCoordinates(i));
	}


	//for (int i = 0; i < Volumes.size(); i++)
	//	Volumes.at(i) /= (double)config.NumParticle();

	return v_1 / max_cellVol;
}

/*particle packings -> layer thickness in the spatial order.
	particle => phase 1
	matrix => phase 2
	thickness 1		0.0		0.0		phase
	thickness 2		0.0		0.0		phase
	...
*/
/** \brief Convert a 1D sphere packing into the ordered sequence of phase-1 (particle) and phase-2 (matrix) layer thicknesses. (CLAUDE)
 *  \param c  Input 1D sphere packing.
 *  \param x  Output list of layers; each entry stores a thickness and a phase label (1=particle, 2=matrix). */
void SpherePacking2Layers(const SpherePacking & c, std::vector<GeometryVector> & x) {
	//box size
	double xi = 0.0;	double xf = sqrt(c.GetBasisVector(0).Modulus2());

	//PeriodicCellList<double> temp(c);
	x.clear();	x.reserve(2 * c.NumParticle());
	
	std::vector<GeometryVector> pos;	double x0;
	for (int i = 0; i < c.NumParticle(); i++) {
		x0 = c.GetCartesianCoordinates(i).x[0];
		pos.emplace_back(GeometryVector(x0 - c.GetCharacteristics(i), x0 + c.GetCharacteristics(i)));
	}

	//At the top boundary
	int i = 0;
	if (pos[i].x[0] < xi) //the first layer is the particle phase
		x.emplace_back(GeometryVector(pos[i].x[1] - pos[i].x[0], 0.0, 0.0, 1));
	else {	//the first layer is the matrix phase
		x.emplace_back(GeometryVector(pos[i].x[0], 0.0, 0.0, 2));
		x.emplace_back(GeometryVector(pos[i].x[1] - pos[i].x[0], 0.0, 0.0, 1));//particle phase
	}
	i++;

	for (i = 1; i < pos.size()-1; i++) {
		x.emplace_back(GeometryVector(pos[i].x[0] - pos[i - 1].x[1], 0.0, 0.0, 2));//matrix phase
		x.emplace_back(GeometryVector(pos[i].x[1] - pos[i].x[0], 0.0, 0.0, 1));//particle phase
	}

	//At the bottom boundary
	double compensate = (pos[0].x[0] < xi)? 0.0 - pos[0].x[0]: 0.0;
	
	if (pos[i].x[1] < xf) {//the last boundary is the matrix phase 
		x.emplace_back(GeometryVector(pos[i].x[0] - pos[i - 1].x[1], 0.0, 0.0, 2));//matrix phase
		x.emplace_back(GeometryVector(pos[i].x[1] - pos[i].x[0], 0.0, 0.0, 1));//particle phase
		x.emplace_back(GeometryVector(xf - pos[i].x[1] - compensate, 0.0, 0.0, 2)); //matrix phase
	}
	else if(pos[i].x[1] > xf){//the bottom boundary is the particle phase
		x.emplace_back(GeometryVector(pos[i].x[0] - pos[i - 1].x[1], 0.0, 0.0, 2));//matrix phase
		x.emplace_back(GeometryVector(pos[i].x[1] - pos[i].x[0], 0.0, 0.0, 1));//particle phase
		x[0].x[0] = x[0].x[0] - (pos[i].x[1] - xf);
	}

	double total = 0.0;
	for (int i = 0; i < x.size(); i++)
		total += x[i].x[0];

	std::cout << "Total size = " << total << "\t";
}

/** \brief Reorder layers so the two phases strictly alternate (chirp the layered pattern). (CLAUDE)
 *  \param input   Input list of layers tagged by phase.
 *  \param result  Output list of layers with the two phases interleaved. */
void chirping(const std::vector<GeometryVector> & input, std::vector<GeometryVector> & result ) {
	result.clear();	result.reserve(input.size());

	std::vector<GeometryVector> phase1;
	std::vector<GeometryVector> phase2;
	phase1.reserve(input.size());
	phase2.reserve(input.size());

	for (int i = 0; i < input.size(); i++) {
		if (input[i].x[3] < 1.5)
			phase1.push_back(input[i]);
		else
			phase2.push_back(input[i]);
	}
	MergeSort(phase1, 0);
	MergeSort(phase2, 0);
	if (phase1.size() > phase2.size()) {
		result.push_back(phase1[0]);
		for (int i = 0; i < phase2.size(); i++) {
			result.push_back(phase2[i]);
			result.push_back(phase1[i+1]);
		}
	}
	else{
		int i = 0;
		for (i = 0; i < phase1.size(); i++) {
			result.push_back(phase2[i]);
			result.push_back(phase1[i]);
		}
		if (phase2.size() > phase1.size())
			result.push_back(phase2[i]);
	}
	

	double total = 0.0;
	for (int i = 0; i < result.size(); i++)
		total += result[i].x[0];

	std::cout << "Total size = " << total << "\t";
}

/*
	phase 2 -> particle
*/
/** \brief Rebuild a 1D sphere packing from a sequence of layer thicknesses. (CLAUDE)
 *  \param thickness  Input list of layers; each entry stores a thickness and a phase label.
 *  \param c          Output 1D sphere packing reconstructed from the particle-phase layers. */
void Layer2SpherePacking(const std::vector<GeometryVector> & thickness, SpherePacking & c) {
	double total = 0;
	for (int i = 0; i < thickness.size(); i++)
		total += thickness[i].x[0];

	std::cout << "Total size = " << total <<"\t";

	std::vector<GeometryVector> basis;	basis.emplace_back(GeometryVector(static_cast<double>(total)));
	c = SpherePacking(1, &basis[0], 1.0);

	GeometryVector x(static_cast<int> (1));	double r = 0, length = 0.0;
	for (int i = 0; i < thickness.size(); i++) {
		if (thickness[i].x[3] < 1.5) {//particle phase!!
			x.x[0] = length + 0.5*thickness[i].x[0]; r = 0.5*thickness[i].x[0];
			c.Insert(r, c.CartesianCoord2RelativeCoord(x));
		}
		length += thickness[i].x[0];		
	}
	
}

#include "DisplaySpheres.h"
#include "HardSphereFluids.h"
#include "GenerateConfigs.h"
#include "ParticleShapes.h"
#include "Dispersions.h"

/** \brief Convert a sphere packing into a dispersion of cubes, sizing each cube from its sphere radius. (CLAUDE)
 *  \param initConfig       Input sphere packing.
 *  \param resultConfig     Output cube dispersion.
 *  \param packingFraction  Target packing fraction; negative selects inscribed cubes. */
void SpherePacking2DispersionCube(const SpherePacking & initConfig, Dispersion_cube & resultConfig, double packingFraction=-1) {
	DimensionType d = initConfig.GetDimension();
	int numPrts = initConfig.NumParticle();
	Cube temp(d, 0.1);
	resultConfig = Dispersion_cube(initConfig, temp);
	
	double phi_init = initConfig.PackingFraction();
	double Radius2SideLength = 0; 
	
	if (packingFraction < 0) {/*inscribed cubes*/
		Radius2SideLength = 2.0 / sqrt(d);
	}
	else {
		Radius2SideLength = pow(HyperSphere_Volume(d,1.0)*  packingFraction / phi_init, 1.0 / (double)d);
	}
	
	double L = 0;
	for (int i = 0; i < numPrts; i++) {
		L = Radius2SideLength * initConfig.GetCharacteristics(i);
		resultConfig.UpdateCharacteristics(Cube(d, L), i);
	}
	std::cout << "phi = " << resultConfig.PackingFraction()<<std::endl;
	if(Verbosity > 4)
		std::cout << "Transformation is done!\n";
}

int CollectiveCoordinateMultiRun(Configuration * pConfig, Potential * pPotential, RandomGenerator & gen, double chi, std::string Prefix, size_t SampleNumber, time_t TimeLimit, const std::string & AlgorithmName);

/** \brief A class to compute the collective-coordinates potential + truncated Lennard-Jones potential
 *  **/
class RepulsiveCCPotential : public Potential 
{
protected:
	const Configuration * pConfig;	//< A configuration of interest
	double rc;						//< The critical radius over which the potential becomes 0.
	const double e_critical = 4.0*(pow(2, -12/6) - pow(2, -6/6));//< LJ potential at the critical radius.
public:
	double e0,	//< ennergy scale for LJ 
		sigma;	//< length scale for LJ
	ShiftedCCPotential * CCPotential;

	/** \brief Construct the potential with a stealthy collective-coordinate kernel and a soft-core repulsion. (CLAUDE)
	 *  \param dim      Spatial dimension of the configuration.
	 *  \param epsilon  Energy scale of the soft-core repulsion.
	 *  \param sigma    Length scale (contact diameter) of the soft-core repulsion. */
	RepulsiveCCPotential(DimensionType dim, double epsilon, double sigma) : Potential (dim) {
		CCPotential = new ShiftedCCPotential(dim);
		e0 = epsilon;	this->sigma = sigma;	this->rc = sigma;/*this->rc = pow(2, 1.0/6.0)*sigma;*/
	}

	/** \brief Destructor. (CLAUDE) */
	virtual ~RepulsiveCCPotential() {}

	/** \brief Compute the truncated LJ potential and collective-coordiante potentials for a given configuraiton. */
	virtual double Energy() {
		double LJ = 0.0L, CC = 0.0L;
		this->pConfig->PrepareIterateThroughNeighbors(this->rc);
#pragma omp parallel for schedule (guided, 1) num_threads(this->ParallelNumThread)
		for (int i = 0; i < pConfig->NumParticle(); i++) {
			double E = 0.0L;
				pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->rc, 
					[&i, &E, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, size_t SourceParticle)->void {
					/*Prevent double-counting */
					if (i < SourceParticle) { 
						/*double r_2 = this->sigma * this->sigma / shift.Modulus2();
 * 						double r_6 = r_2 * r_2 * r_2, r_12 = r_6 * r_6;
 * 												LJ += this->e0*(4.0*(r_12 - r_6) - this->e_critical);*/
						double r = sqrt(shift.Modulus2()) / this->sigma;
						E += this->e0*(1.0 - r)*(1.0 - r);
					}
					}
				);
#pragma omp atomic
				LJ += E;
		}
		CC = this->CCPotential->Energy();

		return LJ + CC;
	}

	/** \brief Compute the force exerted at particle i
	 *  \param[out] force	Net force in a GeometryVector.
	 *  \param[in] i			Particle index. */
	virtual void Force(GeometryVector & force, size_t i) {
		GeometryVector force_LJ(static_cast<int> (this->Dimension));

		pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->rc,
			[&i, &force_LJ, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, size_t SourceParticle)->void {
			/*Prevent self-interaction */
			if (i != SourceParticle) {
				/*double r_2 = this->sigma * this->sigma / shift.Modulus2();
 * 				  double r_6 = r_2 * r_2 * r_2;
 * 				  				  force_LJ = force_LJ - 24.0*this->e0*r_6*(1 - 2.0*r_6)*r_2 / (this->sigma*this->sigma) * shift*/
				force_LJ = force_LJ - 2.0*this->e0 / this->sigma * (1.0 / sqrt(shift.Modulus2()) - 1.0 / this->sigma) * shift;
			}
		});

		GeometryVector force_cc;
		this->CCPotential->Force(force_cc, i);
		force = force_LJ + force_cc;
	}

	/** \brief Compute the force exerted at all particles.
	 *  \param[out] result	result[i] is the force exerted at particle i. */
	virtual void AllForce(std::vector<GeometryVector> & result) {
		this->CCPotential->AllForce(result);
		this->pConfig->PrepareIterateThroughNeighbors(this->rc);
#pragma omp parallel for schedule (guided) num_threads(this->ParallelNumThread)
		for (int i = 0; i < pConfig->NumParticle(); i++) {
			GeometryVector force;
			this->Force(force, i);
			result[i] = result[i] + force;
		}
	}

	/** \brief Define a configuration to compute energy and forces.
	 *  \param c	A configuration of interest. */
	virtual void SetConfiguration(const Configuration & c) {
		this->pConfig = &c;
		this->CCPotential->SetConfiguration(c);
	}

};

/** \brief Relax a configuration toward a stealthy ground state, logging energy and saving snapshots along the way. (CLAUDE)
 *  \param pConfig  Configuration to evolve (modified in place).
 *  \param chi      Target stealthiness parameter setting the number of constrained wavevectors. */
void ObserveEvl(Configuration &pConfig, double chi) {
	DimensionType Dimension = pConfig.GetDimension();
	size_t NumP = pConfig.NumParticle();

	ShiftedCCPotential * pPotential = new ShiftedCCPotential(Dimension);
	size_t min_k_points = std::ceil(chi * NumP * Dimension);
	double max_k;

	max_k = (std::sqrt((double)(Dimension)) / 2) * std::pow((2 * (double)(Dimension)* min_k_points), (1.0 / (double)(Dimension)))*(2 * pi / std::pow(pConfig.PeriodicVolume(), 1.0 / (double)(Dimension)));

	std::vector<GeometryVector> Reciprocal;
	for (DimensionType i = 0; i < Dimension; i++)
		Reciprocal.push_back(pConfig.GetReciprocalBasisVector(i));
	PeriodicCellList<Empty> lReciprocal(Dimension, &Reciprocal[0], std::sqrt(Reciprocal[0].Modulus2())*1.1, false);
	lReciprocal.Insert(Empty(), GeometryVector(Dimension));

	std::vector<GeometryVector> kVecs;
	lReciprocal.IterateThroughNeighbors(GeometryVector(Dimension), max_k, [&Dimension, &kVecs](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
	{
		bool NonTrivial = false;//is not vector zero and is not symmetrical to another vector
		for (DimensionType i = 0; i < Dimension; i++)
		{
			if (PeriodicShift[i] > 0)
			{
				NonTrivial = true;
				break;
			}
			else if (PeriodicShift[i] < 0)
			{
				NonTrivial = false;
				break;
			}
		}
		if (NonTrivial)
			kVecs.push_back(shift);
	});

	std::sort(kVecs.begin(), kVecs.end(), [](const GeometryVector & left, const GeometryVector & right) ->bool {return left.Modulus2() < right.Modulus2(); });
	size_t CountKPoints = 0;
	double LastNorm2 = 0;
	while (CountKPoints < kVecs.size())
	{
		if (!(CountKPoints < min_k_points || std::abs((LastNorm2 - kVecs[CountKPoints].Modulus2()) / LastNorm2) < 1e-10))
			break;
		LastNorm2 = kVecs[CountKPoints].Modulus2();
		pPotential->constraints.push_back(ShiftedCCPotential::KPointConstraint(kVecs[CountKPoints], 1.0));
		CountKPoints++;
	}

	pPotential->ParallelNumThread = omp_get_num_threads();
	pPotential->SetConfiguration(pConfig);
	chi = (double)(pPotential->constraints.size()) / ((NumP - 1)*Dimension);
	
	std::cout << "chi = " << chi << std::endl;
	
	std::cout << pPotential->Energy() << std::endl;


	double energy = 0;
	int stepsize = 1, i = 0;
	std::vector<GeometryVector> log;
	do {
		RelaxStructure_MINOP_withoutPertubation(pConfig, *pPotential, 0.0, 0, 0.0, stepsize);
		pPotential->SetConfiguration(pConfig);
		energy = pPotential->Energy();
		i+=stepsize;
		WriteConfiguration(pConfig, "Evolution/config-"+std::to_string(i));
		log.emplace_back(GeometryVector(i, energy));
		std::cout<< i <<": "<< energy<<std::endl;
	} while (energy > 1e-19);

	WriteFunction(log, "Evolution/logs_RSA2SHU");
}

/** \brief Scratch driver collecting ad-hoc experiments (hard-sphere fluids, stealthy packings, etc.). (CLAUDE)
 *  \return  0 on completion. */
int DebugSimple()
{

	char tempstring[1000];
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[300] = {};
	RandomGenerator rngGod(0), rng(0);
	
         {
        ofile << "Simple code for generating hard sphere fluids with a given packing frations in canonical ensemble\n";
                std::string dir = masterDir + std::string("/results/HS_Fluids/");
                int dimension = 0, N = 0, numConfig = 0, samInter = 0, numThreads=1, startIdx=0;
                double phi = 0.0, rho = 1.0;

                ofile << "dimension = ";                                ifile >> dimension;     ofile << dimension <<"\n";
                ofile << "the number of particles = ";  ifile >> N;             ofile << N << "\n";
                ofile << "the number density = ";               ifile >> rho;   ofile << rho << "\n";
                ofile << "packing fraction = ";                 ifile >> phi;   ofile << phi << "\n";
                ofile << "sampling interval (MC  cycles) = "; ifile >> samInter; ofile << samInter << "\n";
                ofile << "number of configurations = "; ifile >> numConfig; ofile << numConfig << "\n";
                ofile << "numThreads = ";       ifile >> numThreads; ofile << numThreads << "\n";
                ofile << "start idx = ";        ifile >> startIdx; ofile << startIdx <<"\n";
                int iterPerSeed = numConfig / numThreads;
                int count = 0;
                int start = (int)floor(startIdx / (double) numThreads + 0.1);
#pragma omp parallel for num_threads(numThreads)
                for (int i = 0; i < numThreads; i++) {
                        HSF_Gen Seed(dimension, N, phi, rho);
                        Seed.SetInterval(samInter);
                        Seed.setting.SetRandomSeed(i+std::time(NULL));
                        SpherePacking c;

                        char filename[300];
                        sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i);
                        if(start > 0){
                                ReadConfiguration(c, filename);
                                Seed.setting.ReplaceConfig(c);
                        }

                        for (int j = start; j < iterPerSeed; j++) {
                                Seed.GenerateP(c);
                                sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i + j * numThreads);
                                if(c.CheckOverlap())
                                #pragma omp critical
                                {
                                        ofile << i + j * numThreads <<"th config has overlapping particles).\n";
                                }
                                else{
                                	WriteConfiguration(c, filename);
				#pragma omp critical
                                {
                                        ofile << i + j * numThreads <<std::endl;
                                }
                                }
                        }
                }

                auto getC = [&dir, &dimension, &phi, &N](size_t i)->Configuration {
                        char filename[300] = "";
                        sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i);
                        SpherePacking c;
                        ReadConfiguration(c, filename);
                        Configuration result(c,"a");
                        return result;
                };
                omp_set_num_threads(numThreads);
                double K = 4.0*pi / pow(N, 1.0 / (double)dimension);
                std::vector<GeometryVector> result;
                IsotropicStructureFactor(getC, numConfig, 3*K, 3*K, result, K);
                WriteFunction(result, (std::string("temp_")+std::to_string(phi)).c_str());
                ofile << "S(0): simulation \n";
                ofile << result[0];
                ofile << "S(0): theory \n";
                if(dimension==3)
                        ofile << pow(1.0 - phi, 3) / (1.0 + phi + 0.34406 * phi*phi - 0.12802 * phi*phi*phi) << std::endl;
        }
        return 0;




	{
		double K1, K2, val, chi = 0, L = 10, sigma, phi;
		
		DimensionType dim = 2;
		size_t num = L * L*L, numConfig = 300;
		int num_threads = 4;
		ofile << "Dimension = "; ifile >> dim;
		ofile << "K1 = ";	ifile >> K1;
		ofile << "K2 = ";	ifile >> K2;
		ofile << "val = ";	ifile >> val;
		ofile << "sigma = "; ifile >> sigma;
		ofile << "phi = "; ifile >> phi;
		ofile << "# threads = "; ifile >> num_threads;
		ofile << "num particle = "; ifile >>num; 
		ofile << "num configs = "; ifile >> numConfig;
		L = pow(num, 1.0/(double)dim);	
		
		double a = pow(phi/(HyperSphere_Volume(dim, 1.0)), 1.0/(double)dim);
		sprintf(name, "%s/configs/%dD_Stealthy_Ka-%0.1f-%0.3f__%0.1f-%0.3f__L_%0.1f", masterDir.c_str(), dim, K1, K2, val, sigma, L);
		std::string prefix(name);
		
		ofile << "NumParticle = "<<num<<std::endl;
		ofile << "directory = "<< prefix<<std::endl;
		ofile << "a = "<< a<<std::endl;
		

		Configuration pConfig = GetUnitCubicBox(dim, 0.1);
		pConfig.Rescale(L);
		for (size_t i = 0; i < num; i++)
			pConfig.Insert("a", rngGod);

		pConfig.PrepareIterateThroughNeighbors(sigma);
		K1 /= a; K2 /= a;
		RepulsiveCCPotential * potential = new RepulsiveCCPotential(dim, val, sigma);
		{
			potential->ParallelNumThread = num_threads;
			potential->CCPotential->ParallelNumThread = num_threads;
			double K1_modulus = K1 * K1;
			std::vector<GeometryVector> ks_temp = GetKs(pConfig, K2, K2, 1);
			for (auto k = ks_temp.begin(); k != ks_temp.end(); k++) {
				if (k->Modulus2() > K1_modulus) {
					potential->CCPotential->constraints.emplace_back(*k, 1);
					chi++;
				}
				
			}
			chi /= dim * (num - 1);
		}
		ofile << "chi = " << chi << std::endl;

		CollectiveCoordinateMultiRun(&pConfig, potential, rngGod, chi, prefix, numConfig, 10000000000, "MINOP");
		delete potential;

		ConfigurationPack ConfigSet((prefix + "_Success").c_str());
		double rmin = L;
		{
			for (int i = 0; i < ConfigSet.NumConfig(); i++) {
				Configuration x = ConfigSet.GetConfig(i);
				double rm = L;
				for (size_t j = 0; j < x.NumParticle(); j++) {
					double temp = x.NearestParticleDistance(j);
					rm = (temp < rm) ? temp : rm;
				}
				ofile << i << ":\t" << rm << std::endl;
			}

		}




		auto getC = [&rmin, &ConfigSet](size_t i) ->Configuration {
			Configuration c = ConfigSet.GetConfig(i);
			double R = rmin;
			for (size_t j = 0; j < c.NumParticle(); j++) {
				double temp = c.NearestParticleDistance(j);
				R = (temp < R) ? temp : R;
			}
			rmin = (R < rmin) ? R : rmin;
			return c;
		};

		/*Compute g2 function*/
		std::vector<GeometryVector> g2;
		numConfig = ConfigSet.NumConfig();
		IsotropicTwoPairCorrelation(getC, numConfig, L / 2, g2);
		WriteFunction(g2, (prefix + "_g2").c_str());
		phi = 1.0*HyperSphere_Volume(dim, rmin / 2.0);
		ofile << "max packing fraction is " << phi << std::endl;


		auto getS = [&rmin,&a, &ConfigSet](size_t i) ->SpherePacking {
			Configuration c = ConfigSet.GetConfig(i);
			
			SpherePacking s(c, a);
			return s;
		};

		/*Compute the structure factsor*/
		std::vector<GeometryVector> Sk;
		IsotropicSpectralDensity(getS, numConfig, 5, 30, Sk, 2.0*pi/L);
		for(auto iter = Sk.begin(); iter!=Sk.end(); iter++){
			iter->x[0] *= a;
		}
		WriteFunction(Sk, (prefix + "_chik").c_str());
		
		/*Compute the spectral density function*/
		std::vector<GeometryVector> SpecDen;		
	}
	return 0;

	 {
	ofile << "Simple code for generating hard sphere fluids with a given packing frations in canonical ensemble\n";
                std::string dir = masterDir + std::string("/results/HS_Fluids/");
                int dimension = 0, N = 0, numConfig = 0, samInter = 0, numThreads=1, startIdx=0;
                double phi = 0.0, rho = 1.0;

                ofile << "dimension = ";                                ifile >> dimension;     ofile << dimension <<"\n";
                ofile << "the number of particles = ";  ifile >> N;             ofile << N << "\n";
                ofile << "the number density = ";               ifile >> rho;   ofile << rho << "\n";
                ofile << "packing fraction = ";                 ifile >> phi;   ofile << phi << "\n";
                ofile << "sampling interval (MC  cycles) = "; ifile >> samInter; ofile << samInter << "\n";
                ofile << "number of configurations = "; ifile >> numConfig; ofile << numConfig << "\n";
                ofile << "numThreads = ";       ifile >> numThreads; ofile << numThreads << "\n";
		ofile << "start idx = ";	ifile >> startIdx; ofile << startIdx <<"\n";
                int iterPerSeed = numConfig / numThreads;
                int count = 0;
		int start = (int)floor(startIdx / (double) numThreads + 0.1);
#pragma omp parallel for num_threads(numThreads)
                for (int i = 0; i < numThreads; i++) {
                        HSF_Gen Seed(dimension, N, phi, rho);
                        Seed.SetInterval(samInter);
                        Seed.setting.SetRandomSeed(i+std::time(NULL));
                        SpherePacking c;
		
                        char filename[300];
			sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i);
			if(start > 0){
				ReadConfiguration(c, filename);
				Seed.setting.ReplaceConfig(c);		
			}

                        for (int j = start; j < iterPerSeed; j++) {
                                Seed.GenerateP(c);
                                sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i + j * numThreads);
                                if(c.CheckOverlap())
				#pragma omp critical
				{
					ofile << i + j * numThreads <<"th config has overlapping particles).\n";
				}
				else{
				WriteConfiguration(c, filename);
#pragma omp critical
                                {
                                        ofile << i + j * numThreads <<std::endl;
                                }
				}
                        }
                }

                auto getC = [&dir, &dimension, &phi, &N](size_t i)->Configuration {
                        char filename[300] = "";
                        sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i);
                        SpherePacking c;
                        ReadConfiguration(c, filename);
                        Configuration result(c,"a");
                        return result;
                };
                omp_set_num_threads(numThreads);
                double K = 4.0*pi / pow(N, 1.0 / (double)dimension);
                std::vector<GeometryVector> result;
                IsotropicStructureFactor(getC, numConfig, 3*K, 3*K, result, K);
                WriteFunction(result, (std::string("temp_")+std::to_string(phi)).c_str());
                ofile << "S(0): simulation \n";
                ofile << result[0];
                ofile << "S(0): theory \n";
		if(dimension==3)                
			ofile << pow(1.0 - phi, 3) / (1.0 + phi + 0.34406 * phi*phi - 0.12802 * phi*phi*phi) << std::endl;
        }
        return 0;




	{
		omp_set_num_threads(10);
		RSA_Gen gen(2, 400, PhiC(2), 1);
		Configuration c;
		do {
			gen.GenerateC(c);
			while (c.NumParticle() > 400) {
				c.DeleteParticle(size_t(round(c.NumParticle() * rng.RandomDouble())));
			}
		} while (c.NumParticle() != 400);
		WriteConfiguration(c, "Evolution/config-0");
		ofile << "Initial configuration is successfully generated\n";
		ObserveEvl(c, 0.45);

	}

	return 0;




	{/* a piece of code to compute the stealthy hyperuniform packing ... */
		double K1, K2, val, chi = 0, L = 10, sigma, phi = 0.1;

		DimensionType dim = 3;
		size_t num = L * L*L, numConfig = 10;
		int num_threads = 10;
		
		ofile << "num prts = ";	ifile >> num;
		ofile << "chi = ";	ifile >> chi;

		L = pow(num, 1.0/dim);
		
		double a = pow(3.0*phi/(4.0*pi), 1.0/3.0);
		sprintf(name, "./configs/3D_SHU_0.49_%d_soft", num);
		std::string prefix(name);
		
		ofile << "NumParticle = "<<num<<std::endl;
		ofile << "directory = "<< prefix<<std::endl;
		ofile << "a = "<< a<<std::endl;
		

		Configuration pConfig = GetUnitCubicBox(dim, 0.1);
		pConfig.Rescale(L);
		for (size_t i = 0; i < num; i++)
			pConfig.Insert("a", rngGod);

		pConfig.PrepareIterateThroughNeighbors(sigma);
		RepulsiveCCPotential * potential = new RepulsiveCCPotential(dim, val, sigma);
		{
			potential->ParallelNumThread = num_threads;
			potential->CCPotential->ParallelNumThread = num_threads;

			size_t min_k_points = std::ceil(chi * num * dim );
			double max_k = (std::sqrt((double)(dim)) / 2) * std::pow( (2 * (double)(dim) * min_k_points) , (1.0/(double)(dim)) ) *(2*pi/std::pow(pConfig.PeriodicVolume(), 1.0/(double)(dim) ));
			std::vector<GeometryVector> Reciprocal;
			
			for(DimensionType i=0; i<dim; i++)
				Reciprocal.push_back(pConfig.GetReciprocalBasisVector(i));
		
			PeriodicCellList<Empty> lReciprocal(dim, &Reciprocal[0], std::sqrt(Reciprocal[0].Modulus2())*1.1, false);
			lReciprocal.Insert(Empty(), GeometryVector(dim));
			
			std::vector<GeometryVector> kVecs;
			lReciprocal.IterateThroughNeighbors(GeometryVector(dim), max_k, [&dim, &kVecs] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
			bool NonTrivial=false;//is not vector zero and is not symmetrical to another vector
			for(DimensionType i=0; i<dim; i++)
			{
				if(PeriodicShift[i]>0)
				{
					NonTrivial=true;
					break;
				}
				else if(PeriodicShift[i]<0)
				{
					NonTrivial=false;
					break;
				}
			}
			if(NonTrivial)
				kVecs.push_back(shift);
			});

			std::sort(kVecs.begin(), kVecs.end(), [] (const GeometryVector & left, const GeometryVector & right) ->bool {return left.Modulus2()<right.Modulus2();});
		size_t CountKPoints=0;
		double LastNorm2=0;
		while(CountKPoints<kVecs.size())
		{
			if( ! (CountKPoints<min_k_points || std::abs((LastNorm2-kVecs[CountKPoints].Modulus2())/LastNorm2)<1e-10) )
				break;
			LastNorm2=kVecs[CountKPoints].Modulus2();
			potential->CCPotential->constraints.push_back(ShiftedCCPotential::KPointConstraint(kVecs[CountKPoints], 1.0));
			CountKPoints++;
		}
			ofile << "chi = " << CountKPoints / (double)(dim*(num-1)) << std::endl;
		}
	

		CollectiveCoordinateMultiRun(&pConfig, potential, rngGod, chi, prefix, numConfig, 10000000000, "LBFGS");
		delete potential;

		ConfigurationPack ConfigSet((prefix + "_Success").c_str());
		double rmin = L;
		{
			for (int i = 0; i < ConfigSet.NumConfig(); i++) {
				Configuration x = ConfigSet.GetConfig(i);
				double rm = L;
				for (size_t j = 0; j < x.NumParticle(); j++) {
					double temp = x.NearestParticleDistance(j);
					rm = (temp < rm) ? temp : rm;
				}
				ofile << i << ":\t" << rm << std::endl;
			}

		}



		auto getC = [&rmin, &ConfigSet](size_t i) ->Configuration {
			Configuration c = ConfigSet.GetConfig(i);
			double R = rmin;
			for (size_t j = 0; j < c.NumParticle(); j++) {
				double temp = c.NearestParticleDistance(j);
				R = (temp < R) ? temp : R;
			}
			rmin = (R < rmin) ? R : rmin;
			return c;
		};

		/*Compute g2 function*/
		std::vector<GeometryVector> g2;
		numConfig = ConfigSet.NumConfig();
		IsotropicTwoPairCorrelation(getC, numConfig, L / 2, g2);
		WriteFunction(g2, (prefix + "_g2").c_str());
		phi = 1.0*HyperSphere_Volume(dim, rmin / 2.0);
		ofile << "max packing fraction is " << phi << std::endl;

		/*Compute the structure factsor*/
		std::vector<GeometryVector> Sk;
		IsotropicStructureFactor(getC, numConfig, 5, 10, Sk);
		WriteFunction(Sk, (prefix + "_Sk").c_str());		
	}
	return 0;

	
		


	

	{	/*Local volume-fraction variance in direct space.*/
		double R, dR, rho = 1, phi_1, phi_2, dphi, D;
		int d, numSpacing = 1, numSampling, numConfigs;
		std::string latticeName;

		ofile << "Lattice Name = ";
		ifile >> latticeName;
		ofile << "The number of configurations = ";
		ifile >> numConfigs;
		/*ofile << "rho = ";		ifile >> rho;*/
		ofile << "D = ";	ifile >> D;
		ofile << "smallest packing fraction = ";
		ifile >> phi_1;
		/*ofile << "the number of spacing = ";
 * 		ifile >> numSpacing;*/
		ofile << "dphi = ";
		ifile >> dphi;


		ofile << "---------Local volume-fraction variances --------------\n";
		ofile << "R_max = ";	ifile >> R;
		ofile << "dR = ";		ifile >> dR;
		ofile << "the number of samplings = ";	ifile >> numSampling;

		LatticeGen lattice0(latticeName, 1);
		LatticeGen lattice(latticeName, 10);
		{
			SpherePacking c0;
			lattice0.GenerateP(c0);
			double nearestNeighbor = c0.NearestParticleDistance(0);
			d = c0.GetDimension();
			rho = lattice0.GetRho()*pow(D / nearestNeighbor, d);
			lattice0.SetRho(rho);				lattice.SetRho(rho);
			lattice0.GenerateP(c0);

			D = c0.NearestParticleDistance(0);
			ofile << "UnitCell is rescaled to have rho = " << rho << ",\tD = " << D << "\n";
		}
		SpherePacking c, c0;

		/* obtain a sampling points */
		std::vector<GeometryVector> samplingPt;	samplingPt.reserve(numSampling);
		ofile << "Obtaining sampling points \n";
		lattice0.GenerateP(c0);	d = c0.GetDimension();
		lattice.GenerateP(c);

		{
			GeometryVector temp((int)d);
			for (int i = 0; i < numSampling; i++) {
				for (int j = 0; j < d; j++)
					temp.x[j] = rngGod.RandomDouble();

				temp = c0.RelativeCoord2CartesianCoord(temp);
				temp = c.CartesianCoord2RelativeCoord(temp);
				samplingPt.push_back(temp);
			}
		}

		phi_2 = lattice.maxPhi;
		/*double dphi = (phi_2 - phi_1) / numSpacing;*/
		numSpacing = (int)ceil((phi_2 - phi_1) / dphi);

		for (int i = 0; i < numSpacing + 1; i++) {
			double phi = (i == numSpacing - 1) ? phi_2 : phi_1 + dphi * i;
			ofile << phi;
			lattice.SetPackingFraction(phi);
			lattice.GenerateP(c);
			auto getC = [&c, &rngGod](size_t i)->SpherePacking {
				GeometryVector displacement;
				SpherePacking result(c);
				int d = result.GetDimension();
				displacement.Dimension = d;
				for (int i = 0; i < d; i++)
					displacement.x[i] = rngGod.RandomDouble();

				for (int i = 0; i < result.NumParticle(); i++) {
					result.MoveParticle(i, displacement + result.GetRelativeCoordinates(i));
				}

				return result; };
			std::vector<GeometryVector> result;
			omp_set_num_threads(4);
			MC_vfv(getC, numConfigs, result, R, dR, samplingPt);
			ofile << " **** done\n";
			sprintf(name, "%s_lvfv__phi_%0.4f", latticeName.c_str(), phi);
			WriteFunction(result, name);
		}


	}

	return 0;


	



{
	RandomGenerator rngGod;
		std::vector<GeometryVector> v;
		std::string dir;
		DimensionType dim = 0;
		char prefix[200] = "";
		int num = 0, iter=0, orient;
		ofile << "Dimension = ";
		ifile >> dim;
		ofile << "paticle number = ";
		ifile >> num;
		ofile << "iterations = ";
		ifile >> iter;
		ofile << "RandomOrientation (0=y/1=n)";
		ifile >> orient;
		ofile << "directory = ";
		ifile >> dir;
		ofile << "outprefix = ";
		ifile >> prefix;

		double kmin = 0, kLinear = 0, kCircular = 0;
		ofile << "kmin = ";
		ifile >> kmin;
		ofile << "kLinear = ";
		ifile >> kLinear;
		ofile << "kCircular = ";
		ifile >> kCircular;

		Cube prt(dim, 2.0);
		prt.Rotate(rngGod);
		std::vector<GeometryVector> X;
		std::vector<GeometryVector> result2;

		auto getC = [&dir, &orient](size_t i)->Dispersion_cube {
			SpherePacking c;	Dispersion_cube Test;
			char file[300];
			sprintf(file, dir.c_str(), i);
			ReadConfiguration(c, file);
			SpherePacking2DispersionCube(c, Test, 0.01);
			if (orient == 0)
				Test.RotateAllParticles();
			return Test;
		};

		std::vector<GeometryVector> result;
		omp_set_num_threads(10);

		sprintf(name, "%s%s_chik", masterDir.c_str(), prefix);
		ofile << name <<std::endl;

		IsotropicSpectralDensityCube(getC, iter, kCircular, kLinear, result, kmin);

		WriteFunction(result, name);

}
return 0;








{
		ConfigGen * initC;
		ofile << "Input parameters for the RSA packings\n";

		int dim = 0, numPrts, numConfig, numThreads;
		ofile << "num threads = ";
		ifile >> numThreads; ofile << numThreads << std::endl;
		double packingFraction, numberdensity;
		if (dim == 0) {
			ofile << "Spatial dimension = ";	ifile >> dim;	ofile << dim << "\n";
		}
		ofile << "Particle Number = ";
		ifile >> numPrts; ofile << numPrts << std::endl;
		ofile << "Packing Fraction = ";
		ifile >> packingFraction;	ofile << packingFraction << std::endl;
		ofile << "numberdesnity = ";
		ifile >> numberdensity;	ofile << numberdensity << std::endl;
		ofile << "number of configuration = ";
		ifile >> numConfig;	ofile << numConfig << std::endl;
		initC = new RSA_Gen(dim, numPrts, packingFraction, numberdensity);
		
		initC->SetNumThreads(numThreads);

		for (int i = 0; i < numConfig; i++) {
			SpherePacking c;
			initC->GenerateP(c);
			sprintf(name, "%s/RSA_configs/%dD_RSA_%0.4f__%d_config_init_%4d", masterDir.c_str(), dim, packingFraction, numPrts, i);
			WriteConfiguration(c, name);
			ofile << "*" << std::endl;
		}
	}
	return 0;

	
/*
	{
		RandomGenerator rngGod;
		std::vector<GeometryVector> v;

		Parallelepiped test1(2, 1.0);
		//test1.Rotate(rngGod);

		test1.GetVertices(v);

		for (auto iter = v.begin(); iter != v.end(); iter++)
			ofile << *iter;

		test1.GetDimension();
	}
	return 0;
*/
	{
		int dimension = 2;
		double L = 10;
		std::vector<GeometryVector> basis;	
		basis.emplace_back(GeometryVector(L, 0.0));
		basis.emplace_back(GeometryVector(0.0, L));
		for (int intensity = 1; intensity <= 100; intensity *= 10) {
			ofile << "--------------------------------------" << std::endl;
			ofile << "Intensity = " << intensity << std::endl;
			ofile << "--------------------------------------" << std::endl;
			std::vector<GeometryVector> reserve;
			std::vector<Configuration> c_list;
			for (int i = 1; i <= 5; i++) {
				sprintf(name, "./results/HIP/2D-HIP-intensity-%d-run-%d.dat", intensity, i);
				std::fstream read(name, std::fstream::in);
				ReadFunction(reserve, read, 3);
				
				/*  Mixing indices */
				
				Configuration c(dimension, &basis[0], .1);
				
				for (int j = 0; j < reserve.size(); j++) 
					c.Insert("a", GeometryVector(reserve[j].x[1], reserve[j].x[2]));
				
				c_list.push_back(c);
			}

			for (double r_exc = 0.00; r_exc < 1; r_exc += 2) {
				sprintf(name, "./results/HIP/2D-HIP-%d__r_%0.2f", intensity, r_exc);
				ConfigurationPack save(name);
				ofile << "\nr_exc = " << r_exc << std::endl;
				for (int i = 0; i < c_list.size(); i++) {
					std::vector<int> idx_list;	//idx_list[i] == 0 => no overlap
					for (int j = 0; j < c_list[i].NumParticle(); j++)	idx_list.push_back(0);

					for (int j = 0; j < idx_list.size(); j++) {
						if(idx_list[j] == 0)
						c_list[i].IterateThroughNeighbors(c_list[i].GetRelativeCoordinates(j), r_exc, 
							[&j, &idx_list](const GeometryVector & shift, const GeometryVector &dummy, const signed long * dummy2, size_t pID)->void {
							if (pID != j)
								idx_list[pID]++;
						});
					}
					Configuration c_(c_list[i]);	c_.RemoveParticles();

					for (int j = 0; j < idx_list.size(); j++)
						if (idx_list[j] == 0)
							c_.Insert("a", c_list[i].GetRelativeCoordinates(j));

					ofile << i << ": " << c_.NumParticle() << std::endl;
					save.AddConfig(c_);
				}
			}
		}


	}
	return 0;


	{
		int numConfig;
		PoissonPt get(ifile, ofile);
		ofile << "Num Configs = ";
		ifile >> numConfig;
		int d = 1;
		int numPrts = get.GetNumPrts();
		double rho = get.GetRho();

		ofile << "Exclusion volume \n";
		double dR, Rmax;
		ofile << "Rmax = "; ifile >> Rmax;
		ofile << "dR = ";	ifile >> dR;

		std::vector<GeometryVector> Ev;
		std::vector<GeometryVector> MaxR(numConfig, GeometryVector(0.0,0.0));
		for (double r = dR; r <= Rmax; r += dR)
			Ev.emplace_back(GeometryVector(r, 0.0, 0.0, 0.0));
		
#pragma omp parallel for num_threads(4)
		for (int i = 0; i < numConfig; i++) {
			std::vector<double> NND;
			std::vector<GeometryVector> tempEV(Ev);
			char file[200] = "";
			sprintf(file, "./results/Poisson/%dD_Poisson_init_%0.2f_%d__%d", d, rho, numPrts, i);
			SpherePacking c;
			get.GenerateP(c);

			/*Compute Ev(r)*/
			Merge_Sort(c, 0);
			int NUM = numPrts;
			double max_R = 0.0;
			for (int j = 0; j < NUM; j++) {
				double x_beg = 0.0, x_end = 0.0;
				double x_cen = c.GetCartesianCoordinates(j).x[0];
				int idx_beg = j - 1, idx_end = j + 1;
				if (idx_beg < 0) {
					idx_beg = NUM - 1;
					x_beg = -c.GetBasisVector(0).x[0];
				}
				x_beg += c.GetCartesianCoordinates(idx_beg).x[0];

				if (idx_end == NUM) {
					idx_end = 0;
					x_end = +c.GetBasisVector(0).x[0];
				}
				x_end += c.GetCartesianCoordinates(idx_end).x[0];

				/*ind cell bounaries
				x_beg = 0.5*(x_beg + x_cen);
				x_end = 0.5*(x_end + x_cen);
				*/
				x_end = x_end - x_cen;
				max_R = std::max(max_R, x_end);
				NND.push_back(x_end);
			}
			
			MaxR[i].x[0] = max_R;
#pragma omp critical
			{
				ofile << "*";
			}
		}
		ofile << "\ndone!!\n";

		/*for (int i = 0; i < Ev.size(); i++) {
			Ev[i].x[1] /= (double)numConfig;
			Ev[i].x[2] /= (double)numConfig;
			Ev[i].x[2] = sqrt(fabs(Ev[i].x[2] - Ev[i].x[1] * Ev[i].x[1]) / (double)numConfig);
		}*/
		sprintf(name, "./results/Poisson/%dD_Poisson_init_%0.2f_%d__maxR", d, rho, numPrts);
		WriteFunction(MaxR, name);
	}
	return 0;


	{
		int numConfig;
		PoissonPt get(ifile, ofile);
		ofile << "Num Configs = ";
		ifile >> numConfig;
		int d = 1;
		int numPrts = get.GetNumPrts();
		double rho = get.GetRho();

		ofile << "Exclusion volume \n";
		double dR, Rmax;
		ofile << "Rmax = "; ifile >> Rmax;
		ofile << "dR = ";	ifile >> dR;

		std::vector<GeometryVector> Ev;
		for (double r = dR; r <= Rmax; r += dR)
			Ev.emplace_back(GeometryVector(r, 0.0, 0.0, 0.0));

#pragma omp parallel for num_threads(4)
		for (int i = 0; i < numConfig; i++) {
			std::vector<double> NND;
			std::vector<GeometryVector> tempEV(Ev);
			char file[200] = "";
			sprintf(file, "./results/Poisson/%dD_Poisson_init_%0.2f_%d__%d", d, rho, numPrts, i);
			SpherePacking c;
			get.GenerateP(c);

			//Compute Ev(r)
			Merge_Sort(c, 0);
			int NUM = numPrts;
			for (int j = 0; j < NUM; j++) {
				double x_beg = 0.0, x_end = 0.0;
				double x_cen = c.GetCartesianCoordinates(j).x[0];
				int idx_beg = j - 1, idx_end = j + 1;
				if (idx_beg < 0) {
					idx_beg = NUM - 1;
					x_beg = -c.GetBasisVector(0).x[0];
				}
				x_beg += c.GetCartesianCoordinates(idx_beg).x[0];

				if (idx_end == NUM) {
					idx_end = 0;
					x_end = +c.GetBasisVector(0).x[0];
				}
				x_end += c.GetCartesianCoordinates(idx_end).x[0];

				//find cell bounaries
				//x_beg = 0.5*(x_beg + x_cen);
				//x_end = 0.5*(x_end + x_cen);
				x_end = x_end - x_cen;

				NND.push_back(x_end);
			}
			MergeSort(NND);
			int idx = 0, count = 0;
			for (int j = 0; j < tempEV.size(); j++) {
				while (idx < NUM && NND[idx] < tempEV[j].x[0]) {
					count++;
					idx++;
				}
				tempEV[j].x[1] = 1.0-(double)idx / (double)NUM;
				count = 0;
			}
#pragma omp critical
			{
				for (int j = 0; j < Ev.size(); j++) {
					Ev[j].x[1] += tempEV[j].x[1];
					Ev[j].x[2] += tempEV[j].x[1] * tempEV[j].x[1];
				}
			}

			WriteConfiguration(c, file);
		}


		for (int i = 0; i < Ev.size(); i++) {
			Ev[i].x[1] /= (double)numConfig;
			Ev[i].x[2] /= (double)numConfig;
			Ev[i].x[2] = sqrt(fabs(Ev[i].x[2] - Ev[i].x[1]* Ev[i].x[1])/(double)numConfig);
		}
		sprintf(name, "./results/Poisson/%dD_Poisson_init_%0.2f_%d__Ev", d, rho, numPrts);
		WriteFunction(Ev, name);
	}
	return 0;

	{
		std::string prefix;
		int numconfig = 0;
		ofile << "file name = ";
		ifile >> prefix;
		//ConfigurationPack open(tempstring);
		ofile << "Num config = ";
		ifile >> numconfig;
		Configuration c;
		omp_set_num_threads(4);
		ReadConfiguration(c, prefix);

		Display3DStructureFactor(c, 3);

	}
	return 0;

	/*{
		double r[6] = { 0,0,0,0,0,0 };
		std::vector<GeometryVector> result;
		
		std::fstream output("./g2_vs_t.txt", std::fstream::out);

		for (int i = 0; i < 40; i += 2) {
			std::vector<GeometryVector> temp;
			sprintf(name, "./g2_t_%d", i);
			ReadFunction(temp, name);
			int count = 0;
			double R = 0;
			output << i << "\t";
			for (int j = 0; j < temp.size(); j++) {
				R = temp[j].x[0];
				if (((R > 0.92) && (R < 0.96)) || ((R > 3.5) && (R < 3.56))) {
					output << temp[j].x[1] << "\t" << temp[j].x[3] << "\t";

					r[count] += R;
					count++;
				}
			}
			output << "\n";
		}
		for (int i = 0; i < 6; i++)
			ofile << r[i] / 20.0<<"\t";
		output.close();

	}
	return 0;*/
	{
		std::string prefix;
		int numconfig = 0;
		ofile << "file name = ";
		ifile >> prefix;
		//ConfigurationPack open(tempstring);
		ofile << "Num config = ";
		ifile >> numconfig;
		omp_set_num_threads(4);
		std::vector<Configuration> cList;
		/*for (int i = 0; i < numconfig; i+=2) {
			std::vector<GeometryVector> result;
			char filename[1000] = "";
			sprintf(filename, (prefix).c_str(), i);
			Configuration c;
			{
					SpherePacking temp;
				ReadConfiguration(temp, filename);
				c = Configuration(temp, "a");
			}
			auto getC = [&c](size_t j)->Configuration {return c; };
			IsotropicTwoPairCorrelation(getC, 1, 4.0, result);
			sprintf(filename, "g2_t_%d", i);
			WriteFunction(result, filename);
		}*/
		for (int i = 0; i < numconfig; i++) {
			Configuration c;
			{
				char filename[1000] = "";
				sprintf(filename, (prefix).c_str(), i);
				SpherePacking temp;
				ReadConfiguration(temp, filename);
				c = Configuration(temp, "a");
			}
			cList.push_back(c);
		}
		std::vector<GeometryVector> result;
		auto getC = [&cList](size_t i)->Configuration {return cList[i]; };
		IsotropicStructureFactor(getC, numconfig, 0.03, 1, result, 0.002);
		WriteFunction(result, "Sk_2D_HSF");
	}
	return 0;

	{
		int dimension, N, numconfig = 1, series = 0, nskip=1;
		double phi, numdensity;
		bool cont = true;
		omp_set_num_threads(4);
			ofile << "NumConfig = ";	ifile >> numconfig;
			SpherePacking c;
			HSF_Gen ens(ifile, ofile);
			
			double LinearKmax, CircularKmax, dK;
			ofile << "Circular Kmax = "; ifile >>CircularKmax;
			ofile << "Linear Kmax = "; ifile >> LinearKmax;
			ofile << "dK = "; ifile >> dK;

			char prefix[100];	sprintf(prefix, "./results/HS_Fluids/%dD_HSF_%0.4f_%d", ens.setting.Dimension, ens.setting.PackingFraction, ens.setting.NumPrts);
			
			ens.GenerateP(c);
			WriteConfiguration(c, prefix);

			for (; cont;) {
				std::vector<SpherePacking> list;	list.reserve(numconfig);
				for (int i = 0; i < numconfig; i++) {
					ens.GenerateP(c);	list.push_back(c);
				}
				auto getConfig = [&list](size_t i)->Configuration {return Configuration(list[i], "a");  };
				std::vector<GeometryVector> result;
				//Structure factor
				{
					IsotropicStructureFactor(getConfig, numconfig, CircularKmax, LinearKmax, result, dK);
					std::string fname = prefix + std::string("__Sk_")+std::to_string(series);
					WriteFunction(result, fname);
					series++;
				}
				ofile << "continue (y=1/n=0)?"; ifile >> cont;
				if (cont) {
					ofile << "number of 10 MC cycles = ";
					ifile >> nskip;
					double accept = 0.0;
					for (int i = 0; i < nskip; i++) {
						accept = ens.setting.Evolve(10);
						ofile << i << ": acceptance = " << accept << "\n";
					}
				}
				ofile << "\n\n\n";
			}

			//Pair correlation function
			//{
			//	double Rmax, dK;
			//	ofile << "Rmax = "; ifile >> Rmax;
			//	
			//	IsotropicTwoPairCorrelation(getConfig, numconfig, Rmax, result);
			//	std::string fname = prefix + std::string("__g2");
			//	WriteFunction(result, fname);
			//}
			
		
	}
	return 0;

	{
		ofile << "Input file name = ";
		ifile >> tempstring;

		double Kmax = 10.0;
		double dK = 2.0*pi / 100.0;

		bool spectraldensity = true;
		bool config = true;

		ConfigurationPack open(tempstring);
		int numConfig = open.NumConfig();
		double NND = 10000000;
		{		//Determine the nearest neighbor distance
			for (int i = 0; i < numConfig; i++) {
				double nnd = 1000000;
				Configuration temp = open.GetConfig(i);
				for (int j = 0; j < temp.NumParticle(); j++) {
					nnd = std::min(nnd, temp.NearestParticleDistance(j));
				}
				NND = std::min(NND, nnd);
			}

			NND *= 0.5;
		}
		//----------------------------monodisperse packing
		double phaseFraction = 0.0;
		std::vector<SpherePacking> mono;	mono.reserve(open.NumConfig());
		ofile << "-----monodisperse packing------\n";
		for (int i = 0; i < numConfig; i++) {
			double phii;
			Configuration ci = open.GetConfig(i);
			Merge_Sort(ci, 0);
			mono.emplace_back(SpherePacking(ci, NND));
			phii = mono[i].PackingFraction();
			ofile << i << " : " << phii << std::endl;
			phaseFraction += phii;
		}
		ofile << "=> average : " << phaseFraction / (double)numConfig << "\n\n";

		std::vector<SpherePacking> initConfig;
		std::vector<SpherePacking> chirpedConfig;
		auto getConfig1 = [&initConfig](size_t i)->SpherePacking { return initConfig[i]; };
		auto getConfig2 = [&chirpedConfig](size_t i)->SpherePacking { return chirpedConfig[i]; };

		{
			initConfig = std::vector<SpherePacking>(mono);
			//Chirp the packing
			std::vector<GeometryVector> thickness_init, thickness_fin;
			for (int i = 0; i < numConfig; i++) {
				ofile << i << " : ";
				SpherePacking c;
				SpherePacking2Layers(initConfig[i], thickness_init);
				chirping(thickness_init, thickness_fin);
				Layer2SpherePacking(thickness_fin, c);
				chirpedConfig.push_back(c);
				ofile << "\n";
				if (config) {//Config
							 //initial
					sprintf(name, "%s__mono_________PackingFormat_%d", tempstring, i);
					WriteConfiguration(initConfig[i], name);
					sprintf(name, "%s__mono_________ThicknessFormat_%d", tempstring, i);
					WriteFunction(thickness_init, name);

					sprintf(name, "%s__mono_chirped_PackingFormat_%d", tempstring, i);
					WriteConfiguration(chirpedConfig[i], name);
					sprintf(name, "%s__mono_chirped_ThicknessFormat_%d", tempstring, i);
					WriteFunction(thickness_fin, name);
				}
			}

			//Spectral densities
			if (spectraldensity) {
				std::vector<GeometryVector> chiv;
				sprintf(name, "%s__mono_________chiv", tempstring);
				IsotropicSpectralDensity(getConfig1, numConfig, Kmax, Kmax, chiv, dK);
				//WriteFunction(chiv, name);
				PlotFunction_Grace(chiv, name, "k", "\\xc\\0\\sV\\N(k)");
				sprintf(name, "%s__mono_chirped_chiv", tempstring);
				IsotropicSpectralDensity(getConfig2, numConfig, Kmax, Kmax, chiv, dK);
				//WriteFunction(chiv, name);
				PlotFunction_Grace(chiv, name, "k", "\\xc\\0\\sV\\N(k)");
			}


		}

		//--------------------------------------Layered
		{
			ofile << "========Layered structures\n";
			double phi = 0, phi2 = 0;
			//Generate Configurations!!!!!!!!!!
			//initConfig = std::vector<SpherePacking>(mono);
			//Chirp the packing
			std::vector<GeometryVector> thickness_init, thickness_fin;
			for (int i = 0; i < numConfig; i++) {
				thickness_init.clear();	thickness_fin.clear();
				ofile << i << " : ";
				{// Generate in terms of later thickenss!
					thickness_init.emplace_back(GeometryVector(mono[i].GetCartesianCoordinates(0).x[0], 0.0, 0.0, 2.0));	//matrix phase
					for (int j = 0; j < mono[i].NumParticle() - 1; j++) {
						if (j % 2 == 0) {
							thickness_init.emplace_back(GeometryVector(mono[i].GetCartesianCoordinates(j + 1).x[0] - mono[i].GetCartesianCoordinates(j).x[0], 0.0, 0.0, 1.0)); //particle phase
						}
						else {
							thickness_init.emplace_back(GeometryVector(mono[i].GetCartesianCoordinates(j + 1).x[0] - mono[i].GetCartesianCoordinates(j).x[0], 0.0, 0.0, 2.0)); //matrix phase
						}
					}
					thickness_init.emplace_back(GeometryVector(mono[i].GetBasisVector(0).x[0] - mono[i].GetCartesianCoordinates(mono[i].NumParticle() - 1).x[0], 0.0, 0.0, 2.0));//matrix phase
				}

				SpherePacking c;
				//SpherePacking2Layers(initConfig[i], thickness_init);
				Layer2SpherePacking(thickness_init, c);	initConfig[i] = c;	phi += c.PackingFraction(); phi2 += c.PackingFraction()*c.PackingFraction();
				chirping(thickness_init, thickness_fin);
				Layer2SpherePacking(thickness_fin, c);
				chirpedConfig[i] = c;
				ofile << "\n";
				if (config) {//Config
							 //initial
					sprintf(name, "%s__layered_________PackingFormat_%d", tempstring, i);
					WriteConfiguration(initConfig[i], name);
					sprintf(name, "%s__layered_________ThicknessFormat_%d", tempstring, i);
					WriteFunction(thickness_init, name);

					sprintf(name, "%s__layered_chirped_PackingFormat_%d", tempstring, i);
					WriteConfiguration(chirpedConfig[i], name);
					sprintf(name, "%s__layered_chirped_ThicknessFormat_%d", tempstring, i);
					WriteFunction(thickness_fin, name);
				}
			}

			ofile << "Packing fraction = " << phi / numConfig << "+-" << sqrt(phi2 / numConfig - phi * phi / (numConfig*numConfig)) / sqrt(numConfig) << "\n";
			//Spectral densities
			if (spectraldensity) {
				std::vector<GeometryVector> chiv;
				sprintf(name, "%s__layered_________chiv", tempstring);
				IsotropicSpectralDensity(getConfig1, numConfig, Kmax, Kmax, chiv, dK);
				//WriteFunction(chiv, name);
				PlotFunction_Grace(chiv, name, "k", "\\xc\\0\\sV\\N(k)");
				sprintf(name, "%s__layered_chirped_chiv", tempstring);
				IsotropicSpectralDensity(getConfig2, numConfig, Kmax, Kmax, chiv, dK);
				//WriteFunction(chiv, name);
				PlotFunction_Grace(chiv, name, "k", "\\xc\\0\\sV\\N(k)");
			}
		}

		//------------------------------------Tiling
		{
			ofile << "=========Tiling structures\n";
			double phi0 = 0.5, phi = 0;
			//Generate Configurations!!!!!!!!!!
			//initConfig = std::vector<SpherePacking>(mono);
			//Chirp the packing
			std::vector<GeometryVector> thickness_init, thickness_fin;
			for (int i = 0; i < numConfig; i++) {
				thickness_init.clear();	thickness_fin.clear();
				{// Generate in terms of spherepackings
					SpherePacking ctemp(mono[i]); ctemp.RemoveParticles();
					int idx_left, idx_right;
					for (int j = 0; j < mono[i].NumParticle(); j++) {
						idx_left = j - 1;	idx_right = j + 1;
						GeometryVector x(0.0);
						double r = 0;

						/*
						if (idx_left < 0) {
						idx_left += mono[i].NumParticle();
						x = GeometryVector(-1.0);
						r = sqrt(mono[i].GetBasisVector(0).Modulus2());
						}
						if (idx_right >= mono[i].NumParticle()) {
						idx_right -= mono[i].NumParticle();
						x = GeometryVector(1.0);
						r = sqrt(mono[i].GetBasisVector(0).Modulus2());
						}
						*/
						//x = x + mono[i].GetRelativeCoordinates(idx_right) + mono[i].GetRelativeCoordinates(idx_left);	
						//x = 0.5*x;
						//r = r - sqrt((mono[i].GetCartesianCoordinates(idx_right) - mono[i].GetCartesianCoordinates(idx_left)).Modulus2());
						//r = fabs(r)*phi0*0.5;
						GeometryVector left(0.0), right(0.0);
						if (idx_left < 0) {
							idx_left += mono[i].NumParticle();
							left.x[0] = -1.0;
						}
						left = left + mono[i].GetRelativeCoordinates(idx_left);

						if (idx_right >= mono[i].NumParticle()) {
							idx_right -= mono[i].NumParticle();
							right.x[0] = 1.0;
						}
						right = right + mono[i].GetRelativeCoordinates(idx_right);
						x = 0.25*(left + right + 2.0*mono[i].GetRelativeCoordinates(j));
						r = sqrt(mono[i].RelativeCoord2CartesianCoord(0.25*(right - left)).Modulus2())*phi0;
						ctemp.Insert(r, x);
					}
					Merge_Sort(ctemp, 0);
					initConfig[i] = SpherePacking(ctemp);
				}
				ofile << i << " : ";
				SpherePacking c;
				SpherePacking2Layers(initConfig[i], thickness_init);
				//Layer2SpherePacking(thickness_init, c);	initConfig[i] = c;	
				phi += initConfig[i].PackingFraction();
				chirping(thickness_init, thickness_fin);
				Layer2SpherePacking(thickness_fin, c);
				chirpedConfig[i] = c;
				ofile << "\n";
				if (config) {//Config
							 //initial
					sprintf(name, "%s__tiling_________PackingFormat_%d", tempstring, i);
					WriteConfiguration(initConfig[i], name);
					sprintf(name, "%s__tiling_________ThicknessFormat_%d", tempstring, i);
					WriteFunction(thickness_init, name);

					sprintf(name, "%s__tiling_chirped_PackingFormat_%d", tempstring, i);
					WriteConfiguration(chirpedConfig[i], name);
					sprintf(name, "%s__tiling_chirped_ThicknessFormat_%d", tempstring, i);
					WriteFunction(thickness_fin, name);
				}
			}

			ofile << "Packing fraction = " << phi / numConfig << "\n";
			//Spectral densities
			if (spectraldensity) {
				std::vector<GeometryVector> chiv;
				sprintf(name, "%s__tiling_________chiv", tempstring);
				IsotropicSpectralDensity(getConfig1, numConfig, Kmax, Kmax, chiv, dK);
				//WriteFunction(chiv, name);
				PlotFunction_Grace(chiv, name, "k", "\\xc\\0\\sV\\N(k)");
				sprintf(name, "%s__tiling_chirped_chiv", tempstring);
				IsotropicSpectralDensity(getConfig2, numConfig, Kmax, Kmax, chiv, dK);
				//WriteFunction(chiv, name);
				PlotFunction_Grace(chiv, name, "k", "\\xc\\0\\sV\\N(k)");
			}

		}


		//----------------------------

	}

	return 0;





	{
		int dim = 2;
		std::function<SpherePacking(size_t)> getInitConfig;
		ConfigGen * initC = nullptr;

		ofile << "config type : ";
		ifile >> tempstring;
		if (strcmp(tempstring, "RSA") == 0) {
					int num;
					double packingFraction, numberdensity;
					ofile << "Particle Number = ";
					ifile >> num; ofile << num << std::endl;
					ofile << "Packing Fraction = ";
					ifile >> packingFraction;	ofile << packingFraction << std::endl;
					ofile << "numberdesnity = ";
					ifile >> numberdensity;	ofile << numberdensity << std::endl;

					if (dim == 0) {
						ofile << "Spatial dimension = ";	ifile >> dim;	ofile << dim << "\n";
					}

					initC = new RSA_Gen(dim, num, packingFraction, numberdensity);

					getInitConfig = [&initC](size_t i)->SpherePacking {
						SpherePacking result;	initC->GenerateP(result);
						return result;
					};
				}
		else if (strcmp(tempstring, "ImLat") == 0) {
					int num;
					double packingFraction, numberdensity;
					ofile << "Particle Number = ";
					ifile >> num; ofile << num << std::endl;
					ofile << "Packing Fraction = ";
					ifile >> packingFraction;	ofile << packingFraction << std::endl;
					ofile << "numberdesnity = ";
					ifile >> numberdensity;	ofile << numberdensity << std::endl;

					if (dim == 0) {
						ofile << "Spatial dimension = ";	ifile >> dim;	ofile << dim << "\n";
					}

					initC = new ImLatGen(dim, num, packingFraction, numberdensity);

/*					double phiMax = HyperSphere_Volume(dim, 0.5);
					double p = 1.0 - packingFraction / phiMax;
					int NmaxInEachSide = (int)round(pow(num / (1.0 - p), 1.0 / dim));
					int numRemove = pow(NmaxInEachSide, dim) - num;
					if (dim == 1)
						initC = new LatticeGen("Integer", NmaxInEachSide);
					else if (dim == 2)
						initC = new LatticeGen("Square", NmaxInEachSide);
					else if (dim == 3)
						initC = new LatticeGen("Cubic", NmaxInEachSide);
					else {
						ofile << "Dimension is too high\n";
					}*/
					getInitConfig = [&initC](size_t i)->SpherePacking {
						SpherePacking c;	initC->GenerateP(c);
						return c;
					};
				}
		else if (strcmp(tempstring, "HSF") == 0) {
					int num, samplingInterval = 20;
					double packingFraction, numberdensity;
					ofile << "Particle Number = ";
					ifile >> num; ofile << num << std::endl;
					ofile << "Packing Fraction = ";
					ifile >> packingFraction;	ofile << packingFraction << std::endl;
					ofile << "numberdesnity = ";
					ifile >> numberdensity;	ofile << numberdensity << std::endl;

					if (dim == 0) {
						ofile << "Spatial dimension = ";	ifile >> dim;	ofile << dim << "\n";
					}

					initC = new HSF_Gen(dim, num, packingFraction, numberdensity);
					getInitConfig = [&initC](size_t i)->SpherePacking {
						SpherePacking cc;	initC->GenerateP(cc);
						return cc;
					};
				}
	
	
		std::vector<GeometryVector> result;
		IsotropicSpectralDensity(getInitConfig, 100, 10, 10, result, .1);
		WriteFunction(result, "test");
	}

	return 0;
	{
		std::vector<ConfigGen *>test;
		ofile << "System type? ";
		ifile >> tempstring;
		if (strcmp(tempstring, "RSA") == 0) {
			int d, numPrts;
			double numdensity, packingfraction;

			ofile << "dimension = ";
			ifile >> d;
			ofile << "num particles = ";
			ifile >> numPrts;
			ofile << "num density = ";
			ifile >> numdensity;
			ofile << "packing fraction = ";
			ifile >> packingfraction;
			test.emplace_back(new RSA_Gen(d, numPrts, packingfraction, numdensity));
		}
		else if (strcmp(tempstring, "HSF") == 0) {
			int d, numPrts, samplinginterval;
			double numdensity, packingfraction;
			ofile << "dimension = ";
			ifile >> d;
			ofile << "num particles = ";
			ifile >> numPrts;
			ofile << "num density = ";
			ifile >> numdensity;
			ofile << "packing fraction = ";
			ifile >> packingfraction;

			test.emplace_back(new HSF_Gen(d, numPrts, packingfraction, numdensity));
		}
		else if (strcmp(tempstring, "Lattice") == 0) {
			char latticeName[30];
			ofile << "lattice name? ";
			ifile >> latticeName;
			test.emplace_back(new LatticeGen(latticeName, 10));
		}

		SpherePacking cc;
		test[0]->GenerateP(cc);

		WriteConfiguration(cc, "./test");
		PlotPacking2D("./test", cc);
	}

	return 0;

	//{
	//	int numPrts, numConfigs;	double phi;
	//	ifile >> numPrts;
	//	ifile >> numConfigs;
	//	ifile >> phi;
	//	sprintf(name, "./configs/RSA_1D_%0.1f_%d", phi, numPrts);
	//	ofile << name << std::endl;
	//	ConfigurationPack save(name);

	//	double radius = phi / numPrts*0.5;
	//	
	//	for (int i = 0; i < numConfigs; i++) {
	//		Configuration temp = GetUnitCubicBox(1);
	//		SpherePacking c(temp, radius);
	//		RSA ens(c, radius);
	//		ens.RSA_I(numPrts);
	//		temp = Configuration(*ens.pSpheres, "a");
	//		temp.Rescale(numPrts);
	//		save.AddConfig(temp);
	//		ofile << i << " : " << ens.pSpheres->PackingFraction() << "\n";
	//	}
	//}

	//return 0;


	
	{
		SpherePacking c;
		ReadConfiguration(c, "testConfig");		SpherePacking c_(c); c_.RemoveParticles();
		std::vector<GeometryVector> newbasis;
		for (int i = 0; i < 2; i++)
			newbasis.push_back(c.GetBasisVector(i));
		newbasis.emplace_back(GeometryVector(0.0,0.0,1.0));
		container_leak see(c, newbasis, false, 14, 14, 1);
		
		ofile << "prt ID\tnum vertx\tnumnormal\n";
		voro::voronoicell_neighbor vc;
		voro::c_loop_all_periodic routine(see);
		int idx = 0, count = 0;
		double cellVolume, r;	std::vector<int> vs;	GeometryVector n; std::vector<double> rlist;
		std::vector<double> vertices;	std::vector<GeometryVector> centroids(c.NumParticle(), GeometryVector(0.0,0.0,0.0));
		GeometryVector av(0, 0, 0);
		if (routine.start()) {
			do {
				idx = routine.pid();	ofile << idx << "\t";
				if (idx > -1) {/*Particles are in the inner region*/
					if (see.compute_cell(vc, routine)) {
						vc.vertices(vertices);
						int numV = vertices.size()/3;
						numV /= 2;
						ofile << numV << "\t";
						vc.neighbors(vs);
						//ofile << vs.size() << "\n";
						cellVolume = vc.volume();
						rlist.clear();
						GeometryVector  x0 = c.GetCartesianCoordinates(idx);
						for (int i = 0; i < vs.size(); i++) {
							n = c.distance(vs[i], idx);
								
							if (fabs(n.x[2]) < 0.001 && n.Modulus2() > 0.01) {
										rlist.push_back(sqrt(n.Modulus2())/2.0);
							}
						}
						r = rlist[0];
						for (int i = 1; i < rlist.size(); i++)
							r = std::min(r, rlist[i]);
						ofile << rlist.size() <<"\n";
						if (idx < 20)
							ofile << "\t"<< r << "\n";

						c_.Insert(r, c.GetRelativeCoordinates(idx));
						GeometryVector centr(0.0, 0.0, 0.0);
						vc.centroid(centr.x[0], centr.x[1], centr.x[2]);
						centroids[idx] = centr;
						av = av + cellVolume * centr;
					}
					else
						count++;
				}
			} while (routine.inc());

		}
		ofile << "fail : " << count << std::endl;
		ofile << "sum of dr" << av;
		WriteConfiguration(c_, "testConfig_minmax");
		WriteFunction(centroids, "testResult");
		PlotVoronoi2D("testConfig_minmax", c_);
	}

	return 0;


	{
		Verbosity = 4;
		double phi = 0;	int N;
		ofile << "N = ";	ifile >> N;
		ofile << "phi = "; ifile >> phi;
		HSF test(2, N, phi);	test.SetNumThreads(1);
		test.RescaleConfig(1.0);
		test.Prepare();
		WriteConfiguration(* test.pSpheres, "testConfig");

	}


	return 0;

	

	{
		for (;;) {
			ofile << "What's next?\n";
			ifile >> tempstring;

			if (strcmp(tempstring, "Exit") == 0)
				break;
			else {
				int numConfig;
				ofile << "LoadPrefix: ";
				ifile >> name;
				ConfigurationPack load(name);
				numConfig = load.NumConfig();
				ofile << "There are " << numConfig << " configurations \n";
				ofile << "SavePrefix: ";
				ifile >> tempstring;
				for (int i = 0; i < numConfig; i++) {
					sprintf(name, tempstring, i);
					WriteConfiguration(load.GetConfig(i), name);
				}
				ofile << "Done!\n";
			}
		
		}
	
	}


	return 0;



	{
		int numThreads = 6;
		double maxPhi = 0.9999;
		double vmax = 0.001;
		int d = 2;
		int numConfig = 1;
		char dir[200] = "C:/Users/phy00/Dropbox/Document/progress/Hyperuniform construction/CoatedSpheres/";

		std::vector<double> p;	p.push_back(0);	p.push_back(0); p[0] = vmax;
		char fnName[10] = "Power";
		int ni = 1;
		int nf = 100;

		{//INPUT PART
			ofile << "Dimension = ";
			ifile >> d;
			ofile << "max vol = ";
			ifile >> p[0];
			ofile << "parameter = ";
			ifile >> p[1];
			ofile << "function name (Power/Zeta)= ";
			ifile >> fnName;
			ofile << "num steps_init = ";
			ifile >> ni;
			ofile << "num steps_final = ";
			ifile >> nf;
			ofile << "num configs = ";
			ifile >> numConfig;
			ofile << "num threads = ";
			ifile >> numThreads;
		}
		//Generate Unit Simple Cubic
		std::vector<GeometryVector> b;
		for (int i = 0; i < d; i++) {
			b.emplace_back(GeometryVector());
			b[i].Dimension = d;
			b[i].x[i] = 1.0;
		}
		SpherePacking c(d, &b[0], 0.01);

		omp_set_num_threads(numThreads);
		std::vector<GeometryVector> phis;
//#pragma omp parallel for schedule (guided)
		for (int i = 0; i < numConfig; i++) {
			SpherePacking temp(c);
			int tick = 0;
			char name[300] = {};
			sprintf(name, "%sconfigs/%dD_%s_%0.2f_config%d", dir, d, fnName, p[1], i);
			if (ni > 1) {
				char load[300] = {};
				sprintf(load, "%s__%d", name, ni - 1);
				ReadConfiguration(temp, load);
				if (temp.NumParticle()>0)
//#pragma omp critical
				{
					ofile << i << "th config is successfully loaded\n";
				}
			}

			double error = 0.0;
			error = getSpacefillingPacking2(temp, ni, nf, p[0], p[1], fnName, name, numThreads,std::time(NULL) + i);

//#pragma omp critical
			{
				ofile << i << ": covering fraction " << temp.PackingFraction() << "\n";
				ofile << "\terror: " << error << std::endl;
			}

		}

	}
	return 0;


	/*
	{
		{SpherePackingPack Spack("testList");
			for (int i = 0; i < 10; i++) {
				SpherePacking a = GenerateRSAPacking(3, 1000, false, 400);
				Spack.AddConfig(a);
				sprintf(name, "%d", i);
				WriteConfiguration(a, name);
			}
		}
		SpherePackingPack load("testList");
		for (int i = 0; i < load.NumConfig(); i++) {
			SpherePacking a(load.GetConfig(i));
			sprintf(name, "%d__", i);
			WriteConfiguration(a, name);
		}


	
	}
	return 0;
	*/
	
	{
		SpherePacking c_init, c_fin;
		SpherePacking c;
		char name[300]; //sprintf(name, "C:/Users/phy00/Dropbox/Codes/Project2/Project2/results/ImperfectLattices/New folder/sample_%d", i);
		//sprintf(name, "C:/Users/phy00/Dropbox/Codes/Project2/Project2/results/ImperfectLattices/2dImLat_config_0.250___N_1000000_0000___test");
		//ReadConfiguration(c_init, name);
		
		c_init = SpherePacking(GetUnitCubicBox(2),0);
		
		omp_set_num_threads(4);
#pragma omp parallel
		{
			RandomGenerator g(omp_get_thread_num());
#pragma omp for
			for (int i = 0; i < 1000000; i++) {
				GeometryVector x(g.RandomDouble(), g.RandomDouble());
#pragma omp critical
				{
					c_init.Insert(0, x);
				}
			}

		}
		
		bool doCrops = false;
		//ifile >> doCrops;

		if (doCrops) {
			int d = c_init.GetDimension();
			std::vector<GeometryVector> basis(d, GeometryVector(static_cast<int> (d)));
			double L = 0;
			ifile >> L;

			GeometryVector x0(static_cast<int> (d));
			for (int i = 0; i < d; i++)
				x0.x[i] = 0.0;
			std::vector<size_t> xx = Crop(c_init, [d](GeometryVector x)->double {
				double dist = fabs(x.x[0]);
				for (int i = 1; i < d; i++)
					dist = std::max(dist, fabs(x.x[i]));
				return dist; }, x0, L);

			for (int i = 0; i < d; i++)
				basis[i].x[i] = L;

			c = SpherePacking(d, &basis[0], 1.0);
			for (int i = 0; i < xx.size(); i++)
				c.Insert(c_init.GetCharacteristics(0), c.RelativeCoord2CartesianCoord(c_init.GetCartesianCoordinates(i)));

		}
		else {
			c = c_init;
		}

		ofile << "Num Prts = " << c.NumParticle() << std::endl;

		std::vector<double> vols;
		std::vector<int> numVer;
		MakeMono2Poly_parallel(c, c_fin, 0.01, vols, 4, &numVer);
		{
			std::vector<GeometryVector> result(c.NumParticle(), GeometryVector(0.0, 0.0, 0.0, 0.0));
			for (int i = 0; i < result.size(); i++) {
				result[i].x[0] = i;
				result[i].x[1] = vols[i];
				result[i].x[2] = numVer[i];
			}
			WriteFunction(result, "C:/Users/phy00/Dropbox/Codes/Project2/Project2/results/ImperfectLattices/test_100000");
		}
	}
	return 0;


	{

	}

	return 0;



	
	/*
	{
		int numThreads = 6;
		double maxPhi = 0.9999;
		double vmax = 0.001;
		int d = 2;
		int numConfig = 1;
		char dir[200] = "C:/Users/phy00/Dropbox/Document/progress/Hyperuniform construction/CoatedSpheres/";
		std::vector<GeometryVector> b;
		b.emplace_back(GeometryVector(1.0, 0.0));
		b.emplace_back(GeometryVector(0.0, 1.0));
		SpherePacking c(d, &b[0], 0.01);
		std::vector<double> p;	p.push_back(0);	p.push_back(0);
		p[0] = vmax;
		char fnName[10];
		omp_set_num_threads(numThreads);
		

		ofile << "parameter = ";
		ifile >> p[1];
		ofile << "function name (Power/Zeta)= ";
		ifile >> fnName;
		SpacefillingPacking test(c, p, fnName);

		double factor = HyperSphere_Volume(test.dimension, test.maxR)*test.N;
		double sum = 0.0;
		for (int i = 1; i < 500; i++) {
			sum += test.f(i);
			ofile << i << "\t" << factor * sum << "\n";
		}


		

	}
	return 0;*/






	/*
	{
		int numThreads = 6;
		double maxPhi = 0.9999;
		double vmax = 0.01;
		int d = 2;
		int numConfig = 10;
		char dir[200] = "C:/Users/phy00/Dropbox/Document/progress/Hyperuniform construction/CoatedSpheres/";
		std::vector<GeometryVector> b;
		b.emplace_back(GeometryVector(1.0, 0.0));
		b.emplace_back(GeometryVector(0.0, 1.0));
		SpherePacking c(d, &b[0], 0.01);
		std::vector<double> p;	p.push_back(0);	p.push_back(0);
		p[0] = vmax;
		char fnName[10] = "Zeta";
		omp_set_num_threads(numThreads);
		std::vector<int> index;
		index.push_back(1);
		index.push_back(5);
		index.push_back(10);
		index.push_back(15);
		index.push_back(20);
		index.push_back(25);
		index.push_back(30);
		//for (int i = 0; i <= 3; i++)
		//	index.push_back(i * 100);
		int nn = 500;

		for (double p1 = 1.3; p1 < 1.31; p1 += 0.1) {

			
			p[1] = p1;

			std::vector<GeometryVector> phis;
#pragma omp parallel for schedule (guided)
			for (int i = 0; i < numConfig; i++) {
				SpherePacking temp(c);
				int tick = 0;
				char name[300] = {};
				sprintf(name, "%sconfigs/%dD_%s_%0.2f_config%d", dir, d, fnName, p1, i);
				//sprintf(name, "%sconfigs/%dD_%s%d", dir, d, "Zeta", i );
				double error = 0.0;
				error = getSpacefillingPacking2(temp, nn, p[0], p[1], fnName, name);
				sprintf(name, "%sconfigs/%dD_%s_%0.2f_config%d__%d", dir, d, fnName, p1, i,nn);
				ReadConfiguration(temp, name);
				
				SpacefillingPacking tt(temp);
				tt.PlotConfig2D(name, 0.5, "b", "r");
				sprintf(name, "%sconfigs/%dD_%s_%0.2f_config%d__%d_test", dir, d, fnName, p1, i, nn);
				tt.PlotConfig2D(name);


			#pragma omp critical
				{
					ofile << i << ": covering fraction " << tt.pSpheres->PackingFraction() << std::endl;
				}
				//}
			}




		}

	}
	return 0;
	*/
	

	return 0;
	{
		int d = 2;
		std::vector<GeometryVector> b;
		b.emplace_back(GeometryVector(1.0, 0.0));
		b.emplace_back(GeometryVector(0.0, 1.0));
		SpherePacking c(d, &b[0], 0.01);
		std::vector<double> p;	p.push_back(0);	p.push_back(0);
		//p.push_back(1.0);	p.push_back(1.2);
		ofile << "max particle volume: ";
		ifile >> p[0];
		ofile << "parameter: ";
		ifile >> p[1];
		ofile << "rescaling function: ";
		ifile >> name;

		
		SpacefillingPacking test(c, p, name);
		
		ofile << "max prt vol = " << HyperSphere_Volume(d, test.maxR) << std::endl;
		ofile << "num of same prt = " << test.N << std::endl;

		ofile << "i \t phi\n";
		for (int i = 1; i < 60; i++) {
			double r = test.maxR*pow(test.f(i),1.0/(double) d);
			test.AddSpheres(r, test.N);
			sprintf(name, "./results/%d", i);
			//test.PlotConfig2D(name);

			ofile << i << "\t" << test.getPackingFraction() << "\n";
			ofile << "\t" << test.pSpheres ->NumParticle() << "\n";
		}

	}
	return 0;
	
	{
		int numConfigs = 200;
		double K = 4.0*pi;
		double Kmin = 4.0*pi / 100.0;
		int numThreads = 6;
		double dw = 5.0;
		for (double fraction = 0.25; fraction <= 1.00; fraction += 1.25) {
			std::vector<PeriodicCellList<std::complex<double>>> cList;
			std::vector<GeometryVector> result;
			//Load data!
#pragma omp parallel for schedule (guided) num_threads(numThreads)
			for (int i = 0; i < numConfigs; i++) {
				SpherePacking temp;
				std::vector<GeometryVector> numVertices;
				sprintf(name, "C:/Users/phy00/Dropbox/Document/progress/Hyperuniform construction/RSA_configs/2D/configs/2dRSA_config_%0.3f___N_10000_%04d___", fraction, i);
				ReadConfiguration(temp, name);
				{
					//sprintf(name, "C:/Users/phy00/Dropbox/Document/progress/Hyperuniform construction/RSA_configs/2D/configs/2dRSA_config_%0.3f___N_10000_numVert_%04d___", fraction, i);
					//ReadFunction(numVertices, name, 1);
					//for (int j = 0; j < numVertices.size(); j++)
					//	numVertices[j].x[0] -= 6.0;
				}
				RandomGenerator g(i + std::time(NULL));
				numVertices.reserve(temp.NumParticle());
				for (int j = 0; j < temp.NumParticle(); j++) {
					numVertices.emplace_back(GeometryVector(2.0*dw*(g.RandomDouble()-0.5),0.0));
				}

				if (temp.NumParticle() == numVertices.size()) {
					PeriodicCellList<std::complex<double>> c(temp, std::complex<double>(0.0, 0.0));
					c.RemoveParticles();
					for (int j = 0; j < temp.NumParticle(); j++) {
						c.Insert(std::complex<double>(numVertices[j].x[0], 0.0), temp.GetRelativeCoordinates(j));
					}
#pragma omp critical
					{
						cList.push_back(c);
					}
				}
				else {
					ofile << i<< ":\tmismatch: the number of particles\n";
				}
			}
			omp_set_num_threads(numThreads);
		
			auto getC = [&cList](size_t i)->PeriodicCellList<std::complex<double>> {return cList[i]; };	
			IsotropicStructureFactor_weighted(getC, numConfigs, K, K, result, Kmin);
			sprintf(name, "C:/Users/phy00/Dropbox/Document/progress/Hyperuniform construction/RSA_configs/2D/StructureFactor_weighted/Sk_wieghted_N_10000_f_%0.3f_dw_%0.1f", fraction, dw);
			WriteFunction(result, name);
		}
	
	}
	//PlotPacking2D

	
	return 0;

	{
		char conti = 'y';
		for(;conti!='n';){
		double u = 0;
		ifile >> u;
		std::vector<GeometryVector> basis;
		basis.emplace_back(GeometryVector(1.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 1.0));
		std::vector<GeometryVector> v;
		Configuration c = GetLattice_config(basis, 20);
		RandomGenerator g;
		for (int i = 0; i < c.NumParticle(); i++) {
			GeometryVector dx(g.RandomDouble(), g.RandomDouble());
			dx = u * dx;
			c.MoveParticle(i, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(i) + dx));
		}
		for (int i = 0; i < c.NumParticle(); i++) {
		//	ofile << i << " prts";
			GetVoronoiCell(c, i, v);
			PolygonVolume(v, 2);
			//ofile << "\t" << v.size() << "vertices\t"<<PolygonVolume(v, 2) <<"\n";
			ofile << "*";
		}
		//std::vector<GeometryVector> result;
		//{
		//	VoronoiNumSidesDistrubution(getC, 1, result);
		//	PlotFunction_Grace(result, "./results/termite/voronoi_test");
		//}
		ofile << "\ncontinue (y/n)?	";
		ifile >> conti;
		}
	}
	return 0;
	{
		ConfigurationPack pack("configs/termite_model_1");
		Configuration c(pack.GetConfig(14));

		std::vector<GeometryVector> v;
		for (int i = 0; i < c.NumParticle(); i++) {
			GetVoronoiCell(c, i, v);
			ofile << i << "th:\t" << v.size() << "vertices\n";
		}
	
	}
	return 0;
	{
		RandomGenerator g(std::time(NULL));
		int numParticles = 100;
		int numConfig = 100;
		int d = 2;
		
		double tolerance = 1e-20;

		double X = 5.0, Y = 5.0, Z = 1.0;
		for (int i = 0; i < numConfig; i++) {
		std::vector<GeometryVector> b;
		b.emplace_back(GeometryVector(X, 0.0));
		b.emplace_back(GeometryVector(0.0, Y));
		int n = (int)round(sqrt(numParticles));

		int index = 0; GeometryVector xi(static_cast<int> (d));
		
			std::cout << "*";
			voro::container_periodic unitCell(X, 0, Y, 0, 0, Z, n, n, 1, 8);

			//Input random particle coordinates
			double x, y;
			for (int j = 0; j < numParticles; j++) {
				x = X * g.RandomDouble();
				y = Y * g.RandomDouble();
				unitCell.put(j, x, y, 0.5);
			}
			//unitCell.compute_all_cells();

			//std::cout << "compute_all_cells: good\n";
			voro::voronoicell vc;
			voro::c_loop_all_periodic routine(unitCell);
			if (routine.start()) {
				do if (unitCell.compute_cell(vc, routine)) {
					x = routine.x();
					y = routine.y();
					
					std::vector<double> v;
					vc.vertices(v);

					for (int k = 0; k < vc.p; k++) {
						if (v[3 * k + 2] > 0.1) {
							std::cout << "(" << x+v[3*k] << "," << y+v[3*k+1] << ")\n";
						}
					}
				}while (routine.inc());
			}

		}
	}

	return 0;

	{
		int numParticles = 100;
		int numConfig = 100;
		DimensionType d = 2;

		double tolerance = 1e-20;
		GeometryVector * b = new GeometryVector[d];
		b[0] = GeometryVector(5.0, 0.0);
		b[1] = GeometryVector(0.0, 5.0);

		RandomGenerator g;

		std::vector<Configuration> cList;

		for(int i=0; i<numConfig; i++){
			Configuration c(d, b, 1.0);
			for (int j = 0; j < numParticles; j++) {
				c.Insert("a", g);
			}
			cList.push_back(c);
		}

		{//Voro++ part!!
			auto ti = std::time(NULL);
			for (int i = 0; i < numConfig; i++) {
				ofile << i << ":\n";
				Configuration c_;
				Configuration c(cList[i]);
				bool change = transformBasis(c, c_);
				int n = (int)ceil(pow(c.NumParticle(), 1.0 / (double)d));
				int nz = (d == 2) ? 1 : n;
				double z1, z2, z3;
				if (d == 2) {
					z1 = 0;	z2 = 0;	z3 = 1.0;
				}
				else {
					z1 = b[2].x[0];	z2 = b[2].x[1];	z3 = b[2].x[2];
				}
				/*
				double x_min = std::min(0.0, b[1].x[0]);
				double x_max = std::max(b[0].x[0], b[1].x[0]);
				double y_min = std::min(0.0, b[1].x[1]);
				double y_max = std::max(0.0, b[1].x[1]);
				double z_min = std::min(0.0, z3);
				double z_max = std::max(0.0, z3);*/
				
				//voro::pre_container unitCell(x_min, x_max, y_min, y_max, z_min, z_max, false, false, false);
				ofile << b[0].x[0] << ",\t" << b[1].x[0] << ",\t" << b[1].x[1] << ",\t" << z1 << ",\t" << z2 << ",\t" << z3 << ",\t" << n << ",\t" << n << ",\t" << nz << ",\t" << std::endl;
				voro::container_periodic unitCell(b[0].x[0], b[1].x[0], b[1].x[1], z1, z2, z3, n, n, nz, 8);
				//INPUT particle coordiantes
				GeometryVector vi(static_cast<int>(d));
				z3 = 0.5;
				for (int j = 0; j < c.NumParticle(); j++) {
					if (change) {
						vi = c_.RelativeCoord2CartesianCoord(c.GetRelativeCoordinates(j));
					}
					else {
						vi = c.GetCartesianCoordinates(j);
					}
					z1 = vi.x[0];	z2 = vi.x[1];
					if (d == 3)
						z3 = vi.x[2];

					unitCell.put(j, z1, z2, z3);
				}
				int nx=0, ny=0;
				//unitCell.guess_optimal(nx, ny, nz);

				//voro::container_periodic box(b[0].x[0], b[1].x[0], b[1].x[1], z1, z2, z3, nx, ny, nz, 8);
				//unitCell.setup(box);
				//voro::container_

				
				voro::voronoicell vc;
				voro::c_loop_all_periodic routine(unitCell);
				unitCell.compute_all_cells();
					int index = 0; double x = 0, y = 0, dump1 = 0, dump2 = 0;	GeometryVector xi(static_cast<int> (d));
					//if (routine.start()) do if (unitCell.compute_cell(vc, routine)) {
					//	//routine.pos(index, x, y, dump1, dump2);
					//	index = routine.pid();
					//	//x = routine.x();
					//	//y = routine.y();
					//	if (change)
					//		xi = c.GetCartesianCoordinates(index);
					//	else {
					//		xi.x[0] = routine.x();	xi.x[1] = routine.y();
					//		if (d == 3)
					//			xi.x[2] = routine.z();
					//	}

					//	std::vector<double> vList;
					//	
					//	vc.vertices(vList);
					//	for (int j = 0; j < vc.p; j++) {
					//		if (d == 2 && vList[3 * j + 2]>0.1) {
					//			vi.x[0] = vList[3 * j];	vi.x[1] = vList[3 * j + 1]; //vi.x[2] = vList[3 * j + 2];
					//		}
					//		else if (d == 3) {
					//			vi.x[0] = vList[3 * j];	vi.x[1] = vList[3 * j + 1]; vi.x[2] = vList[3 * j + 2];
					//		}

					//		if (change)
					//			vi = c.RelativeCoord2CartesianCoord(c_.CartesianCoord2RelativeCoord(vi));

					//		vi = vi + xi;
					//	}
					//	vc.volume();
					//	
					//}while (routine.inc());
			}

			ofile << "\nVORO++" << std::time(NULL) - ti << "   seconds\n\n\n";
		}
	
		{//CGAL
			auto ti = std::time(NULL);
			for (int i = 0; i < numConfig; i++) {
				ofile << i << ", ";
				std::vector<GeometryVector> v;
				for (int j = 0; j < cList[i].NumParticle(); j++) {
					GetVoronoiCell(cList[i], j, v);
					PolygonVolume(v, d);
				}
			}

			ofile << "VORO++" << std::time(NULL) - ti << "   seconds\n\n\n";
		}
	}
	return 0;

	{
		RandomGenerator g(std::time(NULL));
		int numParticles = 100;
		double tolerance = 1e-20;
		DimensionType d = 2;
		//ofile << "Dimension = ";
		//ifile >> d;
		GeometryVector * b = new GeometryVector[d];
		b[0] = GeometryVector(5.0, 0.0);
		b[1] = GeometryVector(0.0, 5.0);
	

		Configuration c(d, b, 1.0);
		for (int i = 0; i < numParticles; i++) {
			c.Insert("a", g);
		}

		//for (int i = 20; i < 30; i++) {
		//	ofile << "\n\n\n" << i << "th particles!!\n";
		//	std::vector<GeometryVector> vertices;
		//	GetVoronoiCell__(c, i, vertices, nullptr, nullptr);
		//	ofile << "VORO++:\n";
		//	for (int j = 0; j < vertices.size(); j++)
		//		ofile << vertices[j];

		//	ofile << "\n\nCGAL:\n";
		//	GetVoronoiCell(c, i, vertices);
		//	for (int j = 0; j < vertices.size(); j++)
		//		ofile << vertices[j];
		//}



		if (d == 2) {
			double bx = sqrt(b[0].Modulus2());
			double bxy = b[0].Dot(b[1]) / sqrt(b[0].Modulus2()), byy = sqrt(b[1].Modulus2() - bxy*bxy);
			double z_x= 0.0, z_y = 0.0, z_z = 1.0; //strip...
			int n = (int)round(sqrt(c.NumParticle()));
			GeometryVector * b_;
			bool Change = (std::fabs(bx - b[0].x[0]) > tolerance) && (std::fabs(bxy - b[1].x[0]) > tolerance) && (std::fabs(byy - b[1].x[1]) > tolerance);
			Configuration forTrans;
			if (Change) {
				ofile << "need transformation !!\n";
				b_ = new GeometryVector[d];
				b_[0].Dimension = d;
				b_[0].x[0] = bx; b_[0].x[1] = 0.0;
				b_[1].Dimension = d;
				b_[1].x[0] = bxy; b_[1].x[1] = byy;
				forTrans = Configuration(d, b_, 1.0);
			}

			
			{//debug
				std::cout << "|b0|^2 = " << sqrt(bx * bx) << "\n";
				std::cout << "|b1|^2 = " << sqrt(bxy * bxy + byy * byy) << "\n";
			}
			
			voro::container_periodic unitCell(bx, bxy, byy, z_x, z_y,z_z,n ,n ,1,8 );
			GeometryVector temp;
			//Input particle Cartesian coordinates!!
			for (int i = 0; i < c.NumParticle(); i++) {
				if (Change) {
					temp = forTrans.RelativeCoord2CartesianCoord(c.GetRelativeCoordinates(i));
				}
				else {
					temp = c.GetCartesianCoordinates(i);
				}
				unitCell.put(i, temp.x[0], temp.x[1], 0.5);
			}
			//unitCell.compute_all_cells();
			voro::c_loop_all_periodic routine(unitCell);
			voro::voronoicell vc;
			int index = 0; double x = 0, y = 0, dump1 = 0, dump2 = 0;	GeometryVector xx(static_cast<int> (d));
			if (routine.start()) do if (unitCell.compute_cell(vc, routine)) {
				//routine.pos(index, x, y, dump1, dump2);
				index = routine.pid();
				x = routine.x();
				y = routine.y();

				xx.x[0] = x;	xx.x[1] = y;
				if (Change)
					xx = c.RelativeCoord2CartesianCoord(forTrans.CartesianCoord2RelativeCoord(xx));

				{//Debug
					std::cout << "particle index " << index << " at " << xx;
					std::cout << "\t # vertices= " << vc.p/2 << "\n";
					std::cout << "\t cell volume= " << vc.volume() <<"\n";
					{
						std::vector<double> vList;
						GeometryVector vi(static_cast<int>(d));
						vc.vertices(vList);
						for (int j = 0; j < vc.p; j++) {
							if (vList[3 * j + 2]>tolerance) { //on the upper plane
							vi.x[0] = vList[3 * j];	vi.x[1] = vList[3 * j + 1]; //vi.x[2] = vList[3 * j + 2];
								if (Change)
									vi = c.RelativeCoord2CartesianCoord(forTrans.CartesianCoord2RelativeCoord(vi));

								vi = vi + xx;
								std::cout << "\t" << j << ": " << vi;
							}
						}
					}
				/*	std::cout << "-----\n";
					{
						std::vector<GeometryVector> vv;
						GetVoronoiCell(c, index, vv);
						std::cout << "particle index " << index << " at " << c.GetCartesianCoordinates(index);
						std::cout << "\t # vertices= " << vv.size() << "\n";
						std::cout << "\t cell volume= " << PolygonVolume(vv, d) << "\n";
						for (int j = 0; j < vv.size(); j++) {
							std::cout << "\t" << j << ": " << vv[j];
						}
					}*/
					std::cout << "\n\n\n";
				}
			}while (routine.inc());
		}
	}


	return 0;
	/*{
		SpherePacking c;
		ReadConfiguration(c, "./configs/3dRSA_config_1.000___N_100_0008___");
		Configuration cc(c, "a");
		DisplayConfiguration_3D_WithVoronoi(cc,0.1);

		return 0;
	}*/
	{
		ConfigurationPack file("./configs/Stealthy_2D_0.1_2500");
		for (int i = 0; i < 3; i++) {
			Configuration c = file.GetConfig(i);
			sprintf(name, "%d", i);
			//PlotVoronoi2D(name, c);
			Plot(name, c,"",2.0);
		}

	//	int n = file.NumConfig();
	//	ofile << "There are " << n << " configurations\n";
	//	auto get = [&file](size_t i)->Configuration {return file.GetConfig(i); };

	//	std::vector<GeometryVector> result;
	//	IsotropicTwoPairCorrelation(get, n, 20, result,5000,2.0);
	//	PlotFunction_Grace(result, "./results/termite/test_g2");
	//	return 0;
	//}


	//{
	//	std::vector<SpherePacking> list;

		/*for (int i = 0; i < 4; i++) {
			int n = list[i].NumParticle();
			for (int j = 0; j < n; j++) {
				std::vector<GeometryVector> gg;
				ofile << i << "th config: " << j << "th particles\n";
				GetVoronoiCell(list[i], j, gg);
			}
		}*/
		
	}
	return 0;

	{
		char dir[200] = "./results/termite/";
		char cPackPrefix[100] = "termite_model_1";
		char fname[200] = "";
		//ifile >> cPackPrefix;
		
		sprintf(fname, "configs/%s", cPackPrefix);
		ConfigurationPack ens(fname);
		ofile << "There are " << ens.NumConfig() << " configurations\n";
		auto getC = [&ens](size_t i) ->Configuration {return ens.GetConfig(14); };
	
		
		{//Plot some config
			int numC = 3;
			for (int i = 0; i < numC; i++) {
				sprintf(fname,"%sconfig__%s_%d", dir, cPackPrefix,i);
				//Plot(fname, getC(i), fname);
			}
		}

		/*for (int i = 0; i < ens.NumConfig(); i++) {
			Configuration c = getC(i);
			ofile << i << "\n";
		}*/

		/*std::vector<GeometryVector> f;
		Configuration c = getC(0);
		int n = c.NumParticle();
		for (int i = 363; i < n; i++) {
			f.push_back(c.GetCartesianCoordinates(i));
		}
		WriteFunction(f, "./results/termite/test");*/

		

		return 0;

		{//Isotropic g2

			//sprintf(fname, "%sradial_g2__%s_short", dir, cPackPrefix);
			//std::vector<GeometryVector> result;
			//IsotropicTwoPairCorrelation(getC, 100, 10.0, result,500,2.0);
		
			//WriteFunction(result, fname);
			
		/*	for (int i = 0; i < 10; i++) {
				sprintf(fname, "%sradial_g2__%s_%d", dir, cPackPrefix,i);
				auto get = [&ens, i](size_t j)->Configuration {return ens.GetConfig(i); };
			
				std::vector<GeometryVector> result;
				IsotropicTwoPairCorrelation(get, 1, 20.0, result, 500000, 1.0);
				WriteFunction(result, fname);
			}*/
			

		}

		//{//Isotropic Sk
		//	sprintf(fname, "%sSk__%s", dir, cPackPrefix);
		//	double dk = 2.0*pi / 50.0;
		//	double kMax = dk * 70;
		//	std::vector<GeometryVector> result;
		//	IsotropicStructureFactor(getC, ens.NumConfig(), kMax, 2.0*kMax, result, dk);
		//	WriteFunction(result, fname);
		//}

		//{//Directional Structure factor
		//	int length = 100;
		//	std::vector<GeometryVector> result;
		//	std::vector<GeometryVector> grid = DirectionalStructureFactor(getC, ens.NumConfig(), length, result);
		//	sprintf(fname, "%s/Directional_Sk__%s", dir, cPackPrefix);
		//	WriteToContourf(grid, result, fname, 2 * length + 1);
		//}
		return 0;
	}



	{
		char dir1[100] = "C:/Users/phy00/Downloads/Bio_patterns/termite_model/";
		char dir2[100] = "C:/Users/phy00/Downloads/Bio_patterns/termite_model-with_rules_for_more_regularity/";
		char name[200] = "";// "termite_model_snap_centers_00001_0108001";
		
		ConfigurationPack file("configs/termite_model_1");
		file.Clear();

		std::vector<GeometryVector> basis;
		basis.emplace_back(GeometryVector(50.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 50.0));
		Configuration temp(2, &basis[0], 1.0);

		for (int i = 1; i <= 100; i++) {
			std::vector<GeometryVector> mem;
			sprintf(name, "%stermite_model_snap_centers_%05d_0108001", dir1, i);
			//ofile << name << "\n";
			std::fstream ifile(name, std::fstream::in);
			ReadFunction(mem, ifile);
			Configuration c(temp);
			GeometryVector x;
			ofile << "i = " << i << "\tnumPrts = " << mem.size() << "\n";
			
			for (int j = 0; j < mem.size(); j++) {
				GeometryVector y(mem[j].x[0],mem[j].x[1]);
				x = c.CartesianCoord2RelativeCoord(y);
				c.Insert("a", x);
			}
			file.AddConfig(c);
		}
	
		ConfigurationPack file2("configs/termite_model_2");
		file2.Clear();

		for (int i = 1; i <= 100; i++) {
			std::vector<GeometryVector> mem;
			sprintf(name, "%stermite_model_snap_centers_%05d_0108001", dir2, i);
			std::fstream ifile(name, std::fstream::in);
			ReadFunction(mem, ifile);
			Configuration c(temp);
			GeometryVector x;
			ofile << "i = " << i << "\tnumPrts = " << mem.size() << "\n";
			for (int j = 0; j < mem.size(); j++) {
				GeometryVector y(mem[j].x[0], mem[j].x[1]);
				x = c.CartesianCoord2RelativeCoord(y);
				c.Insert("a", x);
			}
			file2.AddConfig(c);
		}
		
		return 0;
	}

	{
		//VoronoiBased(true);
		//return 0;
	}
	{
		CCInversionProblem1(true);
			
		return 0;
	}
	{
		char dir[100];
		char dir0[100] = "results/Reconstructions/";
		std::stringstream stream;
		char temp[200] = {};
		char name[200] = {};

		//------------------
		DimensionType dim = 1;
		size_t Nc = 100;	//number of trial configs
		int numConfigs = 1000; //number of target configurations
		size_t NKSkip = 0;
		char constType[100] = {};
		int numThreads = 6;

		//int numBins = 40;
		int N1 = 10;//min system size
		int N2 = 40;//max system size
		int systType1;
		int systType2;
		int dN = 1;
		int methodNum1 = 0;
		int methodNum2 = 0;
		double EnergyTolerance = 1e-20;
		double PositionTolerance = 1e-6;
		double delta1;
		double delta2;
		double d_delta;

		std::cout << "NKSkip=";
		std::cin >> NKSkip;
		std::cout << "dim=";
		std::cin >> dim;
		std::cout << "min N=";
		std::cin >> N1;
		std::cout << "max N=";
		std::cin >> N2;
		std::cout << "dN=";
		std::cin >> dN;
		std::cout << "numConfigs=";
		std::cin >> numConfigs;
		std::cout << "Nc=";
		std::cin >> Nc;
		//std::cout << "numBins=";
		//std::cin >> numBins;
		std::cout << "system type1:";
		std::cin >> systType1;
		std::cout << "system type2:";
		std::cin >> systType2;
		std::cout << "initial method num: ";
		std::cin >> methodNum1;
		std::cout << "final method num: ";
		std::cin >> methodNum2;
		ofile << "Energ Tolerance = ";
		ifile >> EnergyTolerance;
		ofile << "delta1 = ";
		ifile >> delta1;
		ofile << "delta2 = ";
		ifile >> delta2;
		ofile << "d_delta = ";
		ifile >> d_delta;


		//std::cout << "constraints type: ";
		//std::cin >> constType;
		//	double Kmax = 4* pi;

		if (NKSkip > 0)
			std::cout << "Skip the central " << NKSkip << " k points.\n";

		std::vector<GeometryVector> basis;
		for (int j = 0; j < dim;j++) {
			basis.push_back(GeometryVector());
			basis[j].SetDimension(dim);
			basis[j].x[j] = 1.0;
		}

		//double delta = 0.1;
		for (double delta = delta1; delta <= delta2; delta += d_delta) {
			for (int method_num = methodNum1; method_num <= methodNum2; method_num++) {

				char method[100] = {};
				switch (method_num)
				{
				case 0:
					sprintf(method, "BFGS+MINOP"); break;
				case 1:
					sprintf(method, "Steepest Descent"); break;
				case 2:
					sprintf(method, "ConjugateGradient"); break;
				case 3:
					sprintf(method, "LocalGradientDescent"); break;
				default:
					break;
				}

				std::cout << "method = " << method << std::endl;
				for (int type = systType1; type <= systType2; type++) {
					char systName[200] = {};
					switch (type)
					{
					case 0:
						sprintf(systName, "perturbedLattice(1prt)__delta_%0.4f", delta);	break;
					case 1:
						sprintf(systName, "perturbedLattice(Nprt)__delta_%0.4f", delta);	break;
					case 2:
						sprintf(systName, "Poisson");	break;
					case 3:
						sprintf(systName, "lattice_vacancy");	break;
					case 4:
						sprintf(systName, "lattice_interstitial");	break;
					case 5:
						sprintf(systName, "non-Bravais"); break;
					default:
						break;
					}

					std::cout << " system type = " << systName << std::endl;
					for (int N = N1; N <= N2; N += dN) {

						char filename[400] = {};
						sprintf(filename, "%s%dD%s__%s__N_%d.txt", dir, dim, systName, method, N);
						std::fstream outfile(filename, std::fstream::out);
						sprintf(filename, "%s%dD%s__%s_summary.txt", dir, dim, systName, method);
						std::ofstream summary(filename);

						//wave vectors
						//double k1 = 0;
						//std::vector<GeometryVector> kList;
						int N_eff = N;
						if (type == 3) { N_eff--; }
						else if (type == 4) { N_eff++; }

						//for (int i = 0; i < N_eff; i += 2) {
						//	 k1 += 2.0*pi / (double)N;
						//	 kList.push_back(GeometryVector((double)(k1)));
						// }
						//std::cout << "type = " << systName << "\tN = " << N << "\t numConstraints = " << kList.size() << "\n";
						auto getTargetConfig = [&type, &N, dim, &basis, delta](size_t i)->Configuration {
							Configuration c;
							GeometryVector temp;    temp.SetDimension(dim);
							RandomGenerator g(std::time(NULL) + i);
							if (type == 0) {
								c = GetLattice_config(basis, N);
								int j = (int)floor(N*g.RandomDouble());
								RandomVectorInCube(dim, g, temp);
								temp = delta*temp;
								c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
								char name[200] = {};
								//sprintf(name, "CC_Inv_Log_1PL_d_%d_NKSkip_%d___%d", dim, NKSkip,N);
								//WriteConfiguration(c, name);
							}
							else if (type == 1) {
								c = GetLattice_config(basis, N);
								for (int j = 0; j < N; j++) {
									//int j = (int)floor(N*g.RandomDouble());
									RandomVectorInCube(dim, g, temp);
									temp = delta*temp;
									c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
								}
								char name[200] = {};
								//sprintf(name, "CC_Inv_Log_PL_d_%d_NKSkip_%d___%d", dim, NKSkip, N);
								//WriteConfiguration(c, name);
							}
							else if (type == 2) {
								c = GeneratePoisson(dim, N, g);
								c.Resize(N);
							}
							else if (type == 3) {//vacancies
								c = GetLattice_config(basis, N);
								int j = (int)floor(N*g.RandomDouble());
								c.DeleteParticle(j);
							}
							else if (type == 4) {//interstitial
								c = GetLattice_config(basis, N);
								c.Insert("a", g);
							}
							else if (type == 5) {// random non-Bravais lattice
								GeometryVector *A = new GeometryVector[1];
								A[0].Dimension = 1;
								A[0].x[0] = N;
								c = Configuration(dim, A, 1.0);
								int numBasisAtom = Find_LeastCommonDivider(N);
								std::vector<GeometryVector> atomBasis;
								GeometryVector X(static_cast<int> (1));
								RandomGenerator g(std::time(NULL) + i);
								for (int i = 0; i < numBasisAtom; i++) {
									X.x[0] = (double)numBasisAtom*g.RandomDouble();
									atomBasis.push_back(X);
								}
								X.x[0] = 0;
								for (int i = 0; i < N; i++) {
									int j = i%numBasisAtom;
									int k = (int)floor(i / numBasisAtom);

									X.x[0] = k*numBasisAtom;
									c.Insert("a", c.CartesianCoord2RelativeCoord(X + atomBasis[j]));
								}
							}
							return c;

						};


						//---------------init config
						sprintf(filename, "%stempInit/%d", dir0, N);
						Configuration templateConfig(getTargetConfig(0));
						RandomGenerator g(std::time(NULL));
						auto getInitConfig = [&templateConfig, &g](size_t i)->Configuration {
							Configuration dd(templateConfig);
							size_t numP = dd.NumParticle();
							for (size_t z = 0; z < numP; z++) {
								dd.Insert("a", g);
							}
							return dd;
						};
						ConfigurationPack initConfig(filename);
						if (initConfig.NumConfig() == 0) {
							for (int z = 0; z < Nc; z++) {
								initConfig.AddConfig(getInitConfig(z));
							}
						}

						std::vector<double> data;	data.reserve(numConfigs);

						double min = 0;	double max = 0;	int ind = 0;
						std::vector<GeometryVector> result_suc_real;
						std::vector<GeometryVector> result_ground_real;
						std::vector<GeometryVector> result_suc_imag;
						std::vector<GeometryVector> result_ground_imag;

						std::cout << "\nN=" << N << "\n";
						std::cout << "M\tType\tmaxSols\tminSols\tProb_DiffSols\tav.NumSols\tstdSols\tMaxE_ground\tMaxdx\n";
						for (size_t NumConst = 1; NumConst <= N_eff + 3; NumConst++) {
							outfile << "-----------------------------------------------------------------" << std::endl;
							outfile << "Testing N=" << N << "\tM=" << NumConst << std::endl << std::endl;
							outfile << "num\tconstr\tconfigs";
							for (int z = 0; z < N; z++)
								outfile << "\t\t";
							outfile << "\tminE\t\tmaxE\t\tmaxDx" << std::endl;

							result_suc_real.push_back(GeometryVector(NumConst, 0, 0, 0));
							result_ground_real.push_back(GeometryVector(NumConst, 0, 0, 0));
							result_suc_imag.push_back(GeometryVector(NumConst, 0, 0, 0));
							result_ground_imag.push_back(GeometryVector(NumConst, 0, 0, 0));

							size_t maxSols = 0, minSols = Nc, numDiffNumSols = 0;
							int count = 0;
							double maxE_ground1 = 0, maxE_ground2 = 0, maxDx1 = 0, maxDx2 = 0;
							//size_t numSols_1 = 0, numSols_12 = 0, numSols_2=0, numSols_22=0;
							GeometryVector numSols(0, 0, 0, 0);
							omp_set_num_threads(numThreads);
							double M = 0;
#pragma omp parallel for schedule (guided) num_threads(numThreads) shared(numSols, numDiffNumSols, maxSols, minSols, count, maxE_ground1, maxE_ground2, maxDx1, maxDx2)
							for (int i = 0; i < numConfigs; i++) {
								Configuration c = getTargetConfig(i);
								//outfile << "NumConst=" << NumConst << ", testing...";
								//outfile.flush();
								std::vector<double> result_temp1;
								//std::vector<double> result_temp2;
								std::vector<Configuration> sols1;
								//std::vector<Configuration> sols2;

								FrankInvCCUniqueTest4(c, Nc, NumConst, result_temp1, sols1, outfile, PositionTolerance, EnergyTolerance, 0, "linear", method, filename);
								//FrankInvCCUniqueTest3(c, Nc, NumConst, result_temp2, sols2, outfile, PositionTolerance, EnergyTolerance, 0, "image", method, filename);
								//FrankInvCCUniqueTest2(c, Nc, NumConst, result_temp, outfile, 1e-6, 0, false, method);


								//if (NumConst > N_eff - 1) {
								//	for (size_t q = 1; q < sols1.size(); q++) {
								//		Configuration ref(sols1[q]);
								//		for (size_t p = 1; p < sols2.size(); p++) {
								//			double maxDr = 0;
								//			for (size_t s = 0; s < ref.NumParticle(); s++) {
								//				GeometryVector r = sols2[p].GetRelativeCoordinates(s);
								//				double dr = ref.NearestParticleDistance(r);
								//				maxDr = std::max(dr, maxDr);
								//			}
								//			if (maxDr < 1e-6) {

								//			}
								//		}
								//	}
								//}


#pragma omp critical
								{
									if (M == 0) {
										M = result_temp1[1];
									}
									//double temp_val = (double)Nc / result_temp[4]; //success rate of finding ground states
									double temp_val = (double)result_temp1[4] / Nc; //success rate of finding ground states
									result_ground_real[ind].x[1] += temp_val;	result_ground_real[ind].x[2] += temp_val*temp_val;

									//temp_val = (double)result_temp[5] / Nc; //success rate of reconstructing the ground states.
									temp_val = (double)result_temp1[5] / result_temp1[4]; //success rate of reconstructing the ground states.
									result_suc_real[ind].x[1] += temp_val;	result_suc_real[ind].x[2] += temp_val*temp_val;

									//temp_val = (double)result_temp2[4] / Nc; //success rate of finding ground states
									result_ground_imag[ind].x[1] += temp_val;	result_ground_imag[ind].x[2] += temp_val*temp_val;

									//temp_val = (double)result_temp[5] / Nc; //success rate of reconstructing the ground states.
									//temp_val = (double)result_temp2[5] / result_temp2[4]; //success rate of reconstructing the ground states.
									result_suc_imag[ind].x[1] += temp_val;	result_suc_imag[ind].x[2] += temp_val*temp_val;

									maxSols = std::max(maxSols, sols1.size());
									minSols = std::min(minSols, sols1.size());

									numSols.x[0] += sols1.size();
									//numSols.x[1] += sols2.size();
									numSols.x[2] += sols1.size()*sols1.size();
									//numSols.x[3] += sols2.size()*sols2.size();

									maxDx1 = std::max(maxDx1, result_temp1[2]);
									maxE_ground1 = std::max(maxE_ground1, result_temp1[3]);
									//maxDx2 = std::max(maxDx2, result_temp2[2]);
									//maxE_ground2 = std::max(maxE_ground2, result_temp2[3]);
									if (NumConst > N_eff - 1) {
										//output some details!!!
										outfile << count << "\tTarget\t";
										{//target configuration
											Configuration tempC(sols1[0]);
											Merge_Sort(tempC, 0);
											for (int z = 0; z < tempC.NumParticle(); z++)
												outfile << std::setprecision(7) << tempC.GetCartesianCoordinates(z).x[0] << ",\t";

											outfile << "\n";
										}

										//Real constraints
										//if (NumConst % 2 == 1) 
										if (type != 5) {//when the number of constraints are odd
											if (sols1.size() == 1) {
												outfile << count << "\tReal\tunique";
												for (int m = 0; m < N; m++)
													outfile << "\t\t";

												outfile << result_temp1[6] << "\t" << result_temp1[7] << "\t" << result_temp1[8] << std::endl;

											}
											else {
												for (int t = 0; t < sols1.size(); t++) {
													outfile << count << "\tReal\t";
													Configuration tempC(sols1[t]);
													Merge_Sort(tempC, 0);
													for (int y = 0; y < tempC.NumParticle(); y++)
														outfile << tempC.GetCartesianCoordinates(y).x[0] << ",\t";

													outfile << result_temp1[3 * t + 6] << "\t" << result_temp1[3 * t + 7] << "\t" << result_temp1[3 * t + 8] << std::endl;
													//if (t == 1) {
													//	outfile << result_temp1[3] << "\t"<<result_temp1[2];
													//}
												}
												//outfile << std::endl;
											}

											//if (sols2.size() == 1) {
											//	outfile << count << "\tImag\tunique";
											//	for (int m = 0; m < N; m++)
											//		outfile << "\t\t";

											//	outfile << result_temp2[6] << "\t" << result_temp2[7] << "\t" << result_temp2[8] << std::endl;
											//}
											//else {
											//	for (int t = 0; t < sols2.size(); t++) {
											//		outfile << count << "\tImag\t";
											//		Configuration tempC(sols2[t]);
											//		Merge_Sort(tempC, 0);
											//		for (int y = 0; y < tempC.NumParticle(); y++)
											//			outfile << tempC.GetCartesianCoordinates(y).x[0] << ",\t";

											//		outfile << result_temp2[3 * t + 6] << "\t" << result_temp2[3 * t + 7] << "\t" << result_temp2[3 * t + 8] << std::endl;
											//		//if (t == 1) {
											//		//	outfile << result_temp1[3] << "\t"<<result_temp1[2];
											//		//}
											//	}
											//	//outfile << std::endl;
											//}
										}
										count++;
										outfile << std::endl;
									}

									/*if (sols1.size() != sols2.size()) {
										maxSols = std::max(maxSols, sols2.size());
										minSols = std::min(minSols, sols2.size());
										numDiffNumSols++;
										if (NumConst > N_eff - 1) {
											char temp_name[200] = {};
											for (int l = 0; l < sols1.size(); l++) {
												sprintf(temp_name, "%s%s__%s__%s_N_%d_NumConst_%d_groundstates_%d_%d.txt", dir, systName, method, "real", N, NumConst, i, l);
												WriteConfiguration(sols1[l], temp_name);
											}
											for (int l = 0; l < sols2.size(); l++) {
												sprintf(temp_name, "%s%s__%s__%s_N_%d_NumConst_%d_groundstates_%d_%d.txt", dir, systName, method, "imag", N, NumConst, i, l);
												WriteConfiguration(sols2[l], temp_name);
											}
										}
									}*/

								}
							}
							result_ground_real[ind].x[0] = NumConst;	result_ground_real[ind].x[1] /= (double)numConfigs;	result_ground_real[ind].x[2] = sqrt(result_ground_real[ind].x[2] / (double)numConfigs - result_ground_real[ind].x[1] * result_ground_real[ind].x[1]) / sqrt(numConfigs);
							result_suc_real[ind].x[0] = NumConst;	result_suc_real[ind].x[1] /= (double)numConfigs;	result_suc_real[ind].x[2] = sqrt(result_suc_real[ind].x[2] / (double)numConfigs - result_suc_real[ind].x[1] * result_suc_real[ind].x[1]) / sqrt(numConfigs);
							result_ground_imag[ind].x[0] = NumConst;	result_ground_imag[ind].x[1] /= (double)numConfigs;	result_ground_imag[ind].x[2] = sqrt(result_ground_imag[ind].x[2] / (double)numConfigs - result_ground_imag[ind].x[1] * result_ground_imag[ind].x[1]) / sqrt(numConfigs);
							result_suc_imag[ind].x[0] = NumConst;	result_suc_imag[ind].x[1] /= (double)numConfigs;	result_suc_imag[ind].x[2] = sqrt(result_suc_imag[ind].x[2] / (double)numConfigs - result_suc_imag[ind].x[1] * result_suc_imag[ind].x[1]) / sqrt(numConfigs);
							ind++;

							double meanSols1 = (double)numSols.x[0] / numConfigs;
							//double meanSols2 = (double)numSols.x[1] / numConfigs;
							double stdSols1 = sqrt(((double)numSols.x[2] / numConfigs - meanSols1*meanSols1) / numConfigs);
							//double stdSols2 = sqrt(((double)numSols.x[3] / numConfigs - meanSols2*meanSols2) / numConfigs);


							std::cout << M << "\tReal\t" << maxSols << "\t" << minSols << "\t" << (double)(numDiffNumSols) / numConfigs << "\t" << meanSols1 << "\t" << stdSols1 << "\t" << maxE_ground1 << "\t" << maxDx1 << "\n";
							//std::cout << NumConst << "\tImag\t" << maxSols << "\t" << minSols << "\t" << (double)(numDiffNumSols) / numConfigs << "\t" << meanSols2 << "\t" << stdSols2 << "\t" << maxE_ground2 << "\t" << maxDx2 << "\n";

						}

						//sprintf(filename, "%s%dD%s__%s__N_%d_groundstates_%s.txt", dir, dim, systName, method, N, "imag");
						//WriteFunction(result_ground_imag, filename);
						//sprintf(filename, "%s%dD%s__%s__N_%d_uniquestates_%s.txt", dir, dim, systName, method, N, "imag");
						//WriteFunction(result_suc_imag, filename);

						sprintf(filename, "%s%dD%s__%s__N_%d_groundstates_%s.txt", dir, dim, systName, method, N, "real");
						WriteFunction(result_ground_real, filename);
						sprintf(filename, "%s%dD%s__%s__N_%d_uniquestates_%s.txt", dir, dim, systName, method, N, "real");
						WriteFunction(result_suc_real, filename);
						//Dist.clear();
						data.clear();
					}
				}
			}
		}
	}

	return 0;




	//{
	//	DimensionType d = 2;
	//	int numPrt;
	//	int numConfig;
	//	int numThreads;

	//	double boxSize;
	//	char systName[100] = {};
	//	char savename[200] = {}; //type_dimension_prtNumber_boxSize_numConfigs_characteristic parameter
	//	ofile << "The number of particles = ";
	//	ifile >> numPrt;
	//	ofile << "The number of configurations = ";
	//	ifile >> numConfig;
	//	ofile << "box size = ";
	//	ifile >> boxSize;
	//	ofile << "System type : ";
	//	ifile >> systName;
	//	ofile << "The number of threads = ";
	//	ifile >> numThreads;

	//	double phi;
	//	int timeStep;
	//	int adjTime;
	//	double r;
	//	//ofile << "phi = ";
	//	//ifile >> phi;
	//	ofile << "r = ";
	//	ifile >> r;
	//	ofile << "Time for adjustment";
	//	ifile >> adjTime;
	//	ofile << "time steps = ";
	//	ifile >> timeStep;

	//	//double r = pow(phi / (numPrt*HyperSphere_Volume(d, 1.0)), 1.0 / d);

	//	std::vector<GeometryVector> basis;
	//	basis.push_back(GeometryVector(1.0, 0));
	//	basis.push_back(GeometryVector(0, 1.0));
	//	int num = (int)ceil(pow(numPrt, 1.0 / d));
	//	Configuration templ = GetLattice_config(basis, num);
	//	while (templ.NumParticle() > numPrt) {
	//		templ.DeleteParticle(templ.NumParticle() - 1);
	//	}
	//	templ.Resize(1.0);//unit periodic box
	//	SpherePacking c(templ, r);

	//	//SpherePacking c(d, &templ.GetBasisVector(0), 1.0);
	//	//for (int i = 0; i < numPrt; i++) {
	//	//	c.Insert(r, templ.GetRelativeCoordinates(i));
	//	//}
	//	templ.RemoveParticles();
	//	c.Rescale(boxSize);
	//	HS_pot test(d, (size_t)0);
	//	test.SetConfiguration(c);
	//	test.AdjustDisplacement(adjTime);

	//	for (int i = 0; i < numConfig; i++) {
	//		for (int t = 0; t < timeStep; t++) {
	//			test.Evolve();
	//		}
	//		//Configuration temptemp(test.Output(), "a");
	//		//Merge_Sort(temptemp, 0);
	//		//configs.push_back(temptemp);
	//		
	//		char dd[100];
	//		sprintf(dd, "./results/test_%d", i);
	//		PlotPacking2D(dd, test.Output());

	//		double temp;
	//		ofile << "continue?";
	//		ifile >> temp;
	//	}


	//}

	//return 0;

//	{	
//		char dir[100];
//		char dir0[100] = "results/Reconstructions/";
//	std::stringstream stream;
//	char temp[200] = {};
//	char name[200] = {};
//
//	//------------------
//	DimensionType dim = 1;
//	size_t Nc = 100;	//number of trial configs
//	int numConfigs = 1000; //number of target configurations
//	size_t NKSkip = 0;
//	char constType[100] = {};
//	int numThreads = 6;
//
//	//int numBins = 40;
//	int N1 = 10;//min system size
//	int N2 = 40;//max system size
//	int systType1;
//	int systType2;
//	int dN = 1;
//	int methodNum1 = 0;
//	int methodNum2 = 0;
//	double EnergyTolerance = 1e-20;
//	double PositionTolerance = 1e-6;
//	double delta1;
//	double delta2;
//	double d_delta;
//	
//	std::cout << "NKSkip=";
//	std::cin >> NKSkip;
//	std::cout << "dim=";
//	std::cin >> dim;
//	std::cout << "min N=";
//	std::cin >> N1;
//	std::cout << "max N=";
//	std::cin >> N2;
//	std::cout << "dN=";
//	std::cin >> dN;
//	std::cout << "numConfigs=";
//	std::cin >> numConfigs;
//	std::cout << "Nc=";
//	std::cin >> Nc;
//	//std::cout << "numBins=";
//	//std::cin >> numBins;
//	std::cout << "system type1:";
//	std::cin >> systType1;
//	std::cout << "system type2:";
//	std::cin >> systType2;
//	std::cout << "initial method num: ";
//	std::cin >> methodNum1;
//	std::cout << "final method num: ";
//	std::cin >> methodNum2;
//	ofile << "Energ Tolerance = ";
//	ifile >> EnergyTolerance;
//	ofile << "delta1 = ";
//	ifile >> delta1;
//	ofile << "delta2 = ";
//	ifile >> delta2;
//	ofile << "d_delta = ";
//	ifile >> d_delta;
//
//
//	//std::cout << "constraints type: ";
//	//std::cin >> constType;
//	//	double Kmax = 4* pi;
//
//	if (NKSkip > 0)
//		std::cout << "Skip the central " << NKSkip << " k points.\n";
//
//	std::vector<GeometryVector> basis;
//	for (int j = 0; j < dim;j++) {
//		basis.push_back(GeometryVector());
//		basis[j].SetDimension(dim);
//		basis[j].x[j] = 1.0;
//	}
//
//	//double delta = 0.1;
//	for (double delta = delta1; delta <= delta2; delta += d_delta) {
//		for (int method_num = methodNum1; method_num <= methodNum2; method_num++) {
//
//			char method[100] = {};
//			switch (method_num)
//			{
//			case 0:
//				sprintf(method, "NLOP+MINOP"); break;
//			case 1:
//				sprintf(method, "Steepest Descent"); break;
//			case 2:
//				sprintf(method, "ConjugateGradient"); break;
//			case 3:
//				sprintf(method, "LocalGradientDescent"); break;
//			default:
//				break;
//			}
//
//			std::cout << "method = " << method << std::endl;
//			for (int type = systType1; type <= systType2; type++) {
//				char systName[200] = {};
//				switch (type)
//				{
//				case 0:
//					sprintf(systName, "perturbedLattice(1prt)__delta_%0.4f", delta);	break;
//				case 1:
//					sprintf(systName, "perturbedLattice(Nprt)__delta_%0.4f", delta);	break;
//				case 2:
//					sprintf(systName, "Poisson");	break;
//				case 3:
//					sprintf(systName, "lattice_vacancy");	break;
//				case 4:
//					sprintf(systName, "lattice_interstitial");	break;
//				case 5:
//					sprintf(systName, "non-Bravais"); break;
//				default:
//					break;
//				}
//
//				std::cout << " system type = " << systName << std::endl;
//				for (int N = N1; N <= N2; N += dN) {
//
//					char filename[400] = {};
//					sprintf(filename, "%s%s__%s__N_%d.txt", dir, systName, method, N);
//					std::fstream outfile(filename, std::fstream::out);
//					sprintf(filename, "%s%s__%s_summary.txt", dir, systName, method);
//					std::ofstream summary(filename);
//
//					//wave vectors
//					//double k1 = 0;
//					//std::vector<GeometryVector> kList;
//					int N_eff = N;
//					if (type == 3) { N_eff--; }
//					else if (type == 4) { N_eff++; }
//
//					//for (int i = 0; i < N_eff; i += 2) {
//					//	 k1 += 2.0*pi / (double)N;
//					//	 kList.push_back(GeometryVector((double)(k1)));
//					// }
//					//std::cout << "type = " << systName << "\tN = " << N << "\t numConstraints = " << kList.size() << "\n";
//					auto getTargetConfig = [&type, &N, dim, &basis, delta](size_t i)->Configuration {
//						Configuration c;
//						GeometryVector temp;    temp.SetDimension(dim);
//						RandomGenerator g(std::time(NULL) + i);
//						if (type == 0) {
//							c = GetLattice_config(basis, N);
//							int j = (int)floor(N*g.RandomDouble());
//							RandomVector(dim, temp);
//							temp = delta*temp;
//							c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
//							char name[200] = {};
//							//sprintf(name, "CC_Inv_Log_1PL_d_%d_NKSkip_%d___%d", dim, NKSkip,N);
//							//WriteConfiguration(c, name);
//						}
//						else if (type == 1) {
//							c = GetLattice_config(basis, N);
//							for (int j = 0; j < N; j++) {
//								//int j = (int)floor(N*g.RandomDouble());
//								RandomVector(dim, g, temp);
//								temp = delta*temp;
//								c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
//							}
//							char name[200] = {};
//							//sprintf(name, "CC_Inv_Log_PL_d_%d_NKSkip_%d___%d", dim, NKSkip, N);
//							//WriteConfiguration(c, name);
//						}
//						else if (type == 2) {
//							c = GeneratePoisson(dim, N, g);
//							c.Resize(N);
//						}
//						else if (type == 3) {//vacancies
//							c = GetLattice_config(basis, N);
//							int j = (int)floor(N*g.RandomDouble());
//							c.DeleteParticle(j);
//						}
//						else if (type == 4) {//interstitial
//							c = GetLattice_config(basis, N);
//							c.Insert("a", g);
//						}
//						else if (type == 5) {// random non-Bravais lattice
//							GeometryVector *A = new GeometryVector[1];
//							A[0].Dimension = 1;
//							A[0].x[0] = N;
//							c = Configuration(dim, A, 1.0);
//							int numBasisAtom = Find_LeastCommonDivider(N);
//							std::vector<GeometryVector> atomBasis;
//							GeometryVector X(static_cast<int> (1));
//							RandomGenerator g(std::time(NULL) + i);
//							for (int i = 0; i < numBasisAtom; i++) {
//								X.x[0] = (double)numBasisAtom*g.RandomDouble();
//								atomBasis.push_back(X);
//							}
//							X.x[0] = 0;
//							for (int i = 0; i < N; i++) {
//								int j = i%numBasisAtom;
//								int k = (int)floor(i / numBasisAtom);
//
//								X.x[0] = k*numBasisAtom;
//								c.Insert("a", c.CartesianCoord2RelativeCoord(X + atomBasis[j]));
//							}
//						}
//						return c;
//
//					};
//
//
//					//---------------init config
//					sprintf(filename, "%stempInit/%d", dir0, N);
//					Configuration templateConfig(getTargetConfig(0));
//					RandomGenerator g(std::time(NULL));
//					auto getInitConfig = [&templateConfig, &g](size_t i)->Configuration {
//						Configuration dd(templateConfig);
//						size_t numP = dd.NumParticle();
//						for (size_t z = 0; z < numP; z++) {
//							dd.Insert("a", g);
//						}
//						return dd;
//					};
//					ConfigurationPack initConfig(filename);
//					if (initConfig.NumConfig() == 0) {
//						for (int z = 0; z < Nc; z++) {
//							initConfig.AddConfig(getInitConfig(z));
//						}
//					}
//
//					std::vector<double> data;	data.reserve(numConfigs);
//
//					double min = 0;	double max = 0;	int ind = 0;
//					std::vector<GeometryVector> result_suc_real;
//					std::vector<GeometryVector> result_ground_real;
//					std::vector<GeometryVector> result_suc_imag;
//					std::vector<GeometryVector> result_ground_imag;
//
//					std::cout << "\nN=" << N << "\n";
//					std::cout << "M\tType\tmaxSols\tminSols\tProb_DiffSols\tav.NumSols\tstdSols\tMaxE_ground\tMaxdx\n";
//					for (size_t NumConst = N_eff - 1; NumConst <= N_eff + 1; NumConst++) {
//						outfile << "-----------------------------------------------------------------" << std::endl;
//						outfile << "Testing N=" << N << "\tM=" << NumConst << std::endl << std::endl;
//						outfile << "num\tconstr\tconfigs";
//						for (int z = 0; z < N; z++)
//							outfile << "\t\t";
//						outfile << "\tminE\t\tmaxE\t\tmaxDx" << std::endl;
//
//						result_suc_real.push_back(GeometryVector(NumConst, 0, 0, 0));
//						result_ground_real.push_back(GeometryVector(NumConst, 0, 0, 0));
//						result_suc_imag.push_back(GeometryVector(NumConst, 0, 0, 0));
//						result_ground_imag.push_back(GeometryVector(NumConst, 0, 0, 0));
//
//						size_t maxSols = 0, minSols = Nc, numDiffNumSols = 0;
//						int count = 0;
//						double maxE_ground1 = 0, maxE_ground2 = 0, maxDx1 = 0, maxDx2 = 0;
//						//size_t numSols_1 = 0, numSols_12 = 0, numSols_2=0, numSols_22=0;
//						GeometryVector numSols(0, 0, 0, 0);
//						omp_set_num_threads(numThreads);
//#pragma omp parallel for schedule (guided) num_threads(numThreads) shared(numSols, numDiffNumSols, maxSols, minSols, count, maxE_ground1, maxE_ground2, maxDx1, maxDx2)
//						for (int i = 0; i < numConfigs; i++) {
//							Configuration c = getTargetConfig(i);
//							//outfile << "NumConst=" << NumConst << ", testing...";
//							//outfile.flush();
//							std::vector<double> result_temp1;
//							std::vector<double> result_temp2;
//							std::vector<Configuration> sols1;
//							std::vector<Configuration> sols2;
//
//							FrankInvCCUniqueTest3(c, Nc, NumConst, result_temp1, sols1, outfile, PositionTolerance, EnergyTolerance, 0, "real", method, filename);
//							FrankInvCCUniqueTest3(c, Nc, NumConst, result_temp2, sols2, outfile, PositionTolerance, EnergyTolerance, 0, "image", method, filename);
//							//FrankInvCCUniqueTest2(c, Nc, NumConst, result_temp, outfile, 1e-6, 0, false, method);
//
//
//							//if (NumConst > N_eff - 1) {
//							//	for (size_t q = 1; q < sols1.size(); q++) {
//							//		Configuration ref(sols1[q]);
//							//		for (size_t p = 1; p < sols2.size(); p++) {
//							//			double maxDr = 0;
//							//			for (size_t s = 0; s < ref.NumParticle(); s++) {
//							//				GeometryVector r = sols2[p].GetRelativeCoordinates(s);
//							//				double dr = ref.NearestParticleDistance(r);
//							//				maxDr = std::max(dr, maxDr);
//							//			}
//							//			if (maxDr < 1e-6) {
//
//							//			}
//							//		}
//							//	}
//							//}
//
//
//#pragma omp critical
//							{
//
//								//double temp_val = (double)Nc / result_temp[4]; //success rate of finding ground states
//								double temp_val = (double)result_temp1[4] / Nc; //success rate of finding ground states
//								result_ground_real[ind].x[1] += temp_val;	result_ground_real[ind].x[2] += temp_val*temp_val;
//
//								//temp_val = (double)result_temp[5] / Nc; //success rate of reconstructing the ground states.
//								temp_val = (double)result_temp1[5] / result_temp1[4]; //success rate of reconstructing the ground states.
//								result_suc_real[ind].x[1] += temp_val;	result_suc_real[ind].x[2] += temp_val*temp_val;
//
//								temp_val = (double)result_temp2[4] / Nc; //success rate of finding ground states
//								result_ground_imag[ind].x[1] += temp_val;	result_ground_imag[ind].x[2] += temp_val*temp_val;
//
//								//temp_val = (double)result_temp[5] / Nc; //success rate of reconstructing the ground states.
//								temp_val = (double)result_temp2[5] / result_temp2[4]; //success rate of reconstructing the ground states.
//								result_suc_imag[ind].x[1] += temp_val;	result_suc_imag[ind].x[2] += temp_val*temp_val;
//
//								maxSols = std::max(maxSols, sols1.size());
//								minSols = std::min(minSols, sols1.size());
//
//								numSols.x[0] += sols1.size();
//								numSols.x[1] += sols2.size();
//								numSols.x[2] += sols1.size()*sols1.size();
//								numSols.x[3] += sols2.size()*sols2.size();
//
//								maxDx1 = std::max(maxDx1, result_temp1[2]);
//								maxE_ground1 = std::max(maxE_ground1, result_temp1[3]);
//								maxDx2 = std::max(maxDx2, result_temp2[2]);
//								maxE_ground2 = std::max(maxE_ground2, result_temp2[3]);
//								if (NumConst > N_eff - 1) {
//									//output some details!!!
//									outfile << count << "\tTarget\t";
//									{//target configuration
//										Configuration tempC(sols1[0]);
//										Merge_Sort(tempC, 0);
//										for (int z = 0; z < tempC.NumParticle(); z++)
//											outfile << std::setprecision(7) << tempC.GetCartesianCoordinates(z).x[0] << ",\t";
//
//										outfile << "\n";
//									}
//
//									//Real constraints
//									//if (NumConst % 2 == 1) 
//									if (type != 5) {//when the number of constraints are odd
//										if (sols1.size() == 1) {
//											outfile << count << "\tReal\tunique";
//											for (int m = 0; m < N; m++)
//												outfile << "\t\t";
//
//											outfile << result_temp1[6] << "\t" << result_temp1[7] << "\t" << result_temp1[8] << std::endl;
//
//										}
//										else {
//											for (int t = 0; t < sols1.size(); t++) {
//												outfile << count << "\tReal\t";
//												Configuration tempC(sols1[t]);
//												Merge_Sort(tempC, 0);
//												for (int y = 0; y < tempC.NumParticle(); y++)
//													outfile << tempC.GetCartesianCoordinates(y).x[0] << ",\t";
//
//												outfile << result_temp1[3 * t + 6] << "\t" << result_temp1[3 * t + 7] << "\t" << result_temp1[3 * t + 8] << std::endl;
//												//if (t == 1) {
//												//	outfile << result_temp1[3] << "\t"<<result_temp1[2];
//												//}
//											}
//											//outfile << std::endl;
//										}
//
//										if (sols2.size() == 1) {
//											outfile << count << "\tImag\tunique";
//											for (int m = 0; m < N; m++)
//												outfile << "\t\t";
//
//											outfile << result_temp2[6] << "\t" << result_temp2[7] << "\t" << result_temp2[8] << std::endl;
//										}
//										else {
//											for (int t = 0; t < sols2.size(); t++) {
//												outfile << count << "\tImag\t";
//												Configuration tempC(sols2[t]);
//												Merge_Sort(tempC, 0);
//												for (int y = 0; y < tempC.NumParticle(); y++)
//													outfile << tempC.GetCartesianCoordinates(y).x[0] << ",\t";
//
//												outfile << result_temp2[3 * t + 6] << "\t" << result_temp2[3 * t + 7] << "\t" << result_temp2[3 * t + 8] << std::endl;
//												//if (t == 1) {
//												//	outfile << result_temp1[3] << "\t"<<result_temp1[2];
//												//}
//											}
//											//outfile << std::endl;
//										}
//									}
//									count++;
//									outfile << std::endl;
//								}
//
//								if (sols1.size() != sols2.size()) {
//									maxSols = std::max(maxSols, sols2.size());
//									minSols = std::min(minSols, sols2.size());
//									numDiffNumSols++;
//									//if (NumConst > N_eff - 1) {
//									//	char temp_name[200] = {};
//									//	for (int l = 0; l < sols1.size(); l++) {
//									//		sprintf(temp_name, "%s%s__%s__%s_N_%d_NumConst_%d_groundstates_%d_%d.txt", dir, systName, method, "real", N, NumConst, i, l);
//									//		WriteConfiguration(sols1[l], temp_name);
//									//	}
//									//	for (int l = 0; l < sols2.size(); l++) {
//									//		sprintf(temp_name, "%s%s__%s__%s_N_%d_NumConst_%d_groundstates_%d_%d.txt", dir, systName, method, "imag", N, NumConst, i, l);
//									//		WriteConfiguration(sols2[l], temp_name);
//									//	}
//									//}
//								}
//
//							}
//						}
//						result_ground_real[ind].x[0] = NumConst;	result_ground_real[ind].x[1] /= (double)numConfigs;	result_ground_real[ind].x[2] = sqrt(result_ground_real[ind].x[2] / (double)numConfigs - result_ground_real[ind].x[1] * result_ground_real[ind].x[1]) / sqrt(numConfigs);
//						result_suc_real[ind].x[0] = NumConst;	result_suc_real[ind].x[1] /= (double)numConfigs;	result_suc_real[ind].x[2] = sqrt(result_suc_real[ind].x[2] / (double)numConfigs - result_suc_real[ind].x[1] * result_suc_real[ind].x[1]) / sqrt(numConfigs);
//						result_ground_imag[ind].x[0] = NumConst;	result_ground_imag[ind].x[1] /= (double)numConfigs;	result_ground_imag[ind].x[2] = sqrt(result_ground_imag[ind].x[2] / (double)numConfigs - result_ground_imag[ind].x[1] * result_ground_imag[ind].x[1]) / sqrt(numConfigs);
//						result_suc_imag[ind].x[0] = NumConst;	result_suc_imag[ind].x[1] /= (double)numConfigs;	result_suc_imag[ind].x[2] = sqrt(result_suc_imag[ind].x[2] / (double)numConfigs - result_suc_imag[ind].x[1] * result_suc_imag[ind].x[1]) / sqrt(numConfigs);
//						ind++;
//
//						double meanSols1 = (double)numSols.x[0] / numConfigs;
//						double meanSols2 = (double)numSols.x[1] / numConfigs;
//						double stdSols1 = sqrt(((double)numSols.x[2] / numConfigs - meanSols1*meanSols1) / numConfigs);
//						double stdSols2 = sqrt(((double)numSols.x[3] / numConfigs - meanSols2*meanSols2) / numConfigs);
//
//
//						std::cout << NumConst << "\tReal\t" << maxSols << "\t" << minSols << "\t" << (double)(numDiffNumSols) / numConfigs << "\t" << meanSols1 << "\t" << stdSols1 << "\t" << maxE_ground1 << "\t" << maxDx1 << "\n";
//						std::cout << NumConst << "\tImag\t" << maxSols << "\t" << minSols << "\t" << (double)(numDiffNumSols) / numConfigs << "\t" << meanSols2 << "\t" << stdSols2 << "\t" << maxE_ground2 << "\t" << maxDx2 << "\n";
//
//					}
//
//					sprintf(filename, "%s%s__%s__N_%d_groundstates_%s.txt", dir, systName, method, N, "imag");
//					WriteFunction(result_ground_imag, filename);
//					sprintf(filename, "%s%s__%s__N_%d_uniquestates_%s.txt", dir, systName, method, N, "imag");
//					WriteFunction(result_suc_imag, filename);
//
//					sprintf(filename, "%s%s__%s__N_%d_groundstates_%s.txt", dir, systName, method, N, "real");
//					WriteFunction(result_ground_real, filename);
//					sprintf(filename, "%s%s__%s__N_%d_uniquestates_%s.txt", dir, systName, method, N, "real");
//					WriteFunction(result_suc_real, filename);
//					//Dist.clear();
//					data.clear();
//				}
//			}
//		}
//	}
//	}
//
//		 return 0;




//	char dir[100] = "results/MC/";
//	//molecular dynamics of harmonic lattice.
//
//
//	//------------ system parameters ----------------
//	DimensionType d = 1;
//	int L = 20; //systSize;
//	int N = pow(L, d);
//	double M = 1;//particle mass
//	double K = 1;// spring constant
//	double cL = sqrt(K / M);
//	double a = 1;//lattice constant
//
//	double w_max = 2.0*cL / a;
//	double w_min = cL*2.0*pi / L;
//	double dt = 0.01* 2.0*pi / w_max; //time step
//	double T_f = 3.0*2.0*pi / w_max;
//	double t0 = 2.0*pi / w_min;
//	int n_t = (int)round(T_f / dt); //recording step
//
//									//------------ simulation parameters ------------
//									//int numTimestep = 10000;
//	int numTimestep = 10000;
//	int numSamplingStep = 100;
//	int numTimeTest = 1000;
//	int numEquilIntv = 1000;
//	int NumThreads = 1;
//	double accRateL = 0.45;
//	double accRateU = 0.6;
//	double coolingRate = 0.95;
//
//	double dk = 2.0*pi / L;
//	double k_max = 3.0*pi;
//
//	double temp0 = 0.09; //near melting point
//	double temp1 = 0.05; //first temperature to investigate
//	double temp2 = 0.05; //final temperature to investigate
//	double d_temp = 0.00004;
//	double kT = temp0; //temperature of the systems
//
//	double dx = 6.0* sqrt(kT);
//
//
//	//GeometryVector * basis = new GeometryVector[d];
//	std::vector<GeometryVector> basis;
//	for (int i = 0; i < d; i++) {
//		basis.push_back(GeometryVector(static_cast<int>(d)));
//		basis[i].x[i] = a;
//	}
//	Configuration prtConfig = GetLattice_config(basis, L);
//	//-----------prepare an ensemble
//	std::vector<SimpleHarmonic> ensemble;
//	std::vector<RandomGenerator> g_set;
//
//	SimpleHarmonic temp(d, K, a);
//	temp.SetReference(prtConfig);
//	temp.SetConfiguration(prtConfig);
//
//	for (int i = 0; i < NumThreads; i++) {
//		g_set.push_back(RandomGenerator(std::time(NULL) + i));
//		ensemble.push_back(temp);
//		//give initial sinusoidal displacements
//		//for (int l = 0; l < d; l++) { //for each direction
//		//      GeometryVector A;
//		//      GeometryVector phi;
//		//      GeometryVector q(static_cast<int> (d));
//		//      GeometryVector du(static_cast<int> (d));
//		//      for (int j = 0; j < d; j++) {
//		//              A.x[j] = g_set[i].RandomDouble()*sqrt(kT)/(double)N;
//		//              phi.x[j] = 2.0*pi*g_set[i].RandomDouble();
//		//      }
//		//      q.x[l] = 2.0*pi / (double)L;
//
//		//      for (int j = 0; j < N; j++) {
//		//              //compute the displacement
//		//              for (int k = 0; k < d; k++) {
//		//                      du.x[k] = A.x[k] * sin(phi.x[k] + q.Dot(ensemble[i].confi.GetCartesianCoordinates(j)));
//		//              }
//		//              ensemble[i].updateUn(j, du);
//		//      }
//		//}
//	}
//	std::cout << "end initialization \n";
//
//
//
//	omp_set_num_threads(NumThreads);
//	for (double targetkT = temp1; targetkT >= temp2; targetkT -= d_temp) {
//		//-----------cooling from the initial temperature to the first temperature to investigate
//		int numOfCoolingSteps = (int)round(log(targetkT / kT) / log(coolingRate));
//		coolingRate = pow(targetkT / kT, 1.0 / (double)numOfCoolingSteps);
//		std::cout << "target kT = " << targetkT << "\tcooling rate = " << coolingRate << "\t the number of cooling steps = " << numOfCoolingSteps << "\n";
//		for (int i = 0; i <= numOfCoolingSteps; i++) {
//			std::cout << "kT = " << kT << "\n";
//			SetTrialMovement(ensemble, g_set, dx, kT, accRateL, accRateU, numTimeTest);
//			std::cout << "-----finish the adjustment \n";
//
//			double rate = 0;
//			for (int j = 0; j < numEquilIntv; j++) {//test to whether they are in equilibrium
//				rate = MCevolve(ensemble, g_set, kT, dx);
//				//      std::cout << rate << "\t";
//			}
//			std::cout << "\n\n";
//
//			kT *= coolingRate;
//			if (kT < temp1)
//				kT /= coolingRate;
//		}
//		kT = targetkT; //....
//
//					   //
//		std::cout << "more equilibration \n";
//		for (int i = 0; i < numTimestep; i++) {
//			MCevolve(ensemble, g_set, kT, dx);
//		}
//
//#pragma omp parallel for schedule (guided)
//		for (int i = 0; i < ensemble.size(); i++) {
//			ensemble[i].updateConfiguration();
//		}
//
//		//------energy distribution--------
//		double dE = 0.05*kT;
//		double E_max = 10.00*kT;
//		std::vector<double> E_dist((int)ceil(E_max / dE));
//		double count = 0;
//
//		std::ostringstream filen;
//		//std::string filen = "results\\MC\\energy_DIST L_" +  + "kT_ " +  +".txt";
//		filen << dir << "energy_DIST " << d << "D L_" << L << " kT_ " << kT << ".txt";
//		std::cout << filen.str() << "\n";
//		std::ofstream fname(filen.str());
//		filen.flush();
//		char prefix[100] = {};
//
//
//		//              //generate configurations
//		//#pragma omp parallel for reduction(+:count)
//		//              for (int i = 0; i < NumThreads; i++) {
//		//                      SimpleHarmonic test(d, K, a);
//		//                      RandomGenerator g;
//		//                      test.SetReference(prtConfig);
//		//                      test.SetConfiguration(prtConfig);
//		//
//		//                      //----------------initial seeds
//		//                      for (int j = 0; j < test.confi.NumParticle(); j++) {
//		//                              GeometryVector temp(static_cast<int>(d));
//		//                              for (int k = 0; k < d; k++) {
//		//                                      temp.x[k] = 2.0 * dx *(g.RandomDouble() - 0.5);
//		//                              }
//		//                              test.confi.MoveParticle(j, prtConfig.CartesianCoord2RelativeCoord(prtConfig.GetCartesianCoordinates(j) + temp));
//		//                      }
//		//                      //------------------------------------
//		//
//		//
//		//
//		//                      //----------------evolve configurations
//		//
//		//                      for (int t = -20; t < numTimestep / NumThreads; t++) {
//		//                              int acceptRate = 0;
//		//                              for (int k = 0; k < test.confi.NumParticle(); k++) {
//		//                                      //---- generate a random vector
//		//                                      GeometryVector du(static_cast<int>(d));
//		//                                      for (int l = 0; l < d; l++) {
//		//                                              du.x[l] = dx*(g.RandomDouble() - 0.5);
//		//                                      }
//		//                                      double delta_E = test.Delta_energy(k, du);
//		//                                      double pr = g.RandomDouble();
//		//                                      if (delta_E < 0 || exp(-delta_E / kT) > pr) {
//		//                                              test.confi.MoveParticle(k, test.confi.CartesianCoord2RelativeCoord(test.confi.GetCartesianCoordinates(k) + du));
//		//                                              //if (i<100)
//		//                                              //      acceptRatio++;
//		//                                      }
//		//                              }
//		//
//		//                              //compute the acceptance rate
//		//                              if (t < 100) {
//		//                                      //std::cout << acceptRatio / (double)ttest.confi.NumParticle() << "\n";
//		//std::cout << acceptRatio / (double)ttest.confi.NumParticle() << "\n";
//		//                                      rate.push_back(GeometryVector(acceptRate / (double)test.confi.NumParticle()));
//		//                              }
//		//                              else if (i == 100) {
//		//                                      std::ostringstream f;
//		//                                      f << dir << "accRate " << d << "D L_" << test.confi.NumParticle() << " kT_ " << kT << ".txt";
//		//                                      WriteFunction(rate, f.str());
//		//                              }
//		//
//		//                              for (int j = 0; j < test.confi.NumParticle(); j++) {
//		//                                      double xx = test.Energy(j);
//		//                                      int index = (int)floor(xx / dE);
//		//                                      if (index >= 0 && index < E_dist.size()) {
//		//                                              count++;
//		//                                              E_dist.at(index)++;
//		//                                      }
//		//                              }
//		//
//		//                              if (t >= 0) {
//		//                                      char fileName[100];
//		//                                      sprintf_s(fileName, 100, "%s __%d", prefix, t*NumThreads + i);
//		//                                      //WriteConfiguration(test.confi, fileName);
//		//                              }
//		//                      }
//		//              }
//		//
//		//
//
//		auto getConfig = [&ensemble, numSamplingStep, &g_set, kT, &E_dist, d, dx, &count, E_max, dE, &dir, NumThreads](size_t i)->Configuration {
//			if (i % ensemble.size() == 0) {//evolve configurations
//				for (int j = 0; j<numSamplingStep; j++) {
//					MCevolve(ensemble, g_set, kT, dx);
//				}
//#pragma omp parallel for schedule(guided)
//				for (int j = 0; j < ensemble.size(); j++) {
//					ensemble[j].updateConfiguration();
//					//-------------sampling energy distribution-------
//					//for (int k = 0; k < ensemble[j].confi.NumParticle(); k++) {
//					//      double xx = ensemble[j].Energy(k);
//					//      int index = (int)floor(xx/dE);
//					//      if (index > 0 && index < E_dist.size()) {
//					//#pragma omp critical //sampling energy
//					//              {
//					//                      E_dist.at(index)++;
//					//                      count++;
//					//              }
//
//					//      }
//					//}
//				}
//
//
//
//			}
//			return ensemble[i%NumThreads].confi;
//		};
//
//
//
//
//
//		std::vector<GeometryVector> result;
//
//		for (int oo = 0; oo < 20; oo++) {
//			char fileName[100] = {};
//			sprintf(fileName, "%sconfig_1DthermalXry L_20 kT_0.050000 numSamp_%d", dir, oo);
//			WriteConfiguration(getConfig(oo), fileName);
//			
//		}
//	
//
//		//IsotropicStructureFactor(getConfig, numTimestep, k_max, k_max, result, dk);
//
//		//for (int i = 0; i < E_dist.size(); i++) {
//		//      fname << dE*((double)i+0.5) << "\t" << E_dist[i] / (double(count) *dE) << "\n";
//		//}
//		fname.close();
//		
//		//WriteFunction(result, "results\\MC\\Sk_test");
//		//filen << dir << "Sk "<<d<< "D L_" << L << " kT_" << kT << " numSamp_" << numTimestep;
//		//WriteFunction(result, fileName);
//
//	}
//
//
//
//	return 0;





//	{
//		char dir[100] = "results/CantorSet/";
//		std::stringstream stream;
//		char temp[200] = {};
//		char name[200] = {};
//		DimensionType d = 1;
//		size_t N0 = 100;
//		int numConfigs = 10000;
//		int steps = 9;
//		double phi = 0.5;
//		double kMax = 4.0 * 2.0*pi;
//		double dk = 4.0*pi / (double)N0;
//		std::vector<GeometryVector> basis;
//
//		omp_set_num_threads(6);
//		basis.push_back(GeometryVector(1.0));
//
//		std::vector<double> CantorSet; //In the unusal circumstance, (CantorSet[2*k],CantorSet[2*k+1]) is an open interval deleted...
//		std::vector<double> CantorSet_i;
//		std::vector<SpherePacking> ensemble0;
//		//std::vector<SpherePacking> ensemble;
//		sprintf(name, "%sHS_via_CantorSet_log___phi_%0.1f.txt", dir, phi);
//		std::ofstream ofile(name);
//		std::vector<GeometryVector> result;
//
//		double R;
//		for (int i = 0; i <= steps; i++) {
//			if (i == 0) {
//#pragma omp parallel for schedule(guided)
//				for (int j = 0; j < numConfigs; j++) {
//					RSA setting(d, 2 * N0);
//					setting.gen.seed(std::time(NULL) + j);
//					setting.RSA_I(N0);
//					SpherePacking tempConfig(*(setting.pSpheres));
//					tempConfig.Resize(N0);
//					//SpherePacking tempConfig = GenerateRSAPacking(d, 2 * N0, false, N0);
//					std::vector<GeometryVector> prtList;
//					R = tempConfig.GetMaxRadius();
//					for (int i = 0; i < N0; i++) { prtList.push_back(tempConfig.GetRelativeCoordinates(i)); }
//					MergeSort(prtList, 0);
//					SpherePacking Config(tempConfig);
//					Config.RemoveParticles();
//					for (int k = 0; k < N0; k++) { Config.Insert(phi*R, prtList[k]); }
//#pragma omp critical
//					{
//						ensemble0.push_back(Config);
//					}
//				}
//				ofile << i << "\t" << ensemble0[0].PackingFraction() << "\n";
//				auto getConfig = [&ensemble0](size_t i)->SpherePacking {return ensemble0[i];};
//				IsotropicSpectralDensity(getConfig, numConfigs, kMax, kMax, result, dk);
//				sprintf(name, "%schi_V__HS_via_CantorSet__N0_%d__phi_%0.2f__step_%d", dir, N0, phi, i);
//				WriteFunction(result, name);
//			}
//
//			else {
//				GetCantorSet(CantorSet_i, i);
//				for (int j = 0; j < CantorSet_i.size(); j++) {
//					CantorSet.push_back(CantorSet_i[j]);
//				}
//				/*	int size = CantorSet.size() / 2;
//				GeometryVector center(static_cast<int> (d));
//				double radius = 0;
//				for (int j = 0; j < size; j++) {
//				center.x[0] = 0.5*(CantorSet[2 * j] + CantorSet[2 * j + 1]);
//				radius = phi*0.5*(CantorSet[2 * j + 1] - CantorSet[2 * j]);
//				config.Insert(radius, center);
//				}*/
//
//				//ensemble=std::vector<SpherePacking> (ensemble0);
//#pragma omp parallel for schedule(guided)
//				for (int j = 0; j < numConfigs; j++) {
//					char temp_name[200] = {};
//					SpherePacking temp = ApplyCantorSet2Packing(ensemble0[j], CantorSet, phi, R);
//					sprintf(temp_name, "%stemp/__%d", dir, j);
//					WriteConfiguration(temp, temp_name);
//				}
//				//			ofile << i << "\t" << ensemble[0].PackingFraction()<<"\n";
//
//				auto getConfig = [&dir](size_t i)->SpherePacking {
//					char temp_name[200] = {};
//					sprintf(temp_name, "%stemp/__%d", dir, i);
//					SpherePacking temp_config;
//					ReadConfiguration(temp_config, temp_name);
//					sprintf(temp_name, "%s.txt", temp_name);	remove(temp_name);
//					return temp_config;};
//				IsotropicSpectralDensity(getConfig, numConfigs, kMax, kMax, result, dk);
//				sprintf(name, "%schi_V__HS_via_CantorSet__N0_%d__phi_%0.2f__step_%d", dir, N0, phi, i);
//				WriteFunction(result, name);
//
//			}



			//sprintf(name, "%s__%d", dir,i);
			//result.resize(0);
			//GeometryVector temp;
			//if (i == 0) {
			//	for (int j = 0; j < ensemble0[0].NumParticle(); j++) {
			//		temp = ensemble0[0].GetCartesianCoordinates(j);
			//		result.push_back(GeometryVector(temp.x[0], i, ensemble0[0].GetCharacteristics(j)));
			//	}
			//}
			//else {
			//	for (int j = 0; j < ensemble[0].NumParticle(); j++) {
			//		temp = ensemble[0].GetCartesianCoordinates(j);
			//		result.push_back(GeometryVector(temp.x[0], i, ensemble[0].GetCharacteristics(j)));
			//	}
			//}
			//////WriteConfiguration(HSPacking, name);
			//WriteFunction(result, name);



		}

		//ofile.close();

	
	


	//char dir[] = "results/Reconstructions/";
	//const size_t Nc = 100;
	//size_t NKSkip;
	//size_t dim;
	//int N1;//min system size
	//int N2;//max system size
	//std::cout << "NKSkip=";
	//std::cin >> NKSkip;
	//std::cout << "dim=";
	//std::cin >> dim;
	//std::cout << "min N=";
	//std::cin >> N1;
	//std::cout << "max N=";
	//std::cin >> N2;

	//if (NKSkip>0)
	//	std::cout << "Skip the central " << NKSkip << " k points.\n";
	//
	//for (int methodNum = 0; methodNum < 4; methodNum++) {
	//	char method[100];
	//	switch (methodNum)
	//	{
	//	case 0:
	//		sprintf(method, "NLOP+MINOP"); break;
	//	case 1:
	//		sprintf(method, "Steepest Descent"); break;
	//	case 2:
	//		sprintf(method, "ConjugateGradient"); break;
	//	case 3:
	//		sprintf(method, "LocalGradientDescent"); break;
	//	default:
	//		break;
	//	}
	//	
	//	for (int i = 0; i < 3; i++) {
	//		std::stringstream ss;
	//		ss << dir;
	//		if (i == 0)
	//			ss << "CC_Inv_1PL_d_" << dim << "_NKSkip_" << NKSkip;
	//		else if (i == 1)
	//			ss << "CC_Inv_PL_d_" << dim << "_NKSkip_" << NKSkip;
	//		else if (i == 2)
	//			ss << "CC_Inv_d_" << dim << "_NKSkip" << NKSkip;
	//		else if (i == 3)
	//			ss << "CC_Inv_Fibonacci_d_" << dim << "_NKSkip" << NKSkip;
	//		else
	//			std::cout << "sth wrong \n";

	//		ss << "method_" << method;
	//		char systname[200] = {};
	//		char filename[200] = {};
	//		strncpy(systname, ss.str().c_str(), 200);
	//		sprintf(filename, "%s.txt", systname);
	//		std::fstream outfile(filename, std::fstream::out);
	//		sprintf(filename, "%s_summary.txt", systname);
	//		std::ofstream summary(filename);

	//		std::vector<GeometryVector> basis;
	//		for (int j = 0; j < dim;j++) {
	//			basis.push_back(GeometryVector());
	//			basis[j].SetDimension(dim);
	//			basis[j].x[j] = 1.0;
	//		}

	//		for (size_t N = N1; N < N2; N++)
	//		{
	//			RandomGenerator g(1234);
	//			double delta = 0.2;
	//			GeometryVector temp;	temp.SetDimension(dim);
	//			Configuration c, temp_c;
	//			if (i == 0) {
	//				c = GetLattice_config(basis, N);
	//				int j = (int)floor(N*g.RandomDouble());
	//				RandomVector(dim, temp);
	//				temp = delta*temp;
	//				c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
	//				char name[200] = {};
	//				//sprintf(name, "CC_Inv_Log_1PL_d_%d_NKSkip_%d___%d", dim, NKSkip,N);
	//				//WriteConfiguration(c, name);
	//			}
	//			else if (i == 1) {
	//				c = GetLattice_config(basis, N);
	//				for (int j = 0; j < N; j++) {
	//					//int j = (int)floor(N*g.RandomDouble());
	//					RandomVector(dim, g, temp);
	//					temp = delta*temp;
	//					c.MoveParticle(j, c.CartesianCoord2RelativeCoord(c.GetCartesianCoordinates(j) + temp));
	//				}
	//				char name[200] = {};
	//				//sprintf(name, "CC_Inv_Log_PL_d_%d_NKSkip_%d___%d", dim, NKSkip, N);
	//				//WriteConfiguration(c, name);
	//			}
	//			else if (i == 2) {
	//				c = GeneratePoisson(dim, N, RandomGenerator(1234));
	//				c.Resize(N);
	//			}
	//			else if (i == 3) {
	//				temp_c = GetFibonacciChain(1, 100);
	//				temp = temp_c.GetCartesianCoordinates(N-1);	temp.x[0] += 0.5;
	//				c = Configuration(dim, &temp, 1);
	//				for (int n = 0; n < N; n++)
	//					c.Insert("a", c.CartesianCoord2RelativeCoord(temp_c.CartesianCoord2RelativeCoord(n)));
	//			}

	//			outfile << "testing Nk for N=" << N;
	//			outfile << std::endl;
	//			for (size_t NumConst = 1;; NumConst++)
	//			{
	//				outfile << "NumConst=" << NumConst << ", testing...";
	//				outfile.flush();

	//				std::vector<double> results;
	//				bool result = FrankInvCCUniqueTest(c, Nc, NumConst, results, outfile, 1e-4, 0, false, "NLOPT+MINOP");

	//				if (result)
	//				{
	//					outfile << ", Configuration Unique.\n";
	//					std::cout << N << " " << NumConst << std::endl;
	//					for (int it = 0; it < results.size(); it++)
	//						summary << results[it] << "\t";

	//					summary << "\n";
	//					break;
	//				}
	//				else
	//					outfile << ", Configuration not Unique. Increasing Nk.\n";
	//			}
	//		}
	//	}
	//}
	//	
	//

	//	

	//		return 0;
		










	
	






//	char dir[100] = "results//MC//";
//	//molecular dynamics of harmonic lattice.
//
//
//	//------------ system parameters ----------------
//	DimensionType d = 1;
//	int L = 50; //systSize;
//	int N = pow(L, d);
//	double M = 1;//particle mass
//	double K = 1;// spring constant
//	double cL = sqrt(K / M);
//	double a = 1;//lattice constant
//
//	double w_max = 2.0*cL / a;
//	double w_min = cL*2.0*pi / L;
//	double dt = 0.01* 2.0*pi / w_max; //time step
//	double T_f = 3.0*2.0*pi / w_max;
//	double t0 = 2.0*pi / w_min;
//	int n_t = (int)round(T_f / dt); //recording step
//
//	//------------ simulation parameters ------------
//	int numTimestep = 10000;
//	int numSamplingStep = 1000;
//	int numTimeTest = 1000;
//	int numEquilIntv = 1000;
//	int NumThreads = 1;
//	double accRateL = 0.45;
//	double accRateU = 0.6;
//	double coolingRate = 0.95;
//
//	double dk = 2.0*pi / L;
//	double k_max = 3.0*pi;
//
//	double temp0 = 0.5; //near melting point
//	double temp1 = 0.05; //first temperature to investigate
//	double temp2 = 0.04; //final temperature to investigate
//	double d_temp = 0.02;
//	double kT = temp0; //temperature of the systems
//
//	double dx = 6.0* sqrt(kT);
//
//
//	//GeometryVector * basis = new GeometryVector[d];
//	std::vector<GeometryVector> basis;
//	for (int i = 0; i < d; i++) {
//		basis.push_back(GeometryVector(static_cast<int>(d)));
//		basis[i].x[i] = a;
//	}
//	Configuration prtConfig = GetLattice_config(basis, L);
//	//-----------prepare an ensemble
//	std::vector<SimpleHarmonic> ensemble;
//	std::vector<RandomGenerator> g_set;
//
//	SimpleHarmonic temp(d, K, a);
//	temp.SetReference(prtConfig);
//	temp.SetConfiguration(prtConfig);
//
//	for (int i = 0; i < NumThreads; i++) {
//		g_set.push_back(RandomGenerator(std::time(NULL) + i));
//		ensemble.push_back(temp);
//		//give initial sinusoidal displacements
//		//for (int l = 0; l < d; l++) { //for each direction
//		//      GeometryVector A;
//		//      GeometryVector phi;
//		//      GeometryVector q(static_cast<int> (d));
//		//      GeometryVector du(static_cast<int> (d));
//		//      for (int j = 0; j < d; j++) {
//		//              A.x[j] = g_set[i].RandomDouble()*sqrt(kT)/(double)N;
//		//              phi.x[j] = 2.0*pi*g_set[i].RandomDouble();
//		//      }
//		//      q.x[l] = 2.0*pi / (double)L;
//
//		//      for (int j = 0; j < N; j++) {
//		//              //compute the displacement
//		//              for (int k = 0; k < d; k++) {
//		//                      du.x[k] = A.x[k] * sin(phi.x[k] + q.Dot(ensemble[i].confi.GetCartesianCoordinates(j)));
//		//              }
//		//              ensemble[i].updateUn(j, du);
//		//      }
//		//}
//	}
//	std::cout << "end initialization \n";
//
//
//
//	omp_set_num_threads(NumThreads);
//	for (double targetkT = temp1; targetkT >= temp2; targetkT -= d_temp) {
//		//-----------cooling from the initial temperature to the first temperature to investigate
//		int numOfCoolingSteps = (int)round(log(targetkT / kT) / log(coolingRate));
//		coolingRate = pow(targetkT / kT, 1.0 / (double)numOfCoolingSteps);
//		std::cout << "target kT = " << targetkT << "\tcooling rate = " << coolingRate << "\t the number of cooling steps = " << numOfCoolingSteps << "\n";
//		for (int i = 0; i <= numOfCoolingSteps; i++) {
//			std::cout << "kT = " << kT << "\n";
//			SetTrialMovement(ensemble, g_set, dx, kT, accRateL, accRateU, numTimeTest);
//			std::cout << "-----finish the adjustment \n";
//
//			double rate = 0;
//			for (int j = 0; j < numEquilIntv; j++) {//test to whether they are in equilibrium
//				rate = MCevolve(ensemble, g_set, kT, dx);
//				//      std::cout << rate << "\t";
//			}
//			std::cout << "\n\n";
//
//			kT *= coolingRate;
//			if (kT < temp1)
//				kT /= coolingRate;
//		}
//		kT = targetkT; //....
//
//					   //
//		std::cout << "more equilibration \n";
//		for (int i = 0; i < numTimestep; i++) {
//			MCevolve(ensemble, g_set, kT, dx);
//		}
//
//#pragma omp parallel for schedule (guided)
//		for (int i = 0; i < ensemble.size(); i++) {
//			ensemble[i].updateConfiguration();
//		}
//
//		////------energy distribution--------
//		//double dE = 0.05*kT;
//		//double E_max = 10.00*kT;
//		//std::vector<double> E_dist((int)ceil(E_max / dE));
//		//double count = 0;
//
//		//std::ostringstream filen;
//		////std::string filen = "results\\MC\\energy_DIST L_" +  + "kT_ " +  +".txt";
//		//filen << dir << "energy_DIST " << d << "D L_" << L << " kT_ " << kT << ".txt";
//		//std::cout << filen.str() << "\n";
//		//std::ofstream fname(filen.str());
//		//filen.flush();
//		char prefix[100] = {};
//				
//
//		auto getConfig = [&ensemble, numSamplingStep, &g_set, kT, d, dx, &dir, NumThreads](size_t i)->Configuration {
//			if (i % ensemble.size() == 0) {//evolve configurations
//				for (int j = 0; j<numSamplingStep; j++) {
//					MCevolve(ensemble, g_set, kT, dx);
//				}
//#pragma omp parallel for schedule(guided)
//				for (int j = 0; j < ensemble.size(); j++) {
//					ensemble[j].updateConfiguration();
//					//-------------sampling energy distribution-------
//					//for (int k = 0; k < ensemble[j].confi.NumParticle(); k++) {
//					//      double xx = ensemble[j].Energy(k);
//					//      int index = (int)floor(xx/dE);
//					//      if (index > 0 && index < E_dist.size()) {
//					//#pragma omp critical //sampling energy
//					//              {
//					//                      E_dist.at(index)++;
//					//                      count++;
//					//              }
//
//					//      }
//					//}
//				}
//
//
//
//			}
//			return ensemble[i%NumThreads].confi;
//		};
//
//
//		std::vector<GeometryVector> result;
//
//		
//		omp_set_num_threads(6);
//
//		//IsotropicStructureFactor(getConfig, numTimestep, k_max, k_max, result, dk);
//		//WriteFunction(result, "test");
//
//		char fileName[100] = {};
//		for (int z = 0; z < 8; z++) {
//			sprintf(fileName, "%sconfig_%dDthermalXry L_%d kT_%0.6f numSamp_%d", dir, d, L, kT, z);
//			Configuration temp = getConfig(z);
//			WriteConfiguration(temp, fileName);
//
//		}
//
//		//for (int i = 0; i < E_dist.size(); i++) {
//		//      fname << dE*((double)i+0.5) << "\t" << E_dist[i] / (double(count) *dE) << "\n";
//		//}
//
//		//fname.close();
//		//WriteFunction(result, "results\\MC\\Sk_test");
//		//filen << dir << "Sk "<<d<< "D L_" << L << " kT_" << kT << " numSamp_" << numTimestep;
//		//WriteFunction(result, fileName);
//
//	}
//
//	return 0;
//---------------------------------------------------------------
	//for (int z = 1; z < 5; z++) {
	//	double kT = 0.0001*z;

	//	SimpleHarmonicF test(d, K, a);
	//	test.SetReference(prtConfig);
	//	test.UpdateConfiguration(0);
	//	omp_set_num_threads(NumThreads);

	//	//average of u_q distribution
	//	double u_qCount = 0;
	//	double qMax = pi;
	//	double dq = 2.0*pi / (double)L;
	//	std::vector<GeometryVector> ku2_qDist((int)round(qMax / dq), GeometryVector(static_cast<int>(4)));


	//	double dE = 0.05*kT;
	//	double E_max = 10.00*kT;
	//	double dx = 10.0*sqrt(kT / K);
	//	std::vector<double> E_dist((int)ceil(E_max / dE));
	//	double E_count = 0;


	//	std::vector<double> omega;
	//	for (int i = 0; i < test.q.size(); i++) {
	//		double temp = d;
	//		for (int j = 0; j < d; j++) {
	//			temp -= cos(test.q[i].x[j] * a);
	//		}
	//		omega.push_back(sqrt(2.0*K / M *temp));
	//	}

	//	//std::vector<double> abs_q; //omega...
	//	//for (int i = 0; i < test.q.size(); i++) {
	//	//	abs_q.push_back(sqrt(test.q[i].Modulus2()));
	//	//}
	//	GeometryVector du(static_cast<int>(d));
	//	GeometryVector phi_trial(static_cast<int>(d));
	//	//double energy_orig = test.Energy();
	//	//double energy_trial = 0;
	//	//double dE = 0;
	//	//double Pr = 0;
	//	GeometryVector Edata(4);//etcdata.x[0] = energy_orig, etcdata.x[1] = energy_trial, etcdata.x[2]=dE, etcdata.x[3] =Pr
	//	Edata.x[0] = test.Energy();
	//	RandomGenerator g;
	//	//for (int t = 0; t < 20; t++) {
	//	//	//double rate = 0;
	//	//	for (int i_q = 0; i_q < test.q.size(); i_q++) {
	//	//		//generate a trial coordinate
	//	//		for (int j = 0; j < d; j++) {
	//	//			du.x[j] = dx / sqrt(abs_q[i_q])*(g.RandomDouble() - 0.5);
	//	//			phi_trial.x[j] = 2.0*pi*(g.RandomDouble());
	//	//		}
	//	//		test.update_u_trial(i_q, du, phi_trial);
	//	//		Edata.x[1] = test.Energy_trial();
	//	//		Edata.x[2] = Edata.x[1]-Edata.x[0];
	//	//		Edata.x[3] = g.RandomDouble();

	//	//		if (Edata.x[2] < 0 || exp(-Edata.x[2] / kT) > Edata.x[3]) {//Metropolitan 
	//	//			test.update_u_q(i_q, du, phi_trial);
	//	//			test.update_un();
	//	//			Edata.x[0] = Edata.x[1];
	//	//			//rate++;
	//	//		}
	//	//	}
	//	//}
	//	//test.UpdateConfiguration(1);

	//	auto GC = [&test, &omega, &ku2_qDist, &E_dist, &dx, &g, &Edata, &kT, &dE, &dq](size_t i) ->Configuration {
	//		int d = test.confi.GetDimension();
	//		GeometryVector du(static_cast<int>(d));
	//		GeometryVector phi_trial(static_cast<int>(d));
	//		double rate = 0;

	//		for (int i_q = 0; i_q < test.q.size(); i_q++) {
	//			//generate a trial coordinate
	//			for (int j = 0; j < d; j++) {
	//				du.x[j] = (dx / omega[i_q])*(g.RandomDouble() - 0.5);
	//				phi_trial.x[j] = 2.0*pi*(g.RandomDouble());
	//			}
	//			test.update_u_trial(i_q, du, phi_trial);
	//			Edata.x[1] = test.Energy_trial();
	//			Edata.x[2] = Edata.x[1] - Edata.x[0];
	//			Edata.x[3] = g.RandomDouble();
	//			//std::cout << Edata;
	//			if (Edata.x[2] < 0 || exp(-Edata.x[2] / kT) > Edata.x[3]) {//Metropolitan 
	//				test.update_u_q(i_q, du, phi_trial);
	//				test.update_un();
	//				Edata.x[0] = Edata.x[1];
	//				rate++;
	//			}
	//		}

	//		///std::cout << "acceptance rate = " << rate / (double)test.q.size()<<"\n";

	//		if (i >= 0) {
	//			double xx = 0;	int index = 0;
	//			//obtain the energy distribtion of individual particles
	//			for (int n = 0; n < test.confi.NumParticle(); n++) {
	//				xx = test.Energy(n);
	//				index = (int)floor(xx / dE);
	//				if (index >= 0 && index < E_dist.size()) {
	//					E_dist.at(index)++;
	//				}
	//			}

	//			for (int n = 0; n < test.q.size(); n++) {
	//				xx = test.q[n].Dot(test.u_q[n]);
	//				xx = xx*xx;
	//				index = (int)round(sqrt(test.q[n].Modulus2()) / dq);
	//				ku2_qDist[index].x[2]++;
	//				ku2_qDist[index].x[0] = ku2_qDist[index].x[0] + (xx - ku2_qDist[index].x[0]) / ku2_qDist[index].x[2];
	//				ku2_qDist[index].x[1] = ku2_qDist[index].x[1] + (xx*xx - ku2_qDist[index].x[1]) / ku2_qDist[index].x[2];
	//			}
	//		}

	//		


	//		if (i>0 && i % 100 == 0) {
	//			test.UpdateConfiguration(1); //recompute u_n
	//		}
	//		else {
	//			test.UpdateConfiguration(0);
	//		}

	//		return test.confi;

	//	};


	//	for (int i = -500; i < 0; i++) {
	//		GC(i);
	//	}

	//	std::vector<GeometryVector> result;
	//	char name[100] = {};
	//	IsotropicStructureFactor(GC, numTimestep, 0.5*pi, 0.5*pi, result, pi / (double)L);
	//	sprintf_s(name, 100, "%sSk %dDSHL L_%d numSamp_%d", dir, d, (int)L, numTimestep);
	//	WriteFunction(result, name);
	//	result.resize(0);

	//	double count = 0;
	//	for (int i = 0; i < E_dist.size(); i++) {
	//		count += E_dist[i];
	//	}
	//	for (int i = 0; i < E_dist.size(); i++) {
	//		result.push_back(GeometryVector((i + 0.5)*dE, E_dist[i] / (count*dE)));
	//	}
	//	sprintf_s(name, 100, "%sEdist %dDSHL L_%d numSamp_%d", dir, d, (int)L, numTimestep);
	//	WriteFunction(result, name);

	//	result.resize(0);
	//	for (int i = 0; i < ku2_qDist.size(); i++) {
	//		ku2_qDist[i].x[2] = sqrt((ku2_qDist[i].x[1] - ku2_qDist[i].x[0] * ku2_qDist[i].x[0]) / (double)ku2_qDist[i].x[2]);
	//		ku2_qDist[i].x[1] = ku2_qDist[i].x[0];
	//		ku2_qDist[i].x[0] = dq*(double)i;
	//	}
	//	ku2_qDist.erase(ku2_qDist.begin());
	//	sprintf_s(name, 100, "%sqAq2Average %dDSHL L_%d numSamp_%d", dir, d, (int)L, numTimestep);
	//	WriteFunction(ku2_qDist, name);

	//}
	//

	//return 0;
	






//

	//**********************somewhat molecular dynamics part
	
	//auto inverseC_energy = [kT](double prob) -> double {return -2.0*kT*log(1.0 - prob); };

	////little random perturbation
	////for (int i = 0; i < prtConfig.NumParticle(); i++) {
	////	prtConfig.MoveParticle(i, prtConfig.CartesianCoord2RelativeCoord(prtConfig.GetCartesianCoordinates(i) + GeometryVector(0.1*(g.RandomDouble()-0.5) )));
	////}

	//
	//std::vector<GeometryVector> E_record;

	//ttest.SetConfiguration(prtConfig);
	//dt = 0.0001;
	//ParticleMolecularDynamics chain(prtConfig, dt, M);
	//chain.InitNoseHoover(M*40);
	////chain.FindOptimalTimeStep(ttest, dt, 0.1, 10);
	//std::cout << dt << "\n";
	//chain.TimeStep = dt;
	////initial velocities
	//std::vector<GeometryVector> v0;
	//v0.reserve(L);
	//for (int i = 0; i < L; i++) {
	//	double abs_v = sqrt(2.0 / M *g.RandomDouble(inverseC_energy));
	//	v0.push_back(GeometryVector(static_cast<double>(g.RandomDouble_normal(sqrt(kT) ))));
	//}

	//for (int i = 0; i < v0.size(); i++)
	//	std::cout << v0[i].x[0] << "\n";

	//chain.Velocity = std::vector<GeometryVector> (v0); 
	//ttest.SetConfiguration(prtConfig);
	//chain.Position = prtConfig;

	////velocity distribution
	//double dv = 2.0;
	//double v_max = 10.0*sqrt(kT/M);
	//std::vector<double> v_dist((int)ceil(2*v_max/dv));
	//int count = 0;

	//std::ofstream fname("results\\MD\\energy");
	//std::cout << "Energy = " << chain.GetKineticEnergy() << "\n";

	//char name[100] = {};
	//int numTimestep = (int)floor(2.0*t0 / dt);
	//for (int i = 0; i < numTimestep; i++) {
	//	double t = dt*(double)i/t0;
	//	//sprintf_s(name, 100, "results\\MD\\config_t_%d", i);
	//	//WriteConfiguration(chain.Position, name);

	//	//measure velocity distribution********************
	//	for (int j = 0; j < chain.Position.NumParticle(); j++) {
	//		double xx = chain.Velocity[j].x[0];
	//		int index = (int)floor((xx + v_max) / dv);
	//		if (index >= 0 && index < v_dist.size()) {
	//			count++;
	//			v_dist.at(index)++;
	//		}
	//	}
	//	//***************************************

	//	double PE = chain.GetPotentialEnergy(ttest);
	//	double KE = chain.GetKineticEnergy();
	//	//E_record.push_back(GeometryVector(t, KE, PE, KE+PE));
	//	fname << GeometryVector(t, KE, PE, KE + PE);
	//	chain.NoseHooverEvolve(10, ttest, kT);
	//	
	//	//*************************

	//	if (i % 10000 == 0) {
	//	std::cout << i << "\n";

	//		std::vector<GeometryVector> normalized_dist;
	//		for (int j = 0; j < v_dist.size(); j++) {
	//			normalized_dist.push_back(GeometryVector((j + 0.5)*dv - v_max, v_dist[j] / (dv*(double)count)));
	//		}
	//		sprintf_s(name, 100, "results\\MD\\v_distribution kT_%0.1f %0.4f", kT, t);
	//		WriteFunction(normalized_dist, name);
	//		//std::cout << "\n";
	//	}
	//		
	//	



	//}
	//
	////sprintf_s(name, 100, "results\\MD\\energy");
	////WriteFunction(E_record, name);
	//std::cout << "end";

	//return 0;




	////More realistic simulations for the thermally excited lattices
	//int numConfig = 10000;
	//int numExcitedModes = 10000;
	//double c = 1.0; //speed of sound
	//double m = 1.0; //mass of a particle
	//double energyPrefactor = 0.5*m*c*c;//*(4.0*pi*pi / (10000 * 10000));
	//double phaseFactor = 2.0*pi;
	////int numSampling = 10000;

	//double rho = 1.0;
	//std::vector<GeometryVector> b;
	//b.push_back(GeometryVector(static_cast<double> (1.0)));


	//for (int L = 10000; L <= 10000; L *= 10) {
	//	std::vector<double> q;
	//	for (int i = 0; i < numExcitedModes; i++) {
	//		q.push_back(2.0*pi*(i+1.0) / (double)L);
	//	}
	//	
	//	Configuration tem = GetLattice_config(b, L);
	//	char name[100] = {};

	//	double dk = 2.0*pi / (double)L;
	//	double kMax = 2.0*pi;
	//	std::vector<GeometryVector> kT_S0;
	//	for (double kT = 10; kT < 100; kT += 10) {
	//		if (kT == 10 || kT==50) {
	//			kT += 10;
	//		}
	//		auto inverseC_energy = [kT](double prob) -> double {return -log(1.0 - prob) * kT; };

	//		double sum = 0;
	//		double sum2 = 0;
	//		RandomGenerator g(std::time(NULL));
	//		double * E = new double[q.size()*numConfig];
	//		//prepare energies...
	//		int index = 0;
	//		for (int i = 0; i < numConfig; i++) {
	//			for (int j = 0; j < q.size(); j++) {
	//				E[index] = g.RandomDouble(inverseC_energy);
	//				sum += E[index];
	//				sum2 += E[index] * E[index];
	//				index++;
	//			}
	//		}
	//		double average = sum / (double)(numConfig*q.size());
	//		double sE = sqrt((sum2 / (double)(numConfig*q.size()) - average*average) / (double)(numConfig*q.size()));
	//		std::cout << "energy = " << average << " + " << sE << "\t kT = " << kT << "\t 1/energy = " << 1.0 / average << "\n";

	//		auto getConf = [&tem, &g, &q, &E, &inverseC_energy, energyPrefactor, phaseFactor](size_t j)->Configuration {
	//			Configuration confi(tem);
	//			
	//			std::vector<double> A; //amplitude
	//			std::vector<double> phi; //phase constant
	//			
	//			
	//			for (int i = 0; i < q.size(); i++) {
	//				A.push_back(sqrt(E[i+q.size()*j]/(energyPrefactor*(i+1)*(i+1))));
	//				phi.push_back(phaseFactor*(g.RandomDouble())); //for phase shift
	//				//phi.push_back(pi* floor(g.RandomDouble()));
	//			}

	//			//move particles
	//			
	//			for (int i = 0; i < confi.NumParticle(); i++) {
	//				GeometryVector displacement(static_cast<int> (confi.GetDimension()));
	//				for (int j = 0; j < q.size(); j++) {
	//					displacement.x[0] = displacement.x[0]+ A[j] * sin(q[j] * (double)i + phi[j]);
	//				}

	//				confi.MoveParticle(i, confi.CartesianCoord2RelativeCoord(confi.GetCartesianCoordinates(i) + displacement));
	//			}
	//			return confi;};
	//		std::vector<GeometryVector> sk_result;

	//	//	for (int i = 0; i < numConfig; i++) {
	//	//		getConf(i);
	//	//	}
	//		//IsotropicStructureFactor(getConf, numConfig, 40 * dk, 40 * dk, sk_result, dk);
	//		
	//		IsotropicStructureFactor(getConf, numConfig, 1.2*kMax, 1.2*kMax, sk_result, dk);

	//		sprintf_s(name, 100, "results\\Sk_Z1_sinusoidal perturbation lowest_%dk A2_%0.2f N_%d", q.size(), kT, L);
	//		WriteFunction(sk_result, name);
	//		//WriteConfiguration(conf, name);
	//		//kT_S0.push_back(GeometryVector(average, sk_result[0].x[1], sE , sk_result[0].x[3]));
	//		kT_S0.push_back(GeometryVector(2.0*average, sk_result[0].x[1], 2.0*sE, sk_result[0].x[3]));
	//		sprintf_s(name, 100, "results\\A2_S(0) sinusoidla perturbations lowest_%dk N_%d", q.size(),L);
	//	}

	//	WriteFunction(kT_S0, name);
	//}




	//return 0;




//	//return 0;
//

	










//	//parameter sets
//	//-----------displacement size
//	double ratio_i = 0.1;
//	double ratio_f = 1.0;
//	double c_ratio = 0.2;
//	double dx = 0;
//	//------------system size
//	double phi = 0;
//	double phi_c = 0.91;
//	double dphi = 0.05;
//	double phi_i = 0.85;//phi_c - 3.0*dphi;
//	double phi_f = 0.95;//phi_c + 3.0*dphi;
//	int N = 0; //the number of particles
//	int N_i = 10000;
//	int N_f = 10000;
//	int dN = 10;
//	double L = 0; //box size
//	int numConfig = 100;
//	//-------------time
//	int t_max = 2147483647; //the max of integer
//	int dt = 1000; //save_interval
//	//-------------structure factor
//	double dk = 0;
//	double k_max = 0;
//
//	for (double ratio = ratio_i; ratio <= ratio_f; ratio *= c_ratio) {
//		
//		for (N = N_i; N <= N_f; N *= dN) {
//			L = N;
//			for (phi = phi_i; phi <= phi_f; phi += dphi) {
//				double l = 0.5*phi; //radius of a particle
//				double dx = ratio*l;
//				dk = 3.0*2.0*pi / L;
//				k_max = 3.0*pi ;
//
//				GeometryVector basis[1] = {GeometryVector(static_cast<double> (L))};
//				std::cout << "N = " << N << "\t phi = " << phi << "\t delta = "<<ratio<<"\n";
//				std::vector<Configuration> configList;
//				//generate pt configurations
//#pragma omp parallel for schedule(guided,1)
//				for (int i = 0; i < numConfig; i++) {
//					Configuration temp(1, basis, 1, true);
//					RandomGenerator g(std::time(NULL) + i);
//					for(int num = 0; num<N; num++){ temp.Insert("a", g); }
//#pragma omp critical
//					configList.push_back(temp);
//					
//				}
//
//
//				char post_fix[100] = {};
//				sprintf_s(post_fix, 100, "delta_%0.4f phi_%0.6f N_%d", ratio, phi, N);
//				
//				char f_name[100] = {};
//				sprintf_s(f_name, 100, "f_active %s.txt", post_fix);
//				std::ofstream f_active(f_name);
//
//				std::vector<int> f_list(numConfig,1);
//				int numActive = 1;
//				int t = 0;
//				for (t = 0; t < t_max && numActive != 0; t++) {
//					double f_sum = 0;
//					double f_sum2 = 0;
//					numActive = 0;
//#pragma omp parallel for schedule(guided,1 )
//					for (int i = 0; i < numConfig; i++) {
//
//						if (f_list[i] != 0) {
//							std::vector<size_t> activeIndex;
//							IsActive(configList[i], activeIndex, l);
//
//							if (t%dt == 0) {
//								f_sum += ((double)activeIndex.size() / N);
//								f_sum2 += ((double)activeIndex.size()*activeIndex.size() / (N*N));
//							}
//							f_list.at(i) = activeIndex.size();
//							numActive += activeIndex.size();
//							B_model(configList[i], dx, activeIndex, RandomGenerator(std::time(NULL) + i));
//							
//#pragma omp critical
//							{
//								if (activeIndex.size() == 0) { //becomes inactive at the first time
//									sprintf_s(f_name, 100, "config_1DRO\\config_%s_%d", post_fix, i);
//									WriteConfiguration(configList[i], f_name);
//								}
//							}
//							
//						}
//
//	
//					}
//
//					if (t%dt == 0) {//compute f_active
//						f_sum = (double)f_sum / numConfig;
//						f_sum2 = (double)f_sum2 / numConfig - f_sum*f_sum;
//						//std::cout << t<<"\t"<<numActive << "\n";
//						f_active << t << "\t" << f_sum << "\t" << sqrt(f_sum2/numConfig) << "\n";
//						//if (t % (10 * dt) == 0) {
//						//	std::cout<< t << "\t" << f_sum << "\t" << sqrt(f_sum2 / numConfig) << "\n";
//						//}
//					}
//
//				}
//
//				
//				f_active.close();
//				std::cout << "dynamics ends after t="<<t<<"\n";
//				//compute the structure factor
//				int numInactive = 0;
//				for (int i = 0; i < numConfig; i++) {
//					if (f_list[i] == 0) {
//						numInactive++;
//					}
//				}
//				if (numInactive != numConfig && numInactive > 0) {//there are some active systems => exclude active systems
//					configList.resize(0);
//					configList.reserve(numInactive);
//					for (int i = 0; i < numConfig; i++) {
//						std::vector<GeometryVector> tt;
//						sprintf_s(f_name, 100, "config_1DRO\\config_%s_numInact_%d", post_fix, numInactive);
//						ReadFunction(tt, f_name);
//						
//						Configuration temp(1, basis, 10);
//						for (int j = 0; j < tt.size(); j++) {
//							temp.Insert("a", tt[j]);
//						}
//						configList.push_back(temp);
//					}
//				}
//
//				auto getCONFIG = [&configList](size_t i)-> Configuration {return configList[i];};
//				std::vector<GeometryVector> Sk;
//				IsotropicStructureFactor(getCONFIG, numConfig, k_max, k_max, Sk, dk);
//				sprintf_s(f_name, 100, "Sk_%s", post_fix);
//				WriteFunction(Sk, f_name);
//			}
//
//		}
//
//		
//	}
//
//	//generate configration
//
//
//
//
//
//
//	return 0;







	//GeometryVector q = 2.0*pi / systSize * GeometryVector(static_cast<double> (1)); //wavenumber/wavevector for the phonon.
	//double dk = sqrt(temp.GetReciprocalBasisVector(0).Modulus2())*2.0;
	//double kMax = 3.0*2.0*pi/sqrt(basis[0].Modulus2());

	//char name[100] = {};

	//
	//auto f = [](GeometryVector x)->double {return sqrt(x.Modulus2());};
	//auto win = [d](double x)->double {return 2 * x;};
	//
	//std::vector<GeometryVector> samPrt;
	//RandomGenerator g;
	//for (int i = 0; i < numSampling; i++) {
	//	samPrt.push_back(GeometryVector(static_cast<double> (g.RandomDouble())));
	//}


	////amplitude of the phonon
	//double A = 0.5;
	//for (double j = 1; j < 20; j+= 2) {
	//	GeometryVector Q = (double)j*q;
	//	auto GetConfig = [&temp, Q, A](size_t i) -> Configuration {
	//		Configuration result(temp);
	//		RandomGenerator g(std::time(NULL) + i);

	//		double phi = 2.0*pi*g.RandomDouble(); //random phase of the density wave
	//		for (int i = 0; i < result.NumParticle(); i++) {
	//			result.MoveParticle(i, result.CartesianCoord2RelativeCoord(result.GetCartesianCoordinates(i) + A*cos(phi + Q.Dot(result.GetCartesianCoordinates(i)))));
	//		}
	//		return result;
	//	};
	//	/*for (int i = 0; i < 4; i++) {
	//	sprintf_s(name, 100, "config_phonon wavenumber_%3.3f __%d", sqrt(q.Modulus2()), i);
	//	WriteConfiguration(GetConfig(i), name);
	//	}*/
	//	
	//	std::vector<GeometryVector> sigma2;

	//	/*MC_nv(GetConfig, numConfigs, sigma2, f, Rmax, dR, samPrt, 4, win, 1.0);
	//	sprintf_s(name, 100, "MC_nv Z1 wavenumber_%3.3f Amplitude_%1.1f", sqrt(q.Modulus2()), A);
	//	WriteFunction(sigma2, name);
	//	std::cout << "A = " << A << " end\n";*/

	//	if (A >= 0) {
	//		IsotropicStructureFactor(GetConfig, numConfigs, kMax, kMax, sigma2, dk);
	//		sprintf_s(name, 100, "S_k Z1 wavenumber_%3.3f Amplitude_%1.1f", sqrt(Q.Modulus2()), A);
	//		WriteFunction(sigma2, name);
	//		std::cout << "A = " << A << " end\n";
	//	}	
	//}
	


	
	//----------------------------------------------------------------------
//	//Parameters for configs
//	std::vector<GeometryVector> bv;
//	bv.push_back(GeometryVector(1.0, 0));
//	bv.push_back(GeometryVector(cos(pi / 3.0), sin(pi / 3.0)));
//	double R0 = 0.5* sqrt(bv[0].Modulus2());
//	const int iterMax = 10;
//	const int iter_in_vv = 1;
//	double phi = 0.5; //volume fraction
//	double scales = pow(phi, 1.0 / 2.0);
//	double r0 = scales*R0; //the radius of the largest spheres
//	const DimensionType d = 2;
//
//	//Parameters for samplings
//	double Rmax = 500;
//	double dR = 0.5 / 100.0;
//
//	//Poisson sampling points
//	int numSampling = 10000;
//	std::vector<GeometryVector> samplPrt;
//	RandomGenerator g;
//	for (int i = 0; i < numSampling; i++) {
//		samplPrt.push_back(g.RandomDouble()*bv[0] + g.RandomDouble()*bv[1]);
//	}
//
//	//	int numPrt = (int)Lmax * 4;
//
//	//Generating a HS structure
//	std::vector < std::vector<size_t>> index;
//	std::vector<int> growingSize;
//	std::vector<double> Radius;
//	//Insert lattice particles
//	PeriodicCellList<double> temp_HS = GetLattice(bv, 2, R0);
//
//	std::vector<size_t> tempList;
//	for (int i = 0; i < 4; i++) {
//		Radius.push_back(R0);
//		index.push_back(tempList);
//	}
//	tempList.clear();
//	GeometryVector x0(0, 0);
//	double RR = iSoddyCircles(temp_HS.GetCartesianCoordinates(0), temp_HS.GetCartesianCoordinates(1), temp_HS.GetCartesianCoordinates(2), x0);
//	temp_HS.Insert(RR, temp_HS.CartesianCoord2RelativeCoord(x0));
//	tempList.push_back(0);
//	tempList.push_back(1);
//	tempList.push_back(2);
//	index.push_back(tempList);
//
//	tempList.clear();
//	RR = iSoddyCircles(temp_HS.GetCartesianCoordinates(3), temp_HS.GetCartesianCoordinates(1), temp_HS.GetCartesianCoordinates(2), x0);
//	temp_HS.Insert(RR, temp_HS.CartesianCoord2RelativeCoord(x0));
//	tempList.push_back(3);
//	tempList.push_back(1);
//	tempList.push_back(2);
//	index.push_back(tempList);
//
//	growingSize.push_back(4);
//	growingSize.push_back(6);
//
//	for (int i = 2; i < iterMax; i++) {
//		std::vector<size_t> temp;
//
//		for (int j = growingSize[i - 2]; j < growingSize[i - 1]; j++) {
//
//			for (int k = 0; k < 3; k++) {
//				//get indices of particles contact.
//				std::vector<size_t> ix;
//				for (int l = 0; l < 3; l++) {
//					if (k != l)
//						ix.push_back(index[j].at(l));
//					else
//						ix.push_back(j);
//				}
//				index.push_back(ix);
//				RR = iSoddyCircles(temp_HS.GetCartesianCoordinates(ix[0]), temp_HS.GetCartesianCoordinates(ix[1]), temp_HS.GetCartesianCoordinates(ix[2]), x0);
//				temp_HS.Insert(RR, temp_HS.CartesianCoord2RelativeCoord(x0));
//			}
//		}
//
//		growingSize.push_back(temp_HS.NumParticle());
//	}
//
//	//move a HS configuration into a unit cell
//	for (int i = 0; i < growingSize.size(); i++) {
//		growingSize.at(i) -= 3; //since we removed three particles of the first stage!
//	}
//
//	GeometryVector *v;
//	v = new GeometryVector[2];
//	v[0] = bv[0];
//	v[1] = bv[1];
//	SpherePacking HS(2, v, 0.5);
//	for (int i = 0; i < temp_HS.NumParticle(); i++) {
//		if (1>i || i >= 4)
//			HS.Insert(scales*temp_HS.GetCharacteristics(i), HS.CartesianCoord2RelativeCoord(temp_HS.GetCartesianCoordinates(i)));
//	}
//	temp_HS.RemoveParticles();
//		
//
//	//range functions
//	std::vector<GeometryVector> b;	//reciprocal basis vectors
//	for (int i = 0; i < HS.GetDimension(); i++) {
//		b.push_back(HS.GetReciprocalBasisVector(i));
//	}
//
//	//volume fraction taken by circles over 2nd generation
//	double vol = 0;
//	for (int i = 0; i < growingSize[iter_in_vv]; i++) {
//		vol += HyperSphere_Volume(HS.GetDimension(), HS.GetCharacteristics(i));
//	}
//
//	double b11 = b[0].Modulus2();
//	double b12 = b[0].Dot(b[1]);
//	double b22_s = sqrt(b[1].Modulus2());
//	auto n1 = [&b11](double R, const GeometryVector &x0) -> GeometryVector {
//		double mid = R*sqrt(b11);
//		return GeometryVector(ceil(x0.x[0] - mid), floor(x0.x[0] + mid));
//	};
//	auto n2 = [&b11, &b12, &b22_s](double R, int nx, const GeometryVector &x0, GeometryVector & n2) ->bool {
//		double term1 = nx*b12 / b11;
//		double term2 = R*R - nx*nx / b11;
//
//		if (term2 <0) { //out of the range
//			n2 = GeometryVector(0, 0);
//			return false;
//		}
//		else {
//			term2 = sqrt(term2)*b22_s;
//			n2 = GeometryVector(ceil(x0.x[1] + term1 - term2), floor(x0.x[1] + term1 + term2));
//			return true;
//		}
//	};
//	
//	std::vector<GeometryVector> result1;
//	//std::vector<GeometryVector> result2;
//	for (double R = dR; R <= Rmax; R += dR) {
//
//		std::vector<double> data;
//		data.reserve(numSampling);
//		double volumeWin = HyperSphere_Volume(d, R);
//#pragma omp parallel for schedule(guided, 10)
//		for (int i = 0; i < numSampling; i++) {
//			double sum = 0.0;
//			GeometryVector outer1 = n1(R+r0, samplPrt[i]);
//			GeometryVector inner1 = n1(R, samplPrt[i]);
//			GeometryVector outer2;
//			GeometryVector inner2;
//			int n1_min = (int)outer1.x[0];	int n1_max = (int)outer1.x[1];
//			int n1_inner_min = (int)inner1.x[0]; int n1_inner_max = (int)inner1.x[1];
//						
//			for (int n_1 = n1_min; n_1 <= n1_max; n_1++) {
//				n2(R + r0, n_1, samplPrt[i], outer2);
//				n2(R, n_1, samplPrt[i], inner2);
//			
//				for (int n_2 = outer2.x[0]; n_2 <= outer2.x[1]; n_2++) {
//
//					//within the safe zone (whole fundamental cell in included in the observation window)
//					if (n1_inner_min <= n_1 && n_1 < n1_inner_max && inner2.x[0] <= n_2 && n_2 < inner2.x[1]) {
//						sum += vol;
//					}
//					//away from the safe zone
//					else {
//						GeometryVector displacement(n_1*v[0]+n_2*v[1]-samplPrt[i]);// distance from the center of the window to the fundamental cell of interest.
//						//for (int j = 0; j < HS.NumParticle(); j++) {
//						for (int j = 0; j < growingSize[iter_in_vv-1]; j++) {
//							GeometryVector temp = displacement + HS.GetCartesianCoordinates(j);
//							sum += IntersectionVolume(d, HS.GetCharacteristics(j), R, sqrt(temp.Modulus2()));
//						}
//					}
//				}
//			}
//
//			data.push_back(sum/volumeWin - phi);
//
//		}
//		double average = 0;
//		double sum2 = 0;
//		for (int i = 0; i < data.size(); i++) {
//			average += data[i];
//			sum2 += data[i] * data[i];
//		}
//		average /= (double)data.size();
//		result1.push_back(GeometryVector(R, (sum2 / (double)(data.size()) -average*average)));
//		//result2.push_back(GeometryVector());
//	}
//
//	WriteFunction(result1, "MC_vv HS on triangular lattice");
//
//
//	//std::vector<GeometryVector> k;
//	//k = DirectionalSpectralDensity(getConfig, 1, 20, spe);
//	//WriteToContourf(k, spe, "directional spectral density_Hashin-Shtrikman_triangular lattice", 2 * 20 + 1);
//	return 0;






//--------------------------------------------------------------------


//

//double t_i = GetTickCount();

	
	//int numConfig = 40;
	//int numSamplingPrt = 1000;
	//double Rmax = 200;
	//double size = Rmax*2*4.5;
	//double dR = 0.03;
	//std::vector<int> key;	key.push_back(numSamplingPrt); key.push_back(2);
	//std::vector<GeometryVector> basis;
	//basis.push_back(GeometryVector(1,0));
	//basis.push_back(GeometryVector(0, 1));
	//Configuration Original = GetLattice_config(basis, size);
	//
	//
	//char file[100] = {};

	////Test with Poisson ....
	//GeometryVector test_basis[] = {GeometryVector(size, 0), GeometryVector(0,size)};
	//boost::random::mt19937 rnd(std::time(NULL));

	//auto conf_test = [size, &test_basis, rnd](size_t i)->Configuration {
	//	RandomGenerator g(std::time(NULL) + i);
	//	Configuration temp(2, test_basis, 10);
	//	int num = size*size;//PoissonRNG(rnd, size*size);
	//	for (int j = 0; j < num; j++) {
	//		//temp.DeleteParticle((int)floor(temp.NumParticle()*g.RandomDouble()));
	//		temp.Insert("a", g);
	//	}
	//	std::cout << num << std::endl;
	//	return temp;
	//};

	////WriteConfiguration(conf_test(0), "test_conf");
	//std::vector<GeometryVector> sampling_test;// vanDerCorputPts(numSamplingPrt, key);
	//GeometryVector temp;
	//RandomGenerator gen_test;
	//for (int i = 0; i < numSamplingPrt; i++) {
	//	RandomVector(2, gen_test, temp);
	//	//sampling_test.push_back((2.0 / size)*temp + (0.5 - (1.0 / size))*GeometryVector(1.0, 1.0));
	//	sampling_test.push_back(temp);
	//}

	//auto f = [](const GeometryVector x)->double {return sqrt(x.Modulus2());};
	//auto v = [](double x) ->double {return HyperSphere_Volume(2, x);};
	//std::vector<GeometryVector> result_test;
	//
	//MC_nv(conf_test, numConfig, result_test, f, Rmax, dR, sampling_test, 4, v);
	//sprintf_s(file, 100, "MC_nv_test numConfig_%d numSampling_%d  L_%d (grand canonical)", numConfig, sampling_test.size(), (int)size);
	//WriteFunction(result_test, file);
	//
	//return 0;
	

	//for (double L = 28.5; L < 40; L += 1.5) {
	//	double n_max = floor(size / L);
	//	double n = (n_max > 10) ? n_max : 10;
	//	std::cout << n << " x " << n << ", " << L <<" x " << L << std::endl;
	//	//basis vectors of the original crystals
	//	std::vector<GeometryVector> basis;
	//	basis.push_back(GeometryVector(1, 0));
	//	basis.push_back(GeometryVector(0, 1));

	//	//basis.push_back(GeometryVector(1, 0,0));
	//	//basis.push_back(GeometryVector(0, 1,0));
	//	//basis.push_back(GeometryVector(0, 0, 1));
	//	Configuration sample = GetLattice(basis, basis.size()*L);

	//	//size of domains
	//	std::vector<GeometryVector> b;
	//	b.push_back(GeometryVector(L, 0));
	//	b.push_back(GeometryVector(0, L));

	//	//b.push_back(GeometryVector(5, 0,0));
	//	//b.push_back(GeometryVector(0, 5,0));
	//	//b.push_back(GeometryVector(0, 0, 5));
	//	std::vector<Configuration> test;
	//	int numConfig = 0;

	//	if (n*n > 1000) {
	//		numConfig = 3;
	//	}
	//	else {
	//		numConfig = (int)floor(1000/(n*n));
	//	}

	//	for (int i = 0; i < numConfig; i++) {
	//		test.push_back(GetSimplePolyCrystals(sample, b, n));
	//	}

	//	std::vector<GeometryVector> S_k;
	//	//std::vector<GeometryVector> k;
	//	char file[100] = {};
	//	sprintf_s(file, 100, "Config L_%0.4f n_%0.1f", L, n);
	//	WriteConfiguration(test[0], file);
	//	auto getConf = [&test](size_t i)->Configuration {return test[i];};
	//	double k_max = 2 * pi ;
	//	double dk = 2*2 * pi / size;
	//	IsotropicStructureFactor(getConf, numConfig, 1.05*k_max, k_max, S_k, dk);
	//	
	//	sprintf_s(file, 100, "S_k L_%0.4f n_%0.1f", L, n);
	//	WriteFunction(S_k, file);
	//	std::cout << "end \n";
	//}
	
	//WriteConfiguration(test, "test_config");
	
	//PlotDirectionalStructureFactor2D(test, 2*n*L, "test_sk",true, 0.0,(0.1),(1.0));
	//k = DirectionalStructureFactor(test, 2* n*L, S_k);
	//WriteToContourf(k, S_k, "test", 4*n*L +1);
	
//	DimensionType d = 1;
//	size_t SystNum = 8;
//	for (int N = 2; N < 100; N *= 2) {
//		GeometryVector basis[1];
//		basis[0].Dimension = d;
//		basis[0].x[0] = (double)N;
//		std::vector<Configuration> config;
//		for (size_t i = 0; i < SystNum; i++) {
//			RandomGenerator gen;
//			Configuration temp(d, basis, 1);
//			for (size_t j = 0; j < N; j++) {
//				temp.Insert("a",GeometryVector(gen.RandomDouble()));
//			}
//
//			config.push_back(temp);
//		}
//		
//
//		std::cout << "N = " << N << "\n";
//	
//		double R = 0.4;
//		double delta = 0.01;
//		double dR = 0.01;
//		// true = active & false = passive
//		bool currentState = true;
//		bool lastState = true; 
//		while (fabs(dR) > 1e-10 && ~currentState) {
//
//			
//#pragma omp parallel for schedule(guided)
//			for (int i = 0; i < SystNum; i++) {
//				RandomGenerator rg;
//				IsActive(config[i], ,R);
//				B_model(config[i], R, ,rg);
//			}
//
//			if (currentState != lastState) {
//				dR *= -0.5;
//			}
//			R += dR;
//		}
//		
//
//	}
	
//	RandomGenerator G;
//#pragma omp parallel
//	{
//		RandomGenerator g(std::time(NULL)+omp_get_thread_num());
//#pragma omp for
//		for (int i = 0; i < 30; i++) {
//#pragma omp critical
//			std::cout << omp_get_thread_num()<<"\t"<<g.RandomDouble() << "\n";
//		}
//	}
//
//	return 0;
	

	

	
	
	//double i = 7;
	//double lis[] = { 3,4,6,7,8 };
	//for (int m = 3; m < 4; m++) {
	//	double L = lis[m];
	//	double n = Rmax*11.4;
	//	std::cout << "Rmax = " << Rmax <<"\n";
	//	GeometryVector N;
	//	for (int j = 0; j < 2; j++)
	//		N.x[j] = (n);
	//	N.SetDimension(2);
	//	char filename[100] = {};
	//	char file[100] = {};
	//	sprintf_s(file, 100, "angular_MC_nv_SC_superdisk_p__Lambda.txt");
	//	std::ofstream save(file);
	//	//auto getConfig = [&N, L, n](size_t a)->Configuration {return simplePolyCrystal(N, 0.5*n);};
	//	
	//	auto circle = [](const GeometryVector &x)->double {return sqrt(x.Modulus2());};
	//	DimensionType d = 2;
	//	std::vector<GeometryVector> b_test;
	//	b_test.push_back(GeometryVector(1, 0)); b_test.push_back(GeometryVector(0, 1.0));
	//	int NN = 1024*1024/8; //number of sampling points
	//	int numSamp = 1;
	//	
	//				
	//	Configuration sample = GetLattice(b_test, n);//getConfig(0);
	//	auto getConfig = [&sample](size_t i)->Configuration {return sample;};
	//	
	//	std::vector<GeometryVector> standardError;
	//	std::vector<GeometryVector> averageError;
	//	std::vector<GeometryVector> nv_result;
	//	std::vector<GeometryVector> p_Lambda;
	//	std::vector<int> basis;	basis.push_back(NN); basis.push_back(2);
	//	GeometryVector temp;
	//	std::vector<GeometryVector> PeriodicSampling; PeriodicSampling.reserve(NN);
	//	for (int i = 0; i < NN ; i++) {
	//		temp = vanDerCorputSequence(i, basis);
	//		PeriodicSampling.push_back(sample.CartesianCoord2RelativeCoord(temp));
	//	}

	//	

	//	for (double p = 1.5; p < 16; p *= 1.5) {
	//		auto win = [&p, &d](const GeometryVector &x)->double {
	//			double result = 0;
	//			for (int i = 0; i < d; i++) { result += pow(fabs(x.x[i]), 2.0*p); }
	//			return pow(result, 0.5 / p);};
	//		auto v = [&p, &d](double L) ->double {return Superdisk_vol(d, p, L);};

	//		/*auto win = [&d](const GeometryVector &x)->double {
	//			double result = 0;
	//			for (int i = 0; i < d; i++) { result = (result>fabs(x.x[i]))? result:fabs(x.x[i]); }
	//			return result;};
	//		auto v = [&d](double L) ->double {return pow(2.0*L,d);};*/

	//		angular_MC_nv(getConfig, numSamp, nv_result, win, Rmax, dR, PeriodicSampling, 3, v);

	//		sprintf_s(file, 100, "angular_MC_nv_SC_superdisk_%0.4f__%d", p, NN);
	//		//sprintf_s(file, 100, "test1");
	//		WriteFunction(nv_result, file);
	//		//ReadFunction(nv_result, file);
	//		std::vector<double> X;	std::vector<double> Y;
	//		std::vector<GeometryVector> coeff;


	//		for (int i = 0; i < nv_result.size(); i++) {
	//			X.push_back(nv_result[i].x[0]);
	//			Y.push_back(nv_result[i].x[2]);	//moving average
	//		}
	//		double squared_chi = polyfit(X, Y, 1, coeff);
	//		
	//		//p	\t	c1	\t delta c1	\t squared_chi
	//		GeometryVector dd(p, coeff[1].x[0], coeff[1].x[1], squared_chi);
	//		p_Lambda.push_back(dd);

	//		GeometryVector prt(dd);
	//		prt.x[1] /= Superdisk_surface_area(2, p, 1.0);
	//		std::cout << prt<<"\n" ;

	//		save << dd;
	//	
	//	}
	//	save.close();
	//	
	//	
	//}
	//
	
	//return 0;
	


/** \brief Auto-tune the trial-displacement size dx so the MC acceptance rate falls within a target window. (CLAUDE)
 *  \param pot          List of harmonic potentials, one per configuration.
 *  \param g            List of random number generators, one per configuration.
 *  \param dx           In/out maximum trial-displacement size, adjusted until acceptance is in range.
 *  \param T            Temperature of the Metropolis acceptance.
 *  \param lowAccept    Lower bound of the target acceptance rate.
 *  \param upperAccept  Upper bound of the target acceptance rate.
 *  \param numSteps     Number of trial moves per configuration used to estimate the acceptance rate. */
void SetTrialMovement(std::vector<SimpleHarmonic> &pot, std::vector<RandomGenerator>& g, double &dx, double T, double lowAccept, double upperAccept, int numSteps) {
	double acceptRate = 0.0;
	int numConfig = pot.size();
	int numPrt = pot[0].confi.NumParticle();
	DimensionType d = pot[0].confi.GetDimension();
	
	std::vector<double> count(pot.size(), 0);
		do {
			
			//---sampling the acceptance rate
#pragma omp parallel for schedule(guided)
			for (int i = 0; i < numConfig; i++) {
				GeometryVector trial(static_cast<int> (d));
				int index = 0;
				count[i] = 0;
				for (int l = 0; l < numSteps; l++) {
					//chose a trial displacement and the index of a particle.
					for (int j = 0; j < d; j++) { trial.x[j] = dx*(g[i].RandomDouble() - 0.5); }
					index = (int)floor(numPrt*g[i].RandomDouble()); 
					count[i] += MCevolve(pot[i], g[i], T, trial, index);
				}
			}
			for (int i=0; i<numConfig; i++)
			{
				acceptRate += count[i];
			}

			acceptRate = acceptRate / ((double) numConfig*numSteps);
			//---adjust the max-displacement
			std::cout << "rate = " << acceptRate << "\t dx = " << dx << "\n";;
			if (acceptRate > upperAccept) {
				dx *= 1.1;
			}
			else if (acceptRate < lowAccept) {
				dx *= 0.9;
			}
		
		} while (acceptRate > upperAccept || acceptRate < lowAccept);
	
	
}

