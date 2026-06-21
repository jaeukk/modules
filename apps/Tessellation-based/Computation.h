#ifndef TESS_COMPUTATION_H
#define TESS_COMPUTATION_H

/** \file Computation.h
 *  \brief The `Computation` task hierarchy used by the tiling-problem CLI.
 *
 *  Split out of the legacy monolithic TilingProblem.cpp.  `Computation` is an
 *  abstract task ("compute some statistic over a set of sphere packings, then
 *  write it"); the derived classes implement the individual analyses
 *  (spectral density chi_V, Voronoi cell-volume / vertex distributions, maximum
 *  packing fraction, local number / volume-fraction variance, centroids, the
 *  multiscale family, ...).  Also holds `ImLatGen` (imperfect-lattice generator).
 *
 *  The classes are header-defined (every method body is inline, exactly as in
 *  the legacy single-file design) so both the CLI driver (main.cpp) and
 *  Computation.cpp can instantiate them.
 */

#include "SpacefillingPacking.h"
#include "Voronoi.h"

#include <GenerateConfigs.h>
#include <NumberVariance.h>
#include <StructureFactor.h>

//==============================================================================
// Body sliced verbatim from the legacy TilingProblem.cpp (lines 345-1362):
// ImLatGen + the Computation base class and all of its derivatives.
//==============================================================================
//lattice with vacancies
struct ImLatGen :LatticeGen {
private:
	double p;		//mean fraction of vacancies
	int dist_type;	/* distribution of the number of vacancies
					 * 0: Binomial distribution, 1: Poisson distribution	*/
	int N0;			//the number of particles in initial packings
	std::vector<GeometryVector> getCubicBasis(int dimension) {
		std::vector<GeometryVector> b;
		for (int i = 0; i < dimension; i++) {
			b.emplace_back(GeometryVector(static_cast<int>(dimension)));
			b[i].x[i] = 1.0;
		}
		return b;
	}
public:
	ImLatGen(const std::vector<GeometryVector> & basis, int numPrts, double p, double numberdensity) : LatticeGen(basis, 1){
		SetBasis(basis);
		//if (this->maxPhi > targetPackingFraction)
		//	this->currPhi = targetPackingFraction;
		//else
		//	std::cout << "target packing fraction is too high\n";

		this->p = p;// 1.0 - this->currPhi / this->maxPhi;
		
		this->rho = numberdensity;
		//int LengInOrigLat = (int)round(pow(numPrts/(1.0-this->p),1.0/this->d));
		//int NumInOrigLat = (int)round(pow(LengInOrigLat, this->d));
		int LengInOrigLat = (int)round(pow(numPrts , 1.0 / this->d));
		int NumInOrigLat = (int)round(pow(LengInOrigLat, this->d));
		SetNumInEachSides(std::vector<size_t>(this->d, (size_t)LengInOrigLat));
		
		//readjustment
		this->N0 = this->numPrts;// *this->currPhi / this->maxPhi;
		//this->p = 1.0 - this->currPhi / this->maxPhi;
	}
	ImLatGen(int dimension, int numPrts, double p, double numberdensity) : ImLatGen(getCubicBasis(dimension), numPrts, p, numberdensity ){}
	//set the final packing fraction!!
	virtual void SetPackingFraction(double phi) {
		if (phi < this->maxPhi) {
			this->currPhi = phi;
			this->p = 1.0 - this->currPhi / this->maxPhi;
			int LengInOrigLat = (int)round(pow(this->N0 / (1.0 - this->p), 1.0 / this->d));
			int NumInOrigLat = (int)round(pow(LengInOrigLat, this->d));
			SetNumInEachSides(std::vector<size_t>(this->d, (size_t)LengInOrigLat));
		}
		else {
			std::cout << "the target packing fraction is too high!\n";
		}
	}
	virtual void SetRho(double numberdensity) { this->rho = numberdensity; }
	virtual void SetNumPrts(int numPrts) {
		this->N0 = numPrts;
		//int LengInOrigLat = (int)round(pow(numPrts / (1.0 - this->p), 1.0 / this->d));
		//int NumInOrigLat = (int)round(pow(LengInOrigLat, this->d));
		int LengInOrigLat = (int)round(pow(numPrts, 1.0 / this->d));
		int NumInOrigLat = (int)round(pow(LengInOrigLat, this->d));
		SetNumInEachSides(std::vector<size_t>(this->d, (size_t)LengInOrigLat));
	}

	virtual void GenerateC(Configuration &c) {
		Configuration temp = MultiplicateStructure(this->F_Cell, this->N);
		int currN, idx;	RandomGenerator g(std::time(NULL));
		int numVacancies = (int)round(this->N0*this->p);
		for (int i = 0; i < numVacancies; i++) {
			idx = (int)floor(temp.NumParticle() * g.RandomDouble());
			temp.DeleteParticle(idx);
		}
		c = Configuration(temp);
		//double numDensity = c.NumParticle() / c.PeriodicVolume();
		//double factor = pow(numDensity/ this->rho, 1.0/this->d);
		//c.Rescale(factor);
	}
	virtual void GenerateP(SpherePacking & c) {
		Configuration temp; this->GenerateC(temp);
		//double rho_i = this->N0 / (this->numPrts * this->F_Cell.PeriodicVolume());
		//double factor = pow(rho_i/this->rho, 1.0 / this->d);
		c = SpherePacking(temp, this->rMax);
	}

};

/*-------------------------------------------------------
	CLI for the tiling problems
	The contents also can be used as example codes
  -------------------------------------------------------*/
//Class to simplify computations
class Computation {	
protected:
	int numThreads=1;
	double phi_stored = 0.01; //the packing fraction of constructed dispersions
public:
	bool initialPackingFriendly;
	double origPackingFraction = -1;
	virtual void Compute(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfig) = 0;
	virtual void Write(const std::string OutputPrefix) = 0;
	//virtual void Plot(const std::string OutputPrefix, const std::string & Title) = 0;

	//some computations allow additional options
	//any command not understood by PairStatisticsCLI() is passed to here
	virtual void AdditionalOptions(const std::string & option, std::istream & input, std::ostream & output)
	{
		//by default, accept no additional option.
		output << "Unrecognized command!\n";
	}
	void SetNumThreads(int num) { this->numThreads = num; }
	virtual ~Computation()
	{
	}
};

class chivComputation : public Computation{
public:
	std::vector<GeometryVector> result;
	double CircularKMax, LinearKMax, KPrecision;
	bool lowStorage;	//Activate the low-storage mode of Isotropic spectral densities
	chivComputation(std::istream & ifile, std::ostream & ofile) {
		ofile << "circular K_max=";
		ifile >> CircularKMax; ofile << CircularKMax << "\n";
		ofile << "linear K_max=";
		ifile >> LinearKMax; ofile << LinearKMax << "\n";
		ofile << "K precision (binning width)=";
		ifile >> KPrecision;	ofile << KPrecision << "\n";
		ofile << "Low storage (y=1/n=0)? ";
		ifile >> lowStorage;
		//this->lowStorage = false;
		//SampleProbability = 1.0;
		initialPackingFriendly = true;
	}
	virtual void Compute(std::function<const SpherePacking(size_t i)> GetConfig, size_t numConfigs) {
		result.clear();
		omp_set_num_threads(this->numThreads);
		if (! this->lowStorage)
			::IsotropicSpectralDensity(GetConfig, numConfigs, CircularKMax, LinearKMax, result, KPrecision);
		else
			::IsotropicSpectralDensity_Lstorage(GetConfig, numConfigs, CircularKMax, LinearKMax, result, KPrecision);
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_chiv");
		WriteFunction(result, name);
	}
	virtual void AdditionalOptions(const std::string & option, std::istream & input, std::ostream & output) {
		if (option == "LowStorage") {
			//default = false; true => don't use map <,> in the calculation
			//input >> lowStorage;
			this->lowStorage = true;
		}
		else
			output << "Unrecognized command!\n";
	}
};

// voro++ is in progress
class VoronoiCellDistribution : public Computation {
public:
	std::vector<GeometryVector> result;
	double vMax, dv;
	bool IsInitPacking = true;
	VoronoiCellDistribution(std::istream & ifile, std::ostream & ofile) {
		ofile << "Max Cell-volume = ";
		ifile >> vMax;	ofile << vMax << "\n";
		ofile << "bin width = ";
		ifile >> dv;
		ofile << "do we analyze initial packings? (1=y/0=n)  ";
		ifile >> IsInitPacking; ofile << IsInitPacking << "\n";
		initialPackingFriendly = false;
	}
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfig, size_t numConfigs) {
		std::cout << "Computing Voronoi cell-volume distribution:\n";
		if (IsInitPacking) {
			if (::Verbosity > 2)
				std::cout << "Computing Voronoi cell-volume distribution from the initial packings\n";
#ifdef USE_VORO
			//Initialize the bins
			result.clear();
			int numBins = (int)ceil(vMax / dv);
			for (int i = 0; i < numBins; i++) {
				result.emplace_back(GeometryVector(dv*(i - 1), 0.0, 0.0, 0.0));
			}

			for (int i = 0; i < numConfigs; i++) {
				std::vector<double> vcVols;
				SpherePacking C = GetConfig(i);
				int numPrts = C.NumParticle(); int idx_max = result.size();
				VoronoiCellRoutine(C, this->numThreads, &vcVols, nullptr, nullptr, nullptr);

				std::vector<double> dist_i(result.size(), 0.0);
				for (int j = 0; j < numPrts; j++) {
					int idx = (int)floor(vcVols[j] / dv);
					if (idx < idx_max)
						dist_i[idx]++;
				}
				//compute the distribution for each configuration
				for (int j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts / dv;

					result[j].x[1] += dist_i[j];
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}

			//----final step
			for (int i = 0; i < result.size(); i++) {
				result[i].x[1] /= (double)numConfigs;
				result[i].x[2] = std::fabs(result[i].x[2] / (double)numConfigs - result[i].x[1] * result[i].x[1]);
				result[i].x[2] = sqrt(result[i].x[2] / (double)(numConfigs - 1));
			}
#else
			auto GetC = [&GetConfig](size_t i)->Configuration {return Configuration(GetConfig(i), "a"); };
			::VoronoiVolumeDistrubution(GetC, numConfigs, result, 1000);
#endif		
		}
		else {
			if (::Verbosity > 2)
				std::cout << "Computing Voronoi cell-volume distribution from the constructed dispersions\n";
			//Initialize the bins
			result.clear();
			int numBins = (int)ceil(vMax / dv);
			for (int i = 0; i < numBins; i++) {
				result.emplace_back(GeometryVector(dv*(i - 1), 0.0, 0.0, 0.0));
			}

#pragma omp parallel for schedule(guided) num_threads(this->numThreads)
			for (int i = 0; i < numConfigs; i++) {
				std::vector<double> dist_i(result.size(), 0.0);
				SpherePacking c = GetConfig(i);	int numPrts = c.NumParticle(); int idx_max = result.size();
				int d = c.GetDimension();
				for (int j = 0; j < numPrts; j++) {
					int idx = (int)floor(HyperSphere_Volume(d, c.GetCharacteristics(j)) / phi_stored / dv);
					if (idx < idx_max)
						dist_i[idx]++;
				}
				//compute the distribution for each configuration
				for (int j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts / dv;

#pragma omp atomic
					result[j].x[1] += dist_i[j];
#pragma omp atomic
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}

			//----final step
			for (int i = 0; i < result.size(); i++) {
				result[i].x[1] /= (double)numConfigs;
				result[i].x[2] = std::fabs(result[i].x[2]/(double) numConfigs - result[i].x[1]*result[i].x[1]);
				result[i].x[2] = sqrt(result[i].x[2] / (double)(numConfigs - 1));
			}
		}
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_cellvolDist");
		WriteFunction(result, name);
	}
};

class VoronoiCellVertexDistribution : public Computation {
public:
	std::vector<GeometryVector> result;
	int numMax;
	bool IsInitPacking = true;

	VoronoiCellVertexDistribution(std::istream & ifile, std::ostream & ofile) {
		ofile << "Max vertex number = ";
		ifile >> numMax;	ofile << numMax << "\n";
		//ofile << "do we analyze initial packings? (1=y/0=n)  ";
		//ifile >> IsInitPacking;
		initialPackingFriendly = false;
	}
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfig, size_t numConfigs) {
		std::cout << "Computing Voronoi cell vertex number distribution:\n";
		result.clear();
		int numBins = numMax;
		for (int i = 0; i < numBins; i++) {
			result.emplace_back(GeometryVector(i + 1, 0.0, 0.0, 0.0));
		}
		
#ifdef USE_VORO
			for (int i = 0; i < numConfigs; i++) {
				std::vector<int> numVert;
				SpherePacking C = GetConfig(i);
				int numPrts = C.NumParticle(); int idx_max = result.size();
				VoronoiCellRoutine(C, this->numThreads, nullptr, &numVert, nullptr, nullptr);

				std::vector<double> dist_i(result.size(), 0.0);
				for (int j = 0; j < numPrts; j++) {
					int idx = numVert[j] - 1;
					if (idx < idx_max)
						dist_i[idx]++;
				}
				//compute the distribution for each configuration
				for (int j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts;

					result[j].x[1] += dist_i[j];
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}
#else
#pragma omp parallel for schedule(guided) num_threads(this->numThreads)
			for (int i = 0; i < numConfigs; i++) {
				std::vector<double> dist_i(result.size(), 0.0);
				SpherePacking c = GetConfig(i);	int numPrts = c.NumParticle(); int idx_max = numBins;
				int d = c.GetDimension();	std::vector<GeometryVector> vertices;
				for (int j = 0; j < numPrts; j++) {
					GetVoronoiCell(c, j, vertices);
					int idx = vertices.size();
					if (idx < idx_max)
						dist_i[idx]++;
				}
				//compute the distribution for each configuration
				for (int j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts;
#pragma omp atomic
					result[j].x[1] += dist_i[j];
#pragma omp atomic
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}
#endif		
		//----final step
		for (int i = 0; i < result.size(); i++) {
			result[i].x[1] /= (double)numConfigs;
			result[i].x[2] = std::fabs(result[i].x[2] / (double)numConfigs - result[i].x[1] * result[i].x[1]);
			result[i].x[2] = sqrt(result[i].x[2] / (double)(numConfigs - 1));
		}
		
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_VertexDist");
		WriteFunction(result, name);
	}
};


/*calculate the maximal packing fractions
	type 1: largest particle size is fixed
	type 2: particle centeroid is fixed
	type 3: no constraints
*/
class MaxPackingFraction : public Computation {
public:
	std::vector<GeometryVector> maxPhi;
	std::vector<GeometryVector> result;
	
	int type;	bool IsInitPacking;
	MaxPackingFraction(std::istream & ifile, std::ostream & ofile) {
		ofile << "types of maximal packing fraction (1-3) ";
		ifile >> type;
		ofile << "do we analyze initial packings? (1=y/0=n)  ";
		ifile >> IsInitPacking;	ofile << IsInitPacking <<"\n";
		initialPackingFriendly = false;
	}
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		std::cout << "Computing the Maximal Packing Fractions:\n";
		//initialization
		maxPhi = std::vector<GeometryVector>(numConfigs, GeometryVector(static_cast<int> (type)));
		for (int i = 0; i < numConfigs; i++) {
			SpherePacking c = GetConfigs(i);
			int numPrts = c.NumParticle(), d = c.GetDimension();
			double phi1=0.0, phi2=0.0, phi3=0.0;
			std::vector<double> PHI1(numPrts, 0.0);
			std::vector<double> PHI2, PHI3;

			/*
				Initial particle sizes are known
				Voronoi cell volumes are unknown
			*/
			if (IsInitPacking) {//Initial particle sizes are known
#if defined(USE_VORO)
				VoronoiCellRoutine(c, this-> numThreads, &PHI1, nullptr, nullptr, nullptr);
#else
				std::vector<GeometryVector> v;
				for (int j = 0; j < numPrts; j++) {
					GetVoronoiCell(c, j, v);
					PHI1[j] = PolygonVolume(v, d);
				}
#endif
				if (type > 1)
					PHI2 = PHI1;
				if (type > 2)
					PHI3 = PHI1;

				for (int j = 0; j < numPrts; j++) {
					PHI1[j] = HyperSphere_Volume(d, 0.5*c.NearestParticleDistance(j)) / PHI1[j]; 
					if(type > 1)
						PHI2[j] = HyperSphere_Volume(d, c.GetCharacteristics(j)) / PHI2[j];
					if (type > 2) {
						//PHI3[j] = HyperSphere_Volume(d, c.GetCharacteristics(j)) / PHI1[j];
					}
				}
			}
			/*
				Initial particle sizes are unknown => Assume monodispersity
				Voronoi cell volumes are known
			*/
			else {
				double origPhi = origPackingFraction;
				double prtVol = origPhi / (c.NumParticle()/c.PeriodicVolume());
				double vcVolume;
				double currPackingFraction = c.PackingFraction();
#pragma omp parallel for num_threads(numThreads) schedule(guided)
				for (int j = 0; j < numPrts; j++)
					PHI1[j] = HyperSphere_Volume(d, c.GetCharacteristics(j)) / currPackingFraction; //phi_stored;
				
				if (type > 1)
					PHI2 = PHI1;
				if (type > 2)
					PHI3 = PHI1;

				for (int j = 0; j < numPrts; j++) {
					PHI1[j] = HyperSphere_Volume(d, 0.5*c.NearestParticleDistance(j)) / PHI1[j]; 
					if (type > 1)
						PHI2[j] = prtVol / PHI2[j];
					if (type > 2) {
						//PHI3[j] = HyperSphere_Volume(d, c.GetCharacteristics(j)) / PHI1[j];
					}
				}
			}

			/*Determine the minimum in this configuration*/
			for (int j = 0; j < numPrts; j++) {
				if (j == 0) {
					phi1 = PHI1[j];
					if (type > 1)
						phi2 = PHI2[j];
					if (type > 2)
						phi3 = PHI3[j];
				}
				else {
					phi1 = std::min(PHI1[j], phi1);
					if (type > 1)
						phi2 = std::min(PHI2[j], phi2);
					if (type > 2)
						phi3 = std::min(PHI3[j], phi3);
				}
			}
			maxPhi[i].x[0] = i;
			maxPhi[i].x[1] = phi1;
if(type>1)	maxPhi[i].x[2] = phi2;
if(type>2)	maxPhi[i].x[3] = phi3;
		}
		//Conclusion
		//result[0] = averages
		//result[1] = standard deviation
		//result[2] = mins
		result = std::vector<GeometryVector>(3,GeometryVector(static_cast<int> (type+1)));
		for (int i = 0; i < maxPhi.size(); i++) {
			if (i == 0) {
				result[0].x[1] = maxPhi[i].x[1];
				result[1].x[1] = maxPhi[i].x[1] * maxPhi[i].x[1];
				result[2].x[1] = maxPhi[i].x[1];
				if (type > 1) {
					result[0].x[2] = maxPhi[i].x[2];
					result[1].x[2] = maxPhi[i].x[2] * maxPhi[i].x[2];
					result[2].x[2] = maxPhi[i].x[2];
				}
				if (type > 2) {
					result[0].x[3] = maxPhi[i].x[3];
					result[1].x[3] = maxPhi[i].x[3] * maxPhi[i].x[3];
					result[2].x[3] = maxPhi[i].x[3];
				}
			}
			else {
				result[0].x[1] += maxPhi[i].x[1];
				result[1].x[1] += maxPhi[i].x[1] * maxPhi[i].x[1];
				result[2].x[1] = std::min(maxPhi[i].x[1], result[2].x[1]);
				
				if (type > 1) {
					result[0].x[2] += maxPhi[i].x[2];
					result[1].x[2] += maxPhi[i].x[2] * maxPhi[i].x[2];
					result[2].x[2] = std::min(maxPhi[i].x[2], result[2].x[2]);
				}
				if (type > 2) {
					result[0].x[3] += maxPhi[i].x[3];
					result[1].x[3] += maxPhi[i].x[3] * maxPhi[i].x[3];
					result[2].x[3] = std::min(maxPhi[i].x[3], result[2].x[3]);
				}
				
			}
		}

		for (int i = 1; i <= type; i++) {
			result[0].x[i] /= (double)numConfigs;
			result[1].x[i] = result[1].x[i] / (double)numConfigs - result[0].x[i]*result[0].x[i];
			result[1].x[i] = sqrt(std::fabs(result[1].x[i])/((double)numConfigs-1.0));
		}

		std::cout << "Max Packing Fraction:\n";
		std::cout << "phi_1\tphi_2\tphi_3\n";
		for(int i=0; i<3; i++)
			std::cout << result[i].x[1]<<"\t"<<result[i].x[2]<<"\t"<<result[i].x[3]<<"\n";
//		result[2].x[3] = result[2].x[2];
//		result.x[2] = result.x[1];
//		result.x[1] = result.x[0];
		result[0].x[0] = 0;
		result[1].x[0] = 0;
		result[2].x[0] = 0;
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_MaxPackingFraction");
		WriteFunction(maxPhi, name);
		std::ofstream add;	add.open(name + std::string(".txt"), std::fstream::app);
		add << "---------------------Summary-----------------------\n";
		add << result[0];
		add << result[1];
		add << result[2];
		add.close();
	}
};

#include "NumberVariance.h"
//Compute the local number variance of centroids of particles
class LocalNumberVariance : public Computation {
public:
	int numSamplingPts;
	double Rmax, dR;
	std::vector<GeometryVector> result;
	LocalNumberVariance(std::istream & ifile, std::ostream & ofile) {
		ofile << "Rmax = ";
		ifile >> Rmax;	ofile << "\n" << Rmax << "\n";
		ofile << "dR = ";
		ifile >> dR;	ofile << "\n" << dR << "\n";
		ofile << "the number of sampling points = ";
		ifile >> numSamplingPts;	ofile << "\n" << numSamplingPts << "\n";
		initialPackingFriendly = false;
	}
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		auto GetC = [&GetConfigs](size_t i)->Configuration {return Configuration(GetConfigs(i),"a"); };
		omp_set_num_threads(this->numThreads);
		MC_nv(GetC, numConfigs, result, Rmax, dR, numSamplingPts);
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_localNumVar");
		WriteFunction(result, name);
	}
};
//Compute the loca volume-fraction variance of particle phase
class LocalVolumeFractionVariance : public Computation {
public:
	int numSamplingPts;
	double Rmax, dR;
	std::vector<GeometryVector> result;

	LocalVolumeFractionVariance(std::istream & ifile, std::ostream & ofile) {
		ofile << "Rmax = ";
		ifile >> Rmax;	ofile << "\t" << Rmax << "\n";
		ofile << "dR = ";
		ifile >> dR;	ofile << "\t" << dR << "\n";
		ofile << "the number of sampling points = ";
		ifile >> numSamplingPts;	ofile << "\t" << numSamplingPts << "\n";
		initialPackingFriendly = true;
	}
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		omp_set_num_threads(this->numThreads);
		MC_vfv(GetConfigs, numConfigs, result, Rmax, dR, numSamplingPts);
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string name = OutputPrefix + std::string("_localVolVar");
		WriteFunction(result, name);
	}
};

//Compute the deviations in the centroids
class CentroidComputation : public Computation {
public:
	std::vector<std::vector<GeometryVector>> displacement_i;
	std::vector<GeometryVector> result;
	bool output_each_config = false;
	CentroidComputation(std::istream & ifile, std::ostream & ofile) {
		ofile << "Output displacements in each config? (y=1/n=0)\t";
		ifile >> output_each_config; ofile << "  " << output_each_config << "\n";
	}
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		std::cout << "Computing diplacements between particles and cell centroids:\n";
		displacement_i.reserve(numConfigs);	result.reserve(numConfigs);
#if defined(USE_VORO)
		for (int i = 0; i < numConfigs; i++) {
			std::vector<double> cellvols;
			std::vector<GeometryVector> dr;
			SpherePacking C = GetConfigs(i);	int numPrts = C.NumParticle();

			VoronoiCellRoutine(C, this->numThreads, &cellvols, nullptr, nullptr, &dr);
			GeometryVector sum(static_cast<int> (C.GetDimension()));
			for (int j = 0; j < numPrts; j++) {
				sum = sum + cellvols[j] * dr[j];
			}
			result.push_back(sum);	displacement_i.push_back(dr);
		}
#else
		std::cout << "Please include \"voro++\"\n";
#endif
	}

	virtual void Write(const std::string OutputPrefix) {
		std::string fname = OutputPrefix + std::string("_SumDisplacements");
		WriteFunction(result, fname.c_str());
		if (this->output_each_config) {
			for (int i = 0; i < displacement_i.size(); i++) {
				fname = OutputPrefix + std::string("_SumDisplacements_") + std::to_string(i);
				WriteFunction(displacement_i[i], fname);
			}
		}
	}

};


/*	------------------------------------------------
	Computations for the multiscale tilings
	------------------------------------------------*/
class MultiscaleComputation : public Computation {
protected:
	double tolerance = 1e-15;
	std::vector<std::vector<int>> indices;
	//for configuration i;
	// 0			<= idx_m=1 <=indices[i][0]
	// index[i][0]	< idx_m=2 <= indices[i][1] ...
public:
	void OutputIndexList(std::vector<std::vector<int>> & list) {
		list = std::vector<std::vector<int>>(this->indices);
	}
	void SetIndexList(const std::vector<std::vector<int>> & list) {
		this->indices = std::vector<std::vector<int>>(list);
	}
	/* determine stages without any prior knowledge*/
	void GetIndexList(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		this->indices.resize(numConfigs, std::vector<int>());
		/*	R[i][j] = jth largest radii for configuration i
			it needs when the number of stages are not the same	*/
		std::vector<std::vector<double>> R(numConfigs, std::vector<double>());
#pragma omp parallel for num_threads(this->numThreads) schedule(guided)
		for (int i = 0; i < numConfigs; i++) {
			//packing fraction at each stage in THIS configuration
			SpherePacking c = GetConfigs(i);
			int numprts = c.NumParticle();
			int m = 1, idx_prev = -1;	//stage
			double phi0 = 0.0, r_prev = c.GetCharacteristics(0), r_curr = 0.0;
			for (int j = 0; j < numprts; j++) {
				r_curr = c.GetCharacteristics(j);
				if (r_prev - r_curr > this->tolerance) {//entering to the next stage
					m++;
					R[i].push_back(r_prev);
					this->indices[i].push_back(j - 1);
					idx_prev = j - 1;
					r_prev = r_curr;
				}
			}
			this->indices[i].push_back(numprts - 1);
		}
		/*-------
		for (int i = 0; i < 200; i++) {
			std::cout << i;
			for (int j = 0; j < 5; j++)
				if (i < R[j].size())
					std::cout << "\t" << R[j][i];
				else
					std::cout << "\t\t";
			std::cout << "\n";
		}*/


		bool isThereMissedStage = false;
		//check whether the number of stages is identical for all configurations
		for (int i = 1; i < numConfigs; i++)
			isThereMissedStage |= (R[i - 1].size() != R[i].size());
		if (isThereMissedStage) {
			/*	Particles of a certain radius are missed
				We make a complete list of particle radii.
				While doing so, we also remake the list of particle indices
				R[i][idx_stage[i]] will represent particle radius of configuration i that we are looking at
				*/
			std::vector<double> R_new;	
			std::vector<std::vector<int>> idx_list_new (numConfigs, std::vector<int> ());
			std::vector<size_t> idx_stage(numConfigs, 0);	
			
			bool conti;
			do {//repeat until idx_stage[i] == indices[i].size() - 1 for all i
				conti = false;
				bool radiiAreSame = true;
				double * dev = new double[numConfigs];
				double r0 = R[0][idx_stage[0]];	dev[0] = 0.0;
				for (int i = 1; i < numConfigs; i++) {
					dev[i] = R[i][idx_stage[i]];	//particle radius of config i
					radiiAreSame &= (fabs(dev[i] - r0) < this->tolerance);
				}
				
				if (!radiiAreSame) {
					double r_max = r0;	//find the largest radii at this stage
					for (int i = 1; i < numConfigs; i++)
						r_max = std::max(r_max, R[i][idx_stage[i]]);
														
					R_new.push_back(r_max);
					//update index only when particle radius is the same as the smallest one.
					for (int i = 0; i < numConfigs; i++) {
						idx_list_new[i].push_back(this->indices[i][idx_stage[i]]);
						if (r_max - R[i][idx_stage[i]] < this->tolerance && (idx_stage[i] < R[i].size() - 1)) {
							idx_stage[i]++;
						}
					}
				}
				else {
					R_new.push_back(r0);
					for (int i = 0; i < numConfigs; i++) {
						idx_list_new[i].push_back(this->indices[i][idx_stage[i]]);
						if (idx_stage[i] < R[i].size() - 1)
							idx_stage[i]++;
					}
				}
				
				
				//conti = true if not all of indices become its maximum size
				for (int i = 0; i < numConfigs; i++) 
					conti |= (idx_stage[i] < R[i].size() - 1);
			
			} while (conti);
			this->SetIndexList(idx_list_new);
			//this->indices=idx_list_new;
		}
	}
	/*The indices list is competed after the computation*/
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) = 0;
	virtual void Write(const std::string OutputPrefix) = 0;

};

class SpectralDensity_Multiscale : public MultiscaleComputation {
public:
	int interval;
	double CircularKmax, LinearKmax, KPrecision;
	bool aroundOriginOnly;
	std::vector < std::vector<GeometryVector>> result;
	std::vector<int> outputIndex;
	SpectralDensity_Multiscale(std::istream & ifile, std::ostream & ofile) {
		ofile << "Do we plot chiv(0) vs m only? (y=1/n=0)\t";
		ifile >> aroundOriginOnly; ofile << "\t" << aroundOriginOnly << "\n";
		ofile << "CircularKmax = ";
		ifile >> CircularKmax; ofile << "\t" << CircularKmax << "\n";
		ofile << "LinearKmax = ";
		ifile >> LinearKmax; ofile << "\t" << LinearKmax << "\n";
		ofile << "KPrecision = ";
		ifile >> KPrecision; ofile <<"\t"<<KPrecision<<"\n";
		ofile << "Interval = ";
		ifile >> interval; ofile << "\t" << interval << "\n";
		this->SetIndexList(std::vector<std::vector<int>>());
	}
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		std::cout << "Computing the spectral densities for multiscale coated-sphere model at every " << this->interval << ":\n";
		bool fillTheList = (this->indices.size() == 0);
		//find the index list
		if (fillTheList)
			this->GetIndexList(GetConfigs, numConfigs);

		for (int i = 1; i <= this->indices[0].size(); i++) {
			if (i==1 || i % this->interval == 0 || i == *(this->indices[0].end()))
				outputIndex.push_back(i-1);
		}
	
		struct SkBin
		{
			double Sum1, SumK1, SumK2, SumS, SumS2;
			SkBin() : Sum1(0), SumK1(0), SumK2(0), SumS(0), SumS2(0)
			{}
		};
		size_t NumBin = (size_t)std::floor(LinearKmax / KPrecision) + 1;

		//k-points at which we compute the spectral densities
		std::vector<GeometryVector> ks = GetKs(GetConfigs(0), CircularKmax, LinearKmax);
		//ensemble average of the spectral densities at each stage
		std::vector<std::vector<GeometryVector>> chiV_k(outputIndex.size(), std::vector<GeometryVector>(NumBin, GeometryVector(0.0, 0.0, 0.0, 0.0)));
		std::vector<double> sum1(NumBin, 0.0);
		
#pragma omp parallel for num_threads(this->numThreads) schedule(guided)
		for (int i = 0; i < numConfigs; i++) {
			//CONFIGURATIONS
			//Itilde[l] = chi_v(ks[l]): F.T. indicator fn of config i.
			std::vector<std::complex<double>> Itilde(ks.size(), std::complex<double>(0.0, 0.0));
			//chiV_k_m[m][idx_k] = spectral density at (idx_k) th bin in m-th stage
			std::vector<std::vector<SkBin>> chiV_k_m(outputIndex.size(), std::vector<SkBin>(NumBin, SkBin()));
			SpherePacking C = GetConfigs(i);
			double V = C.PeriodicVolume();

			int idx_init, idx_fin = -1, d=C.GetDimension();
			//STEPS: from stage outputIndex[m] to stage outputIndex[m+1]
			for (int m = 0; m < outputIndex.size(); m++) {
				//Set the particle indices from the stage "outputIndex[m-1]" to "outputIndex[m]"
				idx_init = idx_fin + 1;
				idx_fin = indices[i][outputIndex[m]];
				//Compute Itilde[l]'s
				
				for (int l = 0; l < ks.size(); l++) {//At a given wavevector
					std::map < double, double> mtildes;
					double k2 = ks[l].Modulus2();
					double k = std::sqrt(k2);
					double r = 0.0;
					double mt = 0.0;
					for (int n = idx_init; n <= idx_fin; n++) {//for prts
						r = C.GetCharacteristics(n);
						auto iter = mtildes.find(r);
						if (iter != mtildes.end())
							mt = iter->second;
						else
						{
							mt = Mtilde(k, r, d);
							mtildes.insert(std::make_pair(r, mt));
						}
						Itilde[l] += mt * std::exp(std::complex<double>(0.0, ks[l].Dot(C.GetCartesianCoordinates(n))));
					}

					size_t idx_k = (size_t)std::floor(k / KPrecision);
					double s = (Itilde[l].real() * Itilde[l].real() + Itilde[l].imag() * Itilde[l].imag()) / V;
					chiV_k_m[m][idx_k].Sum1 += 1.0;
					chiV_k_m[m][idx_k].SumK1 += k;
					chiV_k_m[m][idx_k].SumK2 += k2;
					chiV_k_m[m][idx_k].SumS += s;
					chiV_k_m[m][idx_k].SumS2 += s * s;
				}
			}

			/*To chiV_k*/
			for (int m = 0; m < outputIndex.size(); m++) {
				for (int j = 0; j < NumBin; j++) {
					if (chiV_k_m[m][j].Sum1 > 0.0) {
						if (i == 0) {
							//Wavevectors
							double k = chiV_k_m[m][j].SumK1 / chiV_k_m[m][j].Sum1;
							//var in Wavevectors
							double Dk = sqrt(fabs(chiV_k_m[m][j].SumK2 / chiV_k_m[m][j].Sum1 - k * k) / chiV_k_m[m][j].Sum1);
#pragma omp atomic
							chiV_k[m][j].x[0] += k;
#pragma omp atomic
							chiV_k[m][j].x[2] += Dk;
						}

						if (i == 0 && m == 0) {
#pragma omp atomic
							sum1[j] += chiV_k_m[m][j].Sum1;
						}

						double s = chiV_k_m[m][j].SumS / chiV_k_m[m][j].Sum1;
						double s2 = chiV_k_m[m][j].SumS2 / chiV_k_m[m][j].Sum1;
#pragma omp atomic
						chiV_k[m][j].x[1] += s;
#pragma omp atomic
						chiV_k[m][j].x[3] += s2;
					}
				}
			}
		}


		/*	Summarize and Plot

		*/
		for (int m = 0; m < outputIndex.size(); m++) {
			for (int k = 0; k < NumBin; k++) {
				if (sum1[k] > 0) {
					double s = chiV_k[m][k].x[1] / (double)numConfigs;
					double var_s = std::sqrt(std::fabs(chiV_k[m][k].x[3] / (double)numConfigs - s * s) / sum1[k]/(double)numConfigs);

					chiV_k[m][k].x[1] = s;
					chiV_k[m][k].x[3] = var_s;
				}
				else {
					chiV_k[m][k] = GeometryVector(0, 0, 0, 0);
				}
			}
		}
		this->result = chiV_k;
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string filename;
		if (aroundOriginOnly) {//Write chiv(0) vs m
			std::vector<GeometryVector> chiv0;
			for (int i = 0; i < outputIndex.size(); i++) {
				chiv0.emplace_back(GeometryVector(outputIndex[i]+1, result[i][0].x[1], result[i][0].x[3]));
			}
			filename = OutputPrefix + std::string("_chiv0");
			WriteFunction(chiv0, filename);
		}
		else {
			for (int i = 0; i < outputIndex.size(); i++) {
				filename = OutputPrefix + std::string("_chiv_") + std::to_string(outputIndex[i]+1);
				WriteFunction(result[i], filename);
			}
		}

	}
};

class VolumeFraction_Multiscale : public MultiscaleComputation {
public:
	int interval;// , maxstep;
	std::vector<GeometryVector> result;

	VolumeFraction_Multiscale(std::istream & ifile, std::ostream & ofile) {
		ofile << "Interval in stages = \t";
		ifile >> interval;	ofile << "\t" << interval << "\n";
		//ofile << "max in stages = \n";
		//ifile >> maxstep;	ofile << "\t" << maxstep << "\n";
		this->SetIndexList(std::vector<std::vector<int>>());
	}

	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		std::cout << "Computing the volume fraction at every " << this->interval << " stages:\n";
		
		bool fillTheList = (this->indices.size() == 0);
		//find the index list
		if (fillTheList) 
			this->GetIndexList(GetConfigs, numConfigs);

		int binsize = (int)round((this->indices[0].size()) / interval);
		result.clear();
		{
			GeometryVector temp(0, 0, 0, 0);	result.resize(binsize, temp);
		}
		//Compute based on the index list
#pragma omp parallel for num_threads(this->numThreads) schedule(guided)
		for (int i = 0; i < numConfigs; i++) {
			//packing fraction at each stage in THIS configuration
			SpherePacking c = GetConfigs(i);
			int numprts = c.NumParticle(), d = c.GetDimension();
			int m = 1, idx_prev = -1, tick=0;	//stage
			double V = c.PeriodicVolume();
			double phi0 = 0.0;
			for (int j = 0; j < this->indices[i].size(); j++) {
				phi0 += HyperSphere_Volume(d, c.GetCharacteristics(this->indices[i][j]))*(this->indices[i][j] - idx_prev) / V;
				idx_prev = this->indices[i][j];
				if (m%interval == 0) {
#pragma omp atomic
					result[tick].x[1] += phi0;
#pragma omp atomic
					result[tick].x[2] += phi0 * phi0;

					tick++;
				}
				else if (m == this->indices[i].size()) {
#pragma omp atomic
					result[tick].x[1] += phi0;
#pragma omp atomic
					result[tick].x[2] += phi0 * phi0;
				}
				m++;
			}
		}
		//Summary
		for (int i = 0; i < result.size(); i++) {
			if (i < result.size() - 1)
				result[i].x[0] = (i + 1) * this->interval;
			else
				result[i].x[0] = (this->indices[0].size());
			result[i].x[1] /= (double)numConfigs;
			result[i].x[2] = fabs(result[i].x[2] / (double)numConfigs - result[i].x[1] * result[i].x[1]);
			result[i].x[2] = sqrt(result[i].x[2] / ((double)numConfigs - 1));
		}
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string filename = OutputPrefix + std::string("_PackingFractions");
		WriteFunction(result, filename);
	}
};

class VolumeFractionVariance_Multiscale : public MultiscaleComputation {
public:
	int interval;

	double Rmax, dR;
	int numSamplePts;
	std::vector<std::vector<GeometryVector>> result;
	std::vector<int> outputIndex;
	VolumeFractionVariance_Multiscale(std::istream & ifile, std::ostream & ofile) {
		ofile << "Interval in stage = \t";	ifile >> interval; ofile << "\t" << interval << "\n";
		ofile << "R max = \t";	ifile >> Rmax;	ofile << "\t" << Rmax<<"\n";
		ofile << "dR = \t";		ifile >> dR;	ofile << "\t" << dR << "\n";
		ofile << "the number of sample points = \t"; ifile >> numSamplePts; ofile << "\t" << numSamplePts << "\n";
	}

	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		std::cout << "Computing the local volume fraction variance at every " << this->interval << " stages:\n";
		bool fillTheList = (this->indices.size() == 0);
		//find the index list
		if (fillTheList)
			this->GetIndexList(GetConfigs, numConfigs);

		for (int i = 1; i <= this->indices[0].size(); i++) {
			if (i % this->interval == 0 || i == *(this->indices[0].end()))
				outputIndex.push_back(i);
		}

		for (auto iter = outputIndex.begin(); iter != outputIndex.end(); iter++) {
			std::vector<int> idx_prt; //idx_prt[i] = particle idx in the configuration i at the *iter stage
			for (int i = 0; i < numConfigs; i++)
				idx_prt.push_back(this->indices[i][*iter]);
			//Remove particles at stages higher than m
			auto getC = [&GetConfigs, &idx_prt](size_t j)->SpherePacking {
				SpherePacking c = GetConfigs(j);
				int num = c.NumParticle() - 1;
				for (int k = num; k > idx_prt[j]; k--) {
					c.DeleteParticle(j);
				}
				return c;
			};

			std::vector<GeometryVector> RESULT;
			MC_vfv(getC, numConfigs, RESULT, this->Rmax, this->dR, this->numSamplePts);
			this->result.push_back(RESULT);
		}
	}
	virtual void Write(const std::string OutputPrefix) {
		std::string filename;
		for (int i = 0; i < outputIndex.size(); i++) {
			filename = OutputPrefix + std::string("_LocalVolVar_") + std::to_string(outputIndex[i]);
			WriteFunction(result[i], filename);
		}
	}

};

#endif
