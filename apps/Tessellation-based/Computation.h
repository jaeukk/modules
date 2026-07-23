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
#include <VolumeFractionVariance.h>	// MC_vfv (split out of NumberVariance)
#include <StructureFactor.h>
#include <SpectralDensity.h>		// IsotropicSpectralDensity (split out of StructureFactor)

//==============================================================================
// Body sliced verbatim from the legacy TilingProblem.cpp (lines 345-1362):
// ImLatGen + the Computation base class and all of its derivatives.
//==============================================================================
//lattice with vacancies
//
// Rewritten against the current LatticeGen API (models/GenerateConfigs.h). The
// legacy members it relied on -- LatticeGen(basis,int), SetBasis, currPhi,
// SetNumInEachSides, N, rMax -- no longer exist. Instead we build a 1-particle
// cubic unit cell (number density 1; the LatticeGen ctor sets each particle
// radius to half the nearest-neighbour distance = 0.5 here), keep our own
// tile-count and particle radius, then tile the cell and delete random
// vacancies. N0 is the vacancy-free site count = (sites/side)^d, so p is the
// mean vacancy fraction.
struct ImLatGen : LatticeGen {
private:
	double p;						//mean fraction of vacancies
	int N0;							//number of sites in the vacancy-free tiling
	std::vector<size_t> numInSide;	//tile count per direction (base NumInSide is private)
	double rMax;					//particle radius taken from the unit cell

	/** \brief Build the standard cubic (axis-aligned unit) basis in the given dimension.
	 *  \param dimension  Spatial dimension.
	 *  \return           Set of orthonormal cubic basis vectors. */
	static std::vector<GeometryVector> getCubicBasis(int dimension) {
		std::vector<GeometryVector> b;
		for (int i = 0; i < dimension; i++) {
			b.emplace_back(GeometryVector(static_cast<int>(dimension)));
			b[i].x[i] = 1.0;
		}
		return b;
	}

	/** \brief Build a one-particle cubic unit-cell packing from basis vectors.
	 * \param basis Basis vectors of the unit cell.
	 * \return Sphere packing containing one zero-characteristic particle at the origin. */
	static SpherePacking buildUnitCell(const std::vector<GeometryVector> & basis) {
		std::vector<GeometryVector> b(basis);
		SpherePacking cell((DimensionType)b.size(), &b[0], 1.0);
		cell.Insert(0.0, GeometryVector(static_cast<int>(b.size())));
		return cell;
	}
public:
	/** \brief Construct an imperfect-lattice generator from explicit basis vectors.
	 *  \param basis          Basis vectors of the underlying lattice.
	 *  \param numPrts        Number of particles in the initial (vacancy-free) lattice.
	 *  \param p              Mean fraction of vacancies.
	 *  \param numberdensity  Target number density. */
	ImLatGen(const std::vector<GeometryVector> & basis, int numPrts, double p, double numberdensity)
		: LatticeGen(buildUnitCell(basis)) {
		this->p = p;
		this->rho = numberdensity;
		int LengInOrigLat = (int)round(pow((double)numPrts, 1.0 / this->d));
		this->numInSide = std::vector<size_t>(this->d, (size_t)LengInOrigLat);
		this->SetNumCells(this->numInSide);
		this->N0 = (int)round(pow((double)LengInOrigLat, (double)this->d));
		this->rMax = this->F_Cell.GetCharacteristics(0);
	}

	/** \brief Construct an imperfect-lattice generator on a cubic lattice.
	 *  \param dimension      Spatial dimension.
	 *  \param numPrts        Number of particles in the initial (vacancy-free) lattice.
	 *  \param p              Mean fraction of vacancies.
	 *  \param numberdensity  Target number density. */
	ImLatGen(int dimension, int numPrts, double p, double numberdensity) : ImLatGen(getCubicBasis(dimension), numPrts, p, numberdensity) {}
	/** \brief set the final packing fraction!! (and derive the vacancy fraction)
	 *  \param phi  Target packing fraction (must be below maxPhi). */
	virtual void SetPackingFraction(double phi) {
		if (phi < this->maxPhi) {
			this->p = 1.0 - phi / this->maxPhi;
			int Leng = (int)round(pow(this->N0 / (1.0 - this->p), 1.0 / this->d));
			this->numInSide = std::vector<size_t>(this->d, (size_t)Leng);
			this->SetNumCells(this->numInSide);
		}
		else {
			std::cout << "the target packing fraction is too high!\n";
		}
	}

	/** \brief Set the target number density.
	 *  \param numberdensity  Target number density. */
	virtual void SetRho(double numberdensity) { this->rho = numberdensity; }
	/** \brief Set the nominal vacancy-free site count and derive an equal per-axis tile count.
	 * \param numPrts Nominal vacancy-free site count. */
	virtual void SetNumPrts(int numPrts) {
		this->N0 = numPrts;
		int Leng = (int)round(pow((double)numPrts, 1.0 / this->d));
		this->numInSide = std::vector<size_t>(this->d, (size_t)Leng);
		this->SetNumCells(this->numInSide);
	}

	/** \brief Generate a point configuration: tile the cell, then delete random vacancies.
	 *  \param[out] c  Output configuration with vacancies removed. */
	virtual void GenerateC(Configuration &c) {
		Configuration temp = MultiplicateStructure(Configuration(this->F_Cell, "a"), this->numInSide);
		RandomGenerator g(std::time(NULL));
		int numVacancies = (int)round(this->N0*this->p);
		for (int i = 0; i < numVacancies; i++) {
			int idx = (int)floor(temp.NumParticle() * g.RandomDouble());
			temp.DeleteParticle(idx);
		}
		c = temp;
	}

	/** \brief Generate a sphere packing from the vacancy-bearing lattice.
	 *  \param[out] c  Output sphere packing built from the generated configuration. */
	virtual void GenerateP(SpherePacking & c) {
		Configuration temp; this->GenerateC(temp);
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
	/** \brief Compute the statistic over an ensemble of sphere packings.
	 *  \param GetConfigsFunction  Callable returning configuration i.
	 *  \param NumConfig           Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfig) = 0;

	/** \brief Write the computed result to disk.
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix) = 0;
	//virtual void Plot(const std::string OutputPrefix, const std::string & Title) = 0;

	/** \brief some computations allow additional options;
	 *  any command not understood by PairStatisticsCLI() is passed to here.
	 *  \param option  Command keyword not understood by the main CLI.
	 *  \param input   Input stream for reading any option arguments.
	 *  \param output  Output stream for messages. */
	virtual void AdditionalOptions(const std::string & option, std::istream & input, std::ostream & output)
	{
		//by default, accept no additional option.
		output << "Unrecognized command!\n";
	}

	/** \brief Set the number of threads used by Compute.
	 *  \param num  Number of threads. */
	void SetNumThreads(int num) { this->numThreads = num; }
	/** \brief Destroy the Computation. */
	virtual ~Computation()
	{
	}
};

class chivComputation : public Computation{
public:
	std::vector<GeometryVector> result;
	double CircularKMax, LinearKMax, KPrecision;
	bool lowStorage;	//Activate the low-storage mode of Isotropic spectral densities
	/** \brief Construct by reading the circular and linear k limits, binning precision, and low-storage flag.
	 * \param ifile Input stream supplying the parameters.
	 * \param ofile Output stream echoing the prompts. */
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

	/** \brief Compute the isotropic spectral density chi_V(k) over the ensemble.
	 *  \param GetConfig   Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t i)> GetConfig, size_t numConfigs);

	/** \brief Write chi_V(k) to "<prefix>_chiv".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);

	/** \brief Enable low-storage spectral-density mode when option is "LowStorage".
	 * \param option Command keyword.
	 * \param input Unused input stream.
	 * \param output Output stream used for an unrecognized-command message. */
	virtual void AdditionalOptions(const std::string & option, std::istream & input, std::ostream & output);
};

// voro++ is in progress
class VoronoiCellDistribution : public Computation {
public:
	std::vector<GeometryVector> result;
	double vMax, dv;
	bool IsInitPacking = true;
	/** \brief Construct by reading the cell-volume histogram parameters from a stream.
	 *  \param ifile  Input stream supplying vMax, bin width and the init-packing flag.
	 *  \param ofile  Output stream echoing the prompts. */
	VoronoiCellDistribution(std::istream & ifile, std::ostream & ofile) {
		ofile << "Max Cell-volume = ";
		ifile >> vMax;	ofile << vMax << "\n";
		ofile << "bin width = ";
		ifile >> dv;
		ofile << "do we analyze initial packings? (1=y/0=n)  ";
		ifile >> IsInitPacking; ofile << IsInitPacking << "\n";
		initialPackingFriendly = false;
	}

	/** \brief Compute the Voronoi cell-volume distribution over the ensemble.
	 *  \param GetConfig   Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
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

			for (size_t i = 0; i < numConfigs; i++) {
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
				for (size_t j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts / dv;

					result[j].x[1] += dist_i[j];
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}

			//----final step
			for (size_t i = 0; i < result.size(); i++) {
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
			for (size_t i = 0; i < numConfigs; i++) {
				std::vector<double> dist_i(result.size(), 0.0);
				SpherePacking c = GetConfig(i);	int numPrts = c.NumParticle(); int idx_max = result.size();
				int d = c.GetDimension();
				for (int j = 0; j < numPrts; j++) {
					int idx = (int)floor(HyperSphere_Volume(d, c.GetCharacteristics(j)) / phi_stored / dv);
					if (idx < idx_max)
						dist_i[idx]++;
				}
				//compute the distribution for each configuration
				for (size_t j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts / dv;

#pragma omp atomic
					result[j].x[1] += dist_i[j];
#pragma omp atomic
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}

			//----final step
			for (size_t i = 0; i < result.size(); i++) {
				result[i].x[1] /= (double)numConfigs;
				result[i].x[2] = std::fabs(result[i].x[2]/(double) numConfigs - result[i].x[1]*result[i].x[1]);
				result[i].x[2] = sqrt(result[i].x[2] / (double)(numConfigs - 1));
			}
		}
	}

	/** \brief Write the cell-volume distribution to "<prefix>_cellvolDist".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);
};

class VoronoiCellVertexDistribution : public Computation {
public:
	std::vector<GeometryVector> result;
	int numMax;
	bool IsInitPacking = true;

	/** \brief Construct by reading the maximum vertex number from a stream.
	 *  \param ifile  Input stream supplying numMax.
	 *  \param ofile  Output stream echoing the prompts. */
	VoronoiCellVertexDistribution(std::istream & ifile, std::ostream & ofile) {
		ofile << "Max vertex number = ";
		ifile >> numMax;	ofile << numMax << "\n";
		//ofile << "do we analyze initial packings? (1=y/0=n)  ";
		//ifile >> IsInitPacking;
		initialPackingFriendly = false;
	}

	/** \brief Compute the distribution of Voronoi cell vertex counts over the ensemble.
	 *  \param GetConfig   Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfig, size_t numConfigs) {
		std::cout << "Computing Voronoi cell vertex number distribution:\n";
		result.clear();
		int numBins = numMax;
		for (int i = 0; i < numBins; i++) {
			result.emplace_back(GeometryVector(i + 1, 0.0, 0.0, 0.0));
		}
		
#ifdef USE_VORO
			for (size_t i = 0; i < numConfigs; i++) {
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
				for (size_t j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts;

					result[j].x[1] += dist_i[j];
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}
#else
#pragma omp parallel for schedule(guided) num_threads(this->numThreads)
			for (size_t i = 0; i < numConfigs; i++) {
				std::vector<double> dist_i(result.size(), 0.0);
				SpherePacking c = GetConfig(i);	int numPrts = c.NumParticle(); int idx_max = numBins;
				std::vector<GeometryVector> vertices;
				for (int j = 0; j < numPrts; j++) {
					GetVoronoiCell(c, j, vertices);
					int idx = (int)vertices.size() - 1;
					if (idx < idx_max)
						dist_i[idx]++;
				}
				//compute the distribution for each configuration
				for (size_t j = 0; j < dist_i.size(); j++) {
					dist_i[j] = dist_i[j] / numPrts;
#pragma omp atomic
					result[j].x[1] += dist_i[j];
#pragma omp atomic
					result[j].x[2] += dist_i[j] * dist_i[j];
				}
			}
#endif		
		//----final step
		for (size_t i = 0; i < result.size(); i++) {
			result[i].x[1] /= (double)numConfigs;
			result[i].x[2] = std::fabs(result[i].x[2] / (double)numConfigs - result[i].x[1] * result[i].x[1]);
			result[i].x[2] = sqrt(result[i].x[2] / (double)(numConfigs - 1));
		}
		
	}

	/** \brief Write the vertex-count distribution to "<prefix>_VertexDist".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);
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
	/** \brief Construct by reading the max-packing-fraction type and init-packing flag.
	 *  \param ifile  Input stream supplying the type (1-3) and the init-packing flag.
	 *  \param ofile  Output stream echoing the prompts. */
	MaxPackingFraction(std::istream & ifile, std::ostream & ofile) {
		ofile << "types of maximal packing fraction (1-3) ";
		ifile >> type;
		ofile << "do we analyze initial packings? (1=y/0=n)  ";
		ifile >> IsInitPacking;	ofile << IsInitPacking <<"\n";
		initialPackingFriendly = false;
	}

	/** \brief Compute the maximal (achievable) packing fractions over the ensemble.
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		std::cout << "Computing the Maximal Packing Fractions:\n";
		//initialization
		maxPhi = std::vector<GeometryVector>(numConfigs, GeometryVector(static_cast<int> (type)));
		for (size_t i = 0; i < numConfigs; i++) {
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
		for (size_t i = 0; i < maxPhi.size(); i++) {
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

	/** \brief Write the per-config max packing fractions plus a summary to "<prefix>_MaxPackingFraction".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);
};

#include "NumberVariance.h"
//Compute the local number variance of centroids of particles
class LocalNumberVariance : public Computation {
public:
	int numSamplingPts;
	double Rmax, dR;
	std::vector<GeometryVector> result;
	/** \brief Construct by reading Rmax, dR and the sampling-point count from a stream.
	 *  \param ifile  Input stream supplying Rmax, dR and numSamplingPts.
	 *  \param ofile  Output stream echoing the prompts. */
	LocalNumberVariance(std::istream & ifile, std::ostream & ofile) {
		ofile << "Rmax = ";
		ifile >> Rmax;	ofile << "\n" << Rmax << "\n";
		ofile << "dR = ";
		ifile >> dR;	ofile << "\n" << dR << "\n";
		ofile << "the number of sampling points = ";
		ifile >> numSamplingPts;	ofile << "\n" << numSamplingPts << "\n";
		initialPackingFriendly = false;
	}

	/** \brief Compute the local number variance of particle centroids over the ensemble.
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs);

	/** \brief Write the local number variance to "<prefix>_localNumVar".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);
};
//Compute the loca volume-fraction variance of particle phase
class LocalVolumeFractionVariance : public Computation {
public:
	int numSamplingPts;
	double Rmax, dR;
	std::vector<GeometryVector> result;

	/** \brief Construct by reading Rmax, dR and the sampling-point count from a stream.
	 *  \param ifile  Input stream supplying Rmax, dR and numSamplingPts.
	 *  \param ofile  Output stream echoing the prompts. */
	LocalVolumeFractionVariance(std::istream & ifile, std::ostream & ofile) {
		ofile << "Rmax = ";
		ifile >> Rmax;	ofile << "\t" << Rmax << "\n";
		ofile << "dR = ";
		ifile >> dR;	ofile << "\t" << dR << "\n";
		ofile << "the number of sampling points = ";
		ifile >> numSamplingPts;	ofile << "\t" << numSamplingPts << "\n";
		initialPackingFriendly = true;
	}

	/** \brief Compute the local volume-fraction variance of the particle phase over the ensemble.
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs);

	/** \brief Write the local volume-fraction variance to "<prefix>_localVolVar".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);
};

//Compute the deviations in the centroids
class CentroidComputation : public Computation {
public:
	std::vector<std::vector<GeometryVector>> displacement_i;
	std::vector<GeometryVector> result;
	bool output_each_config = false;
	/** \brief Construct by reading whether per-config displacements should be written.
	 *  \param ifile  Input stream supplying the output_each_config flag.
	 *  \param ofile  Output stream echoing the prompts. */
	CentroidComputation(std::istream & ifile, std::ostream & ofile) {
		ofile << "Output displacements in each config? (y=1/n=0)\t";
		ifile >> output_each_config; ofile << "  " << output_each_config << "\n";
	}

	/** \brief Compute the displacements between particles and their Voronoi-cell centroids.
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) {
		std::cout << "Computing diplacements between particles and cell centroids:\n";
		displacement_i.reserve(numConfigs);	result.reserve(numConfigs);
#if defined(USE_VORO)
		for (size_t i = 0; i < numConfigs; i++) {
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

	/** \brief Write the summed displacements (and optionally per-config ones) to "<prefix>_SumDisplacements".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);

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
	/** \brief Copy the current per-config stage-index list out.
	 *  \param[out] list  Receives a copy of the stage-index list. */
	void OutputIndexList(std::vector<std::vector<int>> & list) {
		list = std::vector<std::vector<int>>(this->indices);
	}

	/** \brief Set the per-config stage-index list.
	 *  \param list  Stage-index list to copy in. */
	void SetIndexList(const std::vector<std::vector<int>> & list) {
		this->indices = std::vector<std::vector<int>>(list);
	}

	/** \brief determine stages without any prior knowledge
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	void GetIndexList(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs);

	/** \brief The indices list is competed after the computation
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs) = 0;

	/** \brief Write the multiscale result to disk.
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix) = 0;

};

class SpectralDensity_Multiscale : public MultiscaleComputation {
public:
	int interval;
	double CircularKmax, LinearKmax, KPrecision;
	bool aroundOriginOnly;
	std::vector < std::vector<GeometryVector>> result;
	std::vector<int> outputIndex;
	/** \brief Construct by reading the chi_V(0)-only flag, circular and linear k limits, precision, and stage interval.
	 * \param ifile Input stream supplying the parameters.
	 * \param ofile Output stream echoing the prompts. */
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

	/** \brief Compute the per-stage spectral density chi_V(k) for the multiscale model.
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs);

	/** \brief Write the per-stage chi_V(k), or chi_V(0) vs stage, under OutputPrefix.
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);
};

class VolumeFraction_Multiscale : public MultiscaleComputation {
public:
	int interval;// , maxstep;
	std::vector<GeometryVector> result;

	/** \brief Construct by reading the stage interval from a stream.
	 *  \param ifile  Input stream supplying the stage interval.
	 *  \param ofile  Output stream echoing the prompts. */
	VolumeFraction_Multiscale(std::istream & ifile, std::ostream & ofile) {
		ofile << "Interval in stages = \t";
		ifile >> interval;	ofile << "\t" << interval << "\n";
		//ofile << "max in stages = \n";
		//ifile >> maxstep;	ofile << "\t" << maxstep << "\n";
		this->SetIndexList(std::vector<std::vector<int>>());
	}

	/** \brief Compute the packing fraction at every "interval" stages over the ensemble.
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs);

	/** \brief Write the per-stage packing fractions to "<prefix>_PackingFractions".
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);
};

class VolumeFractionVariance_Multiscale : public MultiscaleComputation {
public:
	int interval;

	double Rmax, dR;
	int numSamplePts;
	std::vector<std::vector<GeometryVector>> result;
	std::vector<int> outputIndex;
	/** \brief Construct by reading the stage interval, Rmax, dR and sample-point count.
	 *  \param ifile  Input stream supplying interval, Rmax, dR and numSamplePts.
	 *  \param ofile  Output stream echoing the prompts. */
	VolumeFractionVariance_Multiscale(std::istream & ifile, std::ostream & ofile) {
		ofile << "Interval in stage = \t";	ifile >> interval; ofile << "\t" << interval << "\n";
		ofile << "R max = \t";	ifile >> Rmax;	ofile << "\t" << Rmax<<"\n";
		ofile << "dR = \t";		ifile >> dR;	ofile << "\t" << dR << "\n";
		ofile << "the number of sample points = \t"; ifile >> numSamplePts; ofile << "\t" << numSamplePts << "\n";
	}

	/** \brief Compute the local volume-fraction variance at every "interval" stages.
	 *  \param GetConfigs  Callable returning configuration i.
	 *  \param numConfigs  Number of configurations in the ensemble. */
	virtual void Compute(std::function<const SpherePacking(size_t)> GetConfigs, size_t numConfigs);

	/** \brief Write the per-stage local volume-fraction variances under OutputPrefix.
	 *  \param OutputPrefix  Output filename prefix. */
	virtual void Write(const std::string OutputPrefix);

};

#endif
