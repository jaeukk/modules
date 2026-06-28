#ifndef HSFluids_INCLUDED
#define HSFluids_INCLUDED

#include <omp.h>
#include <vector>

#include "etc.h"
#include "GeometryVector.h"
#include "PeriodicCellList.h"
#include "RandomGenerator.h"

const double tolerance = 1e-3;
const double error = 8e-3;// -> 1e-2
//Monte Carlo simulations of Hard-sphere fluids
//monodisperse hard spheress
class HSF {

private:
	int ensembleType;
	/*0: canonical ensemble (default) 
	  1: isobaric ensemble */
	std::vector<RandomGenerator> rng;
	double targetAcceptance = 0.3;
	double Diameter;
	bool adjusted, equilibrated;
	bool isMonodisperse = true;

	std::vector<double> pSphereRadii;
public:
	int Dimension;
	size_t NumPrts, NumThreads;
	double PackingFraction, MaxDispl;
	SpherePacking * pSpheres;

	/** \brief Default constructor: build an empty hard-sphere fluid. (CLAUDE) */
	HSF();
	/** \brief Build an initial hard-sphere configuration from a dense lattice packing. (CLAUDE)
	 *  Generate a "Dimension"-dimensional hypercubic simulation box,
	 *  place "NumPrts" particles in the arrangement of densest lattice packings,
	 *  and shrink particle size up to "packingFraction".
	 *  \param Dimension        Spatial dimension of the simulation box.
	 *  \param NumPrts          Number of particles to place.
	 *  \param packingFraction  Target packing fraction after shrinking the particles. */
	HSF(int Dimension, int NumPrts, double packingFraction);
	
	/** \brief Seed every per-thread random number generator. (CLAUDE)
	 *  \param seed  Base seed; generator i is seeded with seed + i*time(NULL). */
	void SetRandomSeed(int seed) {
		for (size_t i = 0; i < this->rng.size(); i++)
			this->rng[i].seed(seed + i*std::time(NULL));
	}
	/** \brief Whether the fluid has reached equilibrium. (CLAUDE)
	 *  \return true if equilibrated. */
	bool IsEquilibrated() { return this->equilibrated; };
	/** \brief Whether the maximum displacement has been adjusted. (CLAUDE)
	 *  \return true if adjusted. */
	bool IsAdjusted() { return this->adjusted; };

	/** \brief Set the number of threads and re-create the per-thread generators. (CLAUDE)
	 *  \param n  Number of threads (must be positive). */
	void SetNumThreads(int n) {
		if (n > 0) {
			this->NumThreads = n;
			this->rng.clear();
			for (int i = 0; i < n; i++) 
				this->rng.emplace_back(RandomGenerator(i));
		}
		else {
			std::cout << "SetNumThreads::invalid input" << std::endl;
		}
	}
	
	/** \brief Replace the current configuration with the given sphere packing and reset derived quantities. (CLAUDE)
	 *  \param newPacking  Sphere packing to copy in as the new configuration. */
	void ReplaceConfig(const SpherePacking & newPacking);// { this->pSpheres = new SpherePacking(newPacking); }
	/** \brief Select the statistical ensemble (NVT or NPT). (CLAUDE)
	 *  \param n  Ensemble type: 0 = canonical (NVT), 1 = isobaric (NPT, not yet implemented). */
	void SetEnsembleType(int n) {
		if (n < 2) {
			ensembleType = n;
			if (n == 0) { std::cout << "NVT (canonical) ensemble" << std::endl; }
			else { std::cout << "NPT (isobaric) ensemble \n not ready yet. Please input 0" << std::endl; }
		}
		else {
			std::cout << "SetEnsembleType::Invalid input" << std::endl;
		}
	}
	/** \brief change the max displacement.
	 *  \param delta  New maximum trial displacement. */
	void SetMaxDisp(double delta) { this->MaxDispl = delta; };
	/** \brief rescale the configuration to be at number density rho.
	 *  \param rho  Target number density. */
	void RescaleConfig(double rho);
	/** \brief run for a single MC cycle (=1 trial move/particle).
	 *  \return acceptance rate. */
	double Evolve();
	/** \brief Run a single MC cycle serially (one trial move per particle). (CLAUDE)
	 *  \return acceptance rate of the trial moves. */
	double Evolve_Serial();
	/** \brief Run a single MC cycle using parallel threads. (CLAUDE)
	 *  \return acceptance rate of the trial moves. */
	double Evolve_Parallel();
	/** \brief run for [cycles] MC cycles (=1 trial move/particle).
	 *  \param cycles  Number of MC cycles to run.
	 *  \return acceptance rate. */
	inline double Evolve(int cycles) {
		double accep = 0.0;
		for (int i = 0; i < cycles; i++)
			accep += Evolve();
		return accep / (double)cycles;
	}
	/** \brief Adjust the max displacement to the target acceptance ratio.
	 *  \return the optimized maximum displacement. */
	double AdjustMaxDisp();
	/** \brief Equilibrate the hard-sphere fluids from the initial configuration.
	 *  \return the number of MC cycles. */
	int Equilibrate();
	/** \brief Calculate the pressure from pair correlation functions.
	 *  \return the reduced pressure of the hard-sphere fluid. */
	double Pressure();
	/** \brief Prepare the system before the production period.
	 *  Initialize the class before starting this function!! */
	void Prepare();
	/** \brief Run this method if the packing is polydisperse. */
	void Switch2Polydisperse(){
		this->isMonodisperse = false;
		for (size_t i = 0; i<this->pSpheres->NumParticle(); i++){
			this->pSphereRadii.push_back(this->pSpheres->GetCharacteristics(i));
		}
		std::cout << "Now, this->Diameter is the maxium particle diameter!!\n";
	};
protected:
	/** \brief Draw a uniform random trial displacement bounded by the maximum displacement. (CLAUDE)
	 *  \param[out] dx  Displacement vector to fill (its dimension is set to the box dimension).
	 *  \param rg       Random number generator used to draw the displacement. */
	inline void TrialMove(GeometryVector &dx, RandomGenerator & rg) {
		dx.Dimension = this->Dimension;
		for (int i = 0; i < Dimension; i++) {
			dx.x[i] = MaxDispl * (rg.RandomDouble() - 0.5);
		}
	}
	/** \brief Calculate g2(D+).
	 *  \return the value of the pair correlation function at contact. */
	double g2_at_D();
	/** \brief Build a lattice packing then thin and shrink it to match the requested parameters. (CLAUDE)
	 *  Generate lattice packing of side length L = ceil(pow(numPrt, 1/d)),
	 *  then some particles are removed and particle radius are reduced to satisfy the given parameters.
	 *  If it's not possible to satisfy the input parameters, it returns an empty packing, and output corresponding message.
	 *  \param[in] numPrts          # of particles in the packing.
	 *  \param[in] packingFraction  target packing fraction.
	 *  \param[in] lattice          lattice vectors.
	 *  \param[out] config          resulting sphere packing. */
	void LatticePacking(int numPrts, double packingFraction, std::vector<GeometryVector> &lattice, SpherePacking & config);
	/** \brief Build a lattice packing from a named lattice, then thin and shrink it. (CLAUDE)
	 *  \param[in] numPrts          # of particles in the packing.
	 *  \param[in] packingFraction  target packing fraction.
	 *  \param[in] lattice          lattice name (Capitalize the first letter).
	 *  \param[out] config          resulting sphere packing. */
	void LatticePacking(int numPrts, double packingFraction, const std::string & lattice, SpherePacking & config);
};

/** \brief Generate configurations of Equilibrium hard-sphere fluids of unit number density.
 *  Pressure is used to check the system is reached into thermal equilibrium.
 *  \param[in] d            spatial dimension.
 *  \param[in] N            the number of particles.
 *  \param[in] phi          packing fraction.
 *  \param[in] saveInterval after equilibration, one configuration is saved in every [saveInterval] MC cycles.
 *  \param[in] numOutput    save [numOuput] configuration at the end.
 *  \param[out] cList       a list of saved configurations!
 *  \param numThreads       Number of threads used during equilibration and sampling. */
inline void HSF_NVT(int d, int N, double phi, int saveInterval, int numOutput, std::vector<SpherePacking> & cList, int numThreads = 1) {
	cList.clear();	cList.reserve(numOutput);
	HSF NVT(d, N, phi);	NVT.SetNumThreads(numThreads);	
	NVT.RescaleConfig(1.0);	//rescale the configuration as the unit number density
	NVT.Prepare();
	
	for (int i = 0; i < numOutput; i++) {
		cList.emplace_back(SpherePacking(* NVT.pSpheres));
		NVT.Evolve(saveInterval);
	}
}

#endif