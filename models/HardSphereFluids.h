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

	/*  */
	HSF();
	/*
		Generate a "Dimension"-dimensional hypercubic simulation box
		Place "NumPrts" particles in the arrangement of densest lattice packings
		Shrink particle size up to "packingFraction"
	*/
	HSF(int Dimension, int NumPrts, double packingFraction);
	
	void SetRandomSeed(int seed) {	
		for (size_t i = 0; i < this->rng.size(); i++)
			this->rng[i].seed(seed + i*std::time(NULL));
	}
	bool IsEquilibrated() { return this->equilibrated; };
	bool IsAdjusted() { return this->adjusted; };

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
	
	void ReplaceConfig(const SpherePacking & newPacking);// { this->pSpheres = new SpherePacking(newPacking); }
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
	//change the max displacement
	void SetMaxDisp(double delta) { this->MaxDispl = delta; };
	//rescale the configuration to be at number density rho 
	void RescaleConfig(double rho);
	/*	run for a single MC cycle (=1 trial move/particle)
		return acceptance rate;	
	*/
	double Evolve();
	double Evolve_Serial();
	double Evolve_Parallel();
	/*	run for [cycles] MC cycles (=1 trial move/particle)
		return acceptance rate;
	*/
	inline double Evolve(int cycles) {
		double accep = 0.0;
		for (int i = 0; i < cycles; i++)
			accep += Evolve();
		return accep / (double)cycles;
	}
	/*	Adjust the max displacement to the target acceptance ratio*/
	double AdjustMaxDisp();
	/*	Equilibrate the hard-sphere fluids from the initial configuration
		return the number of MC cycles
	*/
	int Equilibrate();
	/*	Calculate the pressure from pair correlation functions */
	double Pressure();
	/*  Prepare the system before the production period
		Initialize the class before starting this function!!
	*/
	void Prepare();
	/* Run this method if the packing is polydisperse. */
	void Switch2Polydisperse(){
		this->isMonodisperse = false;
		for (int i = 0; i<this->pSpheres->NumParticle(); i++){
			this->pSphereRadii.push_back(this->pSpheres->GetCharacteristics(i));
		}
		std::cout << "Now, this->Diameter is the maxium particle diameter!!\n";
	};
protected:
	inline void TrialMove(GeometryVector &dx, RandomGenerator & rg) {
		dx.Dimension = this->Dimension;
		for (int i = 0; i < Dimension; i++) {
			dx.x[i] = MaxDispl * (rg.RandomDouble() - 0.5);
		}
	}
	/*	Calculate g2(D+)*/
	double g2_at_D();
	/*	param[in] numPrts = # of particles in the packing
		param[in] packingFraction = target packing fraction
		param[in] lattice = lattice vectors
			Generate lattice packing of side length L = ceil(pow(numPrt, 1/d))
			Then, some particles are removed and particle radius are reduced to satisfy the given parameters
			If it's not possible to satisfy the input parameters, it returns an empty packing, and output corresponding message.
	*/
	void LatticePacking(int numPrts, double packingFraction, std::vector<GeometryVector> &lattice, SpherePacking & config);
	/* param[in] numPrts = # of particles in the packing
		param[in] packingFraction = target packing fraction
		param[in] lattice = lattice name (Capitalize the first letter)
		*/
	void LatticePacking(int numPrts, double packingFraction, const std::string & lattice, SpherePacking & config);
};

/*	Generate configurations of Equilibrium hard-sphere fluids of unit number density 
	Pressure is used to check the system is reached into thermal equilibrium
	param[in] d				spatial dimension
	param[in] phi			packing fraction
	param[in] N				the number of particles
	param[in] saveInterval	after equilibration, one configuration is saved in every [saveInterval] MC cycles
	param[in] numOutput		save [numOuput] configuration at the end
	param[out] cList		a list of saved configurations!
*/
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