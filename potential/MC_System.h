/**	Author	: Ge Zhang
 *	Email	: 
 *	Date	: 2015 
 *  Annotated by Jaeuk Kim (2019). */

/** \file	MC_System.h
	\brief	Header file to perform the simulated annealing */

#ifndef MCSYSTEM_INCLUDED
#define MCSYSTEM_INCLUDED

#include <functional>
#include "MonteCarlo.h"
#include "PeriodicCellList.h"
#include "Move.h"


//construct this class, and set every data member if you want, then call Anneal().
//after that, get the configuration by COPYING (*pList)
/** \brief	A class for performing the simulated annealing.
 *	How to use: construct this class, and set every data member if you want, then call Anneal().
	after that, get the configuration by COPYING (*pList). */
class ParticleSimulatedAnnealing
{
private:
	MonteCarlo<Configuration, Potential *> * pSys;	

public:
	double alpha;	//< a base for a power-law cooling schedule; T = T0 * alpha^n.
	double Tcut;	//< finial temperature. 10^{-6} * pow(0.1, Dimension) by default.		
	double Tinit;	//< initial temperature. 10^5 by default.
	double CellMoveProbability;
	double Pressure;
	RandomGenerator gen;	//< A random number generator.
	Configuration * pList;	//< A point configuration of interest.	After the annealing, this will be a ground-state.
	std::ostream * Output;
	double MinDistance;	//< Particle exclusion diamter. abandon this move immediately if particles are within this distance
	static bool StopSignal;
	std::function<void(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance)> RelaxStructureFunction;

	/** \brief A constructor from a random initial condition.
	 *	\param Dimension		Spatial dimension of the initial point configuration.
	 *	\param NumParticle		Number of particles in the initial point configuration.
	 *	\param InitialDensity	Number density of the initial point configuration; a random packing with a radius "MinDistance" is generated and stored in this->pList.
	 *	\param RandomSeed	A seed for a random generator.
	 *	\param Pressure		Pressure under which an isobaric simulation is performed. Set it to 0 for a constant volume.
	 *	\param pot			A reference of a Potential object, which defines the interaction in this system.
	 *	\param output		An outputstream
	 *	\param CoolingScheduleRescale	>1: slower cooling, <1: faster cooling.
	 *	\param MinDistance		Exclusion diameter in the initial configuraiton.
	 *	\param AlgorithmLevel	0:  Structure Optimization.
								1:  Monte Carlo simulation.
								>1: Parallel tempering with this number of copy.
	 *	\param UseSortedList	TODO(CLAUDE): unclear - review
	 */
	ParticleSimulatedAnnealing(DimensionType Dimension, size_t NumParticle, double InitialDensity, 
		int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, 
		double CoolingScheduleRescale = 1, double MinDistance=0, int AlgorithmLevel = 1, bool UseSortedList=false);
	
	//AlgorithmLevel : 0:Structure Optimization, 1:Monte Carlo, >1:Parallel Tempering with n copies
	
	/** \brief A constrcutor from an initial configuration.
	 *	\param pInitConfig	An initial point configuraiton. It is copied into this->pList.
	 *	\param RandomSeed	A seed for a random generator.
	 *	\param Pressure		Pressure under which an isobaric simulation is performed. Set it to 0 for a constant volume. (CLAUDE)
	 *	\param pot			A reference of a Potential object, which defines the interaction in this system.
	 *	\param output		An outputstream
	 *	\param CoolingScheduleRescale	>1: slower cooling, <1: faster cooling.
	 *	\param MinDistance		Exclusion diameter in the initial configuraiton.
	 *	\param AlgorithmLevel	0:  Structure Optimization.
								1:  Monte Carlo simulation.
								>1: Parallel tempering with this number of copy.
	 *	\param UseSortedList	TODO(CLAUDE): unclear - review
	 */
	ParticleSimulatedAnnealing(const Configuration & pInitConfig, 
		int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, 
		double CoolingScheduleRescale = 1, double MinDistance=0, int AlgorithmLevel = 1, bool UseSortedList=false);
	
	
	/** \brief Run the simulated-annealing (or optimization / parallel-tempering) procedure. (CLAUDE)
	 *	\param pot  Potential defining the interaction; drives the energy/force evaluations. */
	void Anneal(Potential & pot);

	/** \brief A destructor. */
	~ParticleSimulatedAnnealing()
	{
		if(this->pList!=nullptr)
			delete this->pList;
		if(this->pSys!=nullptr)
			delete this->pSys;
	}
};



#endif
