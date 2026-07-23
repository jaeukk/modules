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
#include <PeriodicCellList.h>


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

	/** \brief Construct simulated annealing from a random initial configuration.
	 *	\param Dimension Spatial dimension.
	 *	\param NumParticle Number of particles.
	 *	\param InitialDensity Initial number density.
	 *	\param RandomSeed Random-number seed.
	 *	\param Pressure Pressure used in the enthalpy; set CellMoveProbability to zero after construction for constant volume.
	 *	\param pot Interaction potential.
	 *	\param output Optional logging stream.
	 *	\param CoolingScheduleRescale Greater than one cools more slowly; less than one cools faster.
	 *	\param MinDistance Minimum initial and trial interparticle distance.
	 *	\param AlgorithmLevel Zero selects final structure optimization only, one selects Monte Carlo, and values greater than one select parallel tempering with that many replicas.
	 *	\param UseSortedList Passed to the newly constructed Configuration to select its sorted-list representation.
	 */
	ParticleSimulatedAnnealing(DimensionType Dimension, size_t NumParticle, double InitialDensity, 
		int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, 
		double CoolingScheduleRescale = 1, double MinDistance=0, int AlgorithmLevel = 1, bool UseSortedList=false);
	
	//AlgorithmLevel : 0:Structure Optimization, 1:Monte Carlo, >1:Parallel Tempering with n copies
	
	/** \brief Construct simulated annealing from an existing configuration.
	 *	\param pInitConfig Initial configuration, copied internally.
	 *	\param RandomSeed Random-number seed.
	 *	\param Pressure Pressure used in the enthalpy; set CellMoveProbability to zero after construction for constant volume.
	 *	\param pot Interaction potential.
	 *	\param output Optional logging stream.
	 *	\param CoolingScheduleRescale Greater than one cools more slowly; less than one cools faster.
	 *	\param MinDistance Minimum permitted trial interparticle distance.
	 *	\param AlgorithmLevel Zero selects final structure optimization only, one selects Monte Carlo, and values greater than one select parallel tempering with that many replicas.
	 *	\param UseSortedList Compatibility parameter; ignored because pInitConfig is copied with its existing representation.
	 */
	ParticleSimulatedAnnealing(const Configuration & pInitConfig, 
		int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, 
		double CoolingScheduleRescale = 1, double MinDistance=0, int AlgorithmLevel = 1, bool UseSortedList=false);
	
	
	/** \brief Run the simulated-annealing (or optimization / parallel-tempering) procedure.
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
